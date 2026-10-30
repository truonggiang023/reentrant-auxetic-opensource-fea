"""Experimental OpenCourant export in the native Radioss input format.

The generated model is a second solver formulation, not a translation of the
custom CalculiX UMAT. In particular INIEVO's native characteristic length,
integration at initiation and rate-independent evolution differ from that UMAT.
No calibration or quasi-static qualification follows from deck generation.
Card layouts were checked against OpenRadioss source at
bad6688602d770e019b75773d2d5a48b4f9bd918 (hm_cfg_files/config/CFG).
"""
from pathlib import Path
import hashlib
import json

import numpy as np

from .plastic import validate_material


def _ints(*values):
    return ''.join(f'{int(x):10d}' if x is not None else ' '*10 for x in values)


def _reals(*values):
    return ''.join(f'{float(x):20.12g}' for x in values)


def _group(lines, ident, name, ids):
    lines.extend([f'/GRNOD/NODE/{ident}', name])
    ids = list(ids)
    if not ids:
        raise ValueError('Empty node group: ' + name)
    for start in range(0, len(ids), 10):
        lines.append(_ints(*ids[start:start+10]))


def _bcs(lines, ident, name, group, translation):
    lines.extend([f'/BCS/{ident}', name,
                  f'{translation:>6} 000' + _ints(0, group)])


def _motion(lines, ident, group, amplitude):
    lines.extend([f'/IMPDISP/{ident}', 'Prescribed Y displacement',
                  _ints(20)+'         Y'+_ints(0, 0, group, 0, 0),
                  _reals(1, amplitude, 0, 0)])


def _contact(lines, ident, secondary, main, friction):
    # All-face sets retain internal solid facets for erosion. Idel=1 removes
    # a segment once every solid adjoining it has failed. Edge contact covers
    # core self-impact. Broad flat platens use face contact: including their
    # coplanar edge pairs duplicates the initial end-face contact and generates
    # artificial edge penetrations. Runtime qualification is still required.
    edge_mode = 22 if main == 0 else 1000
    lines.extend([f'/INTER/TYPE25/{ident}', f'Contact {ident}',
                  _ints(secondary, main, 5, 0, 3, 3, None, 1, edge_mode),
                  _ints(0, None)+_reals(0, 0, 0, 0),
                  _reals(0, 0)+_ints(1, 2)+_reals(0),
                  _reals(1, friction, 0, 0, 0),
                  '       000'+_ints(None, 0, 5)+_reals(0)+_ints(0, None)+_reals(0),
                  _ints(0, 0)+_reals(0)+_ints(None, 0, None, None, None, 0)])


def _material(lines, config, plastic, damage):
    lines.extend(['/MAT/PLAS_TAB/1', 'AlSi10Mg tabulated rate-independent plasticity',
                  _reals(config.get('density', 2.67e-9)),
                  _reals(config['young'], config['poisson'], 1e30, 1e30, 1e30),
                  _ints(1, 0)+_reals(0, 0, 1e30), _ints(0)+_reals(1),
                  _ints(10), _reals(1), _reals(0),
                  '/FUNCT/10', 'True flow stress versus plastic strain'])
    lines.extend(_reals(strain, stress) for stress, strain in plastic)
    # FUNCT extrapolates its last segment. The CCX/reference material clamps
    # stress beyond its final supplied plastic strain, so append a constant
    # boundary segment to encode that convention. This is not a measured
    # hardening point and does not change interpolation within the input table.
    last_stress, last_strain = plastic[-1]
    lines.append(_reals(last_strain + max(1., abs(last_strain)), last_stress))
    if damage is None:
        return
    locus = np.asarray(damage.get('fracture_locus', []), float)
    if (locus.ndim != 2 or locus.shape[1] != 2 or len(locus) < 2
            or not np.isfinite(locus).all() or np.any(np.diff(locus[:, 0]) <= 0)
            or np.any(locus[:, 1] <= 0)):
        raise ValueError('Require increasing finite triaxiality locus with positive fracture strains')
    uf = float(damage.get('failure_displacement', 0))
    if not np.isfinite(uf) or uf <= 0:
        raise ValueError('failure_displacement must be finite and positive')
    lines.extend(['/TABLE/1/11', 'Fracture strain versus triaxiality (rate independent)',
                  _ints(1)])
    # TABLE also extrapolates its final/first segments. Preserve the explicitly
    # assumed CCX endpoint clamping so high eta cannot create negative fracture
    # strains. These constant extensions encode a boundary convention, not
    # additional located fracture measurements from the paper.
    lines.append(_reals(locus[0, 0]-max(1., abs(locus[0, 0])), locus[0, 1]))
    lines.extend(_reals(eta, eps) for eta, eps in locus)
    lines.append(_reals(locus[-1, 0]+max(1., abs(locus[-1, 0])), locus[-1, 1]))
    # NINIEVO=1; ILEN=1 initial critical-time-step characteristic length;
    # FAILIP=1. This differs from the CCX UMAT's cube-root-volume length.
    lines.extend(['/FAIL/INIEVO/1', _ints(1, 0, 1, None, None, None, None, 1)+_reals(0),
                  _ints(1, 1, 1, 1), _ints(11)+_reals(1, 1, 0),
                  _ints(0)+_reals(1, 1), _reals(uf, 1, 0)])


def write_radioss_deck(path, nodes, elements, config, plastic, *,
                       duration, drop_ccx_viscosity=False, tensile_coupon=False,
                       ams_time_step=None):
    """Write a native explicit deck and metadata; return the metadata.

    ``duration`` is in seconds and must be chosen explicitly. No mass scaling
    or additional Rayleigh damping is inserted. Native solid artificial
    viscosity and contact-damping defaults remain active and are recorded.
    Nonzero CCX pseudo-time viscosity is
    refused unless its deliberate omission is acknowledged and recorded.
    ``tensile_coupon`` bonds the ends of the supplied mesh for qualification;
    the default uses separate elastic steel platens and sliding TYPE25 contact.
    Only 8-node bricks are accepted; wedge-to-brick degeneration is not hidden.
    ``ams_time_step`` explicitly enables native Advanced Mass Scaling. It
    modifies the high-frequency inertial operator and requires independent
    time-step, loading-speed and energy sensitivity checks. Preserved total
    mass alone does not prove that an AMS result is quasi-static.
    """
    nodes = np.asarray(nodes, float)
    elements = np.asarray(elements)
    if (nodes.ndim != 2 or nodes.shape[1] != 3 or not np.isfinite(nodes).all()
            or len(nodes) < 8 or elements.ndim != 2 or elements.shape[1] != 8
            or not np.issubdtype(elements.dtype, np.integer) or len(elements) < 1
            or elements.min() < 0 or elements.max() >= len(nodes)
            or any(len(set(el)) != 8 for el in elements)):
        raise ValueError('Require finite 3D nodes and nondegenerate 8-node brick connectivity')
    duration = float(duration)
    strain = float(config['compression_strain'])
    if not np.isfinite(duration) or duration <= 0 or not np.isfinite(strain) or not 0 < strain < 1:
        raise ValueError('Require positive duration and strain between zero and one')
    anchor_requested = config.get('core_x_anchor', False)
    if not isinstance(anchor_requested, bool):
        raise ValueError('core_x_anchor must be an explicit boolean')
    anchor_node = None
    if ams_time_step is not None:
        ams_time_step = float(ams_time_step)
        if not np.isfinite(ams_time_step) or not 0 < ams_time_step <= duration/100:
            raise ValueError('AMS time step must be finite, positive and at most duration/100')
    plastic = validate_material(config, plastic)
    damage = config.get('ductile_damage')
    viscosity = float(damage.get('viscosity_time', 0)) if damage else 0
    if not np.isfinite(viscosity) or viscosity < 0:
        raise ValueError('Invalid viscosity_time')
    if viscosity and not drop_ccx_viscosity:
        raise ValueError('CCX pseudo-time viscosity cannot be mapped to physical seconds; explicitly drop it')
    friction = float(config.get('friction', .35))
    self_friction = float(config.get('self_friction', friction))
    density = float(config.get('density', 2.67e-9))
    if not np.isfinite([friction, self_friction, density]).all() or min(friction, self_friction) < 0 or density <= 0:
        raise ValueError('Invalid friction or density')
    xmin, ymin, zmin = nodes.min(axis=0)
    xmax, ymax, zmax = nodes.max(axis=0)
    height = float(ymax-ymin)
    if min(xmax-xmin, height, zmax-zmin) <= 0:
        raise ValueError('Specimen must have nonzero dimensions')
    path = Path(path)
    if path.name.endswith('_0000.rad'):
        prefix = path.name[:-9]
    elif path.suffix == '.rad':
        prefix = path.stem
    else:
        prefix = path.name
    if not prefix or len(prefix) > 60 or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in prefix):
        raise ValueError('Use a short ASCII alphanumeric model prefix')
    parent = path.parent
    parent.mkdir(parents=True, exist_ok=True)
    starter = parent/(prefix+'_0000.rad')
    engine = parent/(prefix+'_0001.rad')
    meta_path = parent/(prefix+'.radioss.json')
    if any(p.exists() for p in (starter, engine, meta_path)):
        raise FileExistsError('Use a new prefix; existing decks are preserved')
    lines = ['#RADIOSS STARTER', '/BEGIN', prefix, _ints(2022, 0),
             f'{"Mg":>20}{"mm":>20}{"s":>20}',
             f'{"Mg":>20}{"mm":>20}{"s":>20}',
             '/TITLE', 'Uncalibrated research candidate: native erosion and TYPE25 contact']
    if ams_time_step is not None:
        # /AMS group=0 selects all parts. Engine /DT/AMS/2 uses the native
        # non-diagonal mass operator and PCG, without adding /DT/NODA/CST.
        lines.extend(['/AMS', _ints(0)])
    _material(lines, config, plastic, damage)
    allnodes = nodes.tolist()
    allelements = [(i+1, 1, tuple(int(v)+1 for v in el)) for i, el in enumerate(elements)]
    ncore = len(nodes)
    if not tensile_coupon:
        margin = float(config.get('platen_margin', 1))
        thick = float(config.get('platen_thickness', 10))
        gap = float(config.get('platen_gap', 0))
        if not np.isfinite([margin, thick, gap]).all() or margin < 0 or thick <= 0 or gap < 0:
            raise ValueError('Invalid platen dimensions')
        for part, low, high in [(2, ymin-gap-thick, ymin-gap), (3, ymax+gap, ymax+gap+thick)]:
            start = len(allnodes)+1
            allnodes.extend([(xmin-margin, low, zmin), (xmax+margin, low, zmin),
                             (xmax+margin, high, zmin), (xmin-margin, high, zmin),
                             (xmin-margin, low, zmax), (xmax+margin, low, zmax),
                             (xmax+margin, high, zmax), (xmin-margin, high, zmax)])
            allelements.append((len(allelements)+1, part, tuple(range(start, start+8))))
        lines.extend(['/MAT/ELAST/2', 'Steel', _reals(config.get('steel_density', 7.85e-9)),
                      _reals(config.get('steel_young', 210000), config.get('steel_poisson', .3))])
    coupled = bool(config.get('platen_reference_points', False))
    if coupled and tensile_coupon:
        raise ValueError('Platen RP coupling requires steel platens')
    if coupled:
        bottom_rp, top_rp = len(allnodes)+1, len(allnodes)+2
        allnodes.extend([(float((xmin+xmax)/2), ymin-gap-thick, float((zmin+zmax)/2)),
                         (float((xmin+xmax)/2), ymax+gap+thick, float((zmin+zmax)/2))])
    lines.append('/NODE')
    lines.extend(_ints(i)+_reals(*xyz) for i, xyz in enumerate(allnodes, 1))
    # HEPH (Isolid=24) has physical hourglass stabilization. It is not claimed
    # to reproduce Abaqus C3D8R or CCX hourglass forces exactly.
    lines.extend(['/PROP/SOLID/1', 'HEPH eight-node solid',
                  _ints(24, 0, None, 0, None, 0, 0, 0)+_reals(0),
                  _reals(0, 0, 0), _reals(0)+_ints(1)])
    for part in range(1, 2 if tensile_coupon else 4):
        lines.extend([f'/PART/{part}', f'Part {part}', _ints(1, 1 if part == 1 else 2)])
        lines.append(f'/BRICK/{part}')
        lines.extend(_ints(eid, *el) for eid, p, el in allelements if p == part)
    tol = max(1e-8, 1e-8*height)
    if tensile_coupon:
        lower = [i+1 for i, n in enumerate(nodes) if abs(n[1]-ymin) < tol]
        upper = [i+1 for i, n in enumerate(nodes) if abs(n[1]-ymax) < tol]
    else:
        lower = [ncore+i for i in (1, 2, 5, 6)]
        upper = [ncore+i for i in (11, 12, 15, 16)]
    _group(lines, 1, 'Bottom loading face', lower)
    _group(lines, 2, 'Top loading face', upper)
    if coupled:
        for ident, rp in ((1, bottom_rp), (2, top_rp)):
            lines.extend([f'/RBE2/{ident}', 'Outer platen face to reference point',
                          _ints(rp)+'   111 111'+_ints(0, ident, 0)])
        _group(lines, 11, 'Bottom reference point', [bottom_rp])
        _group(lines, 12, 'Top reference point', [top_rp])
        lower, upper = [bottom_rp], [top_rp]
        lines.extend(['/BCS/11', 'Bottom RP fixed', '   111 111'+_ints(0,11),
                      '/BCS/12', 'Top RP guided', '   101 111'+_ints(0,12)])
    else:
        _bcs(lines, 1, 'Fixed bottom loading face', 1, '111')
        _bcs(lines, 2, 'Top face XZ restraint', 2, '101')
    if not tensile_coupon:
        _group(lines, 3, 'Core plane-strain Z restraint', range(1, ncore+1))
        _bcs(lines, 3, 'Core Z restraint', 3, '001')
        # Explicit dynamics needs no stiffness-matrix rigid-body gauge. Leave
        # in-plane motion governed by contact unless a legacy sensitivity case
        # explicitly requests the former interior X pin.
        if anchor_requested:
            anchor_node = min(range(len(nodes)), key=lambda i: abs(nodes[i, 0]-(xmin+xmax)/2)) + 1
            _group(lines, 4, 'Legacy interior X roller', [anchor_node])
            _bcs(lines, 4, 'Legacy in-plane pin (not a paper boundary)', 4, '100')
    # Quintic smoothstep displacement: zero velocity/acceleration at both ends.
    lines.extend(['/FUNCT/20', 'Smooth displacement fraction versus seconds'])
    amplitude_samples = [[float(duration*fraction), float(fraction**3*(10-15*fraction+6*fraction**2))]
                         for fraction in np.linspace(0, 1, 101)]
    for time, fraction in amplitude_samples:
        lines.append(_reals(time, fraction))
    _motion(lines, 1, 12 if coupled else 2, height*strain*(1 if tensile_coupon else -1))
    if not tensile_coupon:
        platen_surface = config.get('platen_contact_surface', 'all')
        if platen_surface not in ('all', 'inner_y'):
            raise ValueError('platen_contact_surface must be all or inner_y')
        for part in (1, 2, 3):
            if part == 1 or platen_surface == 'all':
                lines.extend([f'/SURF/PART/ALL/{part}', f'All faces of part {part}', _ints(part)])
            else:
                offsets = (3, 4, 8, 7) if part == 2 else (9, 10, 14, 13)
                lines.extend([f'/SURF/SEG/{part}', f'Inner Y platen face {part}',
                              _ints(part, *(ncore+i for i in offsets))])
        _contact(lines, 1, 1, 0, self_friction)
        _contact(lines, 2, 1, 2, friction)
        _contact(lines, 3, 1, 3, friction)
        lines.extend(['/TH/INTER/1', 'Contact force and energy', 'DEF', _ints(1, 2, 3)])
    lines.extend(['/TH/PART/2', 'Part kinetic internal and hourglass energy', 'DEF',
                  _ints(*range(1, 2 if tensile_coupon else 4)),
                  '/TH/NODE/3', 'Loading face nodal histories',
                  ''.join(f'{variable:<10}' for variable in ('DX', 'DY', 'DZ', 'REACX', 'REACY', 'REACZ'))])
    support_history = bool(config.get('support_reaction_history', False)) or coupled
    history_nodes = upper+lower if support_history else upper
    lines.extend(_ints(node, 0)+f'load-node-{node}' for node in history_nodes)
    lines.append('/END')
    starter.write_text('\n'.join(lines)+'\n')
    engine_lines = ['#RADIOSS ENGINE', f'/RUN/{prefix}/1', _reals(duration),
                    '/PRINT/-1000', '/TFILE/4', _reals(duration/400),
                    '/ANIM/DT', _reals(0, duration/100),
                    '/ANIM/BRICK/TENS/STRESS', '/ANIM/BRICK/EPSP',
                    '/ANIM/BRICK/DAMA', '/ANIM/VECT/DISP']
    if ams_time_step is not None:
        engine_lines.extend(['/DT/AMS/2', _reals(.67, ams_time_step),
                             _reals(1e-4), _ints(1000, 0)])
    engine.write_text('\n'.join(engine_lines)+'\n')
    metadata = dict(backend='OpenCourant / Radioss format experimental', prefix=prefix,
                    starter=str(starter.resolve()), engine=str(engine.resolve()),
                    units=dict(mass='Mg', length='mm', time='s', stress='MPa', force='N'),
                    duration_seconds=duration, height_mm=height, wall_mm=config.get('wall'),
                    displacement_fraction_samples=amplitude_samples,
                    time_history_interval_seconds=duration/400,
                    loading_face_nodes=upper, bottom_support_nodes=lower,
                    platen_reference_points=coupled, support_reaction_history=support_history,
                    loading_node_variables=['DX', 'DY', 'DZ', 'REACX', 'REACY', 'REACZ'],
                    reaction_history_semantics='REAC and interface force channels are cumulative impulses; differentiate in time to obtain interval-averaged forces',
                    interface_force_history_units='N*s',
                    core_nodes=ncore, core_elements=len(elements), target_nominal_strain=strain,
                    reference_area_mm2=config.get('reference_area', (xmax-xmin)*(zmax-zmin)),
                    physical_speed_mm_per_min=2., physical_duration_seconds=height*strain/2*60,
                    peak_imposed_speed_mm_per_s=1.875*height*strain/duration,
                    acceleration_requires_energy_and_speed_sensitivity=True,
                    mass_scaling='native Advanced Mass Scaling' if ams_time_step is not None else False,
                    ams=dict(enabled=ams_time_step is not None,
                             requested_time_step_seconds=ams_time_step,
                             time_step_scale_factor=.67,
                             pcg_relative_tolerance=1e-4,
                             pcg_maximum_iterations=1000,
                             timestep_sensitivity_required=ams_time_step is not None,
                             modifies_inertial_operator=ams_time_step is not None),
                    numerical_damping=dict(additional_rayleigh_requested=False,
                                           native_solid_defaults_active=True,
                                           starter_reported_hourglass_coefficient=.1,
                                           starter_reported_solid_qa=1.1,
                                           starter_reported_solid_qb=.05,
                                           starter_reported_contact_critical_damping=.05,
                                           starter_reported_contact_quadratic_damping=1),
                    completion_verified=False,
                    calibrated=False, experimental_agreement=False,
                    native_erosion_requested=bool(damage),
                    all_face_self_contact_requested=not tensile_coupon,
                    core_self_contact_edges='TYPE25 Iedge=22 (all exterior solid edges)',
                    platen_contact_edges='TYPE25 Iedge=1000 (face contact on broad flat platens)',
                    contact_qualification='not yet established', element_formulation='HEPH Isolid=24',
                    native_characteristic_length='INIEVO ILEN=1 initial critical-time-step formulation',
                    ccx_viscosity_omitted=viscosity, fracture_locus=damage.get('fracture_locus') if damage else None,
                    material_provenance=damage.get('provenance') if damage else None,
                    limitations=['Different solver formulation; requires independent qualification',
                                 'Fracture-locus endpoint coordinates remain assumptions',
                                 'INIEVO and the CCX UMAT differ in initiation increment integration',
                                 'No CCX compression cutoff or viscosity silently transferred'])
    metadata['sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (starter, engine)}
    metadata['core_x_anchor'] = {
        'enabled': anchor_node is not None, 'node': anchor_node,
        'position_mm': nodes[anchor_node-1].tolist() if anchor_node is not None else None,
        'qualification': 'optional legacy boundary sensitivity; the paper does not specify an interior X pin',
    }
    metadata['platen_contact_surface'] = config.get('platen_contact_surface', 'all')
    metadata['hardening_extrapolation'] = {
        'convention': 'constant flow stress after the last supplied plastic strain, matching the CCX law',
        'original_point_count': len(plastic),
        'last_supplied_point': {'plastic_strain': float(plastic[-1][1]),
                                'true_flow_stress_MPa': float(plastic[-1][0])},
        'boundary_segment_point': {'plastic_strain': float(plastic[-1][1] + max(1., abs(plastic[-1][1]))),
                                   'true_flow_stress_MPa': float(plastic[-1][0])},
        'boundary_segment_is_measured_data': False,
        'reason': 'OpenCourant FUNCT otherwise linearly extrapolates the final supplied slope',
    }
    if damage is not None:
        locus = np.asarray(damage['fracture_locus'], dtype=float)
        metadata['fracture_locus_extrapolation'] = {
            'convention': 'constant fracture strain outside the supplied triaxiality interval, matching the assumed CCX endpoint clamp',
            'original_point_count': len(locus),
            'numerical_boundary_points': [
                [float(locus[0, 0]-max(1., abs(locus[0, 0]))), float(locus[0, 1])],
                [float(locus[-1, 0]+max(1., abs(locus[-1, 0]))), float(locus[-1, 1])]],
            'boundary_points_are_measured_data': False,
            'reason': 'OpenCourant TABLE otherwise extrapolates endpoint slopes and can produce nonpositive fracture strains',
        }
    meta_path.write_text(json.dumps(metadata, indent=2)+'\n')
    return metadata
