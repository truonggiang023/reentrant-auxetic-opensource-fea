"""Standalone Abaqus/Explicit comparison decks; no Abaqus import or execution.

Reuse a retained core mesh to distinguish solver effects from remeshing. Native
Abaqus damage is an optional hypothesis, not the custom CalculiX material.
"""
import hashlib
import json
from pathlib import Path
import tarfile

import numpy as np

from .plastic import validate_material


HEX_FACES = ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))


def retained_core(archive):
    """Read the one native starter in an archive without extracting any paths."""
    archive = Path(archive)
    with tarfile.open(archive, 'r:gz') as stream:
        members = [m for m in stream.getmembers() if m.name.endswith('_0000.rad') and m.isfile()]
        if len(members) != 1:
            raise ValueError('Require exactly one retained native starter in the archive')
        member = members[0]
        raw = stream.extractfile(member).read()
    nodes, elements, mode = {}, {}, ''
    for line in raw.decode('utf-8').splitlines():
        if line.startswith('/'):
            mode = line
        elif mode == '/NODE' and line.strip() and not line.startswith('#'):
            fields = line.split()
            ident = int(fields[0])
            if ident in nodes or len(fields) != 4:
                raise ValueError('Invalid or duplicate native node')
            nodes[ident] = list(map(float, fields[1:]))
        elif mode == '/BRICK/1' and line.strip() and not line.startswith('#'):
            fields = list(map(int, line.split()))
            if len(fields) != 9 or fields[0] in elements:
                raise ValueError('Invalid or duplicate core brick')
            elements[fields[0]] = fields[1:]
    used = sorted({n for row in elements.values() for n in row})
    if not elements or any(n not in nodes for n in used):
        raise ValueError('Missing core mesh')
    nodes = {n: nodes[n] for n in used}
    provenance = dict(archive=archive.name, archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
                      starter_member=member.name, starter_sha256=hashlib.sha256(raw).hexdigest())
    return nodes, elements, provenance


def mesh_audit(nodes, elements):
    """Check brick orientation at all eight corners and all eight Gauss points."""
    signs = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                      [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
    xyz = np.asarray([nodes[n] for n in sorted(nodes)], float)
    if not np.isfinite(xyz).all() or min(nodes) <= 0 or min(elements) <= 0:
        raise ValueError('Require finite coordinates and positive IDs')
    derivatives = []
    for q in np.concatenate((signs, signs / np.sqrt(3))):
        d = np.empty((8, 3))
        for k in range(3):
            other = [j for j in range(3) if j != k]
            d[:, k] = signs[:, k] * np.prod(1 + signs[:, other] * q[other], axis=1) / 8
        derivatives.append(d)
    derivatives = np.asarray(derivatives)
    minimum_jacobian, volume, faces = float('inf'), 0., {}
    for row in elements.values():
        if len(row) != 8 or len(set(row)) != 8:
            raise ValueError('Require eight distinct nodes per brick')
        p = np.asarray([nodes[n] for n in row])
        determinants = np.linalg.det(np.einsum('ij,qjk->qik', p.T, derivatives))
        minimum_jacobian = min(minimum_jacobian, float(determinants.min()))
        volume += float(determinants[8:].sum())
        for face in HEX_FACES:
            key = tuple(sorted(row[i] for i in face))
            faces[key] = faces.get(key, 0) + 1
    if minimum_jacobian <= 0 or max(faces.values()) > 2:
        raise ValueError('Inverted/degenerate brick or nonmanifold mesh')
    return dict(nodes=len(nodes), elements=len(elements), bounds_mm=[xyz.min(axis=0).tolist(), xyz.max(axis=0).tolist()],
                dimensions_mm=np.ptp(xyz, axis=0).tolist(), volume_mm3=float(volume),
                minimum_sampled_jacobian_mm3=minimum_jacobian,
                exterior_facets=sum(v == 1 for v in faces.values()),
                interior_facets=sum(v == 2 for v in faces.values()))


def _nset(lines, name, ids):
    ids = sorted(ids)
    if not ids:
        raise ValueError('Empty node set: ' + name)
    lines.append('*NSET, NSET=' + name)
    lines.extend(', '.join(map(str, ids[i:i+12])) for i in range(0, len(ids), 12))


def _plate(nodes, elements, lower, upper, size):
    """Conforming structured C3D8R steel plate with a retained load face."""
    axes = [np.linspace(lo, hi, int(np.ceil((hi-lo)/size))+1) for lo, hi in zip(lower, upper)]
    first_n, first_e = max(nodes)+1, max(elements)+1
    ids = {}
    plate_nodes, plate_elements, ymin_faces, ymax_faces = [], [], [], []
    for k, z in enumerate(axes[2]):
        for j, y in enumerate(axes[1]):
            for i, x in enumerate(axes[0]):
                n = first_n + len(ids)
                ids[i, j, k] = n
                nodes[n] = [float(x), float(y), float(z)]
                plate_nodes.append(n)
    nx, ny, nz = [len(a)-1 for a in axes]
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                e = first_e + len(plate_elements)
                elements[e] = [ids[i, j, k], ids[i+1, j, k], ids[i+1, j+1, k], ids[i, j+1, k],
                               ids[i, j, k+1], ids[i+1, j, k+1], ids[i+1, j+1, k+1], ids[i, j+1, k+1]]
                plate_elements.append(e)
                if j == 0:
                    ymin_faces.append((e, 'S3'))
                if j == ny-1:
                    ymax_faces.append((e, 'S5'))
    return dict(nodes=plate_nodes, elements=plate_elements, ymin_faces=ymin_faces, ymax_faces=ymax_faces,
                ymin_nodes=[n for (i, j, k), n in ids.items() if j == 0],
                ymax_nodes=[n for (i, j, k), n in ids.items() if j == ny])


def write_comparison(path, core_nodes, core_elements, config, plastic, *, duration=.01,
                     damage=False, loading='smooth', platen_size=2.5):
    """Write one self-contained input and a JSON record; preserve existing files."""
    plastic = validate_material(config, plastic)
    if (not np.isfinite([duration, platen_size]).all() or min(duration, platen_size) <= 0
            or loading not in ('smooth', 'linear')):
        raise ValueError('Require positive finite duration/platen size and smooth or linear loading')
    path = Path(path)
    if path.suffix.lower() != '.inp':
        raise ValueError('Use an .inp output path')
    sidecar = path.with_suffix('.json')
    if path.exists() or sidecar.exists():
        raise FileExistsError('Use a fresh output path; existing inputs are preserved')
    core = mesh_audit(core_nodes, core_elements)
    nodes = {int(k): list(v) for k, v in core_nodes.items()}
    elements = {int(k): list(v) for k, v in core_elements.items()}
    lo, hi = np.asarray(core['bounds_mm'])
    height = float(hi[1]-lo[1])
    displacement = height * float(config['compression_strain'])
    margin, thickness = float(config.get('platen_margin', 1)), float(config.get('platen_thickness', 10))
    gap = float(config.get('platen_gap', 0))
    friction, self_friction = float(config.get('friction', .35)), float(config.get('self_friction', config.get('friction', .35)))
    if (not np.isfinite([margin, thickness, gap, friction, self_friction]).all()
            or min(margin, gap, friction, self_friction) < 0 or thickness <= 0):
        raise ValueError('Invalid plate geometry or friction')
    bottom = _plate(nodes, elements, [lo[0]-margin, lo[1]-gap-thickness, lo[2]],
                    [hi[0]+margin, lo[1]-gap, hi[2]], platen_size)
    top = _plate(nodes, elements, [lo[0]-margin, hi[1]+gap, lo[2]],
                 [hi[0]+margin, hi[1]+gap+thickness, hi[2]], platen_size)
    assembly = mesh_audit(nodes, elements)
    lines = ['*HEADING', 'Dong2019 t=0.3mm comparison candidate; NOT calibrated or solver-validated',
             '** Units: mm, s, Mg (tonne), N, MPa. Abaqus/Explicit; no mass scaling.',
             '** Exact retained open-source core mesh; structured elastic steel platens.',
             '** Core U3=0, no in-plane pin; contact only on inward platen faces.',
             '** Damage endpoint coordinates are hypotheses, not recovered paper data.' if damage else
             '** No damage: material/solver comparison baseline; not a fracture model.', '*NODE']
    lines += [f'{n}, {p[0]:.15g}, {p[1]:.15g}, {p[2]:.15g}' for n, p in sorted(nodes.items())]
    for name, ids in [('CORE', sorted(core_elements)), ('BOTTOM_PLATE', bottom['elements']), ('TOP_PLATE', top['elements'])]:
        lines.append('*ELEMENT, TYPE=C3D8R, ELSET=' + name)
        lines += [', '.join(map(str, [e] + elements[e])) for e in ids]
    _nset(lines, 'CORE_NODES', core_nodes)
    _nset(lines, 'BOTTOM_LOAD', bottom['ymin_nodes'])
    _nset(lines, 'TOP_LOAD', top['ymax_nodes'])
    lines += ['*SURFACE, TYPE=ELEMENT, NAME=CORE_EXTERIOR', 'CORE,',
              '*SURFACE, TYPE=ELEMENT, NAME=CORE_INTERIOR', 'CORE, INTERIOR',
              '*SURFACE, NAME=CORE_CONTACT, COMBINE=UNION', 'CORE_EXTERIOR, CORE_INTERIOR']
    for name, facets in [('BOTTOM_INNER', bottom['ymax_faces']), ('TOP_INNER', top['ymin_faces'])]:
        lines.append('*SURFACE, TYPE=ELEMENT, NAME=' + name)
        lines += [f'{e}, {face}' for e, face in facets]
    lines += ['*MATERIAL, NAME=CONFIRMED_ALSI10MG', '*DENSITY', f'{config["density"]:.12g}',
              '*ELASTIC', f'{config["young"]:.12g}, {config["poisson"]:.12g}', '*PLASTIC, HARDENING=ISOTROPIC']
    lines += [f'{stress:.12g}, {strain:.12g}' for stress, strain in plastic]
    if damage:
        d = config.get('ductile_damage')
        locus = np.asarray(d['fracture_locus'] if d else [], float)
        if (locus.ndim != 2 or locus.shape[1] != 2 or len(locus) < 2
                or not np.isfinite(locus).all() or np.any(np.diff(locus[:, 0]) <= 0)
                or np.any(locus[:, 1] <= 0)):
            raise ValueError('Require an explicit ordered fracture-locus hypothesis')
        failure_displacement = float(d['failure_displacement'])
        if not np.isfinite(failure_displacement) or failure_displacement <= 0:
            raise ValueError('Require positive finite failure displacement')
        lines += ['*DAMAGE INITIATION, CRITERION=DUCTILE']
        lines += [f'{eps:.12g}, {eta:.12g}, 0.' for eta, eps in locus]
        lines += ['*DAMAGE EVOLUTION, TYPE=DISPLACEMENT, SOFTENING=LINEAR', f'{failure_displacement:.12g}']
    lines += ['*MATERIAL, NAME=ELASTIC_STEEL', '*DENSITY', f'{config.get("steel_density", 7.85e-9):.12g}',
              '*ELASTIC', f'{config.get("steel_young", 210000):.12g}, {config.get("steel_poisson", .3):.12g}',
              '*SECTION CONTROLS, NAME=CORE_CTRL, HOURGLASS=ENHANCED, ELEMENT DELETION=' + ('YES' if damage else 'NO'),
              '*SECTION CONTROLS, NAME=STEEL_CTRL, HOURGLASS=ENHANCED',
              '*SOLID SECTION, ELSET=CORE, MATERIAL=CONFIRMED_ALSI10MG, CONTROLS=CORE_CTRL', ',',
              '*SOLID SECTION, ELSET=BOTTOM_PLATE, MATERIAL=ELASTIC_STEEL, CONTROLS=STEEL_CTRL', ',',
              '*SOLID SECTION, ELSET=TOP_PLATE, MATERIAL=ELASTIC_STEEL, CONTROLS=STEEL_CTRL', ',']
    for name, value in [('CORE_SELF', self_friction), ('CORE_STEEL', friction)]:
        lines += ['*SURFACE INTERACTION, NAME=' + name, '*FRICTION', f'{value:.12g}',
                  '*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=HARD']
    lines += ['*BOUNDARY', 'CORE_NODES, 3, 3, 0.', 'BOTTOM_LOAD, 1, 3, 0.',
              'TOP_LOAD, 1, 1, 0.', 'TOP_LOAD, 3, 3, 0.',
              '*AMPLITUDE, NAME=LOAD, TIME=STEP TIME, DEFINITION=' + ('SMOOTH STEP' if loading == 'smooth' else 'TABULAR'),
              f'0., 0., {duration:.12g}, 1.', '*STEP, NAME=COMPRESS, NLGEOM=YES', '*DYNAMIC, EXPLICIT',
              f', {duration:.12g}', '*BULK VISCOSITY', '0.06, 1.2',
              '*CONTACT, OP=NEW', '*CONTACT INCLUSIONS', 'CORE_CONTACT, CORE_CONTACT',
              'CORE_CONTACT, BOTTOM_INNER', 'CORE_CONTACT, TOP_INNER', '*CONTACT PROPERTY ASSIGNMENT',
              'CORE_CONTACT, CORE_CONTACT, CORE_SELF', 'CORE_CONTACT, BOTTOM_INNER, CORE_STEEL',
              'CORE_CONTACT, TOP_INNER, CORE_STEEL', '*BOUNDARY, AMPLITUDE=LOAD',
              f'TOP_LOAD, 2, 2, {-displacement:.12g}',
              '*OUTPUT, FIELD, NUMBER INTERVAL=100', '*NODE OUTPUT', 'U, V, A, RF',
              '*ELEMENT OUTPUT, ELSET=CORE', 'S, LE, PEEQ' + (', SDEG, STATUS' if damage else ''),
              '*CONTACT OUTPUT', 'CFORCE, CSTRESS, CDISP',
              f'*OUTPUT, HISTORY, TIME INTERVAL={duration/1000:.12g}', '*NODE OUTPUT, NSET=TOP_LOAD', 'U2, RF2',
              '*NODE OUTPUT, NSET=BOTTOM_LOAD', 'U2, RF2', '*ENERGY OUTPUT',
              'ALLIE, ALLKE, ALLAE, ALLVD, ALLFD, ALLWK, ETOTAL', '*END STEP']
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(lines)+'\n', encoding='ascii')
    metadata = dict(backend='Abaqus/Explicit input, not executed', abaqus_datacheck_performed=False,
                    calibrated=False, experimental_agreement=False, wall_mm=float(config['wall']),
                    core_mesh=core, complete_mesh=assembly, duration_seconds=float(duration),
                    loading=loading, target_compression_mm=displacement,
                    target_nominal_strain=float(config['compression_strain']), reference_area_mm2=float(config.get('reference_area', 1900)),
                    constant_2_mm_per_min_duration_seconds=displacement*30,
                    imposed_mean_speed_mm_per_s=displacement/duration,
                    imposed_peak_speed_mm_per_s=displacement/duration*(1.875 if loading == 'smooth' else 1),
                    mass_scaling=False, bulk_viscosity=[.06, 1.2], hourglass='ENHANCED; independent Abaqus formulation',
                    material=dict(young_MPa=float(config['young']), poisson=float(config['poisson']),
                                  density_Mg_per_mm3=float(config['density']), plastic_rows=len(plastic),
                                  plastic_interpretation='already true stress/plastic strain; no conversion',
                                  after_last_plastic_point='Abaqus constant isotropic hardening extrapolation'),
                    damage_enabled=damage, fracture_locus=config['ductile_damage']['fracture_locus'] if damage else None,
                    failure_displacement_mm=float(config['ductile_damage']['failure_displacement']) if damage else None,
                    platen_mesh_target_size_mm=float(platen_size), steel_elements=len(bottom['elements'])+len(top['elements']),
                    top_loading_nodes=top['ymax_nodes'], bottom_loading_nodes=bottom['ymin_nodes'],
                    core_steel_friction=friction, core_self_friction=self_friction,
                    contact='Explicit general contact: core exterior+interior faces; steel inward Y faces only',
                    input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                    limitations=['Mesh is the retained research core, not the paper 10068-element mesh.',
                                 'Structured elastic steel platens replace the native one-brick plates.',
                                 'Steel side faces are excluded from contact; the older native case included them.',
                                 'Native AMS is not reproduced; no Abaqus mass scaling requested.',
                                 'Hourglass stabilization and contact differ between solvers.',
                                 'Core-core friction .35 is a hypothesis; the paper only specifies platen friction.',
                                 'Damage eta=0/1 coordinates are hypotheses; no CCX cutoff, residual stiffness or viscosity transferred.',
                                 'Accelerated loading needs duration, energy and mesh sensitivity checks.',
                                 'Loading-face RF2 is instantaneous actuator reaction including platen inertia, not isolated core contact force.'])
    sidecar.write_text(json.dumps(metadata, indent=2)+'\n')
    return metadata
