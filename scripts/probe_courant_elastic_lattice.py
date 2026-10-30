"""Run a matched elastic, contact-free lattice with bound native output.

The shared mesh/BC JSON comes from probe_matched_elastic_lattice.py. This is a
stiffness diagnostic, not a compression-paper reproduction or material fit.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from run_radioss import run_generated
from assess_courant import assess
from reentrant_fea.radioss import write_radioss_deck, _group, _bcs, _ints, _reals
from reentrant_fea.radioss_results import read_history_csv


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def generate(mesh_path, output, duration, ams_step, solid):
    shared = json.loads(mesh_path.read_text())
    nr, er = np.asarray(shared['nodes'], float), np.asarray(shared['elements'], int)
    if (not np.array_equal(nr[:, 0], np.arange(1, len(nr)+1))
            or not np.array_equal(er[:, 0], np.arange(1, len(er)+1))):
        raise ValueError('Require original contiguous node and element IDs')
    nodes, elements = nr[:, 1:], er[:, 1:]-1
    bc = shared['boundary']
    top, bottom = bc['top_y_displaced'], bc['bottom_y_fixed']
    displacement = float(bc['prescribed_top_y_mm'])
    height = float(np.ptp(nodes[:, 1]))
    material = shared['material']
    cfg = dict(young=material['young_MPa'], poisson=material['poisson'],
               reference_area=shared['reference_area_mm2'],
               density=2.67e-9, compression_strain=abs(displacement)/height)
    meta = write_radioss_deck(output/'elastic', nodes, elements, cfg,
                              [(225., 0.), (225., 1.)], duration=duration,
                              ams_time_step=ams_step, tensile_coupon=True)
    if meta['loading_face_nodes'] != top:
        raise ValueError('Native and static top-node selections differ')
    starter = Path(meta['starter'])
    original = starter.read_text().splitlines()
    start, stop = original.index('/MAT/PLAS_TAB/1'), original.index('/NODE')
    # Pure elasticity is essential: no plastic table participates in this test.
    original[start:stop] = ['/MAT/ELAST/1', 'Matched linear elastic aluminium',
                            _reals(2.67e-9), _reals(cfg['young'], cfg['poisson'])]
    prop = original.index('/PROP/SOLID/1')
    original[prop+1] = f'Matched elastic solid {solid}'
    original[prop+2] = _ints(solid, 0) + original[prop+2][20:]
    lines, i = [], 0
    while i < len(original):
        if original[i] == '/BCS/1':
            _bcs(lines, 1, 'Bottom Y only', 1, '010')
            i += 3
        elif original[i] == '/BCS/2':
            i += 3
        elif original[i] == '/IMPDISP/1':
            lines += original[i:i+3]+[_reals(1, displacement, 0, 0)]
            i += 4
        elif original[i] == '/TH/NODE/3':
            # The shared writer now includes supported bottom nodes in this
            # history too. Keep this diagnostic's loading block top-only;
            # analyze() pairs it with the separate bottom block added below.
            lines += original[i:i+3]
            lines += [_ints(n, 0)+f'load-node-{n}' for n in top]
            i += 3
            while i < len(original) and not original[i].startswith('/'):
                i += 1
        else:
            lines.append(original[i])
            i += 1
    if lines.pop() != '/END':
        raise ValueError('Unexpected deck layout')
    _group(lines, 3, 'All core Z', bc['all_z_fixed'])
    _bcs(lines, 3, 'All core Z', 3, '001')
    _group(lines, 4, 'One bottom X gauge', [bc['one_bottom_x_pin']])
    _bcs(lines, 4, 'One bottom X gauge', 4, '100')
    lines += ['/TH/NODE/4', 'Bottom nodal histories',
              ''.join(f'{v:<10}' for v in ('DX', 'DY', 'DZ', 'REACX', 'REACY', 'REACZ'))]
    lines += [_ints(n, 0)+f'bottom-{n}' for n in bottom]
    lines += ['/END']
    starter.write_text('\n'.join(lines)+'\n')
    engine = Path(meta['engine'])
    text = engine.read_text().splitlines()
    # Only endpoint animations; the quantitative diagnostic uses dense histories.
    j = text.index('/ANIM/DT')
    text[j+1] = _reals(0, duration)
    engine.write_text('\n'.join(text)+'\n')
    meta['sha256'] = {p.name: digest(p) for p in (starter, engine)}
    meta['element_formulation'] = f'Isolid={solid}, native finite-geometry default'
    meta.pop('hardening_extrapolation', None)
    meta['core_x_anchor'] = dict(enabled=True, node=bc['one_bottom_x_pin'],
                                 qualification='matched bottom rigid-translation gauge')
    meta['diagnostic'] = dict(kind='matched elastic contact-free lattice',
                              mesh_file=str(mesh_path), mesh_sha256=digest(mesh_path),
                              material='ELAST', contacts=False, plasticity=False,
                              damage=False, top_y_displacement_mm=displacement,
                              bottom_nodes=bottom, solid_formulation=solid,
                              script_sha256=digest(__file__))
    meta_path = starter.parent/(meta['prefix']+'.radioss.json')
    meta_path.write_text(json.dumps(meta, indent=2)+'\n')
    return meta


def analyze(meta_path):
    meta = json.loads(meta_path.read_text())
    folder, prefix = meta_path.parent, meta['prefix']
    csv_path = folder/(prefix+'T01.csv')
    actual = assess(meta_path, csv_path)
    if not actual['completed'] or not actual['metadata_binding_verified']:
        raise ValueError('Missing bound completed native evidence')
    starter = (folder/(prefix+'_0000.out')).read_text()
    engine = (folder/(prefix+'_0001.out')).read_text()
    if (not re.search(r'^\s*0\s+WARNING\(S\)\s*$', starter, re.M)
            or re.search(r'\b(?:nan|infinity|inf)\b', engine, re.I)
            or 'WARNING ZERO OR NEGATIVE VOLUME' in engine):
        raise ValueError('Native warnings or invalid response; inspect retained output')
    h = read_history_csv(csv_path)
    labels = list(h)
    v = np.column_stack(list(h.values()))
    top = [j for j, s in enumerate(labels) if s.startswith('Loading face nodal histories')]
    bottom = [j for j, s in enumerate(labels) if s.startswith('Bottom nodal histories')]
    if (len(top) != 6*len(meta['loading_face_nodes'])
            or len(bottom) != 6*len(meta['diagnostic']['bottom_nodes'])):
        raise ValueError('Nodal history ordering differs from verified six-component request')
    time = h['time']
    dt, midpoint = np.diff(time), .5*(time[:-1]+time[1:])
    force = -np.diff(v[:, top[4::6]].sum(1))/dt
    bottom_force = np.diff(v[:, bottom[4::6]].sum(1))/dt
    uy = v[:, top[1::6]].mean(1)
    displacement = -.5*(uy[:-1]+uy[1:])
    ie = .5*(h['INTERNAL ENERGY'][:-1]+h['INTERNAL ENERGY'][1:])
    ke = .5*(h['KINETIC ENERGY'][:-1]+h['KINETIC ENERGY'][1:])
    late = midpoint >= .9*meta['duration_seconds']
    stiffness = force[late]/displacement[late]
    work = float(np.sum(force*-np.diff(uy)))
    mesh_file = Path(meta['diagnostic']['mesh_file'])
    if digest(mesh_file) != meta['diagnostic']['mesh_sha256']:
        raise ValueError('Shared mesh/BC file changed after generation')
    area = json.loads(mesh_file.read_text())['reference_area_mm2']
    nominal_E = stiffness.mean()*meta['height_mm']/area
    # Summing hundreds of rounded nodal impulses and then differentiating can
    # amplify cancellation. Keep that result visible and cross-check against
    # aggregate work and an integrated spring/inertia regression.
    impulse = -v[:, top[4::6]].sum(1)
    observed_u = -uy
    integrated_u = np.r_[0., np.cumsum(.5*(observed_u[1:]+observed_u[:-1])*dt)]
    velocity = np.gradient(observed_u, time)
    selected = (time >= .25*meta['duration_seconds']) & (time <= .95*meta['duration_seconds'])
    design = np.c_[integrated_u[selected], velocity[selected], np.ones(selected.sum())]
    coefficients = np.linalg.lstsq(design, impulse[selected], rcond=None)[0]
    impulse_residual = design@coefficients-impulse[selected]
    final_u = float(observed_u[-1])
    work_k = 2*float(h['EXTERNAL WORK'][-1])/final_u**2
    force_difference = float(np.max(abs(force[late]-bottom_force[late])))
    mean_force = float(np.mean(abs(force[late])))
    reaction_consistent = force_difference <= .01*max(mean_force, 1e-20)
    report = dict(status='matched elastic diagnostic, not calibrated or experimental validation',
                   metadata=str(meta_path), metadata_sha256=digest(meta_path),
                   mesh_sha256=meta['diagnostic']['mesh_sha256'],
                   completed=True, metadata_execution_binding_verified=True,
                   solid_formulation=meta['diagnostic']['solid_formulation'],
                   duration_seconds=meta['duration_seconds'], ams=meta['ams'],
                   core_nodes=meta['core_nodes'], core_elements=meta['core_elements'],
                   raw_nodal_late_stiffness_N_per_mm=float(stiffness.mean()),
                   raw_nodal_late_stiffness_range_N_per_mm=[float(stiffness.min()), float(stiffness.max())],
                   raw_nodal_nominal_effective_modulus_MPa=float(nominal_E),
                   nominal_reference_area_mm2=area,
                   executed_metadata_reference_area_mm2=meta['reference_area_mm2'],
                   area_normalization_note='Effective moduli use the shared 1900 mm2 convention; executed inputs and metadata are preserved.',
                   work_equivalent_stiffness_N_per_mm=work_k,
                   work_equivalent_modulus_MPa=work_k*meta['height_mm']/area,
                   energy_equivalent_stiffness_N_per_mm=2*float(h['INTERNAL ENERGY'][-1])/final_u**2,
                   integrated_impulse_regression=dict(
                       formula='I(t)=k*integral(displacement dt)+m_eff*velocity+offset',
                       time_fraction_interval=[.25, .95],
                       stiffness_N_per_mm=float(coefficients[0]),
                       effective_mass_Mg=float(coefficients[1]),
                       offset_Ns=float(coefficients[2]),
                       impulse_rmse_Ns=float(np.sqrt(np.mean(impulse_residual**2))),
                       qualification='Diagnostic spring/inertia approximation; not an independent static-stiffness validation'),
                   max_late_kinetic_internal_ratio=float(np.max(ke[late]/ie[late])),
                   final_internal_energy_Nmm=float(h['INTERNAL ENERGY'][-1]),
                   final_external_work_Nmm=float(h['EXTERNAL WORK'][-1]),
                   reconstructed_top_displacement_work_Nmm=work,
                   max_late_top_bottom_force_difference_N=force_difference,
                   nodal_reaction_differentiation_qualified=reaction_consistent,
                   analysis_script_sha256=digest(__file__),
                   generation_script_sha256=meta['diagnostic']['script_sha256'],
                   actual_history_samples=len(time),
                   quasi_static_validated=False, experimental_agreement=False,
                   force_history=[dict(time=float(t), displacement_mm=float(u), force_N=float(f),
                                       bottom_force_N=float(b))
                                  for t,u,f,b in zip(midpoint,displacement,force,bottom_force)],
                   qualification='Raw nodal differentiation fails its balance check when flagged false. Work and integrated-impulse values remain diagnostics pending duration/AMS refinement; elastic diagnostics cannot qualify plastic collapse.')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--duration', type=float, default=.004)
    group = p.add_mutually_exclusive_group()
    group.add_argument('--ams-time-step', type=float, default=1e-6)
    group.add_argument('--no-ams', action='store_true', help='Use the native stable time step without Advanced Mass Scaling')
    p.add_argument('--solid', type=int, choices=(14, 18, 24), default=24)
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--timeout-seconds', type=float, default=1800)
    args = p.parse_args()
    if args.threads < 1:
        p.error('Positive thread count required')
    if args.output.exists() and any(args.output.iterdir()):
        p.error('Use a fresh directory; previous evidence is preserved')
    args.output.mkdir(parents=True, exist_ok=True)
    os.environ['OMP_NUM_THREADS'] = str(args.threads)
    meta = generate(args.mesh.resolve(), args.output.resolve(), args.duration,
                    None if args.no_ams else args.ams_time_step, args.solid)
    run_generated(meta, timeout=args.timeout_seconds)
    path = Path(meta['starter']).parent/(meta['prefix']+'.radioss.json')
    report = analyze(path)
    path.with_name('diagnostic.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'force_history'}, indent=2))


if __name__ == '__main__':
    main()
