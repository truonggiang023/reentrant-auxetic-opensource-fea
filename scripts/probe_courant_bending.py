"""Reproduce an elastic HEPH wall-bending stability/stiffness diagnostic.

Run after sourcing /workspace/opencourant/env.sh. This generates original
input decks, not vendor examples. It deliberately retains failed native output:
NORMAL TERMINATION alone does not make a response usable. The three controls
and the high-aspect/default-step case are observations, not calibration fits.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np

from reentrant_fea.radioss import write_radioss_deck


ROOT = Path(__file__).resolve().parents[1]
E, NU, THICKNESS, LENGTH = 69000., .3, .3, 5.66
DISPLACEMENT, DURATION = .001, .002
NX, NY = 2, 16
CASES = [('near_square_control', .3, 1, .9),
         ('high_aspect_default_step', 25., 1, .9),
         ('high_aspect_four_depth_layers', 25., 4, .9),
         ('high_aspect_reduced_step', 25., 1, .3)]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def generate(directory, depth, layers, step_scale):
    x, y, z = np.linspace(0, THICKNESS, NX+1), np.linspace(0, LENGTH, NY+1), np.linspace(0, depth, layers+1)
    nodes = np.array([[xx, yy, zz] for zz in z for yy in y for xx in x])
    index = lambda i, j, k: k*(NX+1)*(NY+1)+j*(NX+1)+i
    elements = np.array([[index(i,j,k), index(i+1,j,k), index(i+1,j+1,k), index(i,j+1,k),
                          index(i,j,k+1), index(i+1,j,k+1), index(i+1,j+1,k+1), index(i,j+1,k+1)]
                         for k in range(layers) for j in range(NY) for i in range(NX)])
    config = {'young': E, 'poisson': NU, 'compression_strain': DISPLACEMENT/LENGTH}
    metadata = write_radioss_deck(directory/'wall', nodes, elements, config,
                                  [(225.,0.), (300.,1.)], duration=DURATION, tensile_coupon=True)
    starter = Path(metadata['starter'])
    original = starter.read_text().splitlines()
    lines, position = [], 0
    while position < len(original):
        if original[position] == '/BCS/1':
            # Global Z restraint below supplies this DOF without duplication.
            lines += ['/BCS/1', 'Bottom XY clamp', f'{"110":>6} 000'+f'{0:10d}{1:10d}']
            position += 3
        elif original[position] == '/BCS/2':
            position += 3  # Tip Y remains free to rotate its section.
        elif original[position] == '/IMPDISP/1':
            lines += ['/IMPDISP/1', 'Prescribed X displacement',
                      original[position+2].replace('Y', 'X'), original[position+3]]
            position += 4
        else:
            lines.append(original[position])
            position += 1
    assert lines.pop() == '/END'
    lines += ['/GRNOD/NODE/5', 'All nodes plane-strain Z restraint']
    for first in range(1, len(nodes)+1, 10):
        lines.append(''.join(f'{node:10d}' for node in range(first, min(first+10, len(nodes)+1))))
    lines += ['/BCS/5', 'All Z fixed', f'{"001":>6} 000'+f'{0:10d}{5:10d}', '/END']
    starter.write_text('\n'.join(lines)+'\n')
    engine = Path(metadata['engine'])
    engine.write_text(engine.read_text()+f'\n/DT\n{step_scale:20.12g}{0:20.12g}\n')
    metadata['sha256'] = {path.name: digest(path) for path in (starter, engine)}
    metadata['diagnostic'] = {'type': 'elastic cantilever bending', 'loading_direction': 'X',
                              'fixed': 'bottom XY, all Z', 'tip_Y_is_free': True,
                              'native_global_time_step_scale': step_scale}
    (directory/'wall.radioss.json').write_text(json.dumps(metadata, indent=2)+'\n')
    return metadata


def execute(directory, executable, arguments, timeout):
    log = directory/(Path(executable).name+'.log')
    environment = os.environ.copy()
    environment['OMP_NUM_THREADS'] = '1'
    with log.open('w') as stream:
        run = subprocess.run([executable, *arguments], cwd=directory, env=environment,
                             stdout=stream, stderr=subprocess.STDOUT, timeout=timeout)
    return run.returncode


def probe(output, label, depth, layers, scale, executables, timeout):
    directory = output/label
    directory.mkdir()
    metadata = generate(directory, depth, layers, scale)
    starter_code = execute(directory, executables[0], ['-i', 'wall_0000.rad', '-np', '1'], timeout)
    engine_code = execute(directory, executables[1], ['-i', 'wall_0001.rad'], timeout)
    starter = (directory/'wall_0000.out').read_text()
    engine = (directory/'wall_0001.out').read_text()
    errors = [line.strip() for line in engine.splitlines() if '**ERROR' in line.upper()]
    volumes = [line.strip() for line in engine.splitlines() if 'WARNING ZERO OR NEGATIVE VOLUME' in line]
    not_finite = bool(re.search(r'\b(?:nan|infinity|inf)\b', engine, re.IGNORECASE))
    cycle_match = re.search(r'TOTAL NUMBER OF CYCLES\s*:\s*(\d+)', engine)
    banner = (directory/(Path(executables[0]).name+'.log')).read_text()
    commit = re.search(r'CommitID:\s*([0-9a-f]{40})', banner)
    result = dict(case=label, directory=str(directory), depth_mm=depth, depth_layers=layers,
                  core_elements=NX*NY*layers, nominal_brick_dimensions_mm=[THICKNESS/NX, LENGTH/NY, depth/layers],
                  largest_to_smallest_edge_ratio=max(THICKNESS/NX,LENGTH/NY,depth/layers)/min(THICKNESS/NX,LENGTH/NY,depth/layers),
                  native_global_time_step_scale=scale, starter_exit_code=starter_code,
                  engine_exit_code=engine_code, engine_normal_termination='NORMAL TERMINATION' in engine,
                  starter_zero_errors='0 ERROR(S)' in starter, starter_zero_warnings='0 WARNING(S)' in starter,
                  engine_error_line_count=len(errors), sample_engine_error=errors[0] if errors else None,
                  negative_volume_message_count=len(volumes), contains_nonfinite_printed_values=not_finite,
                  completed_cycles=int(cycle_match.group(1)) if cycle_match else None,
                  runtime_binary_commit=commit.group(1) if commit else None)
    usable = (starter_code == engine_code == 0 and result['starter_zero_errors']
              and result['starter_zero_warnings'] and result['engine_normal_termination']
              and not errors and not volumes and not not_finite)
    result['elastic_stiffness_comparison_usable'] = usable
    if usable:
        code = execute(directory, executables[2], ['wallT01'], min(timeout, 60))
        if code:
            raise RuntimeError('History conversion failed for '+label)
        with (directory/'wallT01.csv').open() as stream:
            reader = csv.reader(stream)
            headers = next(reader)
            values = np.array([[float(value) for value in row] for row in reader])
        if not np.isfinite(values).all() or np.any(np.diff(values[:,0]) <= 0):
            raise ValueError('Invalid actual history in '+label)
        columns = [i for i,h in enumerate(headers) if h.startswith('Loading face nodal histories')]
        if len(columns) != 6*len(metadata['loading_face_nodes']):
            raise ValueError('Nodal history does not match requested component ordering')
        impulse = values[:,columns[3::6]].sum(axis=1)
        force = np.diff(impulse)/np.diff(values[:,0])
        time = .5*(values[1:,0]+values[:-1,0])
        amplitude = np.asarray(metadata['displacement_fraction_samples'])
        displacement = np.interp(time, amplitude[:,0], amplitude[:,1])*DISPLACEMENT
        late = time > .9*DURATION
        if not late.any() or values[-1,0] < .995*DURATION:
            raise ValueError('No usable late-ramp coverage for '+label)
        stiffness = force[late]/displacement[late]
        reference = 3*(E/(1-NU**2))*(depth*THICKNESS**3/12)/LENGTH**3
        result.update(expected_plane_strain_EB_stiffness_N_per_mm=reference,
                      mean_late_tip_stiffness_N_per_mm=float(stiffness.mean()),
                      late_stiffness_range_N_per_mm=[float(stiffness.min()),float(stiffness.max())],
                      relative_stiffness_error=float(stiffness.mean()/reference-1),
                      last_observed_history_time_seconds=float(values[-1,0]), history_samples=len(values))
    else:
        result['excluded_reason'] = 'Native warning/error/nonfinite response: normal termination and exit0 alone are insufficient.'
    evidence = [directory/name for name in ('wall_0000.rad','wall_0001.rad','wall_0000.out','wall_0001.out','wallT01','wallT01.csv')]
    result['input_and_output_sha256'] = {path.name: digest(path) for path in evidence if path.is_file()}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='New empty output directory; outputs are retained')
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--timeout-seconds', type=float, default=120)
    args = parser.parse_args()
    if not 0 < args.timeout_seconds < float('inf'):
        parser.error('Require finite positive timeout')
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error('Use a new empty output directory to preserve prior evidence')
    executables = [shutil.which(name) for name in ('starter_linux64_gf','engine_linux64_gf','th_to_csv_linux64_gf')]
    if not all(executables):
        parser.error('Source the qualified OpenCourant env.sh first')
    output.mkdir(parents=True, exist_ok=True)
    cases = [probe(output, *case, executables, args.timeout_seconds) for case in CASES]
    cfg = Path(os.environ['RAD_CFG_PATH'])/'config/CFG/radioss2022/PROP/prop_p14_solid.cfg'
    report = dict(status='elastic component audit; no lattice calibration or experimental agreement claim',
                  reproduction=dict(command='source /workspace/opencourant/env.sh && '+
                                     ' '.join([sys.executable,'scripts/probe_courant_bending.py','--output',str(output),'--report',str(args.report)]),
                                     script='scripts/probe_courant_bending.py', script_sha256=digest(__file__)),
                  geometry=dict(wall_thickness_mm=THICKNESS,wall_length_mm=LENGTH,elements_through_wall=NX,elements_along_length=NY),
                  material=dict(young_MPa=E,poisson=NU,density_Mg_per_mm3=2.67e-9,true_flow_table=[[0,225],[1,300]]),
                  boundary_conditions='Bottom X,Y fixed; all Z fixed; top prescribed X with Y free.',
                  loading=dict(tip_displacement_mm=DISPLACEMENT,duration_seconds=DURATION,mass_scaling=False,
                               native_default_numerical_damping=True,amplitude='101-point piecewise-linear quintic smoothstep',
                               history_interval_seconds=DURATION/400),
                  reference_formula='k=3*[E/(1-nu^2)]*[depth*wall_thickness^3/12]/wall_length^3; late-ramp force from Delta(impulse)/Delta(time)',
                  finite_geometry_flag_audit=dict(native_property_cfg=str(cfg),cfg_sha256=digest(cfg),solid_formulation=24,
                                                  corotational_frame=2,ismstr=4,is_small_strain_approximation=False,
                                                  ismstr4_definition='Full Geometric Nonlinearities (/DT/BRICK/CST has No Effect)',source_lines=[202,210]),
                  cases=cases,
                  limitations=['Elastic component diagnostic, not experimental calibration or a universal HEPH aspect-ratio rule.',
                               'Negative-volume/time-step errors and nonfinite output are excluded even if the engine prints NORMAL TERMINATION.',
                               'The AMS lattice uses a different inertial operator and still needs independent time-step, rate, extrusion and contact sensitivity.',
                               'The observed lattice stress mismatch is not attributed to elastic bending stiffness by this diagnostic.'])
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report,indent=2)+'\n')
    for row in cases:
        print(row['case'], 'usable='+str(row['elastic_stiffness_comparison_usable']),
              'relative_stiffness_error='+str(row.get('relative_stiffness_error')))
    print('Observed diagnostic evidence written to '+str(args.report))


if __name__ == '__main__':
    main()
