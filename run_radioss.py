"""Generate an experimental OpenCourant native erosion/contact model.

The .rad format retains its Radioss name after the OpenCourant migration.
Duration must be explicit; accelerated explicit dynamics requires an energy
audit and duration sensitivity before any quasi-static comparison.
"""
import argparse
import hashlib
import json
import os
import signal
import shutil
import subprocess
from pathlib import Path

from reentrant_fea.geometry import core_polygon
from reentrant_fea.mesh import mesh_bricks
from reentrant_fea.plastic import from_reference
from reentrant_fea.radioss import write_radioss_deck
from reentrant_fea.radioss_results import normal_termination, solver_errors


def run_logged(command, directory, log, timeout):
    """Bound a solver process and its children; retain the underlying failure."""
    with Path(log).open('w') as stream:
        process = subprocess.Popen(command, cwd=directory, stdout=stream,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            raise RuntimeError(f'Solver timeout after {timeout:g}s; partial output retained at {log}')
    if code:
        raise subprocess.CalledProcessError(code, command)


def run_generated(metadata, *, timeout, starter='starter_linux64_gf',
                  engine='engine_linux64_gf', converter='th_to_csv_linux64_gf'):
    """Run the sourced official runtime; process success is not qualification."""
    if not 0 < timeout < float('inf'):
        raise ValueError('Require a finite positive timeout')
    executables = []
    for name in (starter, engine, converter):
        path = shutil.which(name)
        if path is None:
            raise FileNotFoundError(f'{name} missing; source /workspace/opencourant/env.sh first')
        executables.append(path)
    directory = Path(metadata['starter']).parent
    prefix = metadata['prefix']
    metadata_file = directory/(prefix+'.radioss.json')
    if json.loads(metadata_file.read_text()) != metadata:
        raise ValueError('Saved metadata does not match the requested model')
    metadata_digest = hashlib.sha256(metadata_file.read_bytes()).hexdigest()
    for name, expected in metadata['sha256'].items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Generated deck changed before execution: '+name)
    for name in (prefix+'_0000.out', prefix+'_0001.out', prefix+'T01', prefix+'T01.csv', prefix+'.run.json'):
        if (directory/name).exists():
            raise FileExistsError('Use a fresh output prefix; prior solver evidence exists: '+name)
    run_logged([executables[0], '-i', Path(metadata['starter']).name, '-np', '1'],
               directory, directory/(prefix+'-starter.log'), timeout)
    starter_out = directory/(prefix+'_0000.out')
    if not starter_out.is_file() or not normal_termination(starter_out.read_text(errors='replace')):
        raise RuntimeError('Starter did not terminate normally; inspect '+str(starter_out))
    if solver_errors(starter_out.read_text(errors='replace')):
        raise RuntimeError('Starter reports input errors; inspect '+str(starter_out))
    run_logged([executables[1], '-i', Path(metadata['engine']).name],
               directory, directory/(prefix+'-engine.log'), timeout)
    engine_out = directory/(prefix+'_0001.out')
    if not engine_out.is_file() or not normal_termination(engine_out.read_text(errors='replace')):
        raise RuntimeError('Engine did not terminate normally; inspect '+str(engine_out))
    if solver_errors(engine_out.read_text(errors='replace')):
        raise RuntimeError('Engine reports errors despite termination; inspect '+str(engine_out))
    history = directory/(prefix+'T01')
    if not history.is_file():
        raise RuntimeError('Missing actual T01 history')
    run_logged([executables[2], history.name], directory,
               directory/(prefix+'-convert.log'), min(timeout, 60))
    csv = history.with_name(history.name+'.csv')
    if not csv.is_file() or not csv.stat().st_size:
        raise RuntimeError('Converter did not produce a history CSV')
    if hashlib.sha256(metadata_file.read_bytes()).hexdigest() != metadata_digest:
        raise ValueError('Model metadata changed during solver execution')
    record = {'input_sha256': metadata['sha256'],
              'metadata_sha256': metadata_digest,
              'executables_sha256': {str(path): hashlib.sha256(Path(path).read_bytes()).hexdigest()
                                     for path in executables},
              'output_sha256': {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                                for path in (starter_out, engine_out, history, csv)},
              'stages': ['starter -np 1', 'engine', 'th_to_csv'],
              'solver_stages_terminated_normally': True}
    (directory/(prefix+'.run.json')).write_text(json.dumps(record, indent=2)+'\n')
    return csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='configs/dong2019_brick_damage.json')
    parser.add_argument('--plastic', default='reference_setup_numeric_source.py')
    parser.add_argument('--output', default='results/opencourant/thin')
    parser.add_argument('--duration-seconds', type=float, required=True,
                        help='Explicit duration; not assumed equivalent to 2 mm/min')
    parser.add_argument('--drop-ccx-viscosity', action='store_true',
                        help='Record intentional omission of CCX pseudo-time viscosity')
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--ams-time-step', type=float,
                        help='Optional native AMS target step; requires its own response sensitivity checks')
    parser.add_argument('--timeout-seconds', type=float, default=600,
                        help='Maximum runtime of each solver stage; partial outputs are retained')
    args = parser.parse_args()
    config = json.loads(Path(args.config).read_text())
    if config.get('element_type') not in ('C3D8', 'C3D8R'):
        parser.error('Select a brick configuration; wedge conversion is not qualified')
    nodes, elements, _, _ = mesh_bricks(core_polygon(config), config['depth'],
                                       config['mesh_size'], config['depth_layers'],
                                       strategy=config.get('brick_mesh', 'recombine'))
    metadata = write_radioss_deck(args.output, nodes, elements, config,
                                 from_reference(args.plastic), duration=args.duration_seconds,
                                 drop_ccx_viscosity=args.drop_ccx_viscosity,
                                 ams_time_step=args.ams_time_step)
    metadata['source_config_sha256'] = hashlib.sha256(Path(args.config).read_bytes()).hexdigest()
    metadata['hardening_source_sha256'] = hashlib.sha256(Path(args.plastic).read_bytes()).hexdigest()
    metadata_path = Path(metadata['starter']).parent/(metadata['prefix']+'.radioss.json')
    metadata_path.write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps(metadata, indent=2))
    if args.run:
        csv = run_generated(metadata, timeout=args.timeout_seconds)
        from assess_courant import assess
        # This three-part writer's verified TYPE25 DEF Y channel is cumulative
        # impulse. Differencing it gives interval-average force, never raw N.
        result = assess(
            metadata_path, csv,
            force_column='Contact force and energy                 3 Contact 3                                var 64',
            force_sign=-1, force_kind='impulse')
        metadata_path.with_name(metadata['prefix']+'.assessment.json').write_text(json.dumps(result,indent=2)+'\n')
        print('Solver stages terminated normally. Energy assessment:', result['global_energy_screen'])
    else:
        print('Generated only. Run with the qualified OpenCourant starter and engine.')
    print('Native erosion/contact are requested; calibration and experimental agreement remain unverified.')


if __name__ == '__main__':
    main()
