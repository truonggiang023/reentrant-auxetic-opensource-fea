"""Short matched native benchmark of ALL versus physical-Y platen surfaces.

Retains the .01-s loading curve and core ALL self-contact. A short engine stop
is intentional, not completion of the lattice loading or a calibration trial.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from run_radioss import run_generated
from reentrant_fea.radioss import _ints, _reals
from reentrant_fea.radioss_results import read_engine_out, read_history_csv


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, default=ROOT/'results/opencourant/thin_free_x_ams1e6.radioss.json')
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--stop-time', type=float, default=1e-4)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--timeout', type=float, default=120)
    ap.add_argument('--modes', choices=('all', 'y_only'), nargs='+', default=['all', 'y_only'])
    ap.add_argument('--print-cycles', type=int, default=10)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError('Preserve existing evidence; use a fresh output directory')
    args.output.mkdir(parents=True)
    source = json.loads(args.source.read_text())
    if not 0 < args.stop_time < source['duration_seconds']:
        raise ValueError('Require an intentionally partial positive stop time')
    folder, baseline = args.source.parent, source['prefix']
    for name, expected in source['sha256'].items():
        if sha(folder/name) != expected:
            raise ValueError('Source input changed: '+name)
    os.environ['OMP_NUM_THREADS'] = str(args.threads)
    results = []
    if args.print_cycles <= 0:
        raise ValueError('Require positive print-cycle interval')
    for mode in args.modes:
        directory = (args.output/mode).resolve()
        directory.mkdir()
        name = 'bench_'+mode
        lines = (folder/(baseline+'_0000.rad')).read_text().replace(baseline, name).splitlines()
        if mode == 'y_only':
            nc = source['core_nodes']
            for ident, ids in [(2, [nc+i for i in (3, 4, 8, 7)]), (3, [nc+i for i in (9, 10, 14, 13)])]:
                at = lines.index(f'/SURF/PART/ALL/{ident}')
                lines[at:at+3] = [f'/SURF/SEG/{ident}', f'Physical platen Y face {ident}', _ints(ident, *ids)]
        starter = directory/(name+'_0000.rad')
        starter.write_text('\n'.join(lines)+'\n')
        el = (folder/(baseline+'_0001.rad')).read_text().replace(baseline, name).splitlines()
        at = el.index('/RUN/'+name+'/1')
        el[at+1] = _reals(args.stop_time)
        el[el.index('/PRINT/-1000')] = f'/PRINT/-{args.print_cycles}'
        engine = directory/(name+'_0001.rad')
        engine.write_text('\n'.join(el)+'\n')
        meta = dict(source, prefix=name, starter=str(starter), engine=str(engine), duration_seconds=args.stop_time)
        meta['diagnostic'] = dict(intentionally_partial=True, full_loading_duration_seconds=source['duration_seconds'],
                                  requested_stop_seconds=args.stop_time, source_metadata_sha256=sha(args.source),
                                  source_starter_sha256=sha(folder/(baseline+'_0000.rad')), platen_surface=mode,
                                  core_all_face_self_contact_preserved=True)
        meta['sha256'] = {p.name: sha(p) for p in (starter, engine)}
        metadata_file = directory/(name+'.radioss.json')
        metadata_file.write_text(json.dumps(meta, indent=2)+'\n')
        prospective = dict(metadata_sha256=sha(metadata_file), input_sha256=meta['sha256'],
                           script_sha256=sha(__file__), omp_threads=args.threads)
        (directory/'execution_started.json').write_text(json.dumps(prospective, indent=2)+'\n')
        started = time.monotonic()
        failure = None
        try:
            history = run_generated(meta, timeout=args.timeout)
        except Exception as exc:
            failure = f'{type(exc).__name__}: {exc}'
            history = None
        elapsed = time.monotonic()-started
        record = dict(mode=mode, elapsed_wall_seconds=elapsed, failure=failure,
                      partial_stages_completed=history is not None, completed_full_loading=False,
                      calibrated=False, quasi_static_validated=False,
                      prospective_execution_record=prospective,
                      input_output_sha256={p.name: sha(p) for p in directory.iterdir() if p.is_file()})
        if (directory/(name+'_0001.out')).is_file():
            record['engine'] = read_engine_out(directory/(name+'_0001.out'))
        if history is not None:
            h = read_history_csv(history)
            t = h['time']
            normal = h[next(k for k in h if k.endswith('var 64'))]
            tangent = h[next(k for k in h if k.endswith('var 67'))]
            amplitude = np.asarray(meta['displacement_fraction_samples'])
            record['last_recorded_time'] = float(t[-1])
            record['last_recorded_nominal_strain'] = float(np.interp(t[-1], amplitude[:, 0], amplitude[:, 1])*meta['target_nominal_strain'])
            record['final_contact_impulse_Ns'] = float(normal[-1]-tangent[-1])
        results.append(record)
        print(json.dumps({k: record[k] for k in ('mode', 'elapsed_wall_seconds', 'failure', 'partial_stages_completed')}), flush=True)
        (args.output/'report.json').write_text(json.dumps(dict(status='Intentionally partial surface timing diagnostic; not calibration', cases=results), indent=2)+'\n')


if __name__ == '__main__':
    main()
