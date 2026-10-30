"""Execute selected prepared study runs sequentially, preserving failed evidence.

Completion is solver evidence, not calibration or quasi-static qualification.
Source /workspace/opencourant/env.sh before executing native study runs.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from calibration_study import (NATIVE_FORCE_COLUMN, NATIVE_TANGENTIAL_COLUMN, ROOT,
                               legacy_native_force_mapping, native_force_mapping, native_assess_kwargs)
from reentrant_fea.plastic import load_plastic, validate_material


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name+'.next')
    temporary.write_text(json.dumps(data, indent=2)+'\n')
    temporary.replace(path)


def select_runs(manifest, candidates=None, role='training'):
    wanted = set(candidates or [])
    if role not in ('training', 'holdout', 'all'):
        raise ValueError('Select training, holdout or all roles')
    known = {row['candidate'] for row in manifest['runs']}
    if wanted-known:
        raise ValueError('Unknown candidate(s): '+', '.join(sorted(wanted-known)))
    rows = [row for row in manifest['runs'] if (not wanted or row['candidate'] in wanted)
            and (role == 'all' or row['role'] == role)]
    if not rows:
        raise ValueError('No runs match the requested candidates and role')
    keys = [(row['candidate'], row['role']) for row in rows]
    if len(set(keys)) != len(keys):
        raise ValueError('Duplicate candidate/role records are ambiguous')
    return rows


def check_inputs(row):
    config_path, source = Path(row['config']), Path(row['hardening_source'])
    if digest(config_path) != row['config_sha256'] or digest(source) != row['hardening_source_sha256']:
        raise ValueError('Configuration or hardening source changed since study preparation')
    config = json.loads(config_path.read_text())
    validate_material(config, load_plastic(source))
    if config['wall'] != row['wall_mm']:
        raise ValueError('Manifest wall thickness does not match configuration')
    return config


def parsed_options(argv, script, allowed, positional=0):
    """Allow only recorded arguments to known repository scripts, never a shell."""
    if (not isinstance(argv, list) or len(argv) < 2+positional
            or not all(isinstance(value, str) for value in argv)
            or (ROOT/argv[1]).resolve() != (ROOT/script).resolve()):
        raise ValueError('Recorded command must invoke the expected repository script')
    options, index = {}, 2+positional
    while index < len(argv):
        key = argv[index]
        if key not in allowed or key in options:
            raise ValueError('Unsupported or duplicate recorded argument: '+key)
        count = allowed[key]
        if index+count >= len(argv):
            raise ValueError('Missing recorded argument value: '+key)
        options[key] = argv[index+1] if count else True
        index += 1+count
    return options


def commands(row, threads, timeout):
    backend = row.get('backend', 'ccx')
    if backend not in ('ccx', 'opencourant'):
        raise ValueError('Unsupported study backend')
    native = backend == 'opencourant'
    script = 'run_radioss.py' if native else 'run_paper.py'
    allowed = {'--config': 1, '--plastic': 1, '--output': 1, '--run': 0}
    allowed.update({'--duration-seconds': 1, '--ams-time-step': 1,
                    '--drop-ccx-viscosity': 0, '--timeout-seconds': 1} if native else
                   {'--threads': 1, '--target-output-only': 0, '--self-contact': 0, '--imperfection-mm': 1})
    options = parsed_options(row['solve_argv'], script, allowed)
    for option, field in (('--config', 'config'), ('--plastic', 'hardening_source')):
        if (ROOT/str(options.get(option, ''))).resolve() != Path(row[field]).resolve():
            raise ValueError('Recorded command differs from manifest '+field)
    if options.get('--run') is not True:
        raise ValueError('Recorded solve command must explicitly request --run')
    expected_output = (Path(row['metadata']).with_name(Path(row['metadata']).name.removesuffix('.radioss.json'))
                       if native else Path(row['deck']))
    if (ROOT/str(options.get('--output', ''))).resolve() != expected_output.resolve():
        raise ValueError('Recorded output differs from manifest')
    if native:
        if float(options.get('--duration-seconds', 'nan')) != row['duration_seconds']:
            raise ValueError('Recorded native duration differs from manifest')
        recorded_ams = float(options['--ams-time-step']) if '--ams-time-step' in options else None
        if recorded_ams != row['ams_time_step'] or not options.get('--drop-ccx-viscosity'):
            raise ValueError('Recorded AMS or viscosity omission differs from manifest')
        options['--timeout-seconds'] = str(timeout)
    else:
        options['--threads'] = str(threads)
    solve = [sys.executable, str(ROOT/script)]
    for key, value in options.items():
        solve.append(key)
        if value is not True:
            solve.append(str(value))
    assess_script = 'assess_courant.py' if native else 'assess_run.py'
    positional = 2 if native else 1
    allowed_assess = {'--output': 1}
    if native:
        allowed_assess.update({'--force-column': 1, '--force-sign': 1, '--force-kind': 1,
                               '--tangential-column': 1, '--force-source': 1})
    opts = parsed_options(row['assess_argv'], assess_script, allowed_assess, positional)
    expected = [row['metadata'], row['history']] if native else [row['deck']]
    for actual, path in zip(row['assess_argv'][2:2+positional], expected):
        if (ROOT/actual).resolve() != Path(path).resolve():
            raise ValueError('Assessment input differs from manifest')
    if (ROOT/str(opts.get('--output', ''))).resolve() != Path(row['assessment']).resolve():
        raise ValueError('Assessment destination differs from manifest')
    if native:
        legacy = legacy_native_force_mapping(row.get('force_mapping'))
        if row.get('force_mapping') == native_force_mapping(reference_points=True):
            if opts.get('--force-source') != 'rp' or any(k in opts for k in ('--force-column','--force-kind','--force-sign','--tangential-column')):
                raise ValueError('RP study assessment must use only the declared RP source')
        else:
            if (opts.get('--force-source') is not None or opts.get('--force-column') != NATIVE_FORCE_COLUMN
                    or opts.get('--force-sign') != '-1' or opts.get('--force-kind') != 'impulse'):
                raise ValueError('Native assessment must use the verified FY impulse channel')
            tangent = opts.get('--tangential-column')
            if tangent not in (None,NATIVE_TANGENTIAL_COLUMN) or (tangent is None and not legacy):
                raise ValueError('Native assessment must use the verified top tangential impulse channel')
            opts['--tangential-column'] = NATIVE_TANGENTIAL_COLUMN
    assessment = [sys.executable, str(ROOT/assess_script), *row['assess_argv'][2:2+positional]]
    for key, value in opts.items():
        assessment.extend([key, str(value)])
    return solve, assessment


def receipt_path(row):
    path = Path(row['assessment'])
    return path.with_name(path.name+'.execution.json')


def verify_completed(row, *, require_receipt=True):
    """Read actual solver evidence; cached completed flags cannot skip a run."""
    config = check_inputs(row)
    if row.get('backend', 'ccx') == 'opencourant':
        from assess_courant import assess
        legacy_mapping = legacy_native_force_mapping(row.get('force_mapping'))
        metadata = json.loads(Path(row['metadata']).read_text())
        if (metadata.get('source_config_sha256') != row['config_sha256']
                or metadata.get('hardening_source_sha256') != row['hardening_source_sha256']
                or metadata.get('wall_mm') != row['wall_mm']
                or metadata.get('duration_seconds') != row['duration_seconds']
                or metadata.get('ams', {}).get('requested_time_step_seconds') != row['ams_time_step']
                or Path(metadata.get('starter', '')).resolve() != Path(row['deck']).resolve()
                or metadata.get('fracture_locus') != config['ductile_damage']['fracture_locus']
                or metadata.get('target_nominal_strain') != config['compression_strain']):
            raise ValueError('Native generated model differs from requested study inputs')
        result = assess(row['metadata'], row['history'], **native_assess_kwargs(row.get('force_mapping')))
        result['legacy_manifest_force_mapping'] = legacy_mapping
        result['force_mapping_recomputed_from_bound_channels'] = True
        if not result.get('run_binding_verified') or not result.get('metadata_binding_verified'):
            raise ValueError('Native execution and metadata bindings are required for resume')
    else:
        from assess_run import assess
        if require_receipt:
            receipt = json.loads(receipt_path(row).read_text())
            if receipt.get('config_sha256') != row['config_sha256']:
                raise ValueError('Execution receipt belongs to a different configuration')
            for path, expected in receipt['evidence_sha256'].items():
                if digest(path) != expected:
                    raise ValueError('CCX output changed since verified execution: '+path)
            required = evidence_paths(row)
            if set(receipt['evidence_sha256']) != {str(path.resolve()) for path in required}:
                raise ValueError('CCX execution receipt has incomplete evidence')
        result = assess(row['deck'])
        if (result.get('generation_config_sha256') != row['config_sha256']
                or result.get('hardening_source_sha256') != row['hardening_source_sha256']
                or result.get('config_wall_mm') != row['wall_mm']):
            raise ValueError('CCX generation does not match requested configuration/material')
    if not result.get('completed'):
        raise ValueError('Actual solver evidence is incomplete or failed')
    return result


def evidence_paths(row):
    deck = Path(row['deck'])
    return [deck, *(deck.with_suffix(suffix) for suffix in
                   ('.material.dat', '.generation.json', '.sta', '.dat', '.frd', '.log')),
            Path(row['assessment'])]


def existing_outputs(row):
    if row.get('backend', 'ccx') == 'opencourant':
        metadata = Path(row['metadata'])
        prefix = metadata.name.removesuffix('.radioss.json')
        matches = list(metadata.parent.glob(prefix+'*'))
    else:
        deck = Path(row['deck'])
        matches = list(deck.parent.glob(deck.stem+'*'))
    if Path(row['assessment']).exists():
        matches.append(Path(row['assessment']))
    return sorted(set(str(path) for path in matches))


def interrupt_process_group(process):
    """Let a child CLI clean detached solver groups before bounded escalation."""
    for requested_signal, grace in ((signal.SIGINT, 8), (signal.SIGTERM, 3), (signal.SIGKILL, 3)):
        try:
            os.killpg(process.pid, requested_signal)
        except ProcessLookupError:
            pass
        try:
            process.wait(timeout=grace)
        except (subprocess.TimeoutExpired, KeyboardInterrupt):
            continue
        # A CLI can exit before its same-group workers. Signal the remaining
        # group too, while native CLIs have already cleaned detached engines.
        if requested_signal == signal.SIGKILL:
            break


def run_logged(argv, log, env, timeout):
    with Path(log).open('x') as stream:
        process = subprocess.Popen(argv, cwd=ROOT, env=env, stdout=stream,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=timeout)
        except KeyboardInterrupt:
            interrupt_process_group(process)
            raise
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            raise RuntimeError(f'Command timed out after {timeout:g}s; partial evidence retained')
    if code:
        raise subprocess.CalledProcessError(code, argv)


def execute(manifest_path, *, candidates=None, role='training', threads=2,
            timeout=3600, dry_run=False, resume=False, keep_going=False):
    if not isinstance(threads, int) or threads < 1 or not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('Require positive threads and a finite positive timeout')
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    rows = select_runs(manifest, candidates, role)
    plans = []
    for row in rows:
        check_inputs(row)
        solve, assess = commands(row, threads, timeout)
        plans.append(dict(candidate=row['candidate'], role=row['role'], solve_argv=solve,
                          assess_argv=assess, existing_outputs=existing_outputs(row)))
        if row.get('backend', manifest.get('backend', 'ccx')) == 'opencourant':
            plans[-1].update(
                legacy_manifest_force_mapping=legacy_native_force_mapping(row.get('force_mapping')),
                derived_force_mapping=native_force_mapping(reference_points=row.get('force_mapping',{}).get('source')=='rp'))
    if dry_run:
        return dict(dry_run=True, runs=plans, calibrated=False, solver_executed=False)
    journal_path = manifest_path.parent/'execution-journal.json'
    journal = json.loads(journal_path.read_text()) if journal_path.exists() else {'events': [], 'calibrated': False}
    env = os.environ.copy()
    for variable in ('OMP_NUM_THREADS', 'CCX_NPROC_STIFFNESS', 'CCX_NPROC_RESULTS', 'CCX_NPROC_EQUATION_SOLVER'):
        env[variable] = str(threads)
    failed = 0
    for row, plan in zip(rows, plans):
        event = dict(candidate=row['candidate'], role=row['role'],
                     manifest_sha256=digest(manifest_path), started_at=datetime.now(timezone.utc).isoformat(),
                     solve_argv=plan['solve_argv'], assess_argv=plan['assess_argv'])
        if 'derived_force_mapping' in plan:
            event.update(derived_force_mapping=plan['derived_force_mapping'],
                         legacy_manifest_force_mapping=plan['legacy_manifest_force_mapping'])
        journal['events'].append(event)
        try:
            if existing_outputs(row):
                if not resume:
                    raise FileExistsError('Existing prefix is preserved; use --resume to verify complete evidence or prepare a fresh study')
                actual = verify_completed(row)
                # Refresh a missing/stale assessment from verified raw output;
                # preserve the completed solver files and their execution hashes.
                write_json(row['assessment'], actual)
                if row.get('backend', 'ccx') == 'ccx':
                    receipt = json.loads(receipt_path(row).read_text())
                    receipt['evidence_sha256'][str(Path(row['assessment']).resolve())] = digest(row['assessment'])
                    write_json(receipt_path(row), receipt)
                event['assessment_refreshed'] = True
                event['status'] = 'skipped_verified_completed'
            else:
                logs = manifest_path.parent/'execution-logs'
                logs.mkdir(parents=True, exist_ok=True)
                identifier = f'{row["candidate"]}_{row["role"]}_{time.time_ns()}'.replace('/', '_')
                event.update(status='running', solve_log=str(logs/(identifier+'-solve.log')),
                             assess_log=str(logs/(identifier+'-assess.log')))
                write_json(journal_path, journal)
                # Native run_logged enforces timeout per starter/engine stage
                # and kills its detached children itself. Permit its cleanup.
                outer_timeout = 3*timeout+120 if row.get('backend') == 'opencourant' else timeout
                run_logged(plan['solve_argv'], event['solve_log'], env, outer_timeout)
                run_logged(plan['assess_argv'], event['assess_log'], env, min(timeout, 600))
                verify_completed(row, require_receipt=False)
                if row.get('backend', 'ccx') == 'ccx':
                    write_json(receipt_path(row), dict(config_sha256=row['config_sha256'],
                        evidence_sha256={str(path.resolve()): digest(path) for path in evidence_paths(row)}))
                event['status'] = 'completed_verified'
        except KeyboardInterrupt:
            event.update(status='interrupted_preserved',
                         error='User interrupt; partial solver evidence was preserved',
                         finished_at=datetime.now(timezone.utc).isoformat())
            write_json(journal_path, journal)
            raise
        except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
            failed += 1
            event.update(status='failed_preserved', error=str(error))
        event['finished_at'] = datetime.now(timezone.utc).isoformat()
        write_json(journal_path, journal)
        print(json.dumps({key: event[key] for key in ('candidate', 'role', 'status')}), flush=True)
        if event['status'] == 'failed_preserved' and not keep_going:
            break
    return dict(failed_runs=failed, journal=str(journal_path), calibrated=False,
                qualification='Execution only; geometry, calibration and explicit quasi-static qualification guards remain unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('--candidate', action='append', help='Repeat to select candidates; omission selects all candidates')
    parser.add_argument('--role', choices=('training', 'holdout', 'all'), default='training')
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--timeout-seconds', type=float, default=3600,
                        help='CCX process timeout; native timeout per solver stage')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--resume', action='store_true', help='Skip only recomputed, bound, completed solver evidence')
    parser.add_argument('--keep-going', action='store_true', help='Continue remaining selected runs after preserving a failed run')
    args = parser.parse_args()
    try:
        result = execute(args.manifest, candidates=args.candidate, role=args.role, threads=args.threads,
                         timeout=args.timeout_seconds, dry_run=args.dry_run, resume=args.resume,
                         keep_going=args.keep_going)
    except KeyboardInterrupt:
        print('Interrupted; partial evidence and the execution journal were preserved.', file=sys.stderr)
        return 130
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2))
    return bool(result.get('failed_runs'))


if __name__ == '__main__':
    sys.exit(main())
