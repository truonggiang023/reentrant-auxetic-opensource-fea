"""Run/export a bound thick-wall diagnostic; never declare calibration.

Archives contain an exact solver-file allowlist, not logs or runtime/environment
directories. --export-only verifies an existing completed case without solving.
"""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tarfile

from assess_courant import assess
from reentrant_fea.calibration import curve_objective
from reentrant_fea.geometry import core_polygon
from scripts.compare_wall_controls import compare, write_outputs


ROOT = Path(__file__).resolve().parents[1]
SOLVER_SUFFIXES = ('_0000.rad', '_0001.rad', '_0000.out', '_0001.out',
                   '.radioss.json', '.run.json', 'T01', 'T01.csv', '.assessment.json')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def export_completed(prefix, config, destination, label, source_commit):
    """Verify completed receipts and export only exact, regular solver files."""
    prefix, config, destination = map(Path, (prefix, config, destination))
    if not re.fullmatch(r'[A-Za-z0-9_-]{1,60}', prefix.name):
        raise ValueError('Require a safe solver prefix')
    metadata_file = Path(str(prefix)+'.radioss.json')
    metadata = json.loads(metadata_file.read_text())
    if metadata.get('prefix') != prefix.name or metadata.get('wall_mm') != 1.:
        raise ValueError('Require the declared thick-wall prefix')
    if metadata.get('source_config_sha256') != sha(config):
        raise ValueError('Configuration does not match executed metadata')
    # assess rechecks every generated input, output and metadata receipt hash.
    actual = assess(metadata_file, Path(str(prefix)+'T01.csv'), force_source='rp')
    if not (actual.get('completed') and actual.get('run_binding_verified')
            and actual.get('metadata_binding_verified')):
        raise ValueError('Only completed receipt-bound runs can be published')
    assessment = Path(str(prefix)+'.assessment.json')
    assessment.write_text(json.dumps(actual, indent=2)+'\n')
    paths = [Path(str(prefix)+suffix) for suffix in SOLVER_SUFFIXES]
    for path in paths + [config]:
        if path.is_symlink() or not path.is_file():
            raise ValueError('Evidence must be a regular file: '+str(path))
    if destination.exists() and any(destination.iterdir()):
        raise FileExistsError('Use a fresh publication directory')
    destination.mkdir(parents=True, exist_ok=True)
    report, curves, references = compare([(label, prefix)])
    case = report['cases'][0]
    case['full_curve_paper_discrepancy'] = curve_objective(
        curves[label], references['thick'], (.04, .45))
    report.update(source_commit=source_commit,
                  protocol='Accelerated mean crosshead speed 576 mm/s; not physical 2 mm/min.',
                  workflow_run_id=os.environ.get('GITHUB_RUN_ID'),
                  calibrated=False)
    write_outputs(report, curves, references, destination/'comparison')
    # Basenames retain the original prefix, so assess can recheck an extracted
    # archive without changing the executed metadata's original path strings.
    entries = {path.name: path for path in paths}
    entries['source_config.json'] = config
    manifest = dict(source_commit=source_commit, calibrated=False,
                    solver_prefix=prefix.name,
                    qualification='Completed diagnostic only; review energy, speed, AMS and mesh controls.',
                    files={name: sha(path) for name, path in entries.items()})
    payload = (json.dumps(manifest, indent=2)+'\n').encode()
    archive_path = destination/'solver_evidence.tar.gz'
    with tarfile.open(archive_path, 'w:gz') as archive:
        for name, path in entries.items():
            archive.add(path, arcname=name, recursive=False)
        info = tarfile.TarInfo('manifest.json')
        info.size = len(payload)
        archive.addfile(info, io.BytesIO(payload))
    with tarfile.open(archive_path) as archive:
        if set(archive.getnames()) != set(entries) | {'manifest.json'}:
            raise ValueError('Unexpected evidence archive entry')
        for name, expected in manifest['files'].items():
            if hashlib.sha256(archive.extractfile(name).read()).hexdigest() != expected:
                raise ValueError('Archive evidence hash mismatch')
    (destination/'manifest.json').write_bytes(payload)
    summary = dict(case=label, completed=True, calibrated=False,
                   observed_end_strain=actual['last_observed_nominal_strain'],
                   last_stress_MPa=actual['reaction_history'][-1]['nominal_stress_MPa'],
                   core_energy_screen_passed=actual['core_energy_screen']['sampled_energy_screen_passed'],
                   global_energy_screen_passed=actual['global_energy_screen']['sampled_energy_screen_passed'],
                   archive_sha256=sha(archive_path))
    (destination/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print('::notice title=Thick-wall diagnostic::'+json.dumps(summary), flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--prefix', type=Path, required=True)
    parser.add_argument('--publish-dir', type=Path, required=True)
    parser.add_argument('--source-commit', required=True)
    parser.add_argument('--timeout-seconds', type=float, default=7200.)
    parser.add_argument('--export-only', action='store_true')
    args = parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9_-]{1,60}', args.case):
        parser.error('Require a safe case label')
    if not re.fullmatch(r'[0-9a-f]{40}', args.source_commit):
        parser.error('Require the complete source commit identity')
    if not args.export_only:
        config = json.loads(args.config.read_text())
        if config.get('wall') != 1. or config.get('compression_strain') != .5:
            parser.error('This trial requires wall1 mm and target strain0.5')
        bounds = core_polygon(config).bounds
        duration = (bounds[3]-bounds[1])*.5/576.
        argv = [sys.executable, str(ROOT/'run_radioss.py'), '--config', str(args.config),
                '--output', str(args.prefix), '--duration-seconds', str(duration),
                '--ams-time-step', '1e-6', '--drop-ccx-viscosity', '--run',
                '--timeout-seconds', str(args.timeout_seconds)]
        subprocess.run(argv, cwd=ROOT, check=True)
    export_completed(args.prefix, args.config, args.publish_dir, args.case, args.source_commit)


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        message = str(error)[:3000].replace('%', '%25').replace('\r', '%0D').replace('\n', '%0A')
        print('::error title=Thick-wall trial failed::'+message, flush=True)
        raise
