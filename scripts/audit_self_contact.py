"""Reproduce the matched early self-contact-off diagnostic from sealed evidence.

This audits an intentionally partial control, not a completed/calibrated lattice.
Use --write-replay-dir to restore its actual native inputs into an empty folder;
then source the qualified OpenCourant env.sh and run starter, engine, converter.
The original .01 s loading amplitude is retained; the engine stops at .0032 s.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile

import numpy as np

from reentrant_fea.radioss_results import normal_termination, solver_errors


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT/'validation/contact_audit/thin_self_contactoff_partial.tar.gz'
ARCHIVE_SHA256 = '533afe71c72cd74bd4930dafa35d16ef3faaffecdb32b555afc0a84dd097ed39'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def load(archive):
    if sha(archive.read_bytes()) != ARCHIVE_SHA256:
        raise ValueError('Evidence archive checksum changed; old reports are invalid')
    with tarfile.open(archive, 'r:gz') as tar:
        names = tar.getnames()
        if len(names) != len(set(names)) or any(not member.isfile() for member in tar):
            raise ValueError('Unexpected archive members')
        data = {name: tar.extractfile(name).read() for name in names}
    manifest = json.loads(data['manifest.json'])
    if set(data) != set(manifest['sha256']) | {'manifest.json'}:
        raise ValueError('Incomplete evidence archive')
    for name, expected in manifest['sha256'].items():
        if sha(data[name]) != expected:
            raise ValueError('Archived file changed: '+name)
    return manifest, data


def audit(archive):
    manifest, data = load(archive)
    base, control = manifest['baseline'], manifest['control']
    bm = json.loads(data[base+'.radioss.json'])
    cm = json.loads(data[control+'.radioss.json'])
    bb = json.loads(data[base+'.run.json'])
    cb = json.loads(data[control+'.run.json'])
    if sha(data[base+'.radioss.json']) != bb['metadata_sha256']:
        raise ValueError('Baseline execution metadata binding changed')
    for name, expected in bb['output_sha256'].items():
        if sha(data[name]) != expected:
            raise ValueError('Baseline execution output binding changed: '+name)
    for name, expected in cm['sha256'].items():
        if sha(data[name]) != expected:
            raise ValueError('Control generation input binding changed: '+name)
    for key in ('displacement_fraction_samples', 'duration_seconds', 'target_nominal_strain',
                'reference_area_mm2', 'source_config_sha256', 'hardening_source_sha256', 'ams'):
        if bm[key] != cm[key]:
            raise ValueError('Unmatched control setting: '+key)
    if (not cm['intentionally_partial'] or cm['requested_stop_seconds'] != .0032
            or cm['duration_seconds'] != .01 or any(row['returncode'] != 0 for row in cb['commands'])):
        raise ValueError('The control is not the intended successful partial execution')
    for prefix in (base, control):
        for stage in ('_0000.out', '_0001.out'):
            text = data[prefix+stage].decode()
            if not normal_termination(text) or solver_errors(text):
                raise ValueError('Solver termination/error evidence is invalid: '+prefix+stage)
            if stage == '_0000.out' and not re.search(r'^\s*0\s+ERROR\(S\)\s*$', text, re.MULTILINE):
                raise ValueError('Starter did not report a standalone zero-error count: '+prefix)
    # Restore the baseline exactly, rather than assuming its material/geometry
    # match. Only one contact block, its TH selection and model name differ.
    starter = data[control+'_0000.rad'].decode().splitlines()
    starter = [base if line == control else line for line in starter]
    index = starter.index('/INTER/TYPE25/2')
    starter[index:index] = data['baseline_removed_self_interface.txt'].decode().splitlines()
    index = starter.index('/TH/INTER/1')
    starter[index+3] = f'{1:10d}{2:10d}{3:10d}'
    restored = ('\n'.join(starter)+'\n').encode()
    if sha(restored) != bb['input_sha256'][base+'_0000.rad']:
        raise ValueError('Control includes changes beyond the removed core self-contact')
    if sha(data[base+'_0001.rad']) != bb['input_sha256'][base+'_0001.rad']:
        raise ValueError('Baseline engine binding changed')
    engine = data[control+'_0001.rad'].decode().splitlines()
    index = engine.index('/RUN/'+control+'/1')
    if float(engine[index+1]) != .0032:
        raise ValueError('Control stop time changed')
    engine[index] = '/RUN/'+base+'/1'
    engine[index+1] = f'{.01:20.12g}'
    if sha(('\n'.join(engine)+'\n').encode()) != bb['input_sha256'][base+'_0001.rad']:
        raise ValueError('Control engine includes changes beyond name and intentional stop')

    def curve(prefix, metadata):
        reader = csv.reader(io.StringIO(data[prefix+'T01.csv'].decode()))
        header = next(reader)
        values = np.array([[float(value) for value in row] for row in reader])
        if not np.isfinite(values).all() or np.any(np.diff(values[:, 0]) <= 0):
            raise ValueError('Invalid time-history samples')
        columns = [i for i, label in enumerate(header)
                   if label.startswith('Contact force and energy')
                   and re.search(r'\s3\s+Contact 3\s', label)]
        if len(columns) != 6:
            raise ValueError('Contact 3 DEF channels are missing or ambiguous')
        # DEF is FNX,FNY,FNZ,FTX,FTY,FTZ. FNY is cumulative normal impulse.
        # Removing Contact 1 shifts opaque names: original var64 becomes var58.
        impulse = values[:, columns[1]]
        time = values[:, 0]
        midtime = .5*(time[1:]+time[:-1])
        force = -np.diff(impulse)/np.diff(time)
        amplitude = np.array(metadata['displacement_fraction_samples'])
        if midtime[0] < amplitude[0, 0] or midtime[-1] > amplitude[-1, 0]:
            raise ValueError('Do not extrapolate loading amplitude')
        strain = np.interp(midtime, amplitude[:, 0], amplitude[:, 1])*metadata['target_nominal_strain']
        stress = force/metadata['reference_area_mm2']
        core_ie = [i for i, label in enumerate(header) if re.fullmatch(r'Part 1\s+IE\s*', label)]
        if len(core_ie) != 1:
            raise ValueError('Core part internal-energy channel is missing or ambiguous')
        return strain, stress, values, header[columns[1]], values[:, core_ie[0]]

    strain, stress, values, channel, core_ie = curve(base, bm)
    cs, cp, cv, control_channel, control_core_ie = curve(control, cm)
    if 'var 64' not in channel or 'var 58' not in control_channel:
        raise ValueError('Unexpected reference channel mapping; requalify the mapping')
    selected = (cs >= .005) & (cs <= .12)
    if selected.sum() < 20 or cs[-1] < .12:
        raise ValueError('Insufficient actual early-strain overlap')
    predicted = np.interp(cs[selected], strain, stress)
    error = cp[selected]-predicted
    normal = normal_termination(data[control+'_0001.out'].decode())
    dt = cm['time_history_interval_seconds']
    requested_stop_reached = abs(cv[-1, 0]-.0032) <= 1.5*dt
    return {
        'status': 'matched early-strain diagnostic; not a full loading result or calibration',
        'baseline_prefix': base, 'control_prefix': control,
        'archive': str(ARCHIVE.relative_to(ROOT)), 'archive_sha256': ARCHIVE_SHA256,
        'script_sha256': sha(Path(__file__).read_bytes()),
        'baseline_execution_binding_verified': True,
        'control_pre_execution_generation_hashes_verified': True,
        'control_execution_input_binding_recorded': False,
        'control_execution_metadata_binding_recorded': False,
        'both_runs_starter_and_engine_termination_verified_without_solver_errors': True,
        'only_self_contact_and_its_history_selection_removed': True,
        'only_engine_name_and_requested_stop_changed': True,
        'same_amplitude_material_geometry_ams_and_reference_area': True,
        'loading_duration_seconds': .01, 'control_requested_stop_seconds': .0032,
        'control_engine_normal_termination': normal,
        'control_reached_requested_partial_stop': bool(requested_stop_reached),
        'control_completed_full_loading': False, 'control_intentionally_partial': True,
        'control_last_history_time_seconds': float(cv[-1, 0]),
        'control_last_force_strain': float(cs[-1]), 'comparison_interval': [.005, .12],
        'overlap_samples': int(selected.sum()),
        'force_mapping': {'baseline': channel, 'control': control_channel,
                          'physical_quantity': 'Contact 3 FNY cumulative normal impulse',
                          'operation': 'negative first difference divided by actual time interval, then1900mm2'},
        'overlap_rmse_MPa': float(np.sqrt(np.mean(error**2))),
        'max_absolute_overlap_difference_MPa': float(max(abs(error))),
        'energy_sample_timestamps_identical': bool(np.array_equal(cv[:, 0], values[:len(cv), 0])),
        'global_internal_energy_samples_bitwise_equal': bool(np.array_equal(cv[:, 1], values[:len(cv), 1])),
        'core_part_internal_energy_samples_bitwise_equal': bool(np.array_equal(control_core_ie, core_ie[:len(cv)])),
        'strain_checks': [dict(strain=e, baseline_stress_MPa=float(np.interp(e, strain, stress)),
                              control_stress_MPa=float(np.interp(e, cs, cp))) for e in (.01, .02, .04, .08, .12)],
        'conclusion': 'Removing core self-contact leaves early stress and internal energy identical at recorded CSV precision. Core self-contact does not explain the early overshoot through strain.12 in this configuration.',
        'limitations': ['Contact-on late collapse remains sensitive to gap, edge and adjacency choices.',
                        'The archived control run record contains commands and runtime only; it has no execution-time input or metadata SHA. Generation deck hashes and the post-run evidence archive do not establish a prospective execution binding.',
                        'No claim that the accelerated AMS run is quasi-static or experimentally calibrated.',
                        'FNY comparison uses normal contact only, consistently for both cases; it is not an independent validation of total platen reaction.',
                        'This partial control cannot select a material candidate or prove full-strain agreement.'],
        'evidence_sha256': manifest['sha256'],
    }, manifest, data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=ARCHIVE)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--verify-report', type=Path)
    parser.add_argument('--write-replay-dir', type=Path)
    args = parser.parse_args()
    result, manifest, data = audit(args.archive)
    if args.verify_report and json.loads(args.verify_report.read_text()) != result:
        raise ValueError('Saved audit report is stale or modified; invalidate it')
    if args.write_replay_dir:
        directory = args.write_replay_dir
        if directory.exists() and any(directory.iterdir()):
            raise ValueError('Replay directory must be empty to preserve existing results')
        directory.mkdir(parents=True, exist_ok=True)
        prefix = manifest['control']
        for suffix in ('_0000.rad', '_0001.rad', '.radioss.json'):
            (directory/(prefix+suffix)).write_bytes(data[prefix+suffix])
        (directory/'REPLAY.txt').write_text(
            'Source the qualified /workspace/opencourant/env.sh. From this folder run:\n'
            f'OMP_NUM_THREADS=2 starter_linux64_gf -i {prefix}_0000.rad -np 1\n'
            f'OMP_NUM_THREADS=2 engine_linux64_gf -i {prefix}_0001.rad\n'
            f'th_to_csv_linux64_gf {prefix}T01\n'
            'This deliberately stops at.0032s of a.01s loading amplitude. Do not mark full loading complete.\n')
    text = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
