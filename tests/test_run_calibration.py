import csv
import json
from pathlib import Path

import pytest

from calibration_study import NATIVE_FORCE_COLUMN, prepare
from run_calibration import ROOT, commands, digest, execute, select_runs, verify_completed


def test_role_selection_and_dry_run_use_current_python_without_execution(tmp_path, monkeypatch):
    manifest = prepare(ROOT/'configs/dong2019_brick_damage.json', tmp_path/'study')
    selected = select_runs(manifest, ['locus_0_1'], 'all')
    assert {row['role'] for row in selected} == {'training', 'holdout'}
    assert len(select_runs(manifest, ['locus_0_1'], 'training')) == 1
    with pytest.raises(ValueError, match='Unknown candidate'):
        select_runs(manifest, ['does-not-exist'])
    monkeypatch.setattr('run_calibration.run_logged', lambda *args: pytest.fail('dry run launched a process'))
    result = execute(tmp_path/'study/manifest.json', candidates=['locus_0_1'], role='holdout', dry_run=True, threads=3)
    assert not result['solver_executed'] and len(result['runs']) == 1
    solve = result['runs'][0]['solve_argv']
    assert solve[solve.index('--threads')+1] == '3'
    assert not (tmp_path/'study/execution-journal.json').exists()
    selected[0]['solve_argv'][1] = '/tmp/arbitrary.py'
    with pytest.raises(ValueError, match='expected repository script'):
        commands(selected[0], 2, 30)


def native_evidence(tmp_path):
    manifest = prepare(ROOT/'configs/dong2019_brick_damage.json', tmp_path/'study',
                       backend='opencourant', duration=.01, ams_time_step=1e-6)
    row = select_runs(manifest, ['locus_0_1'], 'holdout')[0]
    config = json.loads(Path(row['config']).read_text())
    path = Path(row['metadata'])
    path.parent.mkdir(parents=True)
    prefix = path.name.removesuffix('.radioss.json')
    files = {suffix: path.parent/(prefix+suffix) for suffix in
             ('_0000.rad', '_0001.rad', '_0000.out', '_0001.out', 'T01', 'T01.csv')}
    for target in files.values():
        target.write_text('Synthetic binding test evidence; not a physics result\n')
    files['_0001.out'].write_text('NORMAL TERMINATION\n')
    with files['T01.csv'].open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['time', 'INTERNAL ENERGY', 'KINETIC ENERGY', NATIVE_FORCE_COLUMN])
        for i in range(41):
            t = .01*i/40
            writer.writerow([t, 1, 0, -1900*t])
    inputs = {files[suffix].name: digest(files[suffix]) for suffix in ('_0000.rad', '_0001.rad')}
    metadata = dict(prefix=prefix, backend='OpenCourant', sha256=inputs,
        starter=str(files['_0000.rad']), source_config_sha256=row['config_sha256'],
        hardening_source_sha256=row['hardening_source_sha256'], wall_mm=.3, duration_seconds=.01,
        target_nominal_strain=.65, reference_area_mm2=1900,
        displacement_fraction_samples=[[0, 0], [.01, 1]], time_history_interval_seconds=.00025,
        ams={'enabled': True, 'requested_time_step_seconds': 1e-6},
        fracture_locus=config['ductile_damage']['fracture_locus'])
    path.write_text(json.dumps(metadata))
    binding = dict(input_sha256=inputs, metadata_sha256=digest(path),
        output_sha256={files[suffix].name: digest(files[suffix]) for suffix in
                       ('_0000.out', '_0001.out', 'T01', 'T01.csv')})
    record = path.parent/(prefix+'.run.json')
    record.write_text(json.dumps(binding))
    Path(row['assessment']).write_text(json.dumps({'completed': False}))
    return row, files, record, binding


def test_resume_requires_actual_bound_completed_outputs_not_cached_flags(tmp_path, monkeypatch):
    row, files, record, binding = native_evidence(tmp_path)
    monkeypatch.setattr('run_calibration.run_logged', lambda *args: pytest.fail('verified resume launched a process'))
    actual = verify_completed(row)
    assert actual['completed'] and not actual['quasi_static_validated']
    result = execute(tmp_path/'study/manifest.json', candidates=['locus_0_1'], role='holdout', resume=True)
    assert result['failed_runs'] == 0 and not result['calibrated']
    assert json.loads(Path(row['assessment']).read_text())['completed']
    journal = json.loads(Path(result['journal']).read_text())
    assert journal['events'][-1]['status'] == 'skipped_verified_completed'
    # A stale successful cache cannot override a subsequently failed engine.
    Path(row['assessment']).write_text(json.dumps({'completed': True}))
    files['_0001.out'].write_text('**ERROR : FAILED\nNORMAL TERMINATION\n')
    binding['output_sha256'][files['_0001.out'].name] = digest(files['_0001.out'])
    record.write_text(json.dumps(binding))
    before = {str(path): digest(path) for path in files.values()}
    result = execute(tmp_path/'study/manifest.json', candidates=['locus_0_1'], role='holdout', resume=True)
    assert result['failed_runs'] == 1
    assert before == {str(path): digest(path) for path in files.values()}
    with pytest.raises(ValueError, match='incomplete or failed'):
        verify_completed(row)


def test_replaced_history_cannot_be_skipped_or_overwritten(tmp_path, monkeypatch):
    row, files, _, _ = native_evidence(tmp_path)
    files['T01.csv'].write_text(files['T01.csv'].read_text().replace('-19.0', '-20.0'))
    monkeypatch.setattr('run_calibration.run_logged', lambda *args: pytest.fail('partial run was overwritten'))
    result = execute(tmp_path/'study/manifest.json', candidates=['locus_0_1'], role='holdout', resume=True)
    assert result['failed_runs'] == 1
    journal = json.loads(Path(result['journal']).read_text())
    assert 'changed since execution' in journal['events'][-1]['error']
