import csv
import json
from pathlib import Path

import pytest

from calibration_study import NATIVE_FORCE_COLUMN, NATIVE_TANGENTIAL_COLUMN, prepare
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
        writer.writerow(['time', 'INTERNAL ENERGY', 'KINETIC ENERGY', NATIVE_FORCE_COLUMN,
                         NATIVE_TANGENTIAL_COLUMN])
        for i in range(41):
            t = .01*i/40
            writer.writerow([t, 1, 0, -1900*t, 0])
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


def test_legacy_normal_only_resume_recomputes_total_without_changing_raw_evidence(tmp_path, monkeypatch):
    row, files, record, binding = native_evidence(tmp_path)
    # Construct an historical raw history with nonzero FT and a normal-only
    # cached curve. The test execution binding refers to both actual columns.
    with files['T01.csv'].open(newline='') as stream:
        rows = list(csv.reader(stream))
    for values in rows[1:]:
        values[-1] = str(950*float(values[0]))
    with files['T01.csv'].open('w', newline='') as stream:
        csv.writer(stream).writerows(rows)
    binding['output_sha256'][files['T01.csv'].name] = digest(files['T01.csv'])
    record.write_text(json.dumps(binding))
    manifest_path = tmp_path/'study/manifest.json'
    manifest = json.loads(manifest_path.read_text())
    historical = next(item for item in manifest['runs']
                      if item['candidate'] == 'locus_0_1' and item['role'] == 'holdout')
    historical['force_mapping'] = dict(csv_column=NATIVE_FORCE_COLUMN,
                                      compressive_sign=-1, channel_kind='impulse')
    argv = historical['assess_argv']
    index = argv.index('--tangential-column')
    del argv[index:index+2]
    manifest_path.write_text(json.dumps(manifest))
    Path(row['assessment']).write_text(json.dumps(dict(
        completed=True, force_mapping={'contact_force_quantity': 'normal_component_only'},
        reaction_history=[{'strain': 0, 'nominal_stress_MPa': 1},
                          {'strain': .65, 'nominal_stress_MPa': 1}])))
    before = {str(path): digest(path) for path in [*files.values(), record, manifest_path]}
    monkeypatch.setattr('run_calibration.run_logged', lambda *args: pytest.fail('resume reran physics'))
    result = execute(manifest_path, candidates=['locus_0_1'], role='holdout', resume=True)
    assert result['failed_runs'] == 0
    assert before == {str(path): digest(path) for path in [*files.values(), record, manifest_path]}
    actual = json.loads(Path(row['assessment']).read_text())
    assert actual['legacy_manifest_force_mapping']
    assert actual['force_mapping']['component_operation'] == 'FN-FT'
    assert all(point['nominal_stress_MPa'] == pytest.approx(1.5) for point in actual['reaction_history'])
    event = json.loads(Path(result['journal']).read_text())['events'][-1]
    assert event['legacy_manifest_force_mapping']
    assert event['assess_argv'][event['assess_argv'].index('--tangential-column')+1] == NATIVE_TANGENTIAL_COLUMN


@pytest.mark.parametrize('mistake', ['missing-new-tangential-channel', 'normal-used-as-tangent',
                                    'unverified-legacy-tangent', 'sum-instead-of-difference'])
def test_native_commands_reject_misdeclared_total_contact_force(tmp_path, mistake):
    """A study cannot silently fall back to normal-only or FN+FT scoring."""
    manifest = prepare(ROOT/'configs/dong2019_brick_damage.json', tmp_path/'study',
                       backend='opencourant', duration=.01, ams_time_step=1e-6)
    row = select_runs(manifest, ['locus_0_1'], 'holdout')[0]
    argv = row['assess_argv']
    index = argv.index('--tangential-column')
    if mistake == 'missing-new-tangential-channel':
        del argv[index:index+2]
    elif mistake == 'normal-used-as-tangent':
        argv[index+1] = NATIVE_FORCE_COLUMN
    elif mistake == 'unverified-legacy-tangent':
        row['force_mapping'] = dict(csv_column=NATIVE_FORCE_COLUMN,
                                   compressive_sign=-1, channel_kind='impulse')
        argv[index+1] = 'arbitrary friction column'
    else:
        row['force_mapping']['component_operation'] = 'FN+FT'
    with pytest.raises(ValueError, match='verified'):
        commands(row, 2, 30)


def test_rp_study_keeps_observed_rp_assessment_through_runner(tmp_path):
    import json
    from pathlib import Path
    from calibration_study import prepare,ROOT,native_force_mapping,native_assess_kwargs
    from run_calibration import commands
    cfg=json.loads((ROOT/'configs/dong2019_brick_damage.json').read_text())
    cfg['platen_reference_points']=True
    source=tmp_path/'rp_config.json';source.write_text(json.dumps(cfg))
    study=prepare(source,tmp_path/'study',backend='opencourant',duration=.04,holdout_end_strain=.78)
    for row in study['runs']:
        assert row['force_mapping']==native_force_mapping(reference_points=True)
        assert native_assess_kwargs(row['force_mapping'])=={'force_source':'rp'}
        solve,assessment=commands(row,2,3600)
        assert '--force-source' in assessment and 'rp' in assessment
        assert '--force-column' not in assessment and '--tangential-column' not in assessment
