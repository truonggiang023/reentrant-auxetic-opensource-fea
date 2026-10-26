import numpy as np
import pytest

from reentrant_fea.calibration import (fracture_candidate, curve_objective,
                                      score_assessment, select_candidate, sensitivity_rank)


def test_fracture_fit_preserves_published_point_and_does_not_change_input():
    base = {'ductile_damage': {'fracture_locus': [], 'calibrated': True}}
    candidate = fracture_candidate(base, -.1, .8)
    assert candidate['ductile_damage']['fracture_locus'][1] == [.333, .047]
    assert not candidate['ductile_damage']['calibrated']
    assert candidate['ductile_damage']['failure_displacement'] == .001
    assert base['ductile_damage']['fracture_locus'] == []
    with pytest.raises(ValueError, match='bracket'):
        fracture_candidate(base, .4, .8)


def test_partial_perfect_curve_cannot_beat_completed_candidate():
    experiments = [[[0, 0], [.65, .65]], [[0, .02], [.65, .67]]]
    partial = {'completed': False, 'errors': [], 'last_converged_nominal_strain': .14,
               'reaction_history': [{'strain': 0, 'nominal_stress_MPa': .01},
                                    {'strain': .14, 'nominal_stress_MPa': .15}]}
    score = score_assessment(partial, experiments)
    assert score['interval_scored'] == [.02, .14]
    assert score['coverage_fraction'] == pytest.approx(.12/.63)
    assert not score['eligible_for_selection']
    full = curve_objective([[0, .1], [.65, .75]], experiments)
    assert full['overlap_rmse_MPa'] > score['overlap_rmse_MPa']
    full['eligible_for_selection'] = True
    selection = select_candidate([
        {'candidate': 'partial', 'role': 'training', 'wall_mm': 1., 'score': score},
        {'candidate': 'full', 'role': 'training', 'wall_mm': 1., 'score': full}])
    assert selection['selected_candidate'] == 'full'
    assert not selection['calibrated']
    assert select_candidate([{'candidate': 'partial', 'role': 'training',
                              'wall_mm': 1., 'score': score}])['selected_candidate'] is None
    with pytest.raises(ValueError, match='holdout'):
        select_candidate([{'candidate': 'thin', 'role': 'holdout', 'wall_mm': .3, 'score': full}])


def test_no_extrapolation_or_bad_reference_data():
    with pytest.raises(ValueError, match='no usable overlap'):
        curve_objective([[.4, .3], [.5, .3]], [[[0, 0], [.2, .2]]])
    for bad in ([[0, 0], [0, 1]], [[0, 0], [.5, float('nan')]], [[1, 0], [0, 1]]):
        with pytest.raises(ValueError, match='strictly increasing'):
            curve_objective([[0, 0], [1, 1]], [bad])


def test_sensitivity_detects_non_identifiable_parameter_combinations():
    rank = sensitivity_rank([[1, 2], [2, 4], [3, 6]], [1, .1])
    assert rank['rank'] == 1 and not rank['locally_identifiable']
    rank = sensitivity_rank(np.eye(2), [100, .01], tolerance=1e-5)
    assert rank['rank'] == 2 and rank['locally_identifiable']
    assert not sensitivity_rank(np.zeros((3, 2)), [1, 1])['locally_identifiable']


def test_study_intervals_respect_thick_and_thin_reference_support(tmp_path):
    from pathlib import Path
    from calibration_study import ROLE_INTERVALS, prepare, reference_curves
    root = Path(__file__).resolve().parents[1]
    for role, panel in (('training', 'thick'), ('holdout', 'thin')):
        for values in reference_curves(root/f'validation/paper_reference/dong2019_fig12_{panel}.csv'):
            assert values[0][0] <= ROLE_INTERVALS[role][0]
            assert values[-1][0] >= ROLE_INTERVALS[role][1]
    manifest = prepare(root/'configs/dong2019_brick_damage.json', tmp_path/'study')
    import json
    for row in manifest['runs']:
        config = json.loads(Path(row['config']).read_text())
        assert config['compression_strain'] == (.5 if row['role'] == 'training' else .65)
        assert config['wall'] == (1 if row['role'] == 'training' else .3)
    # A completed thick run at .5 covers its own fit interval; it need not
    # extend to the thin-wall .65 validation target.
    assessment = {'completed': True, 'errors': [], 'last_converged_nominal_strain': .5,
                  'reaction_history': [{'strain': 0, 'nominal_stress_MPa': 0},
                                       {'strain': .5, 'nominal_stress_MPa': 10}]}
    refs = reference_curves(root/'validation/paper_reference/dong2019_fig12_thick.csv')
    assert score_assessment(assessment, refs, ROLE_INTERVALS['training'])['eligible_for_selection']


def test_native_study_records_backend_and_nontransferable_assumptions(tmp_path):
    from pathlib import Path
    from calibration_study import NATIVE_FORCE_COLUMN, prepare
    root = Path(__file__).resolve().parents[1]
    with pytest.raises(ValueError, match='explicit finite positive duration'):
        prepare(root/'configs/dong2019_brick_damage.json', tmp_path/'missing-duration', backend='opencourant')
    manifest = prepare(root/'configs/dong2019_brick_damage.json', tmp_path/'study',
                       backend='opencourant', duration=.01, ams_time_step=1e-6)
    assert len(manifest['runs']) == 22
    assert manifest['fit_intervals']['holdout'] == [.03, .64]
    assert not manifest['native_qualification']['qualified']
    assert len({row['metadata'] for row in manifest['runs']}) == len(manifest['runs'])
    assert len({row['assessment'] for row in manifest['runs']}) == len(manifest['runs'])
    for row in manifest['runs']:
        assert row['assessment'].endswith('.assessment.json')
        assert row['metadata'].endswith('.radioss.json')
        assert '--drop-ccx-viscosity' in row['solve_argv']
        assert row['force_mapping'] == dict(csv_column=NATIVE_FORCE_COLUMN,
                                          compressive_sign=-1, channel_kind='impulse')
        assert not row['candidate'].startswith(('viscosity', 'residual'))
        prefix = Path(row['metadata']).name.removesuffix('.radioss.json')
        assert all(c.isalnum() or c in '_-' for c in prefix)
    assert any(row['duration_seconds'] == .02 for row in manifest['runs'])
    assert any(row['ams_time_step'] == 5e-7 for row in manifest['runs'])
    # Exercise a decimal-valued candidate through the actual deck writer;
    # metadata/assessment paths must preserve its full unique model prefix.
    import json
    from reentrant_fea.radioss import write_radioss_deck
    row = manifest['runs'][0]
    config = json.loads(Path(row['config']).read_text())
    nodes = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],
                      [0,0,1],[1,0,1],[1,1,1],[0,1,1]], float)
    prefix = Path(row['solve_argv'][row['solve_argv'].index('--output')+1])
    metadata = write_radioss_deck(prefix, nodes, np.arange(8).reshape(1,8), config,
                                 [(225.,0.), (300.,1.)], duration=row['duration_seconds'],
                                 ams_time_step=row['ams_time_step'], drop_ccx_viscosity=True)
    assert Path(row['metadata']).is_file()
    assert metadata['starter'] == row['deck']


def test_native_score_recomputes_bound_history_and_keeps_physics_gate(tmp_path):
    """Forged JSON completion/curve cannot replace actual engine/history evidence."""
    import csv
    import hashlib
    import json
    from calibration_study import NATIVE_FORCE_COLUMN, score_native_record
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    config = {'ductile_damage': {'fracture_locus': [[0, .08], [.333, .047], [1, .018]],
                                'viscosity_time': .005}}
    source = tmp_path/'hardening.py'
    source.write_text('PLASTIC_DATA = ((225, 0), (300, 1))\n')
    files = {suffix: tmp_path/('model'+suffix) for suffix in
             ('_0000.rad', '_0001.rad', '_0000.out', '_0001.out', 'T01', 'T01.csv')}
    for path in files.values():
        path.write_text('synthetic evidence for binding test only\n')
    files['_0001.out'].write_text('NORMAL TERMINATION\n')
    with files['T01.csv'].open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['time', 'INTERNAL ENERGY', 'KINETIC ENERGY', NATIVE_FORCE_COLUMN])
        for time in np.linspace(0, 1, 41):
            writer.writerow([time, 1, 0, -1900*time])
    sha = {files[suffix].name: digest(files[suffix]) for suffix in ('_0000.rad', '_0001.rad')}
    metadata = dict(prefix='model', backend='OpenCourant', sha256=sha,
                    starter=str(files['_0000.rad']), source_config_sha256='config-hash',
                    hardening_source_sha256=digest(source), wall_mm=.3, duration_seconds=1,
                    time_history_interval_seconds=.025, target_nominal_strain=.65,
                    reference_area_mm2=1900, displacement_fraction_samples=[[0, 0], [1, 1]],
                    fracture_locus=config['ductile_damage']['fracture_locus'], ccx_viscosity_omitted=.005)
    metadata_path = tmp_path/'model.radioss.json'
    metadata_path.write_text(json.dumps(metadata))
    binding = dict(input_sha256=sha, output_sha256={files[suffix].name: digest(files[suffix])
                   for suffix in ('_0000.out', '_0001.out', 'T01', 'T01.csv')})
    (tmp_path/'model.run.json').write_text(json.dumps(binding))
    row = dict(metadata=str(metadata_path), history=str(files['T01.csv']),
               hardening_source=str(source), hardening_source_sha256=digest(source),
               config_sha256='config-hash', wall_mm=.3, deck=str(files['_0000.rad']),
               duration_seconds=1, ams_time_step=None,
               force_mapping=dict(csv_column=NATIVE_FORCE_COLUMN, compressive_sign=-1, channel_kind='impulse'))
    cached = dict(metadata_sha256=digest(metadata_path), history_sha256=digest(files['T01.csv']),
                  completed=False, reaction_history=[{'strain': 0, 'nominal_stress_MPa': 999}])
    references = [[[0, 1], [.65, 1]], [[0, 1], [.65, 1]]]
    result = score_native_record(row, config, cached, references, [.03, .64])
    assert result['solver_completed'] and result['coverage_and_completion_passed']
    assert result['overlap_rmse_MPa'] < 1e-12
    assert result['completion_basis'] == 'observed explicit-solver history'
    assert not result['eligible_for_selection']
    assert result['qualification_blockers']
    # Updating a cache hash after replacing CSV still cannot defeat the
    # independent execution binding to the original actual outputs.
    files['T01.csv'].write_text(files['T01.csv'].read_text().replace('-1900.0', '-2000.0'))
    cached['history_sha256'] = digest(files['T01.csv'])
    with pytest.raises(ValueError, match='changed since execution'):
        score_native_record(row, config, cached, references, [.03, .64])
