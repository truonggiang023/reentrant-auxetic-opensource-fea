"""Common crosshead speed must survive different specimen envelopes/targets."""
import json
from pathlib import Path
import subprocess
import sys

import pytest

from calibration_study import prepare
from reentrant_fea.geometry import core_polygon


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'configs/dong2019_brick_damage.json'


@pytest.mark.parametrize('end_strain', [.65, .78])
def test_physical_speed_uses_each_actual_specimen_height_and_target(tmp_path, end_strain):
    study = prepare(BASE, tmp_path/'study', backend='opencourant',
                    speed_mm_min=2., holdout_end_strain=end_strain)
    protocol = study['loading_protocol']
    assert protocol['physical_reference_speed_mm_min'] == 2.
    assert protocol['acceleration_factor'] == 1.
    assert protocol['loading_profile'] == 'linear'
    assert not study['geometry_qualification']['qualified']
    assert not study['native_qualification']['qualified']
    assert protocol['roles']['training']['actual_envelope_height_mm'] == pytest.approx(48.7)
    assert protocol['roles']['holdout']['actual_envelope_height_mm'] == pytest.approx(48.)
    candidate_rows = [row for row in study['runs'] if row['candidate'] == 'locus_0_1']
    assert len(candidate_rows) == 2
    durations = []
    for row in candidate_rows:
        config = json.loads(Path(row['config']).read_text())
        bounds = core_polygon(config).bounds
        height = bounds[3]-bounds[1]
        measured_mean_speed = height*config['compression_strain']/row['duration_seconds']
        assert measured_mean_speed == pytest.approx(2./60)
        assert row['loading_protocol']['actual_envelope_height_mm'] == pytest.approx(height)
        assert row['loading_protocol']['actual_mean_crosshead_speed_mm_per_s'] == pytest.approx(measured_mean_speed)
        assert row['loading_protocol']['actual_peak_crosshead_speed_mm_per_s'] == pytest.approx(measured_mean_speed)
        assert row['loading_protocol']['constant_physical_reference_speed_applied']
        assert config['loading_profile'] == 'linear'
        assert float(row['solve_argv'][row['solve_argv'].index('--duration-seconds')+1]) == row['duration_seconds']
        assert '--physical-speed-mm-min' not in row['solve_argv']
        durations.append(row['duration_seconds'])
    assert durations[0] != durations[1]
    thin = next(row for row in candidate_rows if row['role'] == 'holdout')
    assert thin['duration_seconds'] == pytest.approx(48.*end_strain/(2./60))
    controls = [row for row in study['runs'] if row['candidate'] == 'duration_double']
    for control in controls:
        base = next(row for row in candidate_rows if row['role'] == control['role'])
        assert control['duration_seconds'] == pytest.approx(2*base['duration_seconds'])
        assert control['loading_protocol']['duration_multiplier'] == 2.
        assert control['loading_protocol']['actual_mean_crosshead_speed_mm_per_s'] == pytest.approx(1./60)
        assert not control['loading_protocol']['constant_physical_reference_speed_applied']
        assert control['numerical_sensitivity']


def test_accelerated_smoothstep_declares_mean_and_peak_speed(tmp_path):
    config = json.loads(BASE.read_text())
    config['loading_profile'] = 'smoothstep'
    source = tmp_path/'source.json'
    source.write_text(json.dumps(config))
    study = prepare(source, tmp_path/'study', backend='opencourant', speed_mm_min=2.,
                    acceleration_factor=18000., ams_time_step=1e-6, holdout_end_strain=.78)
    assert 'mean' in study['loading_protocol']['profile_assumption']
    assert len(study['runs']) == 22
    for row in study['runs']:
        protocol = row['loading_protocol']
        assert protocol['loading_profile'] == 'smoothstep'
        assert protocol['physical_reference_speed_mm_min'] == 2.
        assert protocol['acceleration_factor'] == 18000.
        assert protocol['actual_mean_crosshead_speed_mm_per_s'] == pytest.approx(600./protocol['duration_multiplier'])
        assert protocol['actual_peak_crosshead_speed_mm_per_s'] == pytest.approx(1.875*protocol['actual_mean_crosshead_speed_mm_per_s'])
        assert protocol['physical_reference_duration_seconds'] == pytest.approx(18000.*protocol['base_duration_seconds'])
        assert not protocol['constant_physical_reference_speed_applied']
        assert json.loads(Path(row['config']).read_text())['loading_profile'] == 'smoothstep'


def test_explicit_duration_preserves_legacy_protocol(tmp_path):
    study = prepare(BASE, tmp_path/'study', backend='opencourant', duration=.01)
    assert 'loading_protocol' not in study
    for row in study['runs']:
        assert row['duration_seconds'] == (.02 if row['candidate'] == 'duration_double' else .01)
        assert 'loading_protocol' not in row


def test_physical_speed_rejects_explicit_nonconstant_profile_before_writing(tmp_path):
    config = json.loads(BASE.read_text())
    config['loading_profile'] = 'smoothstep'
    source = tmp_path/'source.json'
    source.write_text(json.dumps(config))
    output = tmp_path/'study'
    with pytest.raises(ValueError, match='Physical speed mode requires loading_profile=linear'):
        prepare(source, output, backend='opencourant', speed_mm_min=2.)
    assert not output.exists()


@pytest.mark.parametrize('kwargs,message', [
    ({'speed_mm_min': 2., 'duration': .01}, 'either duration or speed'),
    ({'speed_mm_min': 2., 'backend': 'ccx'}, 'Speed mode requires OpenCourant'),
    ({'speed_mm_min': 0.}, 'positive speed_mm_min'),
    ({'speed_mm_min': -1.}, 'positive speed_mm_min'),
    ({'speed_mm_min': float('nan')}, 'positive speed_mm_min'),
    ({'speed_mm_min': float('inf')}, 'positive speed_mm_min'),
    ({'speed_mm_min': 2., 'acceleration_factor': 0.}, 'Acceleration factor'),
    ({'speed_mm_min': 2., 'acceleration_factor': -1.}, 'Acceleration factor'),
    ({'speed_mm_min': 2., 'acceleration_factor': float('nan')}, 'Acceleration factor'),
    ({'speed_mm_min': 2., 'acceleration_factor': float('inf')}, 'Acceleration factor'),
    ({'duration': .01, 'acceleration_factor': 2.}, 'requires speed_mm_min'),
    ({'speed_mm_min': 2., 'acceleration_factor': 100000., 'ams_time_step': 1e-4}, 'at most duration/100'),
])
def test_invalid_speed_protocol_does_not_write_a_partial_study(tmp_path, kwargs, message):
    options = {'backend': 'opencourant', **kwargs}
    output = tmp_path/'study'
    with pytest.raises(ValueError, match=message):
        prepare(BASE, output, **options)
    assert not output.exists()


def test_cli_duration_and_speed_are_mutually_exclusive(tmp_path):
    result = subprocess.run([sys.executable, str(ROOT/'calibration_study.py'), 'prepare',
                             '--backend', 'opencourant', '--output', str(tmp_path/'study'),
                             '--duration-seconds', '.01', '--speed-mm-min', '2'],
                            cwd=ROOT, text=True, capture_output=True)
    assert result.returncode == 2
    assert 'not allowed with argument' in result.stderr
    assert not (tmp_path/'study').exists()
