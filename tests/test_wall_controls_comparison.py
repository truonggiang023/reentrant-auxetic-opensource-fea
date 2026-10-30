"""Control reports must retain observed signed response and the correct panel."""
import csv
import json
from pathlib import Path

import numpy as np
import pytest

from scripts import compare_wall_controls as mod


REFERENCE = Path(__file__).resolve().parents[1]/'validation/paper_reference'


def case(tmp_path, monkeypatch, *, wall=1., completed=True, bound=True):
    prefix = tmp_path/'control'
    metadata = dict(wall_mm=wall, fracture_locus=None, duration_seconds=.04,
                    height_mm=48.7, reference_area_mm2=1900., sha256={},
                    target_nominal_strain=.12,
                    displacement_fraction_samples=[[0., 0.], [.04, 1.]])
    Path(str(prefix)+'.radioss.json').write_text(json.dumps(metadata))
    Path(str(prefix)+'T01.csv').write_text('actual test fixture\n')
    Path(str(prefix)+'.run.json').write_text('{}')
    actual = dict(completed=completed, run_binding_verified=bound,
                  metadata_binding_verified=bound,
                  global_energy_screen={'sampled_energy_screen_passed': False},
                  core_energy_screen={'sampled_energy_screen_passed': True},
                  last_observed_nominal_strain=.1199998,
                  reaction_history=[dict(strain=x, nominal_stress_MPa=y, time=t)
                                    for x, y, t in [(0., -.01, 0.), (.005, 2., .0005),
                                                    (.01, 4., .001),
                                                    (.02, 8., .004), (.06, 12., .009),
                                                    (.119999, 7., .0399)]])
    def assess(metadata_path, history_path, **kwargs):
        assert kwargs == {'force_source': 'rp'}
        return actual
    monkeypatch.setattr(mod, 'assess', assess)
    return prefix


def test_bound_observed_curve_not_amplitude_or_exact_nominal_endpoint(tmp_path, monkeypatch):
    prefix = case(tmp_path, monkeypatch)
    report, curves, references = mod.compare([('thick', prefix)], REFERENCE)
    result = report['cases'][0]
    assert result['panel'] == 'thick' and result['model'] == 'J2 without damage'
    assert result['last_force_sample']['strain'] == .119999
    assert not result['stress_at_exact_012']['available']
    assert result['observed_strain_coverage'] == [0., .119999]
    assert curves['thick'][0, 1] == -.01
    assert result['early_slope']['slope_MPa'] == pytest.approx(400.)
    assert result['peak_through_008']['strain'] == .06
    assert result['paper_discrepancy']['interval_scored'] == pytest.approx([.04, .119])
    assert len(result['paper_discrepancy']['overlap_rmse_by_test_MPa']) == 2
    assert not report['calibrated'] and not report['quasi_static_validated']
    mod.write_outputs(report, curves, references, tmp_path/'comparison')
    with (tmp_path/'comparison.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    assert float(rows[0]['nominal_stress_MPa']) == -.01
    assert (tmp_path/'comparison_zoom.png').is_file()
    assert (tmp_path/'comparison_full_paper.png').is_file()


def test_thin_case_selects_thin_experiments(tmp_path, monkeypatch):
    prefix = case(tmp_path, monkeypatch, wall=.3)
    report, _, _ = mod.compare([('thin', prefix)], REFERENCE)
    result = report['cases'][0]
    assert result['panel'] == 'thin'
    assert result['paper_discrepancy']['interval_scored'] == pytest.approx([.03, .119])
    target = report['paper_references']['thin']['replicates'][0]['stress_at_012']
    assert target['stress_MPa'] == pytest.approx(.28122807)


@pytest.mark.parametrize('completed,bound', [(False, True), (True, False)])
def test_unfinished_or_unbound_run_cannot_make_a_completed_report(tmp_path, monkeypatch,
                                                               completed, bound):
    prefix = case(tmp_path, monkeypatch, completed=completed, bound=bound)
    with pytest.raises(ValueError, match='completed solver and bound metadata'):
        mod.compare([('control', prefix)], REFERENCE)


def test_missing_or_unrecognized_wall_is_not_inferred_from_label(tmp_path, monkeypatch):
    prefix = case(tmp_path, monkeypatch, wall=.5)
    with pytest.raises(ValueError, match='wall_mm'):
        mod.compare([('thick', prefix)], REFERENCE)


def test_duplicate_labels_are_rejected():
    with pytest.raises(ValueError, match='unique'):
        mod.compare([('a', Path('one')), ('a', Path('two'))], REFERENCE)


def test_reference_rms_normalization_is_invariant_to_stress_unit_scale():
    observed = np.array([[0., 0.], [.1, 4.], [.2, 8.]])
    experiments = [np.array([[0., 0.], [.1, 2.], [.2, 4.]]),
                   np.array([[0., 0.], [.1, 3.], [.2, 6.]])]
    first = mod.discrepancy(observed, experiments, (.02, .19))
    scaled = [x*np.array([1., .001]) for x in experiments]
    second = mod.discrepancy(observed*np.array([1., .001]), scaled, (.02, .19))
    assert first['overlap_rmse_MPa'] == pytest.approx(1000*second['overlap_rmse_MPa'])
    assert first['relative_RMSE_to_reference_RMS'] == pytest.approx(second['relative_RMSE_to_reference_RMS'])
    assert first['relative_RMSE_to_reference_RMS'] > 0
