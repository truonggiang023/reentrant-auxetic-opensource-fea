import numpy as np
import pytest
from reentrant_fea.paper_metrics import curve_metrics, replicate_metrics


def test_linear_curve_exact_slope_mean_and_energy():
    result = curve_metrics([[0, 0], [.8, 8]])
    assert result['early_slope']['slope_MPa'] == pytest.approx(10)
    assert result['plateau']['mean_stress_MPa'] == pytest.approx(3.25)
    assert result['energy_absorption']['energy_absorption_MPa'] == pytest.approx(3.2)
    assert [p['stress_MPa'] for p in result['target_stresses']] == pytest.approx([1.2, 3.5, 6.5, 7.5])
    assert result['densification']['boundary_peak']


def test_piecewise_plateau_densification_peak():
    # Flat stress to .6; steep densification thereafter makes W/stress decline.
    result = curve_metrics([[0, 1], [.6, 1], [.8, 21]])
    assert result['densification']['strain'] == pytest.approx(.6)
    assert result['densification']['efficiency'] == pytest.approx(.6)
    assert not result['densification']['boundary_peak']
    assert result['plateau']['mean_stress_MPa'] == pytest.approx(1)
    assert result['energy_absorption']['energy_absorption_MPa'] == pytest.approx(2.8)
    assert result['densification']['paper_defined'] is False


def test_truncated_run_cannot_infer_densification_or_unobserved_targets():
    result = curve_metrics([[0, 0], [.12, 1]])
    assert not result['densification']['available']
    assert not result['plateau']['available']
    assert result['target_stresses'][0]['available']
    assert all(not row['available'] for row in result['target_stresses'][1:])


def test_missing_initial_data_does_not_infer_slope_or_zero_origin_energy():
    result = curve_metrics([[.03, 2], [.8, 2]])
    assert not result['early_slope']['available']
    energy = result['energy_absorption']
    assert energy['energy_absorption_MPa'] == pytest.approx(1.54)
    assert not energy['starts_at_zero_strain']
    assert result['densification']['energy_origin_strain'] == .03


def test_zero_stress_efficiency_is_unavailable_and_negative_stress_rejected():
    result = curve_metrics([[0, 0], [.8, 0]])
    assert not result['densification']['available']
    assert result['energy_absorption']['energy_absorption_MPa'] == 0
    with pytest.raises(ValueError, match='nonnegative'):
        curve_metrics([[0, 0], [.8, -1]])


def test_thick_panel_does_not_inherit_thin_intervals():
    result = curve_metrics([[0, 0], [.5, 10]], panel='thick')
    assert not result['plateau']['available']
    assert not result['densification']['available']
    assert not result['target_stresses'][2]['available']
    custom = curve_metrics([[0, 0], [.5, 10]], panel='thick', plateau_interval=(.1, .4))
    assert custom['plateau']['mean_stress_MPa'] == pytest.approx(5)


def test_replicates_remain_separate_and_coverages_are_honest():
    result = replicate_metrics([[[0, 0], [.8, 1]], [[.03, .1], [.7, .8]]])
    assert len(result['replicates']) == 2
    assert result['replicates'][0]['densification']['available']
    assert not result['replicates'][1]['densification']['available']
    with pytest.raises(ValueError):
        replicate_metrics([])


@pytest.mark.parametrize('values', [[[0, 1], [0, 2]], [[.1, 1], [0, 2]],
                                   [[0, 1], [1, np.nan]], [[-.1, 0], [1, 2]]])
def test_invalid_curves_are_not_repaired(values):
    with pytest.raises(ValueError):
        curve_metrics(values)
