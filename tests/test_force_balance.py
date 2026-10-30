import numpy as np
import pytest

from reentrant_fea.force_balance import interval_force, compare_platen_impulses


def test_constant_force_requires_single_impulse_derivative():
    t = np.array([0., .1, .3, .6])
    middle, force = interval_force(t, 12. * t)
    np.testing.assert_allclose(force, 12.)
    np.testing.assert_allclose(middle, [.05, .2, .45])


def test_static_platen_balance_and_nominal_stress():
    t = np.array([0., .1, .2])
    result = compare_platen_impulses(t, 532. * t, -532. * t,
                                     area_mm2=1900)
    np.testing.assert_allclose(result["reaction_nominal_stress_MPa"],
                               [0.28, 0.28])
    np.testing.assert_allclose(result["balance_residual_N"], 0., atol=1e-12)
    assert result["balance_relative_rms"] < 1e-12


def test_accelerating_platen_reaction_differs_from_contact():
    t = np.array([0., .1, .2])
    # R=20 N, C=-12 N, d(momentum)/dt=8 N.
    result = compare_platen_impulses(t, 20*t, -12*t, momentum=8*t)
    np.testing.assert_allclose(result["balance_residual_N"], 0., atol=1e-12)
    np.testing.assert_allclose(result["reaction_force_N"], 20.)
    np.testing.assert_allclose(result["platen_contact_force_N"], -12.)


def test_detects_wrong_contact_sign():
    t = np.array([0., .1, .2])
    result = compare_platen_impulses(t, 20*t, 12*t)
    np.testing.assert_allclose(result["balance_residual_N"], 32.)
    assert result["balance_relative_rms"] > 0.9


@pytest.mark.parametrize("time", [[0.], [0., 0.], [0., -1.], [0., float("nan")]])
def test_rejects_invalid_time(time):
    with pytest.raises(ValueError):
        interval_force(time, np.zeros(len(time)))


def test_rejects_mismatched_and_nonfinite_impulse():
    with pytest.raises(ValueError):
        interval_force([0., 1., 2.], [0., 1.])
    with pytest.raises(ValueError):
        interval_force([0., 1.], [0., float("nan")])
    with pytest.raises(ValueError):
        compare_platen_impulses([0., 1.], [0., 1.], [0.])
    with pytest.raises(ValueError):
        compare_platen_impulses([0., 1.], [0., 1.], [0., -1.], area_mm2=0)
