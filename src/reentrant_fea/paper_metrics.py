"""Curve diagnostics for paper comparison, without fitting or extrapolation.

Thin-wall intervals below are investigator choices for Fig. 12(b), not limits
or a densification definition reported by Dong et al. No thin-wall interval is
silently reused for the thick-wall panel. Printed curves are observations from
which missing initial strain or clipped high-strain response cannot be recovered.
"""
import numpy as np


def _curve(values):
    result = np.asarray(values, dtype=float)
    if (result.ndim != 2 or result.shape[1] != 2 or len(result) < 2
            or not np.isfinite(result).all()
            or np.any(np.diff(result[:, 0]) <= 0)
            or np.any(result[:, 0] < 0) or np.any(result[:, 1] < 0)):
        raise ValueError('Require finite nonnegative strain/stress and strictly increasing strain')
    return result


def _integral(x, y):
    return float(np.sum(np.diff(x)*(y[1:]+y[:-1])/2))


def _window(curve, interval):
    lo, hi = map(float, interval)
    if not np.isfinite([lo, hi]).all() or lo >= hi:
        raise ValueError('Require finite increasing interval limits')
    if curve[0, 0] > lo or curve[-1, 0] < hi:
        return None
    inside = (curve[:, 0] > lo) & (curve[:, 0] < hi)
    x = np.r_[lo, curve[inside, 0], hi]
    return x, np.interp(x, curve[:, 0], curve[:, 1])


def _interval_metric(curve, interval, kind):
    if interval is None:
        return {'available': False, 'reason': 'No panel-specific interval supplied'}
    window = _window(curve, interval)
    result = {'interval_strain': list(map(float, interval)), 'available': window is not None}
    if window is None:
        result['reason'] = 'Observed curve does not cover the whole interval; no extrapolation'
    elif kind == 'mean':
        x, y = window
        result['mean_stress_MPa'] = _integral(x, y)/(x[-1]-x[0])
    else:
        # Equal strain weighting avoids bias from nonuniform digitized/output sampling.
        x = np.linspace(*interval, 101)
        y = np.interp(x, curve[:, 0], curve[:, 1])
        slope, intercept = np.polyfit(x, y, 1)
        result.update(slope_MPa=float(slope), intercept_MPa=float(intercept))
    return result


def curve_metrics(values, *, panel='thin', plateau_interval=None,
                  slope_interval=(.005, .02), densification_interval=None,
                  target_strains=(.12, .35, .65, .75)):
    """Return diagnostics for one nonnegative compression curve.

    Energy absorption is stress integrated over the *observed* strain range,
    with units MPa (= mJ/mm^3). Missing initial samples are never padded by zero.
    Densification is an investigator diagnostic: maximize W(epsilon)/sigma(epsilon)
    in the stated window, where W starts at the first observed sample. It is
    not a fitted or uniquely identified physical densification strain. A peak
    at a search boundary is explicitly flagged.

    Thin defaults: plateau .1-.55; densification search/required coverage .1-.77.
    Thick has no defaults for these two diagnostics; supply intervals explicitly.
    """
    if panel not in ('thin', 'thick'):
        raise ValueError('panel must be thin or thick')
    c = _curve(values)
    if panel == 'thin':
        plateau_interval = (.1, .55) if plateau_interval is None else plateau_interval
        densification_interval = (.1, .77) if densification_interval is None else densification_interval
    targets = np.asarray(target_strains, dtype=float)
    if targets.ndim != 1 or not np.isfinite(targets).all() or np.any(targets < 0):
        raise ValueError('Target strains must be finite and nonnegative')
    points = []
    for target in targets:
        covered = bool(c[0, 0] <= target <= c[-1, 0])
        row = {'strain': float(target), 'available': covered}
        if covered:
            row['stress_MPa'] = float(np.interp(target, c[:, 0], c[:, 1]))
        else:
            row['reason'] = 'Target is outside observed strain range; no extrapolation'
        points.append(row)
    energy = {'available': True, 'observed_interval_strain': c[[0, -1], 0].tolist(),
              'energy_absorption_MPa': _integral(c[:, 0], c[:, 1]),
              'starts_at_zero_strain': bool(c[0, 0] == 0.),
              'convention': 'Integral of stress over observed strain only; no missing initial energy inferred'}
    densification = {'available': False,
        'definition': 'Investigator diagnostic: argmax of observed-range integral(stress)/stress',
        'paper_defined': False, 'energy_origin_strain': float(c[0, 0])}
    if densification_interval is None:
        densification['reason'] = 'No thick-panel densification interval supplied'
    else:
        window = _window(c, densification_interval)
        densification['required_interval_strain'] = list(map(float, densification_interval))
        if window is None:
            densification['reason'] = 'Insufficient strain coverage for densification; no extrapolation'
        else:
            # Include native breakpoints and a dense common strain grid. Every
            # integral is exact for the piecewise linear observed curve.
            lo, hi = densification_interval
            x = np.unique(np.r_[c[c[:, 0] <= hi, 0], np.linspace(lo, hi, 2001)])
            y = np.interp(x, c[:, 0], c[:, 1])
            work = np.r_[0., np.cumsum(np.diff(x)*(y[1:]+y[:-1])/2)]
            eligible = (x >= lo) & (x <= hi) & (y > 0)
            if not np.any(eligible):
                densification['reason'] = 'No positive stress in search interval; efficiency undefined'
            else:
                indices = np.flatnonzero(eligible)
                index = indices[np.argmax(work[eligible]/y[eligible])]
                densification.update(available=True, strain=float(x[index]),
                    efficiency=float(work[index]/y[index]),
                    boundary_peak=bool(np.isclose(x[index], lo) or np.isclose(x[index], hi)),
                    zero_stress_samples_ignored=int(np.count_nonzero((x >= lo) & (x <= hi) & (y == 0))))
    return {'panel': panel, 'coverage_strain': c[[0, -1], 0].tolist(),
            'early_slope': _interval_metric(c, slope_interval, 'slope'),
            'plateau': _interval_metric(c, plateau_interval, 'mean'),
            'target_stresses': points, 'energy_absorption': energy,
            'densification': densification,
            'qualification': 'Curve diagnostics only; no calibration or experimental agreement established'}


def replicate_metrics(experiments, *, panel='thin', **kwargs):
    """Keep each caller-supplied replicate separate, including missing coverage."""
    curves = list(experiments)
    if not curves:
        raise ValueError('At least one experimental replicate required')
    return {'panel': panel, 'replicates': [curve_metrics(c, panel=panel, **kwargs) for c in curves],
            'qualification': 'Replicate diagnostics; printed-curve uncertainty must be reported separately'}
