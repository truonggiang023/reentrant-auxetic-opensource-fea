"""Constrained candidate fitting with coverage and holdout checks.

The available paper does not identify a unique fracture locus. These helpers
keep the published point fixed, expose the missing coordinates as estimates,
and never promote a partial curve to a completed calibration.
"""
from copy import deepcopy

import numpy as np


ANCHOR = (0.333, 0.047)


def fracture_candidate(config, eta_low, eta_high):
    """Return an uncalibrated hypothesis; endpoint coordinates are not paper data."""
    if not np.isfinite([eta_low, eta_high]).all() or not eta_low < ANCHOR[0] < eta_high:
        raise ValueError('Endpoint coordinates must bracket the published anchor')
    result = deepcopy(config)
    law = result['ductile_damage']
    law.update(fracture_locus=[[float(eta_low), .08], list(ANCHOR), [float(eta_high), .018]],
               failure_displacement=.001, calibrated=False)
    law['provenance'] = (
        'Dong 2019 Section 4.1: .047 at eta=.333 and uf=.001 mm. '
        f'Estimated endpoint coordinates eta={eta_low:g}/{eta_high:g}; '
        'endpoint strains .08/.018 are the paper estimate range, not located measurements. '
        'Candidate for thick-wall fitting and independent thin-wall validation; not calibrated.')
    return result


def curve(values, name='curve'):
    """Validate a measured or simulated strain/stress series without sorting it."""
    result = np.asarray(values, dtype=float)
    if (result.ndim != 2 or result.shape[1] != 2 or len(result) < 2
            or not np.isfinite(result).all() or np.any(np.diff(result[:, 0]) <= 0)):
        raise ValueError(f'{name} requires finite, strictly increasing strain/stress pairs')
    return result


def curve_objective(simulation, experiments, interval=(.02, .65), samples=300):
    """RMS discrepancy to both tests on overlap only; report missing coverage.

    Equal strain spacing prevents dense output around cutbacks from dominating
    the objective. Never extrapolate a partial simulation into unrun strain.
    """
    simulation = curve(simulation, 'simulation')
    experiments = [curve(values, 'experiment') for values in experiments]
    if not experiments:
        raise ValueError('At least one experimental curve is required')
    lo, hi = map(float, interval)
    if not np.isfinite([lo, hi]).all() or lo >= hi or samples < 2:
        raise ValueError('Require a finite increasing interval and at least two samples')
    start = max(lo, simulation[0, 0], *(values[0, 0] for values in experiments))
    stop = min(hi, simulation[-1, 0], *(values[-1, 0] for values in experiments))
    if start >= stop:
        raise ValueError('Simulation and experimental curves have no usable overlap')
    strain = np.linspace(start, stop, samples)
    predicted = np.interp(strain, simulation[:, 0], simulation[:, 1])
    reference = np.array([np.interp(strain, values[:, 0], values[:, 1])
                          for values in experiments])
    residual = predicted - reference
    return {
        'interval_requested': [lo, hi], 'interval_scored': [float(start), float(stop)],
        'coverage_fraction': float((stop-start)/(hi-lo)),
        'full_interval_covered': bool(start <= lo and stop >= hi),
        'overlap_rmse_MPa': float(np.sqrt(np.mean(residual**2))),
        'overlap_rmse_by_test_MPa': np.sqrt(np.mean(residual**2, axis=1)).tolist(),
        'overlap_bias_MPa': float(np.mean(residual)),
        'experimental_half_range_rms_MPa': float(np.sqrt(np.mean((np.ptp(reference, axis=0)/2)**2))),
        'samples_per_test': samples,
        'objective_scope': 'shared strain interval only; incomplete curves cannot win calibration',
    }


def score_assessment(assessment, experiments, interval=(.02, .65)):
    history = assessment['reaction_history']
    values = [[row['strain'], row['nominal_stress_MPa']] for row in history]
    result = curve_objective(values, experiments, interval)
    # Explicit native histories have observed strain, not converged Newton
    # increments. Callers must verify their underlying solver files first.
    last_strain = assessment.get('last_converged_nominal_strain',
                                 assessment.get('last_observed_nominal_strain', 0))
    complete = (assessment.get('completed') is True and not assessment.get('errors')
                and last_strain is not None and last_strain >= interval[1]-1e-8)
    result['solver_completed'] = bool(complete)
    result['eligible_for_selection'] = bool(complete and result['full_interval_covered'])
    result['input_sha256'] = assessment.get('input_sha256')
    result['completion_basis'] = ('observed explicit-solver history' if 'last_observed_nominal_strain' in assessment
                                  else 'accepted implicit-solver increments')
    return result


def select_candidate(scored_candidates, training_wall_mm=1.0):
    """Select a candidate only from completed thick-wall fits, never thin holdout."""
    eligible = []
    for row in scored_candidates:
        if row.get('role') != 'training' or not np.isclose(row.get('wall_mm', np.nan), training_wall_mm):
            raise ValueError('Candidate selection requires thick-wall training data, not the thin-wall holdout')
        if row['score'].get('eligible_for_selection'):
            eligible.append(row)
    if not eligible:
        return {'selected_candidate': None, 'status': 'blocked: no completed training curves cover the fit interval',
                'calibrated': False}
    best = min(eligible, key=lambda row: row['score']['overlap_rmse_MPa'])
    return {'selected_candidate': best['candidate'],
            'training_rmse_MPa': best['score']['overlap_rmse_MPa'],
            'status': 'best evaluated candidate; requires sensitivity and independent thin-wall validation',
            'calibrated': False}


def sensitivity_rank(jacobian, parameter_scales, tolerance=1e-3):
    """Diagnose local identifiability from actual response sensitivities.

    Callers supply d(stress residual)/d(parameter), not just sampled parameters.
    Scaling uses the investigated parameter spans, making units comparable.
    """
    jacobian = np.asarray(jacobian, dtype=float)
    scales = np.asarray(parameter_scales, dtype=float)
    if (jacobian.ndim != 2 or scales.shape != (jacobian.shape[1],)
            or not np.isfinite(jacobian).all() or not np.isfinite(scales).all()
            or np.any(scales <= 0) or jacobian.shape[1] == 0 or not 0 < tolerance < 1):
        raise ValueError('Require a finite sensitivity matrix, positive parameter scales and relative tolerance')
    singular = np.linalg.svd(jacobian*scales, compute_uv=False)
    largest = singular[0] if len(singular) else 0
    rank = int(np.sum(singular > tolerance*largest)) if largest else 0
    return {'rank': rank, 'parameter_count': len(scales),
            'locally_identifiable': rank == len(scales),
            'scaled_singular_values': singular.tolist(), 'relative_rank_tolerance': tolerance}
