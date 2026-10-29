"""Independent reaction/impulse bookkeeping for platen force comparisons.

All impulse histories are cumulative in N*s and sampled on the same time grid.
A platen balance is R + C + P = d(momentum)/dt, where R denotes the sum of
constraint reactions, C the force of contact acting on that platen, and P any
other applied force. For a stationary/constant-velocity platen with no other
loads, R = -C. This module does not infer force signs from column names.
"""
import numpy as np


def interval_force(time, cumulative_impulse):
    """Return interval-averaged forces and interval midpoint times.

    Finite differences are taken ONCE. Do not differentiate an instantaneous
    force channel, and do not interpret cumulative impulses as forces.
    """
    t = np.asarray(time, dtype=float)
    impulse = np.asarray(cumulative_impulse, dtype=float)
    if t.ndim != 1 or len(t) < 2 or not np.all(np.isfinite(t)):
        raise ValueError("time must be finite 1D with at least two samples")
    if np.any(np.diff(t) <= 0):
        raise ValueError("time must be strictly increasing")
    if impulse.shape[0] != len(t) or not np.all(np.isfinite(impulse)):
        raise ValueError("impulse must have a finite value at every time")
    dt = np.diff(t)
    shape = (len(dt),) + (1,) * (impulse.ndim - 1)
    return (t[:-1] + t[1:]) / 2, np.diff(impulse, axis=0) / dt.reshape(shape)


def compare_platen_impulses(time, reaction_impulse, contact_impulse,
                            momentum=None, other_applied_impulse=None,
                            area_mm2=None):
    """Compare forces ON THE SAME PLATEN with a signed momentum balance.

    reaction_impulse: sum of prescribed-DOF reaction impulses on the platen.
    contact_impulse: cumulative impulse of contact acting on the platen;
      NOT the opposing impulse acting on the core.
    momentum: platen signed linear momentum (Mg*mm/s expressed in N*s).
    other_applied_impulse: impulses of loads other than reaction/contact.
    Use no momentum only when it can legitimately be neglected; zero momentum
    is assumed in that case and should be explicitly qualified by the caller.
    """
    t = np.asarray(time, float)
    reaction = np.asarray(reaction_impulse, float)
    contact = np.asarray(contact_impulse, float)
    if reaction.ndim != 1 or contact.shape != reaction.shape:
        raise ValueError("reaction and contact impulse must be matching 1D arrays")
    if momentum is None:
        momentum = np.zeros_like(reaction)
    if other_applied_impulse is None:
        other_applied_impulse = np.zeros_like(reaction)
    mid, r = interval_force(t, reaction)
    _, c = interval_force(t, contact)
    _, p = interval_force(t, other_applied_impulse)
    _, inertia = interval_force(t, momentum)
    residual = r + c + p - inertia
    result = {
        "time_mid_s": mid,
        "reaction_force_N": r,
        "platen_contact_force_N": c,
        "other_force_N": p,
        "inertial_force_N": inertia,
        "balance_residual_N": residual,
        "balance_relative_rms": float(np.sqrt(np.mean(residual**2)) /
                                      max(np.sqrt(np.mean((r + c + p)**2)),
                                          np.sqrt(np.mean(inertia**2)), 1e-12)),
    }
    if area_mm2 is not None:
        area = float(area_mm2)
        if not np.isfinite(area) or area <= 0:
            raise ValueError("area_mm2 must be positive and finite")
        result["reaction_nominal_stress_MPa"] = np.abs(r) / area
        result["contact_nominal_stress_MPa"] = np.abs(c) / area
    return result
