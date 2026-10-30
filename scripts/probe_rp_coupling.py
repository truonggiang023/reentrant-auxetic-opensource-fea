"""Run the native direct-face/RBE2 coupon check through strain 0.12.

Source the OpenCourant env.sh and run from the repository root with
PYTHONPATH=.:src. This numerical coupling check is not a calibration.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from compare_rp_reactions import extract
from reentrant_fea.radioss import write_radioss_deck
from run_radioss import run_generated


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path('results/native_rp'))
    args = parser.parse_args()
    # Reject accidental reuse before writing decks or replacing old evidence.
    args.output_dir.mkdir(parents=True, exist_ok=False)
    nodes = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                      [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
    cases = {}
    curves = {}
    for name, reference_point in [('direct', False), ('rp', True)]:
        config = dict(young=69000., poisson=.3, compression_strain=.12,
                      reference_area=1., platen_thickness=.2, platen_margin=.1,
                      friction=.35, platen_reference_points=reference_point,
                      support_reaction_history=True, platen_contact_surface='inner_y')
        prefix = args.output_dir / name
        metadata = write_radioss_deck(prefix, nodes, np.array([list(range(8))]),
                                     config, [(225., 0.), (300., .2)], duration=.001)
        print(run_generated(metadata, timeout=120), flush=True)
        curves[name], cases[name] = extract(prefix)
        receipt = json.loads(prefix.with_suffix('.run.json').read_text())
        cases[name]['solver_stages_terminated_normally'] = receipt['solver_stages_terminated_normally']

    direct, rp = curves['direct'], curves['rp']
    if direct.shape != rp.shape or not np.allclose(direct[:, 0], rp[:, 0]):
        raise ValueError('Comparison requires matching time grids')
    if not np.isfinite(direct).all() or not np.isfinite(rp).all():
        raise ValueError('Non-finite coupon history')
    force_scale = float(np.linalg.norm(direct[:, 2]))
    if force_scale <= 0:
        raise ValueError('Direct coupon did not develop a reaction force')
    relative_rms = float(np.linalg.norm(direct[:, 2] - rp[:, 2]) / force_scale)
    checks = {
        'both_solvers_completed': all(c['solver_stages_terminated_normally'] for c in cases.values()),
        'both_reached_strain012': all(np.isclose(c['final_sample_strain'], .12, rtol=.005)
                                     for c in cases.values()),
        'relative_force_rms_below_0p1_percent': relative_rms < .001,
    }
    report = dict(cases=cases, checks=checks, passed=all(checks.values()),
                  relative_force_rms=relative_rms,
                  max_force_difference_N=float(np.max(abs(direct[:, 2] - rp[:, 2]))),
                  calibrated=False,
                  scope='Native coupon coupling check; not quasi-static lattice or experimental calibration.')
    (args.output_dir / 'comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    for name, curve in curves.items():
        ax.plot(curve[:, 1], curve[:, 4], label=name)
    ax.set(xlabel='Nominal strain', ylabel='Nominal stress (MPa)',
           title='OpenCourant coupon: direct face and RBE2 reference point')
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output_dir / 'comparison.png', dpi=160)
    plt.close(fig)
    print(json.dumps(report, indent=2), flush=True)
    if not report['passed']:
        raise SystemExit('Native coupling regression failed; inspect comparison.json')


if __name__ == '__main__':
    main()
