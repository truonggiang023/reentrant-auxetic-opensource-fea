"""Compare completed, receipt-bound native RP controls with the correct wall panel.

Only observed relative RP motion and signed interval reactions are used. A
completed 0.12 control is not a complete collapse curve or a material fit.
Run from the repository root with PYTHONPATH=.:src.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

import numpy as np

from assess_courant import assess
from calibration_study import reference_curves
from reentrant_fea.calibration import curve, curve_objective
from reentrant_fea.radioss_results import read_engine_out, energy_screen


ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def panel_for_wall(wall):
    if wall is not None and np.isclose(wall, 1., rtol=0, atol=1e-9):
        return 'thick'
    if wall is not None and np.isclose(wall, .3, rtol=0, atol=1e-9):
        return 'thin'
    raise ValueError('Require verified metadata wall_mm=1 or 0.3')


def slope(values, interval=(.005, .02)):
    lo, hi = interval
    result = dict(interval_strain=list(interval), available=False)
    if values[0, 0] <= lo and values[-1, 0] >= hi:
        x = np.linspace(lo, hi, 101)
        coefficient = np.polyfit(x, np.interp(x, values[:, 0], values[:, 1]), 1)
        result.update(available=True, slope_MPa=float(coefficient[0]),
                      intercept_MPa=float(coefficient[1]))
    else:
        result['reason'] = 'Observed interval not covered; no extrapolation'
    return result


def peak(values, end=.08):
    selected = values[(values[:, 0] >= 0) & (values[:, 0] <= end)]
    result = dict(interval_strain=[0., end], available=bool(len(selected)),
                  whole_interval_covered=bool(values[0, 0] <= .005 and values[-1, 0] >= end),
                  definition='Largest observed stress sample through 0.08; not necessarily a local fracture peak')
    if len(selected):
        i = int(np.argmax(selected[:, 1]))
        result.update(strain=float(selected[i, 0]), stress_MPa=float(selected[i, 1]),
                      at_sampled_window_boundary=bool(i in (0, len(selected)-1)))
    return result


def reference_metrics(experiments):
    result = []
    for values in experiments:
        values = curve(values)
        point = dict(strain=.12, available=bool(values[0, 0] <= .12 <= values[-1, 0]))
        if point['available']:
            point['stress_MPa'] = float(np.interp(.12, values[:, 0], values[:, 1]))
        result.append(dict(observed_strain_range=values[[0, -1], 0].tolist(),
                           stress_at_012=point, early_slope=slope(values),
                           peak_through_008=peak(values)))
    return result


def discrepancy(values, experiments, interval):
    result = curve_objective(values, experiments, interval)
    x = np.linspace(*result['interval_scored'], result['samples_per_test'])
    reference = np.asarray([np.interp(x, np.asarray(c)[:, 0], np.asarray(c)[:, 1])
                            for c in experiments])
    reference_rms = float(np.sqrt(np.mean(reference**2)))
    result.update(reference_rms_MPa=reference_rms,
                  relative_RMSE_to_reference_RMS=(result['overlap_rmse_MPa']/reference_rms
                                                  if reference_rms > 0 else None),
                  normalization='Investigator diagnostic: pooled residual RMS divided by pooled experimental stress RMS on the same equal-strain overlap grid; not a paper acceptance criterion')
    return result


def numerical_warnings(text):
    patterns = {
        'zero_or_negative_volume': r'\b(?:ZERO\s+OR\s+NEGATIVE|NEGATIVE)\s+VOLUME\b',
        'nonfinite_value_reported': r'\b(?:nan|inf(?:inity)?|non[- ]?finite)\b',
    }
    flags = {name: bool(re.search(pattern, text, re.I))
             for name, pattern in patterns.items()}
    matching = [line.strip() for line in text.splitlines()
                if any(re.search(pattern, line, re.I) for pattern in patterns.values())]
    flags.update(review_required=any(flags.values()),
                 matching_engine_line_count=len(matching),
                 matching_engine_lines=matching[:20],
                 qualification='Warnings remain visible independently of completion and dense kinetic-energy ratios; inspect numerical validity before accepting a fit')
    return flags


def compare(cases, reference_dir=ROOT/'validation/paper_reference'):
    """Return a report and unmodified observed curves for LABEL/PREFIX pairs."""
    reference_dir = Path(reference_dir)
    if not cases or len({label for label, _ in cases}) != len(cases):
        raise ValueError('Require cases with unique nonempty labels')
    references = {panel: reference_curves(reference_dir/f'dong2019_fig12_{panel}.csv')
                  for panel in ('thick', 'thin')}
    report = dict(status='Completed bound RP control comparison; calibration not established',
                  calibrated=False, quasi_static_validated=False,
                  qualification='Use curve and energy diagnostics together; passing solver/energy screens is not full quasi-static or fracture qualification.',
                  comparison_script_sha256=sha(__file__),
                  strain_measure='Observed relative top/bottom RP DY divided by actual core height',
                  stress_measure='Signed interval reaction from one time derivative of REACY, divided by declared nominal area',
                  paper_references={}, cases=[])
    for panel, experiments in references.items():
        path = reference_dir/f'dong2019_fig12_{panel}.csv'
        report['paper_references'][panel] = dict(path=str(path), sha256=sha(path),
            replicates=reference_metrics(experiments),
            qualification='Digitized experimental traces, not raw tests; paper FE trace is not an experimental fitting target')
    provenance = reference_dir/'dong2019_fig12_provenance.json'
    if provenance.is_file():
        report['paper_reference_provenance_sha256'] = sha(provenance)
    curves = {}
    for label, prefix in cases:
        if not label:
            raise ValueError('Require a nonempty case label')
        prefix = Path(prefix)
        metadata_path = Path(str(prefix)+'.radioss.json')
        history_path = Path(str(prefix)+'T01.csv')
        receipt_path = Path(str(prefix)+'.run.json')
        metadata = json.loads(metadata_path.read_text())
        panel = panel_for_wall(metadata.get('wall_mm'))
        actual = assess(metadata_path, history_path, force_source='rp')
        if not (actual['completed'] and actual.get('run_binding_verified')
                and actual.get('metadata_binding_verified')):
            raise ValueError(f'{label}: require completed solver and bound metadata/outputs')
        # assess() has just verified this output against its execution receipt.
        engine_path = Path(str(prefix)+'_0001.out')
        printed_energy = energy_screen(read_engine_out(engine_path), metadata['duration_seconds'])
        warnings = numerical_warnings(engine_path.read_text(errors='replace'))
        values = curve([[row['strain'], row['nominal_stress_MPa']]
                        for row in actual['reaction_history']], label)
        curves[label] = values
        if 'fracture_locus' not in metadata:
            raise ValueError(f'{label}: metadata does not declare whether damage is enabled')
        damage = metadata['fracture_locus'] is not None
        model = 'J2 with native damage' if damage else 'J2 without damage'
        scoring_interval = (.04, .119) if panel == 'thick' else (.03, .119)
        try:
            difference = discrepancy(values, references[panel], scoring_interval)
        except ValueError as error:
            difference = dict(available=False, reason=str(error))
        last = actual['reaction_history'][-1]
        target = dict(strain=.12, available=bool(values[0, 0] <= .12 <= values[-1, 0]),
                      convention='No extrapolation from the last interval sample to nominal target')
        if target['available']:
            target['stress_MPa'] = float(np.interp(.12, values[:, 0], values[:, 1]))
        report['cases'].append(dict(label=label, prefix=str(prefix), panel=panel,
            wall_mm=float(metadata['wall_mm']), model=model, damage_enabled=damage,
            fracture_locus=metadata['fracture_locus'],
            duration_seconds=metadata['duration_seconds'],
            height_mm=metadata['height_mm'], reference_area_mm2=metadata['reference_area_mm2'],
            ams=metadata.get('ams'),
            observed_strain_coverage=values[[0, -1], 0].tolist(),
            last_observed_nominal_strain=actual.get('last_observed_nominal_strain'),
            target_nominal_strain=metadata.get('target_nominal_strain'),
            scope='Only the observed control strain range; no missing collapse response inferred',
            last_force_sample=last, stress_at_exact_012=target,
            early_slope=slope(values), peak_through_008=peak(values),
            paper_discrepancy=difference,
            global_energy_screen=actual['global_energy_screen'],
            core_energy_screen=actual['core_energy_screen'],
            printed_energy_screen=printed_energy, numerical_warning_flags=warnings,
            core_elements=metadata.get('core_elements'), core_nodes=metadata.get('core_nodes'),
            source_config_sha256=metadata.get('source_config_sha256'),
            hardening_source_sha256=metadata.get('hardening_source_sha256'),
            platen_contact_surface=metadata.get('platen_contact_surface'),
            negative_signed_stress_samples=int(np.count_nonzero(values[:, 1] < 0)),
            input_sha256=metadata['sha256'], metadata_sha256=sha(metadata_path),
            history_sha256=sha(history_path), execution_receipt_sha256=sha(receipt_path),
            calibrated=False))
    return report, curves, references


def write_outputs(report, curves, references, output):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
    with output.with_suffix('.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['case', 'wall_mm', 'observed_nominal_strain', 'nominal_stress_MPa'])
        for case in report['cases']:
            for strain, stress in curves[case['label']]:
                writer.writerow([case['label'], case['wall_mm'], strain, stress])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for kind in ('zoom', 'full_paper'):
        fig, axes = plt.subplots(1, 2, figsize=(13, 4.7))
        for ax, panel, wall, limit in zip(axes, ('thick', 'thin'), (1., .3), (.5, .8)):
            xmax = .12 if kind == 'zoom' else limit
            reference_colors = ('black', 'red') if panel == 'thick' else ('red', 'black')
            for i, reference in enumerate(references[panel], 1):
                values = np.asarray(reference)
                mask = values[:, 0] <= xmax
                ax.plot(values[mask, 0], values[mask, 1], '--', lw=1.5, color=reference_colors[i-1],
                        label=f'Paper experiment {i}')
            native_index = 0
            native_colors = ('tab:blue', 'tab:orange', 'tab:purple', 'tab:cyan')
            for case in report['cases']:
                if case['panel'] != panel:
                    continue
                values = curves[case['label']]
                mask = values[:, 0] <= xmax
                label = f"{case['label']} ({case['model']}; ε≤{values[-1, 0]:.3f})"
                ax.plot(values[mask, 0], values[mask, 1], label=label, lw=1.2,
                        color=native_colors[native_index % len(native_colors)])
                native_index += 1
            ax.set(xlim=(0, xmax), xlabel='Observed nominal strain',
                   ylabel='Nominal stress (MPa)', title=f'Wall t={wall:g} mm')
            ax.grid(alpha=.25)
            ax.legend(fontsize=7)
        fig.suptitle('Bound OpenCourant RP controls; current curves end at observed coverage — not full calibration')
        fig.tight_layout()
        fig.savefig(output.with_name(output.name+'_'+kind+'.png'), dpi=180)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', action='append', required=True, metavar='LABEL=PREFIX')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference-dir', type=Path, default=ROOT/'validation/paper_reference')
    args = parser.parse_args()
    cases = []
    for item in args.case:
        label, separator, prefix = item.partition('=')
        if not separator or not label or not prefix:
            parser.error('--case requires LABEL=PREFIX')
        cases.append((label, Path(prefix)))
    report, curves, references = compare(cases, args.reference_dir)
    write_outputs(report, curves, references, args.output)
    print(json.dumps({k: report[k] for k in ('status', 'calibrated')}, indent=2))
    print('Report: '+str(args.output.with_suffix('.json')))


if __name__ == '__main__':
    main()
