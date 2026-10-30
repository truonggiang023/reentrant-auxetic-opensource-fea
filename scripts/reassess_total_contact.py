#!/usr/bin/env python3
"""Reassess retained native thin-wall histories using the verified FN-FT force.

This only writes derived evidence. Bound decks, metadata, T01 and history CSVs
are read and verified by the assessor; no native run or material is changed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from assess_courant import assess
from plot_paper_comparison import read_reference
from reentrant_fea.calibration import curve_objective


NORMAL = 'Contact force and energy                 3 Contact 3                                var 64'
TANGENT = 'Contact force and energy                 3 Contact 3                                var 67'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def recorded_path(path):
    try:
        return str(Path(path).resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def history_array(assessment):
    return np.asarray([[row['strain'], row['nominal_stress_MPa']]
                       for row in assessment['reaction_history']], float)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, default=ROOT/'results/opencourant/thin_free_x_ams1e6.radioss.json')
    parser.add_argument('--half-ams', type=Path, default=ROOT/'results/opencourant/thin_free_x_ams5e7.radioss.json')
    parser.add_argument('--output-dir', type=Path, default=ROOT/'validation')
    parser.add_argument('--force-audit', type=Path, default=ROOT/'validation/opencourant_total_contact_force_audit.json')
    args = parser.parse_args()
    reference_path = ROOT/'validation/paper_reference/dong2019_fig12_thin.csv'
    references = read_reference(reference_path)
    experiments = [references[key][:, :2] for key in ('experiment_1', 'experiment_2')]
    proof = json.loads(args.force_audit.read_text())
    if proof.get('definition') != 'Contact impulse on SECONDARY = (FNX,FNY,FNZ) - (FTX,FTY,FTZ). Opposite resultant acts on MAIN.':
        raise ValueError('Use the controlled TYPE25 secondary-force convention audit')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    curves, cases = [], []
    for label, path, filename in (
            ('baseline_ams_1e-6', args.baseline, 'opencourant-thin-total-force-assessment.json'),
            ('half_ams_5e-7', args.half_ams, 'opencourant-thin-total-force-half-ams-assessment.json')):
        metadata = json.loads(path.read_text())
        csv_path = path.parent/(metadata['prefix']+'T01.csv')
        result = assess(path, csv_path, force_column=NORMAL, tangential_column=TANGENT,
                        force_sign=-1, force_kind='impulse')
        normal = assess(path, csv_path, force_column=NORMAL, force_sign=-1, force_kind='impulse')
        if not all(result[key] is True for key in ('completed', 'run_binding_verified', 'metadata_binding_verified')):
            raise ValueError('Only complete, input/output/metadata-bound actual runs can enter this comparison')
        assessment_path = args.output_dir/filename
        assessment_path.write_text(json.dumps(result, indent=2)+'\n')
        data, normal_data = history_array(result), history_array(normal)
        curves.append(data)
        grid = np.linspace(.03, .64, 500)
        total_stress = np.interp(grid, data[:, 0], data[:, 1])
        normal_stress = np.interp(grid, normal_data[:, 0], normal_data[:, 1])
        if not data[0, 0] <= grid[0] < grid[-1] <= data[-1, 0]:
            raise ValueError('Do not extrapolate the common comparison interval')
        samples = {}
        for strain in (.12, .35, .64):
            samples[str(strain)] = dict(
                total_contact_MPa=float(np.interp(strain, data[:, 0], data[:, 1])),
                legacy_normal_only_MPa=float(np.interp(strain, normal_data[:, 0], normal_data[:, 1])),
                paper_experiment_MPa=[float(np.interp(strain, curve[:, 0], curve[:, 1])) for curve in experiments])
        cases.append(dict(
            case=label, duration_seconds=metadata['duration_seconds'],
            ams_requested_time_step_seconds=metadata.get('ams', {}).get('requested_time_step_seconds'),
            assessment=filename, assessment_sha256=digest(assessment_path),
            metadata_path=recorded_path(path), metadata_sha256=digest(path),
            raw_history_path=recorded_path(csv_path), raw_history_sha256=digest(csv_path),
            execution_record_path=recorded_path(path.parent/(metadata['prefix']+'.run.json')),
            execution_record_sha256=digest(path.parent/(metadata['prefix']+'.run.json')),
            native_input_sha256=metadata['sha256'],
            hardening_source_sha256=metadata.get('hardening_source_sha256'),
            force_mapping=result['force_mapping'],
            total_force_curve_objective=curve_objective(data, experiments, interval=(.03, .64)),
            historical_normal_only_objective=curve_objective(normal_data, experiments, interval=(.03, .64)),
            normal_only_omission_rms_MPa=float(np.sqrt(np.mean((total_stress-normal_stress)**2))),
            stress_samples=samples, last_force_strain=float(data[-1, 0]),
            last_force_nominal_stress_MPa=float(data[-1, 1]),
            observed_history_strain=result['last_observed_nominal_strain'],
            eroded_elements=len(result['deleted_solid_elements']),
            core_energy_screen=result['core_energy_screen'],
            global_energy_screen=result['global_energy_screen'],
            quasi_static_validated=False, eligible_for_calibration_selection=False))
    grid = np.linspace(.03, .64, 500)
    y0, y1 = (np.interp(grid, data[:, 0], data[:, 1]) for data in curves)
    difference = y1-y0
    integrate = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
    work0, work1 = (float(integrate(y, grid)) for y in (y0, y1))
    sensitivity = dict(
        strain_interval=[.03, .64], samples=500,
        total_force_rmse_difference_MPa=float(np.sqrt(np.mean(difference**2))),
        relative_total_force_curve_L2_difference=float(np.sqrt(np.mean(difference**2)/np.mean(y0**2))),
        nominal_work_per_volume_MPa_baseline=work0, nominal_work_per_volume_MPa_half_ams=work1,
        relative_work_difference_half_minus_baseline=float((work1-work0)/work0),
        quasi_static_validated=False, timestep_convergence_established=False,
        qualification='Two actual AMS steps only; no extrapolation, mesh convergence or quasi-static proof.')
    report = dict(
        status='Corrected postprocessing of two existing completed native runs; calibration remains blocked',
        script='scripts/reassess_total_contact.py', script_sha256=digest(__file__),
        force_quantity='Compressive total Y force on core secondary = -d(FNY-FTY)/dt; interval-averaged N',
        force_convention_audit=dict(path=recorded_path(args.force_audit), sha256=digest(args.force_audit)),
        reference_csv='validation/paper_reference/dong2019_fig12_thin.csv', reference_sha256=digest(reference_path),
        cases=cases, ams_sensitivity=sensitivity,
        calibrated=False, experimental_agreement=False, eligible_for_calibration_selection=False,
        limitations=[
            'Thin-wall results are diagnostics/holdout, not training-data candidate selection.',
            'Native core/platen all-face geometry and AMS dynamics remain unqualified.',
            'Historical normal-only plots/assessments are preserved and superseded for total contact force.',
            'These existing runs retain their original higher-precision hardening inputs; no rounded user table is retrofitted.',
            'CSV finite differences yield interval-averaged forces; T01 REALR4/default-converter precision limits pointwise detail.',
            'Paper experiments are digitized raster lines, not raw laboratory measurements.',
            'No missing fracture-locus coordinates are inferred or calibrated by correcting force bookkeeping.'])
    target = args.output_dir/'opencourant-thin-total-force-qualification.json'
    target.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(report=str(target), cases=[dict(case=row['case'],
                          total_force_RMSE_MPa=row['total_force_curve_objective']['overlap_rmse_MPa'],
                          stress_at_strain_012_MPa=row['stress_samples']['0.12']['total_contact_MPa'])
                          for row in cases], ams_sensitivity=sensitivity), indent=2))


if __name__ == '__main__':
    main()
