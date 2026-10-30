"""Audit Fig.12 targets and actual assessments before selecting material inputs.

Produces a reproducible JSON report and full-range comparison figure. This does
not fit the thin holdout or turn a curve match into calibration qualification.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from plot_paper_comparison import read_reference
from reentrant_fea.calibration import curve_objective
from reentrant_fea.paper_metrics import curve_metrics

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build_report(assessments, reference_dir):
    reference_dir = Path(reference_dir)
    thin_path = reference_dir/'dong2019_fig12_thin.csv'
    thick_path = reference_dir/'dong2019_fig12_thick.csv'
    thin, thick = read_reference(thin_path), read_reference(thick_path)
    refs = {}
    for panel, curves in (('thin', thin), ('thick', thick)):
        refs[panel] = {name: curve_metrics(values[:, :2], panel=panel)
                       for name, values in curves.items() if name.startswith('experiment_')}
    evaluated = []
    for path in assessments:
        path = Path(path)
        a = json.loads(path.read_text())
        values = np.asarray([[row['strain'], row['nominal_stress_MPa']]
                             for row in a['reaction_history']], float)
        area = a.get('reference_area_mm2', a.get('stress_area_convention_mm2'))
        if area != 1900.:
            raise ValueError('Fig12 thin comparison requires declared reference area 1900 mm^2')
        # All supplied assessment paths are explicitly thin-wall comparisons;
        # callers must not pass a thick training assessment as a thin holdout.
        declared = a.get('generation', {}).get('config_wall_mm', a.get('config_wall_mm'))
        if declared is not None and not np.isclose(declared, .3):
            raise ValueError('Expected thin-wall assessment')
        metric = curve_metrics(values, panel='thin')
        objectives = {}
        for label, interval in (('early', (.03, .12)), ('plateau', (.1, .55)),
                                ('full_visible', (.03, .77))):
            try:
                objectives[label] = curve_objective(values,
                    [thin[k][:, :2] for k in ('experiment_1', 'experiment_2')], interval)
            except ValueError as error:
                objectives[label] = {'available': False, 'reason': str(error)}
        evaluated.append({'assessment': str(path), 'sha256': sha(path),
            'backend': a.get('backend', 'CalculiX'),
            'config_wall_mm': declared,
            'wall_thickness_evidence': ('declared in source assessment' if declared is not None else 'caller supplied as thin-wall diagnostic; thickness not independently verified here'),
            'matched_solver_comparison': False,
            'solver_completed': a.get('completed') is True,
            'reference_area_mm2': area,
            'force_mapping': a.get('force_mapping', a.get('curve_source', 'see source assessment')),
            'source_quasi_static_validated': a.get('quasi_static_validated') is True,
            'comparison_scope': 'caller-supplied thin-wall assessment; not an independently rerun solver',
            'metrics': metric, 'discrepancy_to_both_experiments': objectives,
            'calibration_eligible': False})
    return {'status': 'calibration incomplete; numerical and geometry prerequisites unresolved',
        'calibrated': False, 'selected_material': None,
        'paper': {'doi': '10.1016/j.matdes.2019.108036', 'supplied_panel': 'Fig12(b), t=0.3 mm',
                  'training': 'Fig12(a), t=1 mm', 'holdout': 'Fig12(b), t=0.3 mm',
                  'paper_fe_is_experimental_target': False,
                  'reported_thin_plateau_MPa': .28, 'reported_thin_densification_strain': .58,
                  'morphology_target_strains': [.12, .35, .65]},
        'reference_sha256': {p.name: sha(p) for p in (thin_path, thick_path,
                                    reference_dir/'dong2019_fig12_provenance.json')},
        'reference_diagnostics': refs, 'assessments': evaluated,
        'fixed_material_inputs': {'young_MPa': 69000., 'poisson': .3, 'density_Mg_mm3': 2.67e-9,
             'plastic_table': 'already true stress / plastic strain; no reconversion',
             'fracture_anchor': [.333, .047], 'failure_displacement_mm': .001},
        'selection_requirements': [
            'Qualify thick and thin geometry, relative density and meshes separately.',
            'Resolve undamaged elastic stiffness discrepancy with matched BC/formulation controls.',
            'Complete surface-to-RP and direct-face controls with observed relative displacement.',
            'Bind actual solver completion and input/output receipts for thick training candidates.',
            'Qualify core kinetic/internal energy, energy/work balance, speed and numerical sensitivity.',
            'Fit missing locus coordinates only on thick-wall experiments; transfer unchanged to thin holdout.',
            'Cover late thin response through .77 and compare .12/.35/.65 stored deformation frames.',
            'Assess identifiability and uncertainty before claiming a calibrated material.'],
        'metric_limits': 'Observed energy origins can differ between curves; these are not directly comparable total-energy measurements. Intervals and efficiency densification are investigator diagnostics, not paper definitions. Missing coverage is not extrapolated; printed curves have digitization uncertainty.'}


def plot_report(report, reference_dir, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ref = read_reference(Path(reference_dir)/'dong2019_fig12_thin.csv')
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    colors = {'experiment_1': 'red', 'experiment_2': 'black', 'paper_fe': 'blue'}
    for ax in axes:
        for name, values in ref.items():
            ax.plot(values[:, 0], values[:, 1], color=colors[name],
                    ls='-.' if name=='paper_fe' else '-', label=name.replace('_', ' '))
        for entry in report['assessments']:
            data = json.loads(Path(entry['assessment']).read_text())['reaction_history']
            ax.plot([r['strain'] for r in data], [r['nominal_stress_MPa'] for r in data],
                    label=(('OpenCourant contact FN-FT' if entry['backend'].startswith('OpenCourant')
                            else 'CalculiX RP, separate model')+' (unqualified)'))
        ax.set_xlabel('Nominal compressive strain'); ax.set_ylabel('Nominal stress (MPa)')
        ax.grid(alpha=.25)
    axes[0].set_xlim(0, .8); axes[0].set_title('Full Fig12(b) range; no extrapolation')
    axes[1].set_xlim(0, .12); axes[1].set_title('Early response through strain 0.12')
    axes[0].legend(fontsize=8)
    fig.suptitle('Paper revision diagnostics — calibration remains incomplete')
    fig.tight_layout(); fig.savefig(output.with_suffix('.png'), dpi=180)
    fig.savefig(output.with_suffix('.pdf')); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assessment', type=Path, action='append', required=True)
    p.add_argument('--reference-dir', type=Path, default=ROOT/'validation/paper_reference')
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    report = build_report(args.assessment, args.reference_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
    plot_report(report, args.reference_dir, args.output)
    print(args.output.with_suffix('.json'))


if __name__ == '__main__':
    main()
