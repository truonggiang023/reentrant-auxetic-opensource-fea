#!/usr/bin/env python3
"""Compare recorded thin-wall solver history with approximate published traces.

This plots evidence without extending a stopped run or labelling it calibrated.
The reference CSVs are digitized printed lines, not raw experimental records.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def solver_label(assessment):
    """Keep backend identity and dynamics qualification visible on the figure."""
    declared = assessment.get('backend', 'CalculiX')
    name = 'OpenCourant' if declared.startswith('OpenCourant') else declared
    complete = 'completed' if assessment.get('completed') is True else 'incomplete'
    if name == 'CalculiX':
        return name, f'Current CCX (accepted history, {complete})'
    qualification = '' if assessment.get('quasi_static_validated') is True else '; quasi-static unverified'
    if name == 'OpenCourant':
        mapping = assessment.get('force_mapping', {})
        quantity = mapping.get('contact_force_quantity')
        if quantity == 'total_force_on_secondary' and mapping.get('component_operation') == 'FN-FT':
            return name, f'{name} total FN-FT\n({complete}{qualification})'
        column = mapping.get('csv_column', '')
        if quantity == 'normal_component_only' or column.startswith('Contact force and energy'):
            return name, f'{name} legacy normal component only\n({complete}{qualification})'
    return name, f'{name} ({complete}{qualification})'


def read_reference(path: Path) -> dict[str, np.ndarray]:
    groups: dict[str, list] = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            groups.setdefault(row["curve"], []).append([
                float(row["strain"]), float(row["stress_MPa"]), float(row["y_span_MPa"])])
    if set(groups) != {"experiment_1", "experiment_2", "paper_fe"}:
        raise ValueError("Reference must contain two experiments and the paper FE curve")
    result = {key: np.asarray(values) for key, values in groups.items()}
    for key, data in result.items():
        if len(data) < 2 or not np.isfinite(data).all() or np.any(np.diff(data[:, 0]) <= 0):
            raise ValueError(f"Invalid or unsorted reference curve: {key}")
    return result


def read_history(path: Path):
    assessment = json.loads(path.read_text())
    rows = assessment.get("reaction_history", [])
    if len(rows) < 2:
        raise ValueError("At least two accepted reaction-history samples are required")
    data = np.asarray([[row["strain"], row["nominal_stress_MPa"]] for row in rows], dtype=float)
    if not np.isfinite(data).all() or np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError("Reaction history must be finite and strictly increasing in strain")
    # Existing assess_run writes this history only from committed frames. Its
    # own status/qualification are carried into the comparison summary.
    return assessment, data


def resolve_wall_mm(assessment: dict, explicit_wall_mm: float | None = None,
                    root: Path = ROOT) -> tuple[float, str]:
    """Require geometry evidence or an explicit thickness assertion from caller.

    Legacy native-final lacks generation metadata. Recognize it only when the
    complete assessment equals the committed qualification record and that
    record's configuration hash matches the actual configuration file.
    """
    values = []
    generation = assessment.get("generation", {})
    mesh = assessment.get("mesh", {})
    for source, value in (
        ("assessment.config_wall_mm", assessment.get("config_wall_mm")),
        ("generation.config_wall_mm", generation.get("config_wall_mm") if isinstance(generation, dict) else None),
        ("mesh.config_wall_mm", mesh.get("config_wall_mm") if isinstance(mesh, dict) else None),
    ):
        if value is not None:
            if isinstance(value, bool):
                raise ValueError(f"Invalid wall thickness in {source}")
            value = float(value)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"Invalid wall thickness in {source}")
            values.append((source, value))
    if values:
        source, value = values[0]
        if any(not np.isclose(value, other, rtol=0, atol=1e-9) for _, other in values):
            raise ValueError("Conflicting wall thickness metadata in assessment")
    else:
        source = ""
        value = None
        qualification_path = root/"validation/native_qualification.json"
        config_path = root/"configs/dong2019_brick_damage.json"
        if qualification_path.is_file() and config_path.is_file():
            qualification = json.loads(qualification_path.read_text())
            recorded = qualification.get("primary_run", {}).get("results")
            expected_hash = qualification.get("source_sha256", {}).get("configs/dong2019_brick_damage.json")
            actual_hash = hashlib.sha256(config_path.read_bytes()).hexdigest()
            if assessment == recorded and expected_hash == actual_hash:
                value = float(json.loads(config_path.read_text())["wall"])
                source = "matched committed qualification assessment and SHA256-verified configuration"
    if explicit_wall_mm is not None:
        explicit_wall_mm = float(explicit_wall_mm)
        if not np.isfinite(explicit_wall_mm) or explicit_wall_mm <= 0:
            raise ValueError("--wall-mm must be finite and positive")
        if value is not None and not np.isclose(value, explicit_wall_mm, rtol=0, atol=1e-9):
            raise ValueError("--wall-mm conflicts with assessment geometry evidence")
        if value is None:
            value = explicit_wall_mm
            source = "caller assertion via --wall-mm; geometry not independently verified"
    if value is None:
        raise ValueError("Assessment has no verified wall thickness. Supply generation/config_wall_mm metadata or explicitly assert --wall-mm 0.3 for a thin-wall run.")
    if not np.isclose(value, 0.3, rtol=0, atol=1e-9):
        raise ValueError(f"Thin-wall reference requires t=0.3 mm; assessment thickness is {value:g} mm")
    return value, source


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--assessment", type=Path, default=ROOT/"validation/native-final.json")
    ap.add_argument("--reference", type=Path,
                    default=ROOT/"validation/paper_reference/dong2019_fig12_thin.csv")
    ap.add_argument("--wall-mm", type=float,
                    help="Explicit caller assertion of wall thickness for legacy assessments without geometry provenance")
    ap.add_argument("--case-label", help="Visible case description, such as a legacy boundary sensitivity")
    ap.add_argument("--output", type=Path, default=ROOT/"results/paper_comparison/thin_wall_comparison")
    args = ap.parse_args()
    reference = read_reference(args.reference)
    assessment, history = read_history(args.assessment)
    wall_mm, wall_evidence = resolve_wall_mm(assessment, args.wall_mm)
    backend, run_label = solver_label(assessment)
    os.environ.setdefault("MPLCONFIGDIR", "/workspace/mplconfig")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    fig, (ax, early) = plt.subplots(1, 2, figsize=(12, 4.8),
                                   gridspec_kw={"width_ratios": [1.5, 1]}, layout="constrained")
    style = {"experiment_1": ("#d62728", "-", "Paper experiment 1 (digitized)"),
             "experiment_2": ("#222222", "-", "Paper experiment 2 (digitized)"),
             "paper_fe": ("#2455cf", "-.", "Paper FE (digitized)")}
    completed = assessment.get("completed") is True
    for panel in (ax, early):
        for curve, data in reference.items():
            color, ls, label = style[curve]
            panel.plot(data[:, 0], data[:, 1], color=color, linestyle=ls, lw=1.15, label=label)
        panel.plot(history[:, 0], history[:, 1], color="#00836c", lw=2.0,
                   marker="o", markersize=2.4, label=run_label, zorder=5)
        panel.set_xlabel("Nominal compressive strain")
        panel.set_ylabel("Nominal compressive stress (MPa)")
        panel.grid(alpha=0.18)
        panel.spines[["top", "right"]].set_visible(False)
        panel.set_ylim(bottom=0)
    ax.set_xlim(0, 0.8)
    # Keep the published range visible without clipping a mismatching solver
    # peak above it; hiding such a peak would understate the discrepancy.
    stress_ceiling = max(1.4, float(history[:, 1].max())*1.1)
    ax.set_ylim(0, stress_ceiling)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.set_title("Thin wall, t = 0.3 mm — published range")
    end = float(history[-1, 0])
    if end < 0.8:
        ax.axvspan(end, 0.8, color="#00836c", alpha=0.035, zorder=0)
        ax.text((end+0.8)/2, stress_ceiling*.95, f"No {backend} history\nbeyond strain {end:.3f}",
                ha="center", va="top", fontsize=8, color="#286a5d")
    ax.legend(loc="upper center", bbox_to_anchor=(.50, .99), fontsize=7.5,
              frameon=True, framealpha=.92)
    early.set_xlim(0, 0.15)
    early.set_title(f"Early {backend} history")
    early.xaxis.set_major_locator(MultipleLocator(0.05))
    early_values = history[history[:, 0] <= .15, 1]
    early.set_ylim(0, max(0.4, float(early_values.max())*1.15) if len(early_values) else 0.4)
    early.text(0.03, 0.96, "Calibration / experimental agreement\nnot established by this plot",
               transform=early.transAxes, va="top", fontsize=9,
               bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"})
    title = f"Dong et al. (2019), Fig. 12(b), versus {backend}"
    if args.case_label:
        title += '\n'+args.case_label
    fig.suptitle(title, fontsize=13)
    fig.supxlabel("Reference: samples of published raster traces; missing and clipped segments are not extrapolated.", fontsize=8)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".pdf"):
        fig.savefig(args.output.with_suffix(suffix), dpi=220)
    summary = {
        "assessment": str(args.assessment),
        "assessment_sha256": hashlib.sha256(args.assessment.read_bytes()).hexdigest(),
        "reference_csv": str(args.reference),
        "reference_sha256": hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        "wall_mm": wall_mm,
        "wall_thickness_evidence": wall_evidence,
        "current_run_completed": completed,
        "backend": backend,
        "case_label": args.case_label,
        "quasi_static_validated": assessment.get("quasi_static_validated"),
        "force_mapping": assessment.get("force_mapping"),
        "accepted_reaction_samples": len(history),
        "recorded_strain_range": [float(history[0, 0]), float(history[-1, 0])],
        "last_converged_nominal_strain": assessment.get("last_converged_nominal_strain"),
        "stress_area_convention_mm2": assessment.get("stress_area_convention_mm2"),
        "assessment_qualification": assessment.get("qualification"),
        "experimental_agreement": "not established",
        "extrapolated": False,
        "limitations": ["Published reference is digitized, not raw experimental data.",
                        "No solver data is manufactured after the last accepted recorded frame.",
                        "Completion and calibration require independent validation beyond this comparison."],
    }
    args.output.with_suffix(".json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
