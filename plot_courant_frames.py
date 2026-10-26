#!/usr/bin/env python3
"""Render actual OpenCourant animation frames without retaining eroded cells.

The official anim_to_vtk converter supplies deformed coordinates and explicit
erosion status. No interpolation between frames or deformation amplification is
performed. Strains follow the exact imposed-displacement samples in metadata.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np


def _block(text, header, size, dtype=float):
    match = re.search(header, text, re.MULTILINE)
    if not match:
        raise ValueError("Missing VTK block: " + header)
    tail = text[match.end():]
    end = re.search(r"\n(?:POINTS|CELLS|CELL_TYPES|POINT_DATA|CELL_DATA|SCALARS|VECTORS|TENSORS|FIELD)\b", tail)
    if end:
        tail = tail[:end.start()]
    values = np.fromstring(tail, sep=" ", dtype=dtype)
    if values.size != size:
        raise ValueError(f"VTK block has {values.size} values, expected {size}")
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite VTK values")
    return values


def frame_time(text):
    match = re.search(r"^TIME\s+1\s+1\s+\w+\s*\n([^\n]+)", text, re.MULTILINE)
    if not match:
        raise ValueError("Animation VTK lacks TIME field")
    value = float(match.group(1).strip().replace("D", "E"))
    if not np.isfinite(value) or value < 0:
        raise ValueError("Invalid animation timestamp")
    return value


def read_frame(text):
    if "\nASCII\nDATASET UNSTRUCTURED_GRID\n" not in text[:150]:
        raise ValueError("Require official ASCII unstructured-grid VTK output")
    points_header = re.search(r"^POINTS (\d+) \w+\s*$", text, re.MULTILINE)
    cells_header = re.search(r"^CELLS (\d+) (\d+)\s*$", text, re.MULTILINE)
    if not points_header or not cells_header:
        raise ValueError("Missing VTK mesh headers")
    npoints = int(points_header.group(1)); ncells, nitems = map(int, cells_header.groups())
    points = _block(text, r"^POINTS \d+ \w+\s*$", 3*npoints).reshape(-1, 3)
    values = _block(text, r"^CELLS \d+ \d+\s*$", nitems, int)
    if nitems != ncells*9:
        raise ValueError("Only eight-node hexahedral animation cells are supported")
    values = values.reshape(ncells, 9)
    if not np.all(values[:, 0] == 8):
        raise ValueError("Only eight-node hexahedral animation cells are supported")
    cells = values[:, 1:]
    if np.any(cells < 0) or np.any(cells >= npoints):
        raise ValueError("Invalid VTK connectivity")
    types = _block(text, r"^CELL_TYPES \d+\s*$", ncells, int)
    if not np.all(types == 12):
        raise ValueError("Require VTK_HEXAHEDRON cell type")
    displacement = _block(text, r"^VECTORS Displacement \w+\s*$", 3*npoints).reshape(-1, 3)
    scalar = lambda name: _block(text, rf"^SCALARS {name} \w+ 1\s*\nLOOKUP_TABLE \w+\s*$", ncells)
    erosion = scalar("EROSION_STATUS")
    if not np.isin(erosion, [0, 1]).all():
        raise ValueError("Unsupported erosion status; expected 0=eroded, 1=active")
    return {"time": frame_time(text), "points": points, "initial_points": points-displacement,
            "cells": cells, "active": erosion == 1, "part_id": scalar("PART_ID"),
            "element_id": scalar("ELEMENT_ID"), "damage": scalar("3DELEM_MAX_DAMAGE_ELEMENT")}


def front_polygons(frame, part=1, tol=1e-3):
    """Select actual initial front faces of active hexes, keeping current coords."""
    selected = (frame["part_id"] == part) & frame["active"]
    indices = np.flatnonzero(selected)
    if not len(indices):
        return [], np.empty(0)
    cells = frame["cells"]
    # VTK's two extrusion caps; depth layers beyond the front are excluded by
    # their initial Z coordinates, not their potentially folded current Z.
    zmin = float(frame["initial_points"][cells[frame["part_id"] == part], 2].min())
    polygons, cell_indices = [], []
    for index in indices:
        for face in (cells[index, :4], cells[index, 4:]):
            if np.all(np.abs(frame["initial_points"][face, 2]-zmin) <= tol):
                polygons.append(frame["points"][face, :2])
                cell_indices.append(index)
    return polygons, np.asarray(cell_indices, dtype=int)


def strain_at_time(time, metadata):
    samples = np.asarray(metadata["displacement_fraction_samples"], dtype=float)
    if (samples.ndim != 2 or samples.shape[1] != 2 or len(samples) < 2
            or not np.isfinite(samples).all() or np.any(np.diff(samples[:, 0]) <= 0)
            or np.any(np.diff(samples[:, 1]) < 0)):
        raise ValueError("Invalid imposed displacement samples")
    if time < samples[0, 0]-1e-10 or time > samples[-1, 0]+1e-10:
        raise ValueError("Frame timestamp falls outside imposed loading history")
    return float(np.interp(time, samples[:, 0], samples[:, 1])*metadata["target_nominal_strain"])


def projected_diagnostics(polygons, element_ids, tolerance=1e-8):
    """Screen active front projections only; this is not a 3-D contact proof."""
    from shapely.geometry import Polygon
    from shapely.strtree import STRtree
    valid, ids, invalid = [], [], []
    for points, eid in zip(polygons, element_ids):
        area = .5*np.sum(points[:, 0]*np.roll(points[:, 1], -1)
                         - points[:, 1]*np.roll(points[:, 0], -1))
        poly = Polygon(points)
        if area <= 0 or not poly.is_valid:
            invalid.append(int(eid))
        else:
            valid.append(poly); ids.append(int(eid))
    tree = STRtree(valid)
    count = 0; maximum = 0.; total = 0.
    for i, poly in enumerate(valid):
        for j in tree.query(poly):
            if j <= i:
                continue
            area = poly.intersection(valid[int(j)]).area
            if area > tolerance:
                count += 1; maximum = max(maximum, area); total += area
    return {"active_front_faces": len(polygons), "invalid_projected_element_ids": invalid,
            "overlap_pairs": count, "max_overlap_mm2": maximum, "sum_pair_overlap_mm2": total,
            "qualification": "front projection of active cells only; does not establish complete 3-D contact validity"}


def convert(path, converter):
    before = path.stat()
    result = subprocess.run([str(converter), str(path.resolve())], capture_output=True, text=True)
    if result.returncode:
        raise ValueError(f"Animation conversion failed for {path.name}: {result.stderr.strip()}")
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError("Animation changed during conversion: " + path.name)
    return result.stdout


def paper_morphology_comparison(pdf, targets, summaries, frames, output):
    """Use user-supplied Fig. 11(b) crops; never vendor the paper's images."""
    from PIL import Image
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    # Positions in the source's original 1472x2680 embedded JPEG. Printed row
    # captions were independently checked on rendered PDF page 9.
    crops = {
        .12: [(80, 1720, 590, 2030), (920, 1720, 1450, 2045)],
        .35: [(100, 2145, 610, 2425), (930, 2145, 1465, 2430)],
        .65: [(100, 2540, 610, 2680), (930, 2540, 1465, 2680)],
    }
    if any(target not in crops for target in targets):
        raise ValueError("Paper morphology crops are available only for targets .12, .35, .65")
    with tempfile.TemporaryDirectory(prefix="dong2019-fig11-") as temp:
        prefix = Path(temp)/"fig11"
        subprocess.run(["pdfimages", "-f", "9", "-l", "9", "-j", str(pdf), str(prefix)], check=True)
        image = Image.open(Path(str(prefix)+"-000.jpg")).convert("RGB")
        if image.size != (1472, 2680):
            raise ValueError("Expected original Fig. 11 raster (1472x2680); wrong paper/version")
        reference = [[np.asarray(image.crop(box)) for box in crops[target]] for target in targets]
    fig, axes = plt.subplots(len(targets), 3, figsize=(12, 3.1*len(targets)), squeeze=False)
    fig.subplots_adjust(top=.82, bottom=.08, left=.08, right=.97, hspace=.28, wspace=.12)
    for i, (target, summary, frame) in enumerate(zip(targets, summaries, frames)):
        for j in (0, 1):
            axes[i, j].imshow(reference[i][j])
            axes[i, j].set_axis_off()
        axes[i, 0].text(-.03, .5, f"Paper strain\n{target:.2f}", transform=axes[i, 0].transAxes,
                        ha="right", va="center", fontsize=10)
        current = axes[i, 2]
        if frame is None:
            current.text(.5, .5, "No actual frame available", transform=current.transAxes, ha="center")
        else:
            polygons, _ = front_polygons(frame)
            current.add_collection(PolyCollection(polygons, facecolors="#2368a2", edgecolors="none"))
            # Orphan nodes of eroded cells may move far outside the surviving
            # specimen. They neither render nor determine the active-mesh view.
            used = np.unique(frame["cells"][(frame["part_id"] == 1) & frame["active"]])
            xy = frame["points"][used, :2]
            current.set_xlim(xy[:, 0].min()-.5, xy[:, 0].max()+.5)
            current.set_ylim(xy[:, 1].min()-.5, xy[:, 1].max()+.5)
            current.set_aspect("equal")
            current.set_title(f"Actual strain {summary['actual_strain']:.5f}; {summary['eroded_core_cells']} eroded", fontsize=9)
        current.set_axis_off()
    for ax, label in zip(axes[0], ["Paper experiment (Fig. 11b)", "Paper FE (Fig. 11b)", "Current OpenCourant"]):
        position = ax.get_position()
        fig.text(position.x0+position.width/2, .91, label, ha="center", va="center", fontsize=11)
    fig.suptitle("Thin wall t=0.3 mm — morphology comparison", fontsize=15, y=.98)
    fig.text(.5, .025, "Source: Dong et al. (2019), user-supplied PDF. Current frames use true deformation, omit eroded cells, and are not interpolated.\nImages have different camera/rendering settings; Bi-V/Z/collapse agreement is not established by this visual comparison.",
             ha="center", fontsize=8)
    for suffix in (".png", ".pdf"):
        fig.savefig(output.with_name(output.name+"_paper_comparison").with_suffix(suffix), dpi=210)
    return {"pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(), "source": "Dong et al. (2019), Fig. 11(b), PDF page 9",
            "method": "Crops of original embedded raster at verified printed strain rows; separate actual current meshes; no paper image copied into tracked repository files"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("metadata", type=Path)
    ap.add_argument("--converter", type=Path,
                    default=Path("/workspace/opencourant/runtime/OpenCourant/exec/anim_to_vtk_linux64_gf"))
    ap.add_argument("--targets", nargs="+", type=float, default=[.12, .35, .65])
    ap.add_argument("--maximum-strain-error", type=float, default=.005)
    ap.add_argument("--output", type=Path, default=Path("results/opencourant/deformation_frames"))
    ap.add_argument("--paper-pdf", type=Path,
                    help="Optional user-supplied Dong PDF for a separate Fig. 11(b) morphology comparison")
    args = ap.parse_args()
    metadata = json.loads(args.metadata.read_text())
    if not np.isclose(metadata.get("wall_mm", float("nan")), .3):
        raise ValueError("This plot requires verified thin-wall metadata (wall_mm=0.3)")
    if not args.targets or not np.isfinite(args.targets).all() or min(args.targets) < 0:
        raise ValueError("Require finite nonnegative target strains")
    if not np.isfinite(args.maximum_strain_error) or args.maximum_strain_error <= 0:
        raise ValueError("Require positive maximum strain error")
    prefix = metadata["prefix"]
    files = [p for p in args.metadata.parent.glob(prefix+"A*")
             if re.fullmatch(re.escape(prefix)+r"A\d+", p.name)]
    if not files:
        raise ValueError("No animation files found for metadata prefix")
    available = []
    failures = []
    for path in sorted(files):
        try:
            text = convert(path, args.converter)
            time = frame_time(text)
            available.append((strain_at_time(time, metadata), time, path))
        except ValueError as exc:
            failures.append(str(exc))
    if not available:
        raise ValueError("No readable animation frames: " + "; ".join(failures))
    os.environ.setdefault("MPLCONFIGDIR", "/workspace/mplconfig")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize

    fig, axes = plt.subplots(1, len(args.targets), figsize=(5*len(args.targets), 5.1), squeeze=False)
    fig.subplots_adjust(left=.045, right=.91, bottom=.17, top=.79, wspace=.22)
    summaries = []
    selected_frames = []
    collections = []
    for target, ax in zip(args.targets, axes[0]):
        strain, time, path = min(available, key=lambda entry: abs(entry[0]-target))
        summary = {"target_strain": target, "available": abs(strain-target) <= args.maximum_strain_error}
        if summary["available"]:
            text = convert(path, args.converter)
            frame = read_frame(text)
            selected_frames.append(frame)
            polygons, indices = front_polygons(frame, part=1)
            colors = frame["damage"][indices]
            collection = PolyCollection(polygons, array=colors, cmap="viridis", norm=Normalize(0, 1),
                                        linewidths=0, edgecolors="none")
            ax.add_collection(collection)
            collections.append(collection)
            for part in (2, 3):
                platens, _ = front_polygons(frame, part=part)
                ax.add_collection(PolyCollection(platens, facecolors="#aab2bc", edgecolors="#4b5563", linewidths=.5))
            # The existing geometry spans the same initial 48 mm height.
            initial = frame["initial_points"]
            used = np.unique(frame["cells"][frame["active"]])
            current = frame["points"][used]
            ax.set_xlim(min(initial[:, 0].min(), current[:, 0].min())-2,
                        max(initial[:, 0].max(), current[:, 0].max())+2)
            ax.set_ylim(initial[:, 1].min()-1.5, initial[:, 1].max()+1.5)
            core = frame["part_id"] == 1
            eroded = int(np.sum(core & ~frame["active"]))
            ax.set_title(f"Target {target:.2f}; actual {strain:.5f}\nt = {time:.6f} s; eroded {eroded}/{core.sum()}", fontsize=10)
            summary.update(animation=str(path), animation_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                           actual_strain=strain, actual_time_seconds=time, strain_error=strain-target,
                           total_core_cells=int(core.sum()), eroded_core_cells=eroded,
                           active_core_front_faces=len(polygons),
                           active_max_element_damage=float(np.max(colors)) if len(colors) else None,
                           projected_geometry=projected_diagnostics(polygons, frame["element_id"][indices]))
        else:
            selected_frames.append(None)
            ax.text(.5, .5, f"Target {target:.2f}\nNo actual frame within\n±{args.maximum_strain_error:g} strain", ha="center", va="center", transform=ax.transAxes)
            summary.update(nearest_actual_strain=strain, nearest_actual_time_seconds=time)
            ax.set_axis_off()
        ax.set_aspect("equal")
        ax.set_xlabel("X (mm)"); ax.set_ylabel("Y (mm)")
        summaries.append(summary)
    if collections:
        cax = fig.add_axes([.935, .25, .012, .43])
        fig.colorbar(collections[0], cax=cax, label="Active-element damage (solver output)")
    fig.suptitle("Thin-wall t=0.3 mm — OpenCourant actual deformed frames", fontsize=14, y=.98)
    fig.supxlabel("Front projection; true displacement scale; eroded elements excluded. Morphology and experimental agreement remain unvalidated.", fontsize=9, y=.03)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".pdf"):
        fig.savefig(args.output.with_suffix(suffix), dpi=210)
    report = {"metadata": str(args.metadata), "metadata_sha256": hashlib.sha256(args.metadata.read_bytes()).hexdigest(),
              "converter": str(args.converter), "converter_sha256": hashlib.sha256(args.converter.read_bytes()).hexdigest(),
              "method": "Official ASCII VTK conversion; deformed POINTS rendered directly; EROSION_STATUS==1 required; front face selected from initial POINTS-Displacement Z.",
              "strain_source": "Piecewise-linear interpolation of exact imposed displacement_fraction_samples at actual VTK TIME, multiplied by target_nominal_strain; no mesh interpolation.",
              "frames": summaries, "conversion_failures": failures,
              "experimental_mode_agreement": "not established", "deformation_amplification": 1,
              "eroded_elements_rendered": False, "mesh_interpolation": False}
    if args.paper_pdf:
        report["paper_morphology_reference"] = paper_morphology_comparison(
            args.paper_pdf, args.targets, summaries, selected_frames, args.output)
    args.output.with_suffix(".json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
