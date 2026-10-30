#!/usr/bin/env python3
"""Digitize Fig. 12's embedded raster curves, retaining source limitations.

Requires poppler's pdfimages, numpy, and Pillow. The paper PDF is supplied by the
user and is not redistributed. Outputs are approximate samples of published
lines, not raw experimental measurements or an identified fracture locus.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from PIL import Image


# Bounds were read from the original 1884x772 image, not a resized screenshot.
# The origin and axis maxima below correspond to labeled tick marks.
PANELS = {
    "thin": {"wall_mm": 0.3, "axis": [1094, 1862, 591, 21],
             "max_strain": 0.8, "max_stress_MPa": 1.4,
             "legend": [1130, 1460, 45, 180],
             "colors": {"experiment_1": "red", "experiment_2": "black", "paper_fe": "blue"}},
    "thick": {"wall_mm": 1.0, "axis": [119, 894, 591, 15],
              "max_strain": 0.5, "max_stress_MPa": 20.0,
              "legend": [170, 515, 45, 170],
              "colors": {"experiment_1": "black", "experiment_2": "red", "paper_fe": "blue"}},
}


def digitize_panel(rgb: np.ndarray, panel: dict) -> list[dict]:
    """Sample colored traces by column, excluding borders, ticks, and legends.

    Median coordinates avoid selecting anti-aliased line edges. A vertical trace
    segment cannot be reduced uniquely to y(x): its complete detected span is
    reported, and callers must not regard the median as an exact measurement.
    Missing columns (dash gaps, crossings, occlusion) remain absent.
    """
    if rgb.shape != (772, 1884, 3):
        raise ValueError(f"Expected original Fig. 12 RGB dimensions 1884x772, got {rgb.shape}")
    r, g, b = np.moveaxis(rgb.astype(int), 2, 0)
    masks = {
        "red": (r > 150) & (g < 130) & (b < 140) & (r-g > 65) & (r-b > 65),
        "blue": (b > 130) & (r < 155) & (g < 155) & (b-r > 45) & (b-g > 40),
        "black": (np.maximum.reduce([r, g, b])-np.minimum.reduce([r, g, b]) < 45)
                 & ((r+g+b) < 330),
    }
    x0, x1, y0, y1 = panel["axis"]
    lx0, lx1, ly0, ly1 = panel["legend"]
    rows = []
    for curve, color in panel["colors"].items():
        mask = masks[color].copy()
        mask[ly0:ly1+1, lx0:lx1+1] = False
        # Three pixels at each border exclude axes and inward-pointing ticks.
        # At the top use six: clipped curves have no recoverable y coordinate.
        # Right-axis inward ticks are black; remove their entire reach for that
        # curve rather than misreporting a tick as a late experimental sample.
        for x in range(x0+4, x1-(18 if color == "black" else 3)):
            ys = np.flatnonzero(mask[y1+6:y0-3, x]) + y1+6
            # Interior top ticks reach below the frame at x tick locations;
            # reject isolated axis fragments close to the top.
            ys = ys[ys >= y1+15]
            if len(ys) < 2:
                continue
            # A single curve can have a vertical segment, but unrelated stray
            # marks should not be joined. Split into contiguous pixel groups.
            groups = np.split(ys, np.flatnonzero(np.diff(ys) > 2)+1)
            groups = [group for group in groups if len(group) >= 2]
            if not groups:
                continue
            # Multiple groups usually mean intersections/JPEG artifacts. Keep
            # the most strongly supported group and report its span.
            group = max(groups, key=len)
            y = float(np.median(group))
            sy = panel["max_stress_MPa"]/(y0-y1)
            rows.append({"curve": curve,
                         "strain": (x-x0)*panel["max_strain"]/(x1-x0),
                         "stress_MPa": (y0-y)*sy,
                         "y_span_MPa": float(group[-1]-group[0])*sy})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("pdf", type=Path)
    ap.add_argument("--output-dir", type=Path,
                    default=Path(__file__).resolve().parents[1]/"validation/paper_reference")
    args = ap.parse_args()
    with tempfile.TemporaryDirectory(prefix="dong2019-digitize-") as tmp:
        prefix = Path(tmp)/"page10"
        subprocess.run(["pdfimages", "-f", "10", "-l", "10", "-j", str(args.pdf), str(prefix)], check=True)
        image_path = Path(str(prefix)+"-000.jpg")
        rgb = np.asarray(Image.open(image_path).convert("RGB"))
        image_sha = hashlib.sha256(image_path.read_bytes()).hexdigest()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries = {}
    for name, panel in PANELS.items():
        rows = digitize_panel(rgb, panel)
        target = args.output_dir/f"dong2019_fig12_{name}.csv"
        with target.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, ["curve", "strain", "stress_MPa", "y_span_MPa"],
                                    lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({key: f"{value:.9g}" if isinstance(value, float) else value
                                 for key, value in row.items()})
        summaries[name] = {**panel, "file": target.name, "counts": {
            curve: sum(row["curve"] == curve for row in rows) for curve in panel["colors"]},
            "pixel_size_strain": panel["max_strain"]/(panel["axis"][1]-panel["axis"][0]),
            "pixel_size_stress_MPa": panel["max_stress_MPa"]/(panel["axis"][2]-panel["axis"][3])}
    provenance = {
        "source": "Dong et al., Materials & Design 182 (2019) 108036, Fig. 12, PDF page 10",
        "doi": "10.1016/j.matdes.2019.108036",
        "pdf_sha256": hashlib.sha256(args.pdf.read_bytes()).hexdigest(),
        "embedded_image_sha256": image_sha,
        "embedded_image": "PDF page 10 first image, RGB JPEG 1884x772; extracted with pdfimages -j",
        "method": "Color thresholding of original embedded image, per-column median of largest contiguous trace group; legends and axes excluded. No interpolation of missing columns during extraction.",
        "uncertainty": "At least +/-2 image pixels per coordinate for line width, antialiasing and axis selection, plus the reported vertical y_span_MPa/2. This is a coordinate-resolution bound, not a statistical confidence interval. Larger errors are possible at intersections/occlusion; curves are samples of printed lines, not raw data.",
        "limitations": ["Original experiments are not available as raw data.",
                        "Trace crossings and dashed-line gaps can be missing or ambiguous.",
                        "Near-vertical segments are represented by a per-column median and span.",
                        "The top of the thin-wall curves is clipped by the published plot; no extrapolation is supplied.",
                        "Compression curves do not uniquely identify the fracture locus.",
                        "Experimental labels follow the legend of each panel; its colors swap between panels."],
        "panels": summaries,
    }
    (args.output_dir/"dong2019_fig12_provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    print(json.dumps({name: value["counts"] for name, value in summaries.items()}))


if __name__ == "__main__":
    main()
