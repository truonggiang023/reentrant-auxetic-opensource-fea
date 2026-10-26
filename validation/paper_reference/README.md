# Published comparison curves

These CSVs are approximate samples of the printed curves in Fig. 12 of Dong et
al., *Materials & Design* 182 (2019) 108036, DOI
10.1016/j.matdes.2019.108036. They are **not raw experimental measurements**, a
complete material calibration, or evidence that the current simulation agrees.

The source PDF embeds Fig. 12 as a 1884×772 JPEG. The extraction uses this original
image, without rescaling. `dong2019_fig12_provenance.json` records the source
hashes, axes, color mapping, coordinate resolution and limitations. No copy of
the paper or its image is redistributed here. The input is the user-supplied PDF.

Columns are `curve`, `strain`, `stress_MPa`, `y_span_MPa`. `curve` is
`experiment_1`, `experiment_2` or `paper_fe`, as labelled in the relevant panel.
The experiments' colors swap between thick and thin panels. `y_span_MPa` is the
vertical extent of detected trace pixels at that column, not a confidence
interval. Missing columns and clipped curve tops are absent. Coordinates have
at least ±2 source-pixel uncertainty plus half the vertical span; crossings and
occlusion can introduce larger error. Lines plotted between samples are visual
interpolation, not additional observations.

To reproduce the CSVs, from the repository root, with `pdfimages` (Poppler),
NumPy and Pillow available:

```sh
/workspace/fea-venv/bin/python scripts/digitize_dong2019.py /path/to/paper.pdf
```

To plot the thin-wall (`t = 0.3 mm`) published curves against the committed,
accepted reaction history:

```sh
MPLCONFIGDIR=/workspace/mplconfig /workspace/fea-venv/bin/python \
  plot_paper_comparison.py --assessment validation/native-final.json \
  --output results/paper_comparison/thin_wall_comparison
```

The command writes PNG, PDF and a JSON evidence summary. The default assessment
contains 45 recorded reaction samples through strain 0.143; its last accepted
increment is 0.14421875. Those are different sampling limits. Neither is extended
to the requested final strain. The paper FE curve is a separate published
reference, never a replacement for the current solver output. When supplying a
different assessment, verify that it is the thin-wall geometry and uses the same
nominal stress-area convention (the existing assessment uses 1900 mm²).
The plotting command checks `config_wall_mm` / `generation.config_wall_mm`
metadata and rejects a thickness other than 0.3 mm. The known legacy assessment
is recognized through its complete qualification record and the SHA256 of its
configuration. Other legacy assessments must explicitly supply `--wall-mm 0.3`;
the JSON summary identifies this as a caller assertion, not verified geometry.

For a calibration study, use the thick-wall curves as the fitting target and
the thin-wall curves as a held-out transfer check, consistent with Section 4.1.
Fitting the thin-wall data and then calling comparison with it independent
validation would be incorrect. Curve matching cannot identify a unique missing
fracture locus without additional constraints and sensitivity analysis.

## Actual explicit-solver deformation frames

`plot_courant_frames.py` converts existing OpenCourant animation files using the
official `anim_to_vtk` utility and renders the initial front faces of active
hexahedra at their actual deformed coordinates. It requires the generated
`.radioss.json` metadata to map the VTK timestamp through the exact imposed
displacement samples. It does not assume displacement is linear in time.

```sh
MPLCONFIGDIR=/workspace/mplconfig /workspace/fea-venv/bin/python \
  plot_courant_frames.py results/opencourant/RUN.radioss.json \
  --output results/opencourant/RUN_frames \
  --paper-pdf /path/to/user-supplied-paper.pdf
```

The default converter path is
`/workspace/opencourant/runtime/OpenCourant/exec/anim_to_vtk_linux64_gf` and can be
changed with `--converter`. The optional PDF produces a separate comparison of
Fig. 11(b)'s experiment, published FE and the current solver at requested strains
0.12, 0.35 and 0.65. Paper images remain in local output artifacts, not tracked
repository files. Their source hash and crop provenance are recorded.

The chosen frame must lie within 0.005 strain of its requested target. Every
panel reports its actual strain/time; an unavailable frame is shown as missing.
No eroded mesh is interpolated to create an exact target frame. The VTK converter
retains failed cells, so the plot explicitly requires `EROSION_STATUS == 1`.
Their orphan nodes also do not control the active-mesh viewing bounds. `POINTS`
already contains deformed positions: the displacement field is used only to
recover the original front plane, never added again. The JSON records original
animation hashes, original element IDs for invalid projected faces, erosion
counts and active-front overlap diagnostics. A front-plane overlap screen and
similar-looking images cannot establish complete 3-D contact validity or
experimental mode agreement.
