# Native Advanced Mass Scaling

OpenCourant's native `/AMS` and `/DT/AMS/2` controls are available through
`write_radioss_deck(..., ams_time_step=...)` and `run_radioss.py
--ams-time-step`. The option is explicit and off by default. The exporter
selects all parts, requests the supplied target step, uses a 0.67 step scale
factor, and sets the native PCG relative tolerance to 1e-4 with a maximum of
1000 iterations. It does not add a `/DT/NODA/CST` control.

For example, to **generate** a research candidate from the repository root:

```sh
source /workspace/opencourant/env.sh
/workspace/fea-venv/bin/python run_radioss.py \
  --config configs/dong2019_brick_damage.json \
  --output results/opencourant/thin_ams_example \
  --duration-seconds .01 --ams-time-step 1e-6 --drop-ccx-viscosity
```

Use a fresh prefix for each run. Append `--run --timeout-seconds 1200` only when
the corresponding solve is intended. Inspect the generated metadata and the
current qualification reports before interpreting its output. An accelerated
0.01-second explicit loading is not the paper's physical 2 mm/min test.

AMS modifies the high-frequency inertial operator. Total mass may stay unchanged
while the force history changes substantially. The current contact-coupon
comparison demonstrates this: increasing the native target step reduced the
cycle count, but changed the common-time platen impulse by approximately
20–48% relative to the unscaled run. Those coupon results do **not** qualify AMS
for fitting experimental data. Peak force alone would miss these differences.

The native INIEVO coupon still erodes its solid element with AMS enabled. That
is a component capability check, not a quasi-static response check. The
performance-only lattice pilot reached 1000 cycles at a 1e-6-second step in
63 seconds; its material table predates the constant-tail correction, so its
trajectory is excluded from physical comparison.

Before using an AMS lattice result for calibration, refine the target step,
repeat at slower loading, and compare full force/strain histories. Audit core
kinetic energy separately from the kinetic energy of prescribed moving steel
platens. Retain the global energy balance, contact energy, hourglass energy,
and added-mass checks. Neither subtracting platen motion nor ignoring startup
alone establishes a quasi-static solution. Mesh, contact, and experimental
holdout comparisons remain required.

Native artificial viscosity and damping defaults remain active. The starter
reports solid QA=1.1, QB=0.05, hourglass coefficient 0.1, and TYPE25 critical
damping 0.05 with its quadratic damping flag enabled. These are recorded in
metadata; omission of an additional Rayleigh damping card does not disable
the native defaults.

Exact observed outcomes and artifact hashes are in
`opencourant_ams_qualification.json`.

The completed thin-wall lattice was also rerun with the target step halved
from 1e-6 to 5e-7 s, preserving the exact mesh, material, contact, boundary and
loading curve. Both runs completed to strain approximately .65 with verified
execution bindings. Their full stress histories differ by 16.80% in relative
L2 and .1515 MPa RMS; integrated nominal work changes by 5.70%. Erosion changes
from 254 to 243 solids. The .1–.5 plateau mean changes by only 1.04%, hiding
substantial shifts in failure events. **AMS convergence is not established.**
See the [quantitative report](opencourant_ams_lattice_sensitivity.json),
[comparison plot](opencourant-thin-ams-sensitivity.png), and
[half-step assessment](opencourant-thin-half-ams-final.json).
