# OpenCourant native fracture/contact workflow

OpenCourant is the community continuation of OpenRadioss. The official
`latest-20261011` Linux archive is pinned and checked against its published
SHA256 by `scripts/setup_opencourant.sh`. Installation and the official smoke
test were executed, not just described. See `opencourant_runtime.json`.

This backend uses native tabulated plasticity, INIEVO initiation and linear
displacement softening, and TYPE25 contact. It is a different material/element
formulation from the custom CalculiX UMAT. The complete fracture locus and
experimental agreement remain unestablished.

## Verified components

- Official TWISBEAM: both executables ran, with ten reference checks passing.
- Constrained tensile coupon: native solid erosion, ERODED changing from zero
  to one, and load loss after rupture. This is not a uniaxial-stress calibration.
- Sliding steel-platen compression coupon: opposing transmitted loads and
  completion to prescribed nominal strain .12, with no starter errors/warnings.
- Official TYPE25 all-face erosion benchmark: 24 deleted solids while contact
  operates. Its large global energy error excludes it as calibration evidence;
  see `opencourant_contact_qualification.json` for exact scope and differences.
- A controlled initially internal-face test: after the covering solid erodes,
  ALL transfers about 6110 N while the exterior-only control transfers zero.
  No contact remains at the removed solid's old position. See
  `opencourant_fresh_face_qualification.json`.
- A controlled 180-degree rotation test: native self-contact detects the
  initially back-facing face near 160 degrees; the no-contact control remains
  unloaded. See `opencourant_rotation_contact_qualification.json`.
- Four free-lateral/constrained material probes verify true flow-stress input,
  the located tensile fracture anchor, and the declared fracture-locus endpoint
  convention at high and negative triaxiality. See
  `opencourant_uniaxial_qualification.json`.
- Thin lattice: 21,870 core nodes and 7,392 HEPH bricks start with zero errors
  and warnings. Core self-contact retains edge contact; broad steel platens
  use face contact to avoid duplicate coplanar edge intersections.

These targeted tests establish fresh-face activation and detection after large
rotation for their specified geometries. They do not establish correct contact
throughout an arbitrary collapsing lattice or a calibrated material model.

Native FUNCT and TABLE initially extrapolated their endpoint slopes, unlike
the CCX boundary conventions. The exporter now adds constant numerical
segments for both hardening and fracture-locus endpoints. They preserve all
original input points and are explicitly labelled as numerical conventions,
not additional measured material data. Runs from before these fixes are
excluded from material-law comparison. Native compression damage still differs
from the CCX example: INIEVO accumulates damage at negative triaxiality, whereas
the CCX configuration has a compression cutoff. That difference is explicit;
the backend is not represented as an exact UMAT translation.

## Commands

From the checkout:

```sh
export PATH=/workspace/bin:/workspace/fea-venv/bin:$PATH
source /workspace/opencourant/env.sh
export CCX_EXECUTABLE=/workspace/bin/ccx
export CCX_DAMAGE_EXECUTABLE=/workspace/bin/ccx-damage
OMP_NUM_THREADS=2 python -m pytest -q -rx
```

Reinstall or refresh the official runtime with:

```sh
PATH=/workspace/fea-venv/bin:$PATH bash scripts/setup_opencourant.sh
```

Generate a model with an explicitly selected duration:

```sh
python run_radioss.py --config configs/dong2019_brick_damage.json \
  --output results/opencourant/new_case --duration-seconds .01 \
  --drop-ccx-viscosity
```

Appending `--run --timeout-seconds 600` runs starter, engine and the bundled
T01-to-CSV converter. A fresh prefix is required to preserve previous results.
Each solver stage is bounded; incomplete outputs remain for diagnosis. A .01 s
lattice pilot estimated about six hours on the present four-thread allocation;
it was deliberately stopped after measuring cost. The duration is accelerated
relative to physical 2 mm/min and needs an energy and speed sensitivity study.

CCX's pseudo-time viscosity cannot be translated to physical seconds. The
explicit omission flag records its removal; it does not assert equivalence.
Native characteristic length and HEPH stabilization also differ from CCX.
No mass scaling is inserted by default. Native bulk viscosity, hourglass and
contact damping defaults remain active and must be audited; zero-valued input
fields must not be mistaken for disabled damping.

Native AMS can be requested explicitly with `--ams-time-step 1e-6`; it is off
by default and modifies the inertial operator even if total mass is unchanged.
Read `opencourant_ams.md` for actual speed and sensitivity results. Those tests
do not yet qualify AMS for parameter fitting.

## Assessing actual output

The bundled converter reads binary T01:

```sh
th_to_csv_linux64_gf new_caseT01
```

Use `/TFILE/4`; the legacy ASCII path was observed to repeat its header without
advancing, and is not used. A successful process or normal termination alone
does not qualify quasi-static response.

Interface FY and nodal REACY are accumulated impulses in this runtime, not
instantaneous forces. Convert differences to interval-average forces using
the actual time intervals. The verified top-platen channel for the current
three-part writer is:

```sh
python assess_courant.py results/opencourant/new_case.radioss.json \
  results/opencourant/new_caseT01.csv \
  --force-column 'Contact force and energy                 3 Contact 3                                var 64' \
  --force-sign -1 --force-kind impulse \
  --output results/opencourant/new_case.assessment.json
```

The assessor checks input and output hashes, uses the actual prescribed
displacement curve to map strain, and does not extrapolate. Runs without a
matching execution manifest remain diagnostic only. Inspect core-part energies
as well as whole-model energies: driven platen inertia can dominate global KE.
Dense energies, energy balance, hourglass/contact work, mesh and loading-speed
sensitivity are still needed even if a sampled energy screen passes.
New execution records also bind the metadata hash before solver launch. The
historical completed baseline predates this extra binding; that limitation is
reported explicitly rather than inventing an execution-time metadata hash.

## Current lattice result

The free-X baseline uses the current boundary default and completed normally:
10,002 cycles, 934.93 seconds on four threads, 254 eroded core elements and no
engine error messages. Last stored force strain is .6499977; actual prescribed
motion at the last history sample is .6499985. Execution-time hashes bind both
decks and metadata to the native output. Its evidence is committed:

- [Force/stress history and assessment](opencourant-thin-free-x-final.json)
- [Paper stress comparison](opencourant-thin-free-x-comparison.png)
- [Actual, erosion-filtered deformation frames](opencourant-thin-free-x-frames.png)
- [Quantitative qualification](opencourant_thin_free_x_qualification.json)
- [Exact generated input archive](opencourant-thin-free-x-inputs.tar.gz)

Against both digitized experiments over strain .03–.64, RMS discrepancy is
.6074 MPa with +.5858 MPa bias. At strain .12 it gives 1.1848 MPa versus
.2812/.2591 MPa experimentally. Removing the interior X pin changes the curve
by about 12% in relative L2 over the same interval and does not resolve this
overprediction. The maximum core KE/IE after the first 5% of loading time is
.609; over strain .12–.35 it is .0156, and over .35–.65 it is .0080. These
selected windows do not establish loading-speed independence, work balance,
or AMS/mesh convergence. This is **not a calibrated reproduction**.

A completed rerun halving the AMS target to 5e-7 s changes the full force
history by 16.80% in relative L2 and integrated nominal work by 5.70%; erosion
changes from 254 to 243 elements. The [AMS sensitivity report](opencourant_ams_lattice_sensitivity.json)
and [plot](opencourant-thin-ams-sensitivity.png) demonstrate that timestep
convergence remains unresolved. Matching a plateau mean alone is insufficient.

Reproduce this baseline with current defaults and a fresh prefix:

```sh
OMP_NUM_THREADS=4 python run_radioss.py \
  --config configs/dong2019_brick_damage.json \
  --output results/opencourant/free_x_replay --duration-seconds .01 \
  --ams-time-step 1e-6 --drop-ccx-viscosity --run --timeout-seconds 1800
```

The archive preserves exact generated inputs, avoiding changes from future
mesher or exporter versions. It contains project-generated decks, not solver
binaries or images from the paper.

## Diagnosing the early discrepancy

Actual early native frames at strains .00412 and .02068 have zero damage and
erosion. The [physical audit](opencourant_physical_discrepancy_audit.json)
finds no dimensional or reference-area error large enough to explain the
stress excess. E=69000 MPa, nu=.3, density=2.67e-9 Mg/mm3 and the already
converted true-stress/plastic-strain interpretation are
[explicit user-confirmed inputs](materials/user_material_confirmation.json).
They are retained; printed coupon-raster slopes do not override them.

A matched control removes only core self-contact and deliberately stops at
.0032 s of the unchanged .01 s loading amplitude. Over strain .005–.12 its
stress and core internal-energy samples are identical to the baseline at CSV
precision. Core self-contact therefore does not explain the early overshoot
in this case. This intentionally partial diagnostic has generation hashes
but lacks a prospective execution input/metadata binding, and cannot be used
as a completed fit. Recheck the sealed evidence and report with:

```sh
python scripts/audit_self_contact.py \
  --verify-report validation/contact_audit/thin_self_contactoff_partial.json
```

The [plastic-bending diagnostic](opencourant_plastic_bending_audit.json)
also shows sensitivity to solid formulation and wall discretization. Both
formulations agree in the elastic range and differ after yielding. Its
prescribed tip face introduces a substantial finite-rotation couple, so the
ideal free-tip beam plastic moment is only a reference. The finer cases also
use a different loading duration. These controls do not establish a solver
defect, clean mesh convergence, or the cause of the lattice discrepancy. They
show why elastic coupon checks alone cannot qualify plastic wall bending.

## Earlier pinned baseline

The corrected-law, legacy interior-X-pin baseline completed to strain .65:
10,002 cycles, 976.37 seconds on four threads, 252 eroded core elements, no
engine error messages and normal termination. Its last stored force point is
strain .6499977. Data and figures are committed:

- [Force/stress history and assessment](opencourant-thin-pinned-final.json)
- [Paper stress comparison](opencourant-thin-pinned-comparison.png)
- [Actual, erosion-filtered deformation frames](opencourant-thin-pinned-frames.png)
- [Quantitative qualification](opencourant_thin_qualification.json)

The strain .03–.64 curve discrepancy is about .573 MPa RMS with +.551 MPa bias.
The global printed energy error reaches 5.2%; core energy transients and AMS,
loading-speed and mesh sensitivities remain unqualified. Morphology at .35
differs from the paper. This completed solve is **not a calibrated reproduction**.

The boundary audit found the old X roller at an interior first-row node. The
paper specifies out-of-plane restraint and platen friction, not this interior
X pin. Explicit dynamics does not require a static stiffness-matrix gauge.
The current exporter therefore leaves core X motion free by default; set
`core_x_anchor=true` only to reproduce the legacy boundary assumption. See
`opencourant_boundary_audit.json`. This correction does not itself establish
the cause of the stress discrepancy.

Exact legacy inputs are preserved in `opencourant-thin-pinned-inputs.tar.gz`.
Extract into a fresh directory and run its `_0000.rad` then `_0001.rad` with
the pinned runtime to replay that baseline. Regenerating with current default
boundary settings creates a different, unpinned model and must be labelled
separately. The input archive contains project-generated decks, not solver
binaries or images from the uploaded paper.

The independent bending audit is reproducible with
`scripts/probe_courant_bending.py`. Stable controls agree with beam theory to
about .63%. A high-aspect, one-depth-layer case fails at timestep scale .9 but
passes at .3. This confirms the need for timestep/discretization qualification;
it does not prove that HEPH bending stiffness causes the lattice mismatch.

## Excluded early pilot

The `.0001 s` dynamic integration smoke was stopped by its 360 s timeout at
observed time `4.375822e-5 s`, approximately strain `.250`. Recorded core
ERODED reached 1,175; later printed rupture messages include additional events.
The maximum core kinetic/internal-energy ratio after startup was approximately
6.9. It is strongly dynamic and is explicitly rejected as quasi-static
calibration or experimental-agreement evidence. It must not be joined to the
paper curve or extrapolated to .65.

This early integration smoke also predates the corrected hardening and
fracture-table endpoint conventions and is excluded from material comparison.
Exact input hashes and rejection reasons are retained in
`opencourant_lattice_pilots.json`.

Paper comparison data and the accepted historical CCX comparison are described
in `paper_reference/README.md`. Thick-wall fitting and independent thin-wall
validation must remain separate; see `calibration_workflow.md`.

## Exact physical velocity versus accelerated diagnostic

The user-provided paper paragraph does not identify Standard/Explicit or the
step formulation. C3D8R alone does not identify the solver. Physical constant
2 mm/min is distinct from the accelerated smoothstep durations used above.

Generate the declared physical loading condition (without launching billions
of explicit cycles) with:

```sh
python run_radioss.py --config configs/dong2019_paper_velocity012.json \
  --physical-speed-mm-min 2 --drop-ccx-viscosity \
  --output results/paper_velocity/new_case
```

Actual mesh height 48 mm gives T=172.8 s at strain .12; the selected linear
amplitude has both mean and peak velocity 2 mm/min. Smoothstep has peak
velocity 1.875 times its average. Neither an accelerated run nor switching
its output time units establishes physical-rate or quasi-static equivalence.
The current hardening law contains one curve and no constitutive rate effect.
See `paper_loading_conditions.json` for specified and unreported quantities.

For supported RP bookkeeping, use generated `platen_reference_points=true` and
assess the completed, bound run with `assess_courant.py META CSV --force-source rp`.
This uses observed relative RP DY and differentiates cumulative REACY once;
nominal stress is minus the signed top reaction divided by area. It is not an
independent proof that native REACY and Abaqus RF2 have identical inertial terms.
Core-only AMS (`ams_scope=core`) is an explicitly labelled numerical control
that avoids including steel/platen constraints in the modified mass operator;
no quasi-static claim follows from successful completion. Optional
`restart_cycle_interval` writes two rotating restart files; creation was tested,
continuation has not been qualified.
