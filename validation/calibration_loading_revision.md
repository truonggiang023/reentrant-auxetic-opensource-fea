# Loading and mesh correction after the supplied paper paragraph

The paper paragraph does not identify Standard/Explicit, time integration,
mass scaling or an acceleration protocol. C3D8R alone does not settle that.
It specifies elastic steel plates, self-contact, constrained honeycomb U3,
core/platen friction 0.35, indenter speed 2mm/min and the stated damage anchor /
displacement evolution. Core self-friction is not specified in that paragraph.
See [the recorded conditions](paper_loading_conditions.json).

## New completed controls

The [bound summary](calibration_loading_revision.json) separates the observations:

- Full original-mesh RP control completed strain 0.119999708 in 0.01s, J2 only,
  core-only AMS 1e-6. Last interval RP stress 1.185053MPa and same-run contact
  stress 1.184505MPa. Both still exceed the paper's two samples near 0.12,
  approximately 0.281/0.259MPa. A change to RP extraction does not remove the gap.
- This run fails the energy screens: maximum core KE/IE 0.352 and global 2.223.
  The signed whole-platen momentum screen also remains unqualified (normalized
  residual 0.733). No ad hoc momentum scaling was applied to make it pass.
- The wall-aligned structured mesh preserves the polygon area/perimeter and
  has 10544 core bricks. Its native elastic work stiffness 590.1698N/mm agrees
  with the same-mesh static C3D8I 585.1354N/mm within 0.86%. Original-mesh error
  was 72.27%. This is component work/energy agreement; raw nodal derivative
  balance remains flagged false, and nonlinear/fracture/contact convergence
  is not established. The mesh does not match paper edge seeds 19/10 exactly.
- Depth scaling 25→0.3mm changes normalized native short-test stiffness only
  -0.0178%; high extrusion aspect alone did not explain the excess.
- Reduced-pressure Icpre1 gives 733.7875N/mm in a slow original-mesh control,
  but softens an otherwise accurate regular beam by 13.45%. It was not adopted
  as a curve-fitting correction. The default pressure treatment is preserved.

The structured strategy is optional (`brick_mesh=structured_walls`, mesh size
0.8mm), with geometry, conformity and positive-Jacobian tests. The library
reproduced the executed diagnostic node coordinates/connectivity exactly.
The existing default mesh has not been silently replaced in old configurations.

## Physical velocity and a controlled rate test

At actual core height 48mm, the physical constant nominal strain rate is
(2/60)/48=6.94444e-4 s^-1. Constant-speed durations are 172.8s at 0.12,
936s at 0.65 and 1123.2s at 0.78. A smoothstep amplitude is not constant velocity:
its peak velocity is 1.875 times its average.

`run_radioss.py --physical-speed-mm-min 2` now derives time from actual mesh
height and uses linear displacement. The physical deck at
`results/paper_conditions/velocity 012` was generated and checked, not executed.
No mass scaling is requested in that generation. A literal physical-time
Explicit analysis can be exceptionally costly; generation is not validation.

The old long no-AMS accelerated control was stopped with partial evidence
preserved. A matched original-mesh/RP/J2/core-AMS rate pair now compares 0.01s
and 0.04s, changing only duration and output/checkpoint settings. Its slow run
is active. `compare_loading_rates.py` automatically assesses the bound result
and writes the report/plot only after actual completion. Live state:
`results/calibration_completion/rate_comparison.status.json`.
The current material law is rate independent; this experiment measures dynamic
loading sensitivity, not a fitted constitutive strain-rate law.

## Abaqus variant

[run_paper_velocity012.py](../abaqus/run_paper_velocity012.py) keeps the user's
material table, selects physical duration, linear amplitude, all-node core U3,
RP histories and no requested mass scaling. Explicit remains an implementation
choice. See [instructions](../abaqus/README_PAPER_VELOCITY.md).
The script was syntax/contract checked but not executed in licensed Abaqus.
The Abaqus and native material/mesh/platen assumptions are still different;
this file does not establish a matched cross-solver reproduction.

Calibration is not declared complete from these controls. Loading-speed and
nonlinear mesh/contact/erosion validation precede fracture-parameter fitting.

The retained [input/output archive](calibration_loading_inputs.tar.gz) includes
the completed RP case's raw T01/CSV, solver outputs and receipt; the structured
elastic metadata/input/report and static reference; and the pressure-control
metadata/input/report. Large elastic T01/CSV and animation files remain in
local results and are not included, so that archive alone does not reverify
all component output hashes. Original executed metadata is preserved.

`structured_walls` can now be selected in the normal generators. The native
and static controls were reproduced with the same node coordinates and
connectivity. The physical-speed structured example is
`configs/dong2019_structured_velocity012.json`; it was generated with the
physical-speed flag and gives 172.8 s at 0.12, but that literal-time explicit
analysis has not been run. It is not labelled calibrated.
