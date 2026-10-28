# Fracture calibration and independent thin-wall validation

The paper does not publish a complete located fracture locus. Its Section 4.1
specifies fracture strain 0.047 at triaxiality 0.333 and failure displacement
0.001 mm. Other fracture strains 0.08–0.018 are estimates without their
triaxialities. The calibration tools preserve the located point and vary the
missing coordinates explicitly. They do not recover unavailable measurements.

Section 4.1 says the authors determined inputs using the thick-wall honeycomb
and reused them unchanged for the thin wall. Accordingly, **t=1 mm is training**
and **t=0.3 mm is an independent holdout**. Fitting the thin-wall curve and then
calling that same curve independent validation would overstate the evidence.

The existing hardening table is labelled user-supplied true stress/plastic
strain in `reference_setup_numeric_source.py`. It is not newly digitized coupon
data. Section 2.2 reports a tensile fracture strain near .047, while Section 4.1
uses .047 in its damage model; the plotted tensile strain must not silently be
treated as a measured equivalent plastic-strain history.

From the repository root, with the prepared Python environment:

```sh
/workspace/fea-venv/bin/python calibration_study.py prepare \
  --plastic validation/materials/user_plastic_table_2026-10-26.inp \
  --output /workspace/fea-research/new-calibration-study
```

Use an empty output directory. This creates a manifest and 24 configurations:
nine three-point-locus hypotheses and three numerical sensitivity variants,
each paired across the two wall thicknesses. Bounds are investigator choices,
not statistical confidence bounds. The central paper point and failure
displacement stay fixed; no generated configuration sets `calibrated=true`.

For native OpenCourant, use an explicitly selected duration and optional AMS:

```sh
source /workspace/opencourant/env.sh
/workspace/fea-venv/bin/python calibration_study.py prepare \
  --backend opencourant --duration .01 --ams-time-step 1e-6 \
  --plastic validation/materials/user_plastic_table_2026-10-26.inp \
  --output /workspace/fea-research/new-native-study
```

This prepares 22 paired configurations: nine locus hypotheses plus independent
duration and AMS-step sensitivities. It does not launch their solves. Native
variants do not pretend that CCX pseudo-time viscosity or residual stiffness
are transferable fit parameters. Each record includes native generate/run and
impulse-aware assessment commands. The native holdout force interval is
.03–.64 because interval-average force samples need not land at exactly .65;
the deformation target remains .65. Source config, geometry, hardening,
metadata, actual execution outputs and force-channel mapping are checked.
Cached completion flags or copied curves cannot qualify a native run.
The selected hardening file is validated before study files are written and
its path and SHA flow through both backend commands and scoring. Omitting
`--plastic` preserves the original higher-precision source; supplied data
already uses true stress/plastic strain and is not converted again.
The current native study remains selection-blocked pending actual core/work
balance, loading-speed, AMS and mesh/contact qualification.

The completed native baseline and its step-halving control now demonstrate
14.97% relative L2 difference in corrected total force in the stress history; see
[the AMS report](opencourant-thin-total-force-qualification.json).
The [physical audit](opencourant_physical_discrepancy_audit.json) also finds
early excess stiffness before damage, and a matched self-contact-off control
leaves the strain .005–.12 response unchanged. Inverse fracture fitting cannot
resolve that initial mismatch. Numerical formulation and discretization
qualification must precede selecting a material candidate. User-confirmed
elastic constants and the already converted true/plastic table are retained;
see [material provenance](materials/README.md).

Changing only `wall` in the inherited centerline geometry changes the specimen's
outer dimensions and relative density. These generated thick-wall configurations
are therefore hypotheses, not qualified thick-wall replicas. Check dimensions,
relative density, member arrangement and mesh against both paper specimens before
material fitting. The manifest records `geometry_qualification.qualified=false`
and blocks candidate selection until that prerequisite is resolved with evidence.

Each manifest record contains `solve_argv` and `assess_argv`, with absolute
configuration/output paths. Run them from the manifest's `working_directory`
using the prepared interpreter and solver environment. Preserve the `.inp`,
`.material.dat`, solver outputs and resulting assessment together. These are
potentially costly lattice simulations; generating a study does not execute
them. Numerical sensitivity variants vary viscosity or residual stiffness and
are excluded from material-candidate selection.

After completed thick-wall simulations and accepted-output assessments exist:

```sh
/workspace/fea-venv/bin/python calibration_study.py score \
  /workspace/fea-research/new-calibration-study/manifest.json \
  --output /workspace/fea-research/new-calibration-study/scores.json
```

The scorer uses both digitized experimental replicates, never the paper FE
curve as an experimental target. It samples at equal strain spacing over the
shared interval and reports RMS discrepancy, bias, experimental half-range and
coverage. It does not extrapolate missing solver strain. Partial runs receive
an overlap diagnostic but cannot win candidate selection even if their early
curve error is small. The thick-wall fit interval is strain .04–.45 because
Fig. 12(a) ends at .5 and one digitized replicate first becomes separable near
.03484. The thin holdout interval is .03–.65. Both reference coverage and an
actually completed solver run are required; new thick configurations end at .5,
while thin configurations end at .65. Legacy manifests with a single interval
use these role-specific defaults. Deck/material hashes, hardening, material
descriptors and the fracture locus are checked against each candidate. The
assessment also needs the matching `generation_config_sha256` recorded by the
runner after actually executing that configuration; an assessment lacking that
binding can provide a diagnostic but cannot win selection.

The best evaluated candidate is still not proof of a unique calibration or
experimental agreement. Transfer its parameters unchanged to the thin-wall
holdout; report the fit and holdout discrepancies separately. Assess local
identifiability from actual response sensitivities with
`reentrant_fea.calibration.sensitivity_rank`, and perform mesh, increment,
penalty and viscosity studies. Incomplete contact/deletion physics cannot be
compensated by parameter fitting. Deformation modes, accepted frames and energy
balance require their own checks.

Reference curves are image-digitized approximations. See
`paper_reference/dong2019_fig12_provenance.json` for extraction details and
uncertainty, and retain the paper/measurement provenance in any fitted report.

## Prospective acceptance checks

Before ranking material parameters, require normal completion, bound input and
output evidence, no unresolved initial contact warnings, and all requested
strain coverage. Review core kinetic/internal and hourglass/internal energy
ratios against 5% screening limits, added/original mass against 1%, and global
energy balance against 5%. These are investigator screening choices, not
measurements or acceptance limits stated by Dong. Erosion, prescribed platen
motion and native numerical damping must be included in the interpretation.
Passing a sampled energy screen is insufficient by itself.

Refine AMS target step and loading speed independently. Compare entire force
histories, not only peaks or conserved total mass; use a 5% normalized RMS
response change as a preliminary convergence target and report the actual
values. Perform a separate mesh/contact sensitivity check. An unqualified
AMS run remains a research trajectory even if its curve resembles the paper.

For an experimental-agreement assessment, report discrepancy to both tests,
their mutual discrepancy, digitization uncertainty, plateau stress,
densification strain and the available deformation modes. A preliminary
assessment may target plateau error within 10% and densification strain within
.03, but those are investigator choices and do not establish agreement alone.
The Bi-V/Z/collapse evidence needs actual stored frames with their real strain
labels; curve fitting cannot substitute for the morphology comparison.
Parameter identifiability and estimated-locus uncertainty must remain separate
from an observed curve match. No automated candidate ranking in this workflow
sets `calibrated=true` or asserts experimental reproduction.

## Current open-source calibration correction

All new native studies and resumed legacy studies reconstruct the total
secondary contact force as FN−FT before differentiating the impulse channels.
The material table and E remain fixed. User observations from a separate
Abaqus model are context only; this workflow does not require Abaqus.

To investigate the tested platen surface control in a new native study:

```sh
python calibration_study.py prepare --help
```

Use `configs/dong2019_inner_y_contact.json` as the base config via the prepare
command's config argument. The inner-Y control is not yet qualified for full
collapse; selection remains blocked pending numerical qualification. Do not
interpret a lower fitting error from an unqualified trajectory as calibration.
