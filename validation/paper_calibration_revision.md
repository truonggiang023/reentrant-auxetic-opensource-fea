# Paper-following calibration revision

Calibration remains incomplete. No material parameters have been selected or
rescaled to match the supplied thin-wall image. The figure is Dong et al.
(2019), Fig.12(b), t=0.3 mm. Its two experimental traces are the thin holdout;
the paper's blue FE curve is not an experimental fitting target. Section 4.1
uses the thick t=1 mm specimen for fitting before transferring inputs unchanged.

## Current quantitative comparison

The [full comparison](paper_calibration_revision.png) and
[machine-readable diagnostics](paper_calibration_revision.json) use the
existing digitized reference CSVs and actual saved solver assessments.
At strain 0.12 the two experiments are 0.28123 and 0.25912 MPa; the corrected
native baseline is about 1.39007 MPa. The plateau means over the explicitly
chosen .1–.55 window are 0.30568/0.27116 MPa versus native 0.90367 MPa.
That averaging interval is an investigator diagnostic, not the paper's
reported plateau definition (paper scalar approximately 0.28 MPa).

The experimental efficiency-max diagnostic occurs at about 0.573/0.556;
the paper reports densification approximately 0.58. Missing initial red-trace
samples mean its energy integral starts at 0.02604. No absent initial energy
is fabricated. Native output ending near 0.65 cannot supply a diagnostic
requiring .1–.77 coverage. Morphology at .12/.35/.65 must still be checked
using actual frames, independently of the stress-curve score.

The optional CalculiX RP J2 overlay is a separate model and has different
mesh, platen and material-law choices. It is not a matched native/CCX solver
comparison or Abaqus validation.

## Executed and implemented revision

- Added pure curve diagnostics with explicit unavailable states for incomplete
  slope/plateau/target/densification coverage and separate experimental replicates.
- Added `--holdout-end-strain .78` to study preparation. Thick training remains
  .5 and historical force fit windows stay unchanged. The local revised native
  study contains 20 unrun hypotheses and remains geometry/physics blocked.
- Fixed CCX selection requiring actual verified execution receipts rather than
  cached completion/provenance fields alone.
- Restored legacy top-only native histories by default. Bottom histories are
  explicitly requested for support balance and automatically for RP coupling.
- Repaired the matched elastic probe to retain separate top/bottom channel
  groups after the new support-history option.
- Executed a fresh static C3D8I reference on the exact retained 7392-element,
  21870-node core mesh: stiffness 584.0973 N/mm, effective modulus 14.7561 MPa.
  This is a component diagnostic with no platen contact, plasticity or damage.
- Completed the matched native HEPH no-AMS control in 275567 cycles (744 s).
  Work-equivalent stiffness is 1006.2094 N/mm versus static C3D8I 584.0973 N/mm,
  approximately 72.26% higher. The native archived AMS result is 1006.2032 N/mm;
  removing AMS changes this work-based value by only about 0.00062%.
  Maximum late KE/IE is 3.87e-5; final internal/work discrepancy is about
  -4.57e-6. Nodal reaction differentiation passes the component balance screen.
  These are elastic component diagnostics; collapse/quasi-static material
  calibration is still unqualified. See [bound comparison](paper_revision_elastic_noams.json)
  and [stiffness plot](paper_revision_elastic_noams.png).

## What must pass before fitting

1. Qualify thick and thin geometry/density/mesh against the two specimens.
2. Resolve undamaged stiffness differences using matched elastic components;
   E=69000 MPa, nu=.3, density 2.67e-9 Mg/mm3 and the already true/plastic table
   remain fixed. Fracture tuning cannot repair an undamaged slope error.
3. Qualify full-lattice loading speed, core energies, global work/contact balance,
   mesh/contact and AMS sensitivity where AMS is used. Existing baseline
   KE/IE and approximately 14.97% AMS sensitivity fail the prospective screens.
4. Execute and bind completed thick training candidates, preserving fracture
   anchor .047 at triaxiality .333 and failure displacement .001 mm. Fit only
   missing locus-coordinate hypotheses; evaluate identifiability.
5. Transfer the selected inputs unchanged to the thin holdout, run to .78,
   and assess both experiments, late densification and real deformation frames.

The long no-AMS RP/direct honeycomb pair remains a separate pending control.
An RP trajectory reaching only 0.12 cannot complete full fracture calibration.
Native qualification blockers have not been bypassed. No bare `qualified=true`
flag or a low curve error is accepted as evidence of experimental agreement.

Reassess the executed elastic controls from the original local output paths:

```sh
MPLCONFIGDIR=/workspace/mplconfig python scripts/report_elastic_revision.py \
  --native results/revision_paper/elastic_noams_complete/elastic.radioss.json \
  --static results/revision_paper/shared/ccx_elastic_report.json \
  --output validation/paper_revision_elastic_noams
```

The report rechecks native execution/history binding, static output hashes
and exact mesh/BC/material equality. The archived AMS value is explicitly
labelled historical, not a new run. `paper_revision_elastic_inputs.tar.gz`
retains input decks, mesh/BC descriptions, execution receipts and source
scripts; large solver histories remain in local `results/revision_paper`.
The archive alone does not independently reproduce the reported outputs.
