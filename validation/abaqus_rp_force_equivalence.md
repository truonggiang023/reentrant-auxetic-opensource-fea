# Abaqus RP reaction versus OpenCourant platen forces

**Purpose:** diagnose the early thin-wall lattice force overshoot before any
fracture-locus calibration. Abaqus RF2 at a displacement-driven, coupled
reference point is a constraint reaction; the existing OpenCourant TYPE25
history is a *contact impulse acting on the secondary/core*. They need not
match when the platen accelerates, other contacts/loads exist, or constraints
and loading differ.

## Test plan (no calibration, no material changes)

1. Freeze the *same* core mesh, platen geometry, thickness, contact surface
   definitions, friction, BCs, displacement vs time, material and area in both
   solvers. Record whether Abaqus uses kinematic or distributing coupling, and
   whether the RP is on a rigid or deformable platen. Do not call a different
   model an equivalent RP comparison.
2. In Abaqus export RP U2/RF2, platen mass/inertia and time; separately sum top
   and bottom reaction forces and, if available, platen contact forces.
3. In OpenCourant request T01 nodal REACY history for every *actually prescribed*
   top-platen node and for bottom supported nodes. REACY is cumulative impulse
   for the verified current runtime: sum the nodes first, differentiate once
   using actual T01 timestamps. Do not silently include core nodes or
   unconstrained directions. Capture top-platen velocity/momentum and all
   applied loads. Export both top-core interface (FN-FT, with action/reaction
   sign conversion) and bottom-core interface separately. A core-contact
   secondary force is not the force acting on the platen.
4. For each platen separately check, with signed Y components:
   sum(R_impulses) + sum(C_impulses) + other_impulses = momentum change.
   Use `reentrant_fea.force_balance.compare_platen_impulses` on a shared time
   grid. Never compare signs after taking absolute values; absolute value is
   appropriate only for the final magnitude stress plot.
5. Use the same reference area (currently 1900 mm² if Abaqus agrees), nominal
   strain defined from the actual top-bottom relative displacement / initial
   48 mm gauge length, and 0–0.12 strain grid. Report raw force, smooth
   interval-mean force (without hiding oscillations), balance residual,
   kinetic/internal energy and platen motion. Take care that T01 finite
   differences represent *interval-averaged* forces, whereas RP RF2 is
   an instantaneous sampled reaction.
6. Repeat with/without AMS at matched duration and with a slower loading
   history before attributing the discrepancy to material plasticity. If the
   platen reaction and contact forces agree (after inertia) but are both too
   high relative to Abaqus, inspect element formulation and boundary/contact
   implementations rather than rescaling E or the fracture locus.

## Acceptance / reporting

- A TYPE25 sign test has already found the current writer's contact force on
  secondary = FN-FT (impulse); contact force on platen is opposite.
- Existing lattice at nominal strain 0.12: corrected total contact nominal
  stress = 1.390 MPa; experimental lines ~0.281 and 0.259 MPa.
  This is NOT yet an RP-equivalent reaction measurement.
- Report maximum and RMS signed balance residual, interval timings, all
  source column names, and which platen/core each force acts on.
- A completed solver run and small residual establish bookkeeping at best;
  neither proves quasi-static conditions nor experimental agreement.
- The Python unit tests are *synthetic only*. No native OpenCourant/Abaqus
  execution or measured RP equivalence is claimed by this document.

### Local synthetic unit tests

```sh
python -m pytest -q tests/test_force_balance.py
```

### Minimal usage once signed cumulative impulses have been extracted

```python
from reentrant_fea.force_balance import compare_platen_impulses
# time_s, top_prescribed_reaction_Ns, on_top_platen_contact_Ns
# must be arrays on the SAME time grid. Oppose core secondary force.
result = compare_platen_impulses(
    time_s, top_prescribed_reaction_Ns, on_top_platen_contact_Ns,
    momentum=top_platen_momentum_Ns, area_mm2=1900.0)
print(result["balance_relative_rms"])
```

**Not implemented yet:** extraction of all lattice prescribed-node REACY
channels and core/platen sign mapping from the actual 6x5 T01 deck; native
solver execution and Abaqus RP export. Do not manufacture a numerical
comparison from existing total-contact-only data.

## Native RBE2 check (2026-10-30)

`platen_reference_points=true` now selects an RBE2 coupling from each outer
platen face to its central reference node. RP rotations are fixed, the bottom
RP is fixed and the top RP receives Y displacement. Inner faces remain elastic.
Both support and loading node histories are requested. Default direct-face
loading remains available as a numerical control.

Executed `scripts/probe_rp_coupling.py` with the pinned official runtime:
zero starter errors/warnings and normal engine termination for both cases.
At sampled strain 0.1199997, comparing interval forces obtained from REACY
impulses gives relative RMS difference 3.749666e-5 (0.00375%) and maximum
absolute difference 0.008866 N. See `native_rp_coupon_comparison.json`.
This is a small solid coupon, not the honeycomb. Its maximum KE/IE=0.798
fails the quasi-static screen. It establishes numerical agreement of the two
loading implementations only; it does not establish Abaqus agreement.

Reproduce from the repository root, with a fresh output directory:

```sh
source /workspace/opencourant/env.sh
PYTHONPATH=.:src /workspace/fea-venv/bin/python scripts/probe_rp_coupling.py
/workspace/fea-venv/bin/python scripts/compare_rp_reactions.py \
  results/native_rp/direct results/native_rp/rp \
  --output results/native_rp/comparison.json
```

The full honeycomb trial at duration 0.01 s, target strain 0.12, native AMS
1e-6 s, J2 plasticity and inner-Y platen contact failed with
`AMS IS LIKELY DIVERGING` after 1721 cycles. The trailing normal-termination
message does not qualify that run. No completed 0.12 lattice RP comparison
or Abaqus validation follows from it. A no-AMS pilot is required to isolate
AMS instability and estimate the full runtime before any material calibration.

The no-AMS full-mesh pilot was bounded to 90 s and timed out; its initial
stable step was 1.4515e-8 s. This is an incomplete run, not a convergence
failure or a completed comparison. A longer no-AMS RP/direct pair has been
launched locally under `results/native_rp_lattice/full_*`; results remain
pending and are not included in the committed qualification evidence.
