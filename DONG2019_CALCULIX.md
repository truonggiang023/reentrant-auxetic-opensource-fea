# Dong et al. (2019): CalculiX replication candidate

**Research status: not a validated reproduction.** The reference is Dong et al.,
Materials & Design 182 (2019), 108036. The uploaded paper is evidence for model
conditions, not an instruction to claim agreement. Current solver outcomes are
recorded in `validation/native_qualification.json`; only valid, completed
increments can support a deformation claim.

## Reference conditions and missing calibration

Section 4, p. 8 gives C3D8R elements (10,068 total), elastic steel platens,
friction 0.35, out-of-plane constraints, self-contact, and compression at
2 mm/min. It specifies fracture plastic strain 0.047 at triaxiality 0.333,
other estimated fracture strains in the range 0.08–0.018, and linear
**displacement** softening with failure displacement 0.001 mm. It does not
publish the triaxialities corresponding to those other values.

The demonstration locus assigns 0.08 at eta=0 and 0.018 at eta=1, with the
published central point retained. Endpoint clamping, a compression cutoff at
eta=-1/3, residual stiffness 1e-4, and viscosity time 0.005 are explicit
assumptions. They are not a recovered or fitted paper calibration. The hardening
table is read from the repository's numeric reference source without executing
Abaqus Python. No claim is made that it is a newly digitized experimental curve.

Target frames are Bi-V at nominal strain 0.12, Z at 0.35, and severe collapse at
0.65. Experimental mean thin-wall plateau stress is approximately 0.28 MPa and
mean densification strain 0.58. These reference quantities are not acceptance
criteria that can be satisfied by tuning a numerical parameter without a
sensitivity study.

## Prepared workflow

From `/workspace/reentrant-auxetic-opensource-fea`:

```sh
export PATH=/workspace/bin:/workspace/fea-venv/bin:$PATH
OMP_NUM_THREADS=2 CCX_EXECUTABLE=/workspace/bin/ccx \
  CCX_DAMAGE_EXECUTABLE=/workspace/bin/ccx-damage pytest -q -rx
python run_paper.py --config configs/dong2019_brick_damage.json \
  --self-contact --target-output-only --output results/native.inp --run --threads 2
MPLCONFIGDIR=/workspace/mplconfig python plot_run.py results/native.inp
```

Rebuild the qualified solver/library with `scripts/build_ccx_damage.sh`. Keep
`native.inp` and `native.material.dat` together. The runner selects the compiled
backend and supplies the sidecar path. The source/material hashes are recorded
in `/workspace/ccx-user/build-sha256.txt`.

`--imperfection-mm .003` applies a documented 1%-of-wall-thickness X perturbation
for mode-sensitivity studies. It is not a measured manufacturing imperfection.
Contact partitioning uses the unperturbed reference topology, so a smooth small
perturbation does not explode the number of surface patches. Actual solver node
coordinates include the perturbation.

`--target-output-only` requests exact target strain frames. With the example's
`output_interval=0.005`, it additionally records a regular time grid for reaction
histories. Native damage is integrated internally at every solver increment;
output sampling does not replace that integration. The custom build also fixes
a CCX 2.20 retry bookkeeping defect that can schedule a zero increment after
accepting an output time; see the native README. Offline initiation diagnostics
from sparse PEEQ histories have additional sampling error.

## Current solver evidence

The final test suite records **41 passed and 3 known mortar XFAILs**. An accepted
full-lattice frame at nominal strain 0.12 has no detected front-plane overlap or
inverted projected elements. Its nominal stress is about 0.660 MPa (reaction
1254 N / 1900 mm2); this is solver evidence, not an experimental agreement claim.
The plotted cells remain mostly regular, and Bi-V agreement is not established.

A 0.003 mm imperfection comparison reached accepted strain 0.14518 with maximum
recorded integration-point damage 0.49995. It was deliberately stopped after
repeated contact cutbacks to prioritize the corrected backend; this is not a
completed run or a solver-declared failure. Its accepted states remain usable.

The original zero-imperfection run hit a zero-step error after an output-point
retry at time 0.185. The patched full-lattice rerun accepted that retry and
continued with positive increments beyond it. See the qualification JSON for
its final available strain and the status of the 0.35/0.65 targets. The corrected
run was stopped at a 20-minute qualification runtime limit after repeated
cutbacks, with accepted nominal strain 0.14422. This is a deliberately stopped
run, not a solver-declared failure. Neither lattice reached 0.35 or 0.65. These results
do not establish the complete deformation sequence, mesh convergence, or a
calibrated constitutive response.

## Material and element implementation

The solver-side C material uses multiplicative finite-strain J2 plasticity,
isotropic Hencky elasticity, triaxiality-dependent initiation, and crack-band
linear displacement softening. It returns PK2 stress and a symmetric numerical
material tangent to CalculiX. Trial state starts from `xstateini`, so retries and
line searches do not accumulate plasticity or damage. See [native/README.md](native/README.md).

The safe interface uses nine descriptor properties and a read-only full-table
sidecar. The complete 178-point hardening table is covered by an actual solver
regression. **Oversized inline tables are unsafe in CalculiX 2.20**, whose local
buffers hold 21 properties. Earlier development lattice trials that used such
tables are invalid and are excluded from evidence.

SDV1–6 hold the inverse plastic metric, SDV7 PEEQ, SDV8 initiation integral,
SDV9 post-initiation displacement, and SDV10 damage. Integration-point SDVs are
written to DAT. FRD nodal extrapolations must not drive failure decisions.

C3D8R generation uses blossom quad recombination and conforming refinement,
including splitting the small number of residual triangles. It gives two cells
through initially one-quad wall widths and fewer skew quads than splitting every
triangle into three. `brick_mesh=triangle_subdivision` retains the earlier
comparison mesh. These meshes differ from the paper's topology and discretization;
through-width, depth, hourglass and global mesh sensitivity remain required.

The scoped compatibility patch uses E*(1-D) consistently for `@REENTRANT`
hourglass force and stiffness paths. It avoids unset stress-only tangent arrays
and initializes the stiffness modulus from the actual material data. It changes
no other material. Hourglass loading and a two-brick free-DOF bending test match
the builtin elastic solver reference. That is an interface check, not proof of
continuum bending accuracy or negligible artificial hourglass energy.

## Contact and boundary conditions

The platens are separate C3D8 steel solids. Only the lower outer support face and
upper outer loading face are constrained. Their contact faces remain free to
deform elastically. Setting `rigid_platen_kinematics=true` recovers the earlier
all-node translation comparison model. A single brick per platen is still a
coarse approximation of its bending/contact response.

All core U3 DOFs are fixed by default, including interior depth layers. A single
bottom centreline X roller removes initial rigid translation before friction
carries load. Its influence on mode selection needs sensitivity checks.

Penalty contact uses positive normal and stick stiffnesses. The revised example
uses 1e4 MPa/mm normal stiffness and 1e3 MPa/mm stick stiffness, with mu=0.35.
Check actual penetration and stiffness sensitivity; convergence alone does not
qualify a penalty choice. Zero friction omits the `*FRICTION` card.

Optional wall contact creates disjoint slave/master patches. It excludes
shared-node neighbors, paths entirely through initial solid, and initially
back-facing pairs. Each unordered pair occurs once. A closure benchmark compares
runs with and without wall contact, checks pre-contact reactions, and measures
reported penetration. This algorithm can miss contacts made possible only by
large rotations of initially back-facing patches. It does not create fracture
surfaces when material fails.

## Verification and scientific limits

Compiled tests cover rigid-rotation objectivity, plastic-volume preservation,
repeated trials, elastic unloading, compression cutoff, viscosity relaxation,
softening energy for three band lengths, full-table solver-side tensile load
loss, hourglass loading and free bending. Contact tests compress a small specimen
to 0.12, 0.35 and 0.65 and exercise wall-gap closure. These are component checks,
not honeycomb deformation-mode validation.

Three coarse frictional mortar cases remain known numerical XFAILs in stock CCX
2.20. They are marked only after the log confirms cutback exhaustion/minimum
increment failure; other errors fail the tests. MORTAR, LINMORTAR and Petrov
variants were investigated, and penalty contact remains the selected path.

The runner requires clean completion and an accepted final status time, and
checks the native material-load signature. `frd.committed_fields` rejects the
unconverged “best solution” that CalculiX can write after failure.
`frd.projected_overlap` detects overlapping front-plane element projections;
it is a diagnostic for missed contact, not a 3-D contact proof.

Important remaining limits:

- Failed material retains residual stiffness and its contact faces. There is no
  automatic deletion or new crack topology; ghost walls can affect late collapse.
- The fracture locus, regularization, imperfections and full-mesh response are
  uncalibrated. Static pseudo-time does not establish equivalence to 2 mm/min or
  to an explicit dynamic calculation.
- The buffered-strut envelope is about 76.50 x 48 x 25 mm, not exactly the nominal
  76 x 48 x 25 mm. Reaction/1900 mm2 is the paper's nominal-stress convention;
  it must be distinguished from stress based on the generated envelope.
- A converged image does not establish Bi-V or Z. Force curves, contact overlap,
  material states, mesh and numerical-parameter sensitivities must also agree
  before any experimental-reproduction claim.
