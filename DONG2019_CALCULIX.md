# Dong et al. (2019): CalculiX development candidate

**Incomplete reproduction.** Python creates meshes/input decks and reads results;
CalculiX performs the finite element solve. The tested Linux executable is Debian
CalculiX 2.20-1, using SPOOLES. The bundled Windows 2.13 executable has not been
validated with these changes. No experimental agreement is claimed.

## Run on the prepared cloud machine

From `/workspace/reentrant-auxetic-opensource-fea`:

```sh
export PATH=/workspace/bin:/workspace/fea-venv/bin:$PATH
OMP_NUM_THREADS=2 CCX_EXECUTABLE=/workspace/bin/ccx pytest -q -rx
python run_paper.py --run --threads 2
```

The second command attempts the full lattice; convergence is **not established**.
It writes a log and requires both a clean solver finish and a completed `.sta`
step at pseudo-time 1. Old outputs are retained in a unique `*-previous-*`
directory. A zero process exit status alone is insufficient.

Default output contains every converged increment, including unaveraged CORE
integration-point stress and PEEQ in `.dat`. `--target-output-only` instead asks
CalculiX for exact nominal strain frames 0.12, 0.35 and 0.65 (if in the loading
range); sparse output is unsuitable for accurate history integration. Static
time is a load fraction, not experimental seconds. Nominal strain is total
prescribed displacement / initial height (48 mm here); platen and contact
compliance mean it is not exactly the core's material strain.

## Changes and verification

- Validate positive elasticity, admissible Poisson ratio, finite monotone true
  stress/plastic strain tables starting at zero plastic strain, compression
  range, friction, contact stiffnesses and increment ordering. Stress softening
  in a `*PLASTIC` table is rejected as a substitute for ductile damage.
- Use explicit `ROT NODE` nodes with translational DOFs 1–3 for rigid-grip
  rotations in the advanced writer. Reference-node DOFs 4–6 are not used as
  CalculiX rigid-body rotation constraints.
- Keep separate meshed steel platens, positive normal penalty and stick slope,
  friction 0.35 and out-of-plane constraints. Add one bottom centreline X roller
  to remove the initial in-plane rigid translation before friction carries load.
  This explicit numerical anchor must be included in boundary-condition
  sensitivity studies; it is not a measured experimental support.
- Use the actual CCX 2.20 minimum increment floor of 1e-6 of step time, rather
  than requesting 1e-10 which CCX silently increases with a warning.
- Expose penalty (`SURFACE TO SURFACE`), `MORTAR`, and `LINMORTAR` contact choices
  for controlled numerical studies. Penalty remains the default. Mortar is not
  a demonstrated cure for nonlinear convergence.
- Reject same-surface self-contact pairs. Reusing one exterior surface as both
  master and slave is not a benchmarked wall self-contact algorithm. Crushing
  requires disjoint nonadjacent wall patches, appropriate search and exclusions,
  and contact-force/penetration benchmarks before enabling it.

Actual solver regressions compress a **two-wedge rectangular specimen**, with
frictional steel platens, through nominal strains 0.12, 0.35 and 0.65. Penalty
contact completes all three and transmits compressive reaction. These are
contact/execution tests, not honeycomb morphology tests. Mortar fails the same
coarse specimen after cutbacks; the tests report these three known numerical
failures as XFAIL only when the log confirms cutback exhaustion or the minimum
increment limit. Any other failure fails the test.

Full lattice attempts use 10,026 core nodes and 6,880 C3D6 elements. The original
penalty deck exhibits enormous rigid translation corrections (~1e9 mm). The
anchored penalty deck still exhausts cutbacks in its first increment. Exploratory
mortar variants make limited progress, but were stopped because they advance
very slowly near nominal strain 0.021; they did not reach any target frame.
See `validation/ccx220_summary.json` for captured solver evidence. No lattice
Bi-V at 0.12, Z at 0.35, or collapse at 0.65 is verified.

## Ductile initiation and unresolved coupled damage

`reentrant_fea.ccx_results.integration_history(dat)` reads unaveraged integration
point stress and PEEQ from CCX 2.20 `*EL PRINT`. It computes triaxiality as mean
stress / von Mises stress, with positive tensile hydrostatic stress. Separate
integration points are never merged into an element-average history.

`initiation_indices(dat, locus)` integrates `dPEEQ / fracture_strain(eta)` on each
integration point. It refuses a single-point locus or stress triaxialities outside
the calibration range; it never assumes the tensile datum 0.047 at eta=0.333 is
valid in compressive or shear states. Synthetic test loci verify plumbing only.
Output sampling and interpolation introduce integration error; refine output
and loading increments for sensitivity studies.

This is an **uncoupled diagnostic**, not a damage constitutive law. Stock CCX
2.20 `*PLASTIC` does not supply the paper's Abaqus ductile initiation plus
0.001 mm displacement softening and automatic deletion. `--ductile-damage` in
the advanced CLI remains an explicit unsupported-feature error. No Python FE
solver, externally imposed element deletion, or unsupported `*DAMAGE` keywords
are substituted for a solver-side material law.

A physically appropriate coupled implementation remains required: a compiled,
finite-strain CalculiX user material with committed integration-point state,
calibrated tension/shear/compression fracture locus, characteristic-length or
fracture-energy regularization, stress degradation and verified tangent. It
needs unloading, rotating-body, mesh-objectivity and energy benchmarks, plus a
supported treatment of failed elements and contact. The repository's standalone
small-strain Python kernel cannot be attached to CCX as that implementation.
The missing calibration and this solver extension remain blockers to a complete
damage reproduction.

## Scientific limitations

The geometry is a buffered-strut approximation; its envelope is approximately
76.50 × 48 × 25 mm, not exactly the nominal 76 × 48 × 25 mm. C3D6 wedges with two
depth layers are not the paper's 10,068 C3D8R elements. End-contact discretization,
mesh aspect ratios, penalty sensitivity, reaction balance and convergence at
large buckling/collapse require further work. Implicit static equilibrium can
fail at snap-through; no dynamics or mass scaling has been introduced to hide
that failure. Frictional work, penetration and mesh convergence are not validated.

Dong et al., Materials & Design 182 (2019), 108036 reports the target Bi-V → Z
sequence, mean plateau stress ~0.28 MPa and densification strain ~0.58. None of
these experimental quantities is validated by the current lattice calculations.
