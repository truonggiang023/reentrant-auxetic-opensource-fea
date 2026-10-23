# CalculiX 2.20 finite-strain ductile material

Build with `scripts/build_ccx_damage.sh` after installing the Python development
environment. This builds CalculiX itself with external-behaviour support, plus
`libREENTRANT.so`. The shared library is **not** a Python FE solver.

`run_paper.py` selects `ccx-damage` for a configuration containing
`ductile_damage`. It passes `REENTRANT_MATERIAL_DATA` to that process only. The
writer produces an input deck and matching `.material.dat` sidecar. Keep both
files together. The nine deck properties are E, nu, failure displacement,
residual stiffness fraction, characteristic length, negative hardening-point
count, fracture-point count, compression cutoff, and viscosity time. The full
stress/plastic-strain and triaxiality/fracture-strain tables are loaded once from
the sidecar before solver threads start. No credentials are involved.

**CalculiX 2.20 has 21-entry local property buffers. Never place a full hardening
table under `*USER MATERIAL`.** Earlier development trials with oversized inline
tables are invalid and must not be interpreted as solver evidence. The regression
suite exercises the full 178-point repository table through the safe sidecar.

The law uses a multiplicative plastic metric and isotropic Hencky elasticity.
Spectral logarithmic elastic trial strains give a radial J2 return, with true
Cauchy flow stresses converted to Kirchhoff stresses using J. The plastic metric
remains isochoric. Initiation integrates dPEEQ/epsilon_f(eta) above a specified
compression cutoff. Post-initiation equivalent plastic displacement is L*dPEEQ,
including only the fraction after initiation in a crossing increment. Linear
softening degrades the complete stress, with a nonzero residual stiffness.
Characteristic length is the cube root of initial element volume, binned to
three significant figures. Highly anisotropic meshes still require a directional
regularization and mesh sensitivity study.

Optional viscosity uses backward Euler: D_new=(tau*D_old+dt*D_target)/(tau+dt).
It is an explicit numerical regularization parameter in solver time units,
not an experimentally measured strain-rate law. Zero disables it. The example
uses tau=0.005 with unit static pseudo-time; compare smaller values before
interpreting results. The symmetric numerical material tangent is an approximation
because CCX's native interface accepts 21 symmetric components. Stress-only calls
also return a valid tangent, avoiding unset arrays in solver predictor paths.

State variables: SDV1–6 inverse plastic metric (11,22,33,12,13,23), SDV7 PEEQ,
SDV8 initiation integral, SDV9 post-initiation displacement, SDV10 damage. State
updates start from `xstateini`; Newton trial calls do not accumulate irreversible
history. `PEEQ` output from the builtin plastic material is not the user material's
PEEQ: use SDV7. Integration-point SDVs are printed in DAT; FRD values can be
extrapolated/averaged and are unsuitable for pointwise failure decisions.

The scoped `patch_hourglass.py` compatibility patch uses E*(1-D) consistently
for hourglass forces and stiffness in `@REENTRANT` C3D8R elements. It avoids
unset stress-only tangent data and negative or inconsistent hourglass moduli.
It changes no other material. Both force and stiffness paths are regression-tested
against builtin elastic hourglass loading and a two-brick bending problem.
Original vendor checksums are checked before the patch is applied or recognized.
The build records solver, library and source hashes outside the checkout.

The checksum-verified `patch_timepoints.py` also corrects CCX 2.20 output-point
bookkeeping after a face-contact retry. A retry can accept a scheduled time after
clearing the output-point flag and rewinding its pointer; the next scheduled
increment then becomes zero. The patch advances past completed points using only
machine-roundoff tolerance, restores the output flag, and leaves future points
and nonlinear convergence tolerances intact. `test_timepoint_retry.py` compiles
the actual helper and checks completed, roundoff-close and future points. The full-lattice rerun
also accepted the previously failing time 0.185 and continued with positive
increments beyond it; qualification outcomes are recorded in validation/.

Supported checks include rigid-rotation objectivity, plastic incompressibility,
repeat trials, elastic unloading, compression cutoff, regularized softening energy
for three band lengths, viscosity relaxation, solver-side tensile load loss, and
free-DOF elastic bending. These checks qualify components, not a Dong reproduction.

Limitations: isothermal isotropic law; no qualified initial-stress treatment;
no automatic element deletion or creation of fracture surfaces. Residual elements
remain in the mesh and on contact surfaces, and can act as ghost contact walls.
The initial fracture-locus endpoints, compression cutoff and viscosity in the
example are assumptions. The paper gives only epsilon_f=0.047 at eta=0.333,
failure displacement 0.001 mm, and a range 0.08–0.018 for other fracture strains;
it does not provide the missing triaxiality coordinates. No complete calibrated
locus or experimental agreement can be inferred from those statements.
