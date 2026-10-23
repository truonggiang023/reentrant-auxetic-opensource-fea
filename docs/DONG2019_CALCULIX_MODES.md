# Dong et al. (2019) CalculiX deformation-mode replication candidate

This branch is a **research prototype**, not a validated reproduction of Fig. 11(b).

## Paper targets

Thin wall 0.30 mm, five columns and six rows, nominal 76 x 48 x 25 mm. Target snapshots nominal compressive strain 0.12 (Bi-V), 0.35 (Z) and 0.65 (deep folding). Experimental thin-wall plateau stress 0.28 MPa, densification at about 0.58. Paper FE baseline: C3D8R with 10,068 elements, steel platens, friction 0.35, out-of-plane constraint, and Abaqus ductile damage initiation 0.047 plastic fracture strain at triaxiality 0.333 with displacement-based linear softening 0.001 mm.

Our current model **still uses C3D6 wedge elements, not C3D8R**. The prior observed relative density ~0.11854 differs from the paper's thin-wall 0.109. Both can affect deformation modes.

## One command: CCX 2.13 elastoplastic with platen contact and snapshots

From repository root with installed Python dependencies:

    python run_paper_modes.py --generate-only
    python run_paper_modes.py --threads 16

The runner uses tools/calculix/bin/ccx/ccx213.exe by default. It regenerates results/reentrant.inp and solves it. Alternatively, generate-only then execute the existing python run.py. Generated model: elastoplastic AlSi10Mg (full numeric reference data), two elastic C3D8 steel platens, two unilateral frictional contact pairs, bottom platen fixed, top driven in Y, core front/back U3 fixed. Three nonlinear static steps end at strains 0.12, 0.35 and 0.65. Open results/reentrant.frd to compare step-end mode shapes.

**Self-contact is experimental** in this model and is not claimed to work for ccx213:

    python run_paper_modes.py --self-contact --threads 16

The same-surface contact definition may be rejected or fail to enforce contact; this is not solver-validated. Large compression without self-contact is not physically representative after cell-wall collision.

## Experimental CCX 2.23 Rice-Tracey element deletion

Native DAMAGE INITIATION was introduced in CCX 2.23. There is **no displacement-based Abaqus DAMAGE EVOLUTION** in this pipeline. The optional config uses illustrative Rice-Tracey card parameters 1.0 and 0.047. **These are not calibrated to Dong et al.** In particular, the tensile fracture strain 0.047 cannot automatically be treated as the Rice-Tracey limit.

    python run_paper_modes.py --damage --ccx "F:\path\to\ccx223.exe" --threads 16

This mode is only a speculative deck-generation trial, not a validated drop-in replacement. The official CCX 2.23 deletion implementation has known singular-matrix failure when removed elements leave isolated nodes (CalculiX Discourse October 2026). The source paper's displacement-regularized damage law would require a custom solver extension or different solver. Never claim damage or mode validation based solely on file generation.

## Implemented vs unimplemented

- Implemented in Python generator: steel platen discretization, friction contact surfaces, out-of-plane constraint, strain targets, full original plastic hardening data, reference geometry config.
- Optional and unverified: lateral wall self-contact; CCX >=2.23 Rice-Tracey deletion with noncalibrated parameters.
- **Not implemented:** the paper's Abaqus displacement-softened damage evolution and equivalently calibrated fracture, C3D8R reference mesh, true contact-verified large compression, guaranteed Bi-V to Z modes.

For comparison: read displacement snapshots from FRD at the end of each step; compute nominal strain = -U2/H0 and nominal compressive stress = |sum RF2 on top platen|/(76*25) MPa. Do not tune fracture damage to force the mode shape before checking geometry, contact and boundary conditions.

Reference: Dong et al., Materials & Design 182 (2019) 108036.
