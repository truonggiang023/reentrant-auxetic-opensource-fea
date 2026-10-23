# Re-entrant auxetic FEA — independent Python / VS Code project

**Status: research prototype; NOT a validated Abaqus-equivalent or Dong-2019-equivalent damage solver.**

## Implemented and tested

- Parametric 5 × 6 re-entrant XY geometry (Shapely union of finite-width struts), including `lambda_m`.
- Native Gmsh 2D triangular mesh, extruded in Python into genuine 6-node wedge (C3D6) solids.
- CalculiX `ccx` nonlinear finite-strain (`NLGEOM`) J2 elastic–plastic baseline input deck.
- Top/bottom nodal rigid-grip kinematic approximation: common imposed vertical displacement at all top nodes, fixed bottom nodes. No steel platen elements in current deck.
- Exact plastic table ingestion from original setup Python file using AST literal parsing; no Abaqus installation/import needed.
- Separate Python material-point *prototype* for triaxiality-based ductile initiation, characteristic-length softening, and failure flag, with unit tests. **Not connected to the FE engine.**

## Not implemented / not validated

1. **Deformable steel platens and bonded face-to-face tie**: current model is a simplified rigid-grip displacement boundary condition, not physical platen meshes.
2. **Self-contact and friction**: not active in generated deck; crushing after cell-wall contact will not be physically valid.
3. **Coupled ductile damage evolution and element deletion**: material-point algorithm exists only as a tested demonstration, and is **not integrated into CalculiX**. `*PLASTIC` does not simulate damage or deletion.
4. Explicit dynamics, mass scaling, finite-strain algorithm matching, solver energy equivalence, mesh convergence, calibration, plateau stress and mode-shape validation are not completed.
5. Original reference uses C3D8R; modified target is C3D6. Geometry uses geometric buffers, not guaranteed to match Abaqus mitered offsets exactly.

Do **not** report any generated result as a validated ductile-damage simulation. An appropriate next stage is to integrate a custom constitutive and deletion-capable open-source solver, and add real self-contact, then validate against experimental Bi-V and Z modes and plateau stress.

## Windows + VS Code

Prerequisites: Python 3.10+, [CalculiX ccx](http://www.calculix.de/), VS Code, optionally ParaView. `ccx` must be on PATH or passed via `--ccx`.

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
pytest -q
python -m reentrant_fea.cli --config configs/dong2019_thin.json --plastic reference_setup_numeric_source.py --output results/reentrant.inp
ccx -i reentrant
```

Run `ccx -i reentrant` from inside the `results` folder; or append `--run` to the Python command. The original input Python file is **numeric source data only**, not executed and not required by Abaqus.

For a clean self-contained project, the reference input is bundled as `reference_setup_numeric_source.py`; the CLI reads only the numeric plastic table.

## Test coverage

`pytest -q` tests geometry validity, morphology transition, basic damage initiation/evolution/failure and generated solver deck syntax. This does not test convergence or validate CalculiX solver behavior.

## Scientific reference

Dong et al., *Materials & Design* 182 (2019) 108036. Thin wall t=0.3 mm, 5 × 6 cells, nominal 76 × 48 × 25 mm, experimental plateau 0.28 MPa and densification strain 0.58. Original published FE: C3D8R (10,068 elements), friction coefficient 0.35, compression 2 mm/min. Do not equate these with the modified wedge / bonded-grip model.

## Units

mm, N, MPa; geometry extrusion along Z and compression along Y. The baseline is implicit static; time is pseudo-time, not experimental seconds.

## Advanced experimental solver deck (2026-10-23)

The project now includes `advanced_solver.py`, `contact.py`, and `damage_history.py`.

- `--advanced` creates **rigid-body bonded grips** using CalculiX `*RIGID BODY` and reference nodes; this is **not a meshed steel platen**. The rigid-grip deck now passes a small CalculiX 2.20 compression regression; the full lattice remains unvalidated.
- `--self-contact` now raises an explicit unsupported-feature error. Same-surface contact pairs are not accepted as a verified crushing self-contact implementation.
- `--ductile-damage` **intentionally raises `NotImplementedError`**. The project does not pretend that an offline damage tracker or `*PLASTIC` provides coupled ductile softening or element deletion.
- `damage_history.evaluate_damage_csv(...)` integrates the prototype damage law on a user-provided integration history CSV. It can flag elements for *analysis only*. It **does not delete elements from an active nonlinear simulation**. Its triaxiality must be computed from stress tensors exported by a solver; `PEEQ` alone is insufficient.

Example:

```powershell
python -m reentrant_fea.cli --config configs/dong2019_thin.json --plastic reference_setup_numeric_source.py --output results/advanced.inp --advanced
```

### What is still required to honestly claim a complete damage solver

A verified constitutive implementation with stress return mapping, triaxiality-dependent initiation, mesh-objective displacement regularization, damage-degraded tangent/stress, and integration-point state management, plus solver-supported **element deletion during the nonlinear solve**. The contact implementation must support self-contact between crushing internal walls without artificial initial penetrations, and needs benchmark checks. These are not provided by the CalculiX baseline in this repository. Do not use the current output for calibrated failure or plateau claims.

## CalculiX 2.20 solver evidence

See [DONG2019_CALCULIX.md](DONG2019_CALCULIX.md) for Linux setup, actual contact benchmarks, integration-point ductile initiation diagnostics, and unresolved lattice convergence, coupled damage and self-contact. The three target lattice modes remain unverified.
