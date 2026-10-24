# Re-entrant honeycomb FEM with CalculiX

Research implementation for the Dong et al. (2019) compression problem.
**Experimental agreement and the complete Bi-V → Z → collapse sequence are not
established.** CalculiX performs the FEM solve; Python handles geometry, meshing,
automation, tests and postprocessing.

The primary workflow now includes a compiled finite-strain J2 ductile-softening
material for CalculiX 2.20, separate steel platens with sliding contact, and
optional disjoint wall-contact patches. See [model and evidence](DONG2019_CALCULIX.md)
and the [native material contract](native/README.md) for qualifications and limits.

## Linux cloud setup

The prepared machine uses `/workspace/fea-venv`, `/workspace/bin/ccx` (packaged
CalculiX 2.20) and `/workspace/bin/ccx-damage` (the qualified external-material
build). No service needs to remain running.

```sh
python -m venv /workspace/fea-venv
/workspace/fea-venv/bin/pip install -e '.[dev,post]'
bash scripts/build_ccx_damage.sh
export PATH=/workspace/bin:/workspace/fea-venv/bin:$PATH
OMP_NUM_THREADS=2 CCX_EXECUTABLE=/workspace/bin/ccx \
  CCX_DAMAGE_EXECUTABLE=/workspace/bin/ccx-damage pytest -q -rx
```

The build uses signed Debian package/source metadata and checks the vendor files
before applying scoped hourglass and output-time retry compatibility patches. It installs locally in
`/workspace`, preserves TLS/package verification, and records artifact hashes.
It currently targets Debian amd64 with GCC 14 and the downloaded Fortran 12
compiler. The bundled Windows 2.13 executable is not qualified for this material.

## Generate and run

```sh
# Builtin J2 plasticity, wedge mesh, platen contact:
python run_paper.py --config configs/dong2019_contact.json --run

# Native coupled softening, all-quad C3D8R mesh, wall-contact candidate:
D

# Plot accepted solver frames; missing targets are shown explicitly:
MPLCONFIGDIR=/workspace/mplconfig python plot_run.py results/dong2019_contact.inp
```

The damage configuration is **explicitly uncalibrated**: the paper publishes one
fracture-locus point and a range for other fracture strains, without the missing
triaxiality coordinates. Its assumptions and viscosity are written in the JSON.
Use `--imperfection-mm .003` only as a labelled imperfection sensitivity study;
that value is not a measured specimen defect.

The writer produces `.inp` and `.material.dat` together. The runner sets
`REENTRANT_MATERIAL_DATA` for that solver process and chooses `ccx-damage` when
a damage law is present. Do not replace the nine descriptor properties with a
large inline `*USER MATERIAL` table: CalculiX 2.20 has 21-entry local buffers.

The runner saves logs, preserves old result files, and verifies a completed status
file and material-load signature. A zero process exit alone is insufficient.
`frd.committed_fields` filters out unconverged “best solution” frames after a
failure. Target-output mode includes exact strains 0.12, 0.35 and 0.65; the example
also requests a regular output grid for reaction histories. Failed or unrun
frames must not be presented as deformation evidence.

## Model scope

Units are mm, N and MPa. Compression is along Y and extrusion along Z. The
paper-oriented deck fixes all core U3 DOFs and removes one in-plane rigid
translation with a centreline roller. Only the outer steel platen faces receive
the prescribed/support displacement; their contact faces remain elastically
free. Static time is pseudo-time, not experimental seconds.

The C3D8R mesher uses Gmsh quad recombination and conforming refinement. Wedges
and fully integrated C3D8 bricks are also available. Geometry and mesh convergence
remain to be demonstrated; an element-type match is not a reproduction.

Damage is coupled to stress and state inside the compiled material. It uses
characteristic-length softening with retained residual stiffness. **It does not
automatically delete elements or create new crack-contact surfaces.** Retained
failed faces can affect collapse/contact. The wall algorithm uses initially
facing disjoint patches and cannot guarantee all new contacts after severe
rotations. The test suite exposes known frictional mortar convergence failures
as XFAILs; penalty contact is the selected formulation.

The legacy `reentrant-fea`/`--advanced` writers and standalone Python material-point
prototypes remain for comparison. The compiled material workflow is
`run_paper.py` with an explicit `ductile_damage` configuration, rather than the
legacy advanced writer's unsupported damage switch.
