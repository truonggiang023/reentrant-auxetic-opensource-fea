# Re-entrant open-source FEA migration

This is a migration/verification starter, NOT a validated Abaqus-equivalent damage solver.

## Reference
- Abaqus script: `ReEntrant_C3D6_BondedPlatens_Validation.py`
- Dong et al. (2019): original FE model uses C3D8R; modified user's model uses C3D6 and tied platen interfaces.

## Install (Windows VS Code)
- Python 3.11+; `pip install numpy pandas matplotlib meshio gmsh`
- Gmsh executable (optional if Python gmsh package works)
- CalculiX `ccx` executable; ParaView
- Add CalculiX executable folder to PATH.

## Required solver fidelity work
1. Export actual mesh and surface sets from Abaqus to INP or build exact equivalent geometry with Gmsh.
2. Map Abaqus C3D6 nodes to CalculiX C3D6, check connectivity and orientation.
3. Replicate end platen tie (or prescribed face displacement for a rigid-platen idealization); avoid double contact/tie.
4. Replicate finite sliding and self-contact and friction coefficients.
5. Check plastic true stress vs true plastic strain, finite strain behavior and rate effects.
6. Implement and verify stress-triaxiality-dependent ductile initiation, displacement-based softening, and element deletion in a supported custom solver/constitutive integration. Native CalculiX *PLASTIC alone is not equivalent.
7. Perform mesh convergence and compare force-displacement, mode shapes, damage and energies against Abaqus.

## Baseline script
`baseline_from_inp.py` reads an exported Abaqus mesh INP containing *NODE and *ELEMENT, TYPE=C3D6, and creates a CalculiX baseline *STATIC deck with imposed face displacements. This is only an elastic-plastic debugging baseline: it deliberately does NOT claim to simulate damage or platen contact. Use it to verify mesh/material before implementing full equivalence.

Example:
`python baseline_from_inp.py --input exported_core_mesh.inp --output core_baseline.inp --height 47.7 --displacement 28.62`

The original user's Abaqus script has the complete plasticity table. For numerical fidelity, extract the exact *PLASTIC data from an Abaqus-generated .inp and pass it using `--plastic-inp`.
