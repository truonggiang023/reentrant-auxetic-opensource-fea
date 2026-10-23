# Coupled damage verification track

A standalone Python small-strain C3D6 FEM verification package was developed separately from the CalculiX baseline. It includes 6-point wedge integration, nonlinear sparse equilibrium assembly, material-point J2 plasticity with triaxiality-dependent damage, and removal of fully failed elements from stiffness/internal-force assembly in subsequent increments.

The implementation has passed five local unit tests (geometry volume, rigid translation, elastic equilibrium, deletion, and post-deletion zero reaction). The Python code package is not yet committed to this branch; it must be reviewed and uploaded separately.

**Not yet completed:** finite-strain constitutive integration, consistent tangent, reliable large-strain crushing convergence, frictional self-contact, mesh-regularized fracture, full honeycomb results, experimental validation against Dong et al. (2019), or a custom CalculiX executable. The existing CalculiX baseline and visualize.py remain unchanged.
