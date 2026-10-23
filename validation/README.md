# Native CalculiX qualification evidence

`native_qualification.json` describes the selected backend, exact source/build
and input hashes, test outcomes, available target frames, calibration assumptions
and numerical limits. `native-final.json` contains the selected run assessment;
`native-final.sta.txt` preserves its accepted increments and failed attempts.
`native-tests.txt` records the full suite. These qualify implementation components;
they do not establish experimental reproduction.

The changes add C3D8R quad meshing, elastic steel contact faces, disjoint wall
contact with topology exclusions, and a solver-side multiplicative J2 material
with triaxiality-dependent initiation and length-regularized displacement
softening. The material uses nine descriptors plus a read-only hardening/locus
sidecar to stay within CCX 2.20's 21-property local buffers. Checksum-verified
patches make native hourglass force/stiffness consistent and prevent a zero
scheduled increment after an accepted output-point retry.

The tests include native material invariance/history/softening-energy checks,
full-table solver tensile load loss, elastic hourglass and free bending
comparisons, platen compression at all three target strains, wall-gap closure,
and accepted-frame filtering. Three coarse frictional mortar cases remain known
XFAILs; penalty counterparts pass. Do not interpret the coupon target strains
as the full honeycomb reaching those targets.

The paper supplies only one located fracture-locus point. Other triaxiality
coordinates, compression cutoff and viscosity remain explicitly uncalibrated.
Damaged elements retain stiffness and contact faces; no deletion, new fracture
surfaces or complete general self-contact is implemented. Global mesh,
hourglass-energy, contact-penalty and viscosity sensitivity remain unqualified.
The full Bi-V/Z/severe-collapse sequence and experimental agreement are not
established. The PR should remain a draft research implementation. The corrected lattice
was deliberately stopped after 20 minutes of qualification, at accepted strain
0.14422 with repeated contact cutbacks; it did not reach 0.35 or 0.65.

`native-imperfect-stopped.*` is a labelled 0.003 mm imperfection comparison,
deliberately stopped to prioritize the corrected backend. Its accepted fields
are valid, but the run did not complete. `native-zero-before-timepoint-fix.*`
records the earlier zero-step failure. Development trials with oversized inline
material tables and incorrect earlier hourglass stiffness are excluded.

The large raw FRD/DAT files are local ignored results, not public artifacts.
Assessment and plotting scripts read only accepted status times and exclude
CalculiX's unconverged best-solution frames. Front-plane overlap is a geometric
diagnostic, not a proof of complete 3-D contact validity. The uploaded paper and
its rendered figures are not included.
