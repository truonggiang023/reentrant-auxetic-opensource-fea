# Dong et al. 2019 — CalculiX replication candidate

Reference: Dong et al., Materials & Design 182 (2019) 108036, Fig. 11(b).

## Scope
This branch adds **CalculiX-only FEM solving** with separately meshed top and bottom steel platen elements, a contact pair at each core/platen interface, friction 0.35, out-of-plane U3 constraint, original plastic hardening table, nominal strain 0.65, and optional explicitly imposed geometric imperfection. Python is strictly pre/post-processing. The repository's original run.py and baseline are untouched.

It **does not** reproduce Abaqus Explicit general contact, ductile damage initiation/evolution, C3D8R meshing, or verified large-strain honeycomb self-contact. The new single-surface CCX self-contact is intentionally **experimental** and may be rejected or behave incorrectly with ccx213.exe. A completed solver exit alone does not demonstrate mode equivalence.

## Run (Windows PowerShell; existing portable CalculiX)
```powershell
.\.venv\Scripts\python.exe run_paper.py --run --threads 16
.\.venv\Scripts\python.exe convert_paper.py
```
If the CCX213 bundle is installed at a different location, pass `--ccx path/to/ccx.exe`.
Run exploratory single-surface self contact with `--self-contact`, but expect that it may be unsupported. To explore mode selection via a **documented** deterministic imperfection (without fitting to Fig. 11), use `--imperfection-mm 0.01`. Never silently call this a measured specimen imperfection.

Contact physics is explicitly different from the baseline's rigid bonded grips. Large deformation and contact-induced densification may still fail convergence or produce other modes. Contact normal stiffness and load increments are numerical parameters requiring sensitivity analysis.

## Reference conditions and targets
- 5 columns x 6 rows; thin wall 0.30 mm; 76 x 48 x 25 mm nominal specimen
- Material AlSi10Mg; original tensile true stress vs plastic strain table
- Reference paper finite element: 10,068 C3D8R, 19 horizontal / 10 inclined edge divisions; self contact; `mu = 0.35` at platen interfaces, out-of-plane constrained, prescribed compression 2 mm/min
- Target mode frames **Bi-V at strain 0.12**, **Z at strain 0.35**, strong collapse at strain 0.65
- Experimental thin-wall plateau mean ~0.28 MPa, densification mean strain ~0.58; paper FE 5x6 plateau ~0.31 MPa

## Important scientific gaps
1. **Damage**: ccx213.exe has no equivalent of the Abaqus ductile initiation table plus displacement-based linear softening (`u_f=0.001 mm`) and in-solver deletion used by the paper. We have not inserted unsupported keywords. A Python post-processing damage flag is not material damage in the solver.
2. **Core self-contact**: this deck can emit a tentative same-surface contact pair but does not establish support or valid evolution of folded walls. Without successful self-contact, compressed walls may interpenetrate.
3. **Elements**: uses C3D6 wedge mesh, not paper's C3D8R brick. Mesh morphology, nominal geometry, relative density, and integration quality must be checked.
4. **Quasi-static equivalence**: the deck is implicit static, not the Abaqus Explicit simulation in the supplied source script. Loading speed is not a physical variable for rate-independent implicit static.
5. **Validation**: the solution has not been run on your Windows ccx213.exe in this session; figure and curve matching are unconfirmed. Run tests with `pytest tests/test_paper_contact_solver.py -q`.

ParaView: open the converted VTU time series, use Warp By Vector on U, compare geometries at the requested nominal strain frame. Use top reaction magnitude / 1900 mm² for nominal stress only when the exact nominal envelope is appropriate. A proper RF-strain extractor and qualitative Bi-V/Z comparison remain necessary before claiming reproduction.

### Literature-specific caution
The uploaded paper says the thin honeycomb showed Bi-V → Z, with a measured mean plateau 0.28 MPa and mean densification 0.58. It also says damage initiation used fracture plastic strain 0.047 at triaxiality 0.333 (additional fitted values), followed by 0.001 mm linear displacement softening. This **cannot** be represented by an ordinary `*PLASTIC` deck.
