"""Dong et al. 2019 thin-wall replication candidate: CalculiX is the ONLY FE solver.

Creates real, separately meshed steel platens and sliding contact instead of
bonding or clamping the honeycomb end nodes. Self-contact is opt-in because
single-surface self-contact requires verification in the installed CCX build.
This is a numerical research candidate, NOT a validated Dong-2019 reproduction.
"""
from pathlib import Path
import numpy as np
from .contact import exterior_faces


def _nset(lines, name, ids):
    ids = [int(x) for x in ids]
    if not ids:
        raise ValueError("Empty NSET: " + name)
    lines.append("*NSET,NSET=" + name)
    for i in range(0, len(ids), 12):
        lines.append(", ".join(map(str, ids[i:i + 12])))


def _surface(lines, name, entries):
    if not entries:
        raise ValueError("Empty contact surface: " + name)
    lines.append("*SURFACE,NAME=" + name + ",TYPE=ELEMENT")
    lines.extend(f"{eid}, {face}" for eid, face in entries)


def write_paper_contact_deck(
    path, nodes, elements, config, plastic, *,
    self_contact=False, imperfections=None,
):
    """Create a CalculiX nonlinear static contact model (mm, N, MPa).

    imperfections: optional per-node XY perturbations in mm, shape (N, 2).
    They are deterministic INPUT geometric imperfections, not fitted modes.
    """
    nodes = np.asarray(nodes, dtype=float).copy()
    elements = np.asarray(elements, dtype=int)
    if len(nodes) < 6 or len(elements) < 1 or not np.isfinite(nodes).all():
        raise ValueError("Invalid core mesh")
    if imperfections is not None:
        perturb = np.asarray(imperfections, dtype=float)
        if perturb.shape != (len(nodes), 2):
            raise ValueError("imperfections must be (N, 2)")
        nodes[:, :2] += perturb
    if len(plastic) < 2:
        raise ValueError("Expected original calibrated true stress / plastic strain table")
    ymin, ymax = np.min(nodes[:, 1]), np.max(nodes[:, 1])
    zmin, zmax = np.min(nodes[:, 2]), np.max(nodes[:, 2])
    xmin, xmax = np.min(nodes[:, 0]), np.max(nodes[:, 0])
    height = ymax - ymin
    if height <= 0:
        raise ValueError("Zero specimen height")
    # Planens span the full specimen in X and depth in Z.
    width = xmax - xmin
    side_margin = max(float(config.get("platen_margin", 1.0)), 0.0)
    x0, x1 = xmin - side_margin, xmax + side_margin
    gap = float(config.get("platen_gap", 0.0))
    thick = float(config.get("platen_thickness", 10.0))
    if thick <= 0 or gap < 0:
        raise ValueError("Invalid platen parameters")

    node_lines = [f"{i}, {x:.12g}, {y:.12g}, {z:.12g}"
                  for i, (x, y, z) in enumerate(nodes, 1)]
    ncore = len(nodes)
    # Numbering: 1-4 front (zmin), 5-8 back (zmax).
    def plate(xlo, xhi, ylo, yhi, start):
        xyz = [(xlo, ylo, zmin), (xhi, ylo, zmin),
               (xhi, yhi, zmin), (xlo, yhi, zmin),
               (xlo, ylo, zmax), (xhi, ylo, zmax),
               (xhi, yhi, zmax), (xlo, yhi, zmax)]
        return ([f"{start+i}, {x:.12g}, {y:.12g}, {z:.12g}"
                 for i, (x,y,z) in enumerate(xyz)], list(range(start,start+8)))

    lower, lower_ids = plate(x0,x1,ymin-gap-thick,ymin-gap,ncore+1)
    upper, upper_ids = plate(x0,x1,ymax+gap,ymax+gap+thick,ncore+9)
    lines = ["*HEADING", "Dong2019 candidate: sliding platen / large deformation",
             "** Units: mm, N, MPa; stress=top reaction/(76*25) only for exact paper envelope",
             "*NODE"] + node_lines + lower + upper
    lines.append("*ELEMENT,TYPE=C3D6,ELSET=CORE")
    for eid, el in enumerate(elements,1):
        lines.append(f"{eid}, " + ", ".join(str(int(i)+1) for i in el))
    lower_el, upper_el = len(elements)+1, len(elements)+2
    lines += ["*ELEMENT,TYPE=C3D8,ELSET=PLATENS",
              f"{lower_el}, " + ", ".join(map(str,lower_ids)),
              f"{upper_el}, " + ", ".join(map(str,upper_ids))]

    faces = exterior_faces(elements)
    tol = max(1e-7, float(config.get("contact_face_tol", 1e-6)))
    core_bottom, core_top, sides = [], [], []
    for eid, face, verts in faces:
        y = nodes[list(verts), 1]
        if np.all(np.abs(y-ymin) <= tol):
            core_bottom.append((eid,face))
        elif np.all(np.abs(y-ymax) <= tol):
            core_top.append((eid,face))
        else:
            sides.append((eid,face))
    _surface(lines, "CORE_BOTTOM", core_bottom)
    _surface(lines, "CORE_TOP", core_top)
    _surface(lines, "BOTTOM_PLATEN_CONTACT", [(lower_el,"S5")])
    _surface(lines, "TOP_PLATEN_CONTACT", [(upper_el,"S3")])
    if self_contact:
        _surface(lines, "CORE_LATERAL", sides)
    _nset(lines,"BOTTOM_PLATEN_NODES",lower_ids)
    _nset(lines,"TOP_PLATEN_NODES",upper_ids)
    z_faces = [i+1 for i,p in enumerate(nodes)
               if abs(p[2]-zmin)<tol or abs(p[2]-zmax)<tol]
    _nset(lines,"CORE_PLANE_STRAIN",z_faces)

    lines += ["*MATERIAL,NAME=ALSI10MG","*ELASTIC",
              f"{config['young']:.12g}, {config['poisson']:.12g}",
              "*PLASTIC"]
    lines.extend(f"{s:.12g}, {p:.12g}" for s,p in plastic)
    lines += ["*MATERIAL,NAME=STEEL","*ELASTIC",
              f"{config.get('steel_young',210000.0):.12g}, {config.get('steel_poisson',0.3):.12g}",
              "*SOLID SECTION,ELSET=CORE,MATERIAL=ALSI10MG",
              "*SOLID SECTION,ELSET=PLATENS,MATERIAL=STEEL",
              "*SURFACE INTERACTION,NAME=SLIDING",
              "*SURFACE BEHAVIOR,PRESSURE-OVERCLOSURE=LINEAR",
              str(config.get("normal_penalty", 1.0e5)),
              "*FRICTION",
              f"{config.get('friction',0.35):.12g}, {config.get('stick_slope',100000.0):.12g}",
              "*CONTACT PAIR,INTERACTION=SLIDING,TYPE=SURFACE TO SURFACE",
              "CORE_BOTTOM,BOTTOM_PLATEN_CONTACT",
              "*CONTACT PAIR,INTERACTION=SLIDING,TYPE=SURFACE TO SURFACE",
              "CORE_TOP,TOP_PLATEN_CONTACT"]
    if self_contact:
        lines += ["** EXPERIMENTAL: single-surface self contact is CCX-build-dependent",
                  "*CONTACT PAIR,INTERACTION=SLIDING,TYPE=SURFACE TO SURFACE",
                  "CORE_LATERAL,CORE_LATERAL"]

    lines += ["*BOUNDARY","BOTTOM_PLATEN_NODES,1,3,0.",
              "TOP_PLATEN_NODES,1,1,0.","TOP_PLATEN_NODES,3,3,0.",
              "CORE_PLANE_STRAIN,3,3,0.",
              "*STEP,NLGEOM=YES,INC=100000",
              "*STATIC",
              f"{config.get('initial_increment',1e-5):.12g}, 1., {config.get('minimum_increment',1e-10):.12g}, {config.get('maximum_increment',.002):.12g}",
              "*BOUNDARY",
              f"TOP_PLATEN_NODES,2,2,{-float(config['compression_strain'])*height:.12g}",
              "*NODE FILE,GLOBAL=YES","U,RF",
              "*EL FILE","S,E,PEEQ",
              "*NODE PRINT,NSET=TOP_PLATEN_NODES,TOTALS=ONLY","RF",
              "*END STEP"]
    path = Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text("\n".join(lines)+"\n",encoding="utf-8")
    return {"path":str(path),"core_nodes":ncore,"core_elements":len(elements),
            "top_faces":len(core_top),"bottom_faces":len(core_bottom),
            "initial_height":height,"nominal_area":float(config.get("reference_area",76.0*25.0)),
            "self_contact_experimental":bool(self_contact)}
