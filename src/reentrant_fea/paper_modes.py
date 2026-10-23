"""Dong 2019 thin-wall replication candidate for CalculiX.

Generates elastic steel platens, unilateral frictional platen/core contact,
optional EXPERIMENTAL core self-contact and three deformation snapshots.
This is NOT a validated damage/deletion implementation.
"""
from pathlib import Path
import numpy as np
from .contact import exterior_faces


def _nset(lines, name, ids):
    lines.append("*NSET,NSET=" + name)
    ids = list(ids)
    for k in range(0, len(ids), 12):
        lines.append(", ".join(map(str, ids[k:k + 12])))


def _surface(lines, name, faces):
    if not faces:
        raise ValueError("Empty contact surface " + name)
    lines.append("*SURFACE,NAME=" + name + ",TYPE=ELEMENT")
    lines.extend(f"{eid}, {side}" for eid, side in faces)


def _platen(x0, x1, y0, y1, depth, nid, eid):
    """One first-order C3D8 hex platen, bottom-face S1 / top-face S2."""
    xyz = [
        (x0, y0, 0), (x1, y0, 0), (x1, y0, depth), (x0, y0, depth),
        (x0, y1, 0), (x1, y1, 0), (x1, y1, depth), (x0, y1, depth),
    ]
    node_ids = list(range(nid, nid + 8))
    return xyz, node_ids, (eid, node_ids)


def write_paper_deck(
    path, nodes, elements, top, bottom, config, plastic,
    *, self_contact=False, damage=False
):
    """Generate one 3-step nonlinear CalculiX model with real platen contact.

    damage=True explicitly fails because CCX 2.13 does not implement
    Abaqus ductile damage evolution + element deletion.
    """
    if damage:
        raise NotImplementedError(
            "CCX 2.13 cannot reproduce Abaqus ductile damage/deletion. "
            "A verified extended CalculiX engine is required."
        )
    if self_contact and not config.get("experimental_self_contact_ack", False):
        raise ValueError("Set experimental_self_contact_ack=true in config")
    nodes = np.asarray(nodes, float)
    elements = np.asarray(elements, int)
    if len(plastic) < 2:
        raise ValueError("Missing plastic table")
    ylo, yhi = float(nodes[:, 1].min()), float(nodes[:, 1].max())
    xlo, xhi = float(nodes[:, 0].min()), float(nodes[:, 0].max())
    zlo, zhi = float(nodes[:, 2].min()), float(nodes[:, 2].max())
    height = yhi - ylo
    width = float(config.get("platen_width", 76.0))
    depth = zhi - zlo
    margin = (width - (xhi - xlo)) / 2
    if margin < -1e-5:
        raise ValueError(
            f"Core width {xhi-xlo:.5g} mm exceeds nominal platen width {width}. "
            "Correct geometry before replication."
        )
    xmid = (xlo + xhi) / 2
    px0, px1 = xmid - width / 2, xmid + width / 2
    th = float(config.get("platen_thickness", 25.0))
    gap = float(config.get("platen_initial_gap", 0.0))
    if th <= 0 or gap < 0:
        raise ValueError("Invalid platen geometry")

    lines = [
        "*HEADING", "Dong2019 candidate: platen contact, thin-wall Bi-V/Z modes",
        "** Core generated with Python+Gmsh; CalculiX handles the nonlinear solve",
        "** NOT a verified damage or self-contact benchmark",
        "*NODE",
    ]
    for i, v in enumerate(nodes, 1):
        lines.append(f"{i}, {v[0]:.10g}, {v[1]:.10g}, {v[2]:.10g}")
    n0 = len(nodes) + 1
    e0 = len(elements) + 1
    lower = _platen(px0, px1, ylo - gap - th, ylo - gap, depth, n0, e0)
    upper = _platen(px0, px1, yhi + gap, yhi + gap + th, depth, n0 + 8, e0 + 1)
    for xyz, ids, _ in (lower, upper):
        for pid, (x, y, z) in zip(ids, xyz):
            lines.append(f"{pid}, {x:.10g}, {y:.10g}, {z:.10g}")

    lines.append("*ELEMENT,TYPE=C3D6,ELSET=CORE")
    for i, el in enumerate(elements, 1):
        lines.append(f"{i}, " + ", ".join(str(int(j) + 1) for j in el))
    lines.append("*ELEMENT,TYPE=C3D8,ELSET=BOTTOM_PLATEN")
    lines.append(f"{e0}, " + ", ".join(map(str, lower[1])))
    lines.append("*ELEMENT,TYPE=C3D8,ELSET=TOP_PLATEN")
    lines.append(f"{e0+1}, " + ", ".join(map(str, upper[1])))

    _nset(lines, "TOP_PLATEN_NODES", upper[1])
    _nset(lines, "BOTTOM_PLATEN_NODES", lower[1])
    tol = max(1e-6, depth * 1e-7)
    flat = np.flatnonzero(
        (np.abs(nodes[:, 2] - zlo) < tol) |
        (np.abs(nodes[:, 2] - zhi) < tol)
    ) + 1
    _nset(lines, "CORE_OUT_OF_PLANE", flat)

    lines += [
        "*MATERIAL,NAME=ALSI10MG", "*ELASTIC",
        f"{float(config['young']):.10g}, {float(config['poisson']):.10g}",
        "*PLASTIC",
    ]
    lines += [f"{s:.12g}, {p:.12g}" for s, p in plastic]
    lines += [
        "*SOLID SECTION,ELSET=CORE,MATERIAL=ALSI10MG",
        "*MATERIAL,NAME=STEEL", "*ELASTIC",
        f"{float(config.get('steel_young',210000)):.10g}, "
        f"{float(config.get('steel_poisson',.3)):.10g}",
        "*SOLID SECTION,ELSET=TOP_PLATEN,MATERIAL=STEEL",
        "*SOLID SECTION,ELSET=BOTTOM_PLATEN,MATERIAL=STEEL",
    ]

    # Core free-face inventory. These are the contact-facing end faces;
    # front/back triangular faces must NOT participate in wall self-contact.
    outer = exterior_faces(elements)
    eps = max(1e-6, float(config.get("mesh_size", 1.0)) * 0.01)
    top_faces = []
    bottom_faces = []
    wall_faces = []
    for eid, side, verts in outer:
        xyz = nodes[list(verts)]
        if np.all(abs(xyz[:, 1] - yhi) < eps):
            top_faces.append((eid, side))
        elif np.all(abs(xyz[:, 1] - ylo) < eps):
            bottom_faces.append((eid, side))
        elif (
            not np.all(abs(xyz[:, 2] - zlo) < eps)
            and not np.all(abs(xyz[:, 2] - zhi) < eps)
        ):
            wall_faces.append((eid, side))
    _surface(lines, "CORE_CONTACT_TOP", top_faces)
    _surface(lines, "CORE_CONTACT_BOTTOM", bottom_faces)
    _surface(lines, "STEEL_CONTACT_TOP", [(e0 + 1, "S1")])
    _surface(lines, "STEEL_CONTACT_BOTTOM", [(e0, "S2")])
    if self_contact:
        _surface(lines, "CORE_WALLS", wall_faces)

    mu = float(config.get("friction", 0.35))
    lines += [
        "*SURFACE INTERACTION,NAME=PLATEN_CONTACT",
        "*SURFACE BEHAVIOR,PRESSURE-OVERCLOSURE=HARD",
        "*FRICTION", f"{mu:.8g}",
        "*CONTACT PAIR,INTERACTION=PLATEN_CONTACT,TYPE=SURFACE TO SURFACE",
        "CORE_CONTACT_TOP,STEEL_CONTACT_TOP",
        "*CONTACT PAIR,INTERACTION=PLATEN_CONTACT,TYPE=SURFACE TO SURFACE",
        "CORE_CONTACT_BOTTOM,STEEL_CONTACT_BOTTOM",
    ]
    if self_contact:
        lines += [
            "** EXPERIMENTAL self-contact: same-surface support is solver-version dependent",
            "*CONTACT PAIR,INTERACTION=PLATEN_CONTACT,TYPE=SURFACE TO SURFACE",
            "CORE_WALLS,CORE_WALLS",
        ]

    lines += [
        "*BOUNDARY",
        "BOTTOM_PLATEN_NODES,1,3,0.",
        "TOP_PLATEN_NODES,1,1,0.",
        "TOP_PLATEN_NODES,3,3,0.",
    ]
    if config.get("constrain_out_of_plane", True):
        lines.append("CORE_OUT_OF_PLANE,3,3,0.")

    targets = tuple(config.get("snapshot_strains", [0.12, 0.35, 0.65]))
    if sorted(targets) != list(targets) or min(targets) <= 0:
        raise ValueError("snapshot_strains must be positive and increasing")
    initial_inc = float(config.get("initial_increment", 1e-4))
    max_inc = float(config.get("max_increment", 0.01))
    for idx, strain in enumerate(targets, 1):
        lines += [
            f"** SNAPSHOT_{idx}: nominal strain {strain:.6g}",
            "*STEP,NLGEOM=YES,INC=100000",
            "*STATIC",
            f"{initial_inc:.10g},1.0,1e-10,{max_inc:.10g}",
            "*BOUNDARY",
            f"TOP_PLATEN_NODES,2,2,{-strain*height:.10g}",
            "*NODE FILE",
            "U,RF",
            "*EL FILE",
            "S,E,PEEQ",
            "*NODE PRINT,NSET=TOP_PLATEN_NODES",
            "U,RF",
            "*END STEP",
        ]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf8")
    return path
