"""Dong et al. 2019 thin-wall replication candidate: CalculiX is the ONLY FE solver.

Creates real, separately meshed steel platens and sliding contact instead of
bonding or clamping the honeycomb end nodes. Optional wall contact uses
disjoint initially facing patches with adjacency and solid-interior exclusions.
This is a numerical research candidate, NOT a validated Dong-2019 reproduction.
"""
from pathlib import Path
import numpy as np
from .contact import exterior_faces, wall_contact_patches
from .plastic import validate_material
from .native_damage import material_cards


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
    plastic = validate_material(config, plastic)
    contact_type = config.get('contact_type', 'SURFACE TO SURFACE')
    if contact_type not in ('SURFACE TO SURFACE', 'MORTAR', 'LINMORTAR'):
        raise ValueError('Unsupported contact_type for CCX 2.20')
    for key, default in [('normal_penalty', 1e5), ('stick_slope', 1e5)]:
        if key == 'stick_slope' and float(config.get('friction', .35)) == 0:
            continue
        value = float(config.get(key, default))
        if not np.isfinite(value) or value <= 0:
            raise ValueError(key + ' must be finite and positive')
    friction = float(config.get('friction', .35))
    if not np.isfinite(friction) or friction < 0:
        raise ValueError('friction must be finite and nonnegative')
    increments = [float(config.get(k, d)) for k, d in
                  [('minimum_increment', 1e-6), ('initial_increment', 1e-5), ('maximum_increment', .002)]]
    if not np.isfinite(increments).all() or not 0 < increments[0] <= increments[1] <= increments[2] <= 1:
        raise ValueError('Require 0 < minimum <= initial <= maximum increment <= 1')
    nodes = np.asarray(nodes, dtype=float).copy()
    elements = np.asarray(elements, dtype=int)
    if nodes.ndim != 2 or nodes.shape[1] != 3 or len(nodes) < 6 or len(elements) < 1 or not np.isfinite(nodes).all():
        raise ValueError("Invalid core mesh")
    reference_nodes=nodes.copy()
    if imperfections is not None:
        perturb = np.asarray(imperfections, dtype=float)
        if perturb.shape != (len(nodes), 2):
            raise ValueError("imperfections must be (N, 2)")
        if not np.isfinite(perturb).all():
            raise ValueError('Nonfinite imperfection')
        nodes[:, :2] += perturb
    if len(plastic) < 2:
        raise ValueError("Expected original calibrated true stress / plastic strain table")
    ymin, ymax = np.min(nodes[:, 1]), np.max(nodes[:, 1])
    zmin, zmax = np.min(nodes[:, 2]), np.max(nodes[:, 2])
    xmin, xmax = np.min(nodes[:, 0]), np.max(nodes[:, 0])
    height = ymax - ymin
    if height <= 0:
        raise ValueError("Zero specimen height")
    if zmax <= zmin:
        raise ValueError('Zero specimen depth')
    if elements.ndim != 2 or elements.shape[1] not in (6,8) or elements.min() < 0 or elements.max() >= len(nodes):
        raise ValueError('Invalid C3D6/C3D8 connectivity')
    element_type = config.get('element_type', 'C3D8R' if elements.shape[1]==8 else 'C3D6')
    if element_type not in ('C3D6','C3D8','C3D8R') or ((element_type=='C3D6') != (elements.shape[1]==6)):
        raise ValueError('element_type and connectivity disagree')
    # Planens span the full specimen in X and depth in Z.
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
    lines.append(f"*ELEMENT,TYPE={element_type},ELSET=CORE")
    for eid, el in enumerate(elements,1):
        lines.append(f"{eid}, " + ", ".join(str(int(i)+1) for i in el))
    lower_el, upper_el = len(elements)+1, len(elements)+2
    lines += ["*ELEMENT,TYPE=C3D8,ELSET=PLATENS",
              f"{lower_el}, " + ", ".join(map(str,lower_ids)),
              f"{upper_el}, " + ", ".join(map(str,upper_ids))]

    faces = exterior_faces(elements)
    tol = max(1e-7, float(config.get("contact_face_tol", 1e-6)))
    core_bottom, core_top = [], []
    for eid, face, verts in faces:
        y = nodes[list(verts), 1]
        if np.all(np.abs(y-ymin) <= tol):
            core_bottom.append((eid,face))
        elif np.all(np.abs(y-ymax) <= tol):
            core_top.append((eid,face))
    _surface(lines, "CORE_BOTTOM", core_bottom)
    _surface(lines, "CORE_TOP", core_top)
    _surface(lines, "BOTTOM_PLATEN_CONTACT", [(lower_el,"S5")])
    _surface(lines, "TOP_PLATEN_CONTACT", [(upper_el,"S3")])
    wall_pairs = wall_contact_patches(elements, reference_nodes) if self_contact else []
    for i,(slave,master) in enumerate(wall_pairs):
        _surface(lines,f'WALL_SLAVE_{i}',[(eid,face) for eid,face,_ in slave])
        _surface(lines,f'WALL_MASTER_{i}',[(eid,face) for eid,face,_ in master])
    rigid_platens=bool(config.get('rigid_platen_kinematics', False))
    lower_load=lower_ids if rigid_platens else [lower_ids[i] for i in (0,1,4,5)]
    upper_load=upper_ids if rigid_platens else [upper_ids[i] for i in (2,3,6,7)]
    _nset(lines,"BOTTOM_PLATEN_NODES",lower_load)
    _nset(lines,"TOP_PLATEN_NODES",upper_load)
    z_faces = (list(range(1,ncore+1)) if config.get('plane_strain',True) else
               [i+1 for i,p in enumerate(nodes) if abs(p[2]-zmin)<tol or abs(p[2]-zmax)<tol])
    _nset(lines,"CORE_PLANE_STRAIN",z_faces)
    # Remove the otherwise free in-plane rigid translation at initially unloaded
    # frictional contact. A centreline roller leaves end faces free to slide.
    bottom_nodes = np.flatnonzero(np.abs(nodes[:, 1] - ymin) <= tol)
    centre = (xmin + xmax) / 2
    anchor = bottom_nodes[np.argmin(np.abs(nodes[bottom_nodes, 0] - centre))]
    _nset(lines, "CORE_X_ANCHOR", [anchor + 1])
    times = [strain / float(config['compression_strain']) for strain in (.12, .35, .65)
             if strain <= float(config['compression_strain'])]
    target_output_only = bool(config.get('target_output_only', False))
    interval=config.get('output_interval')
    if interval is not None:
        interval=float(interval)
        if not np.isfinite(interval) or not 0<interval<=1:
            raise ValueError('output_interval must lie in (0,1]')
        times=sorted(set(times+[float(v) for v in np.arange(interval,1,interval)]+[1.]))
    if times and target_output_only:
        lines += ['*TIME POINTS,NAME=TARGET_FRAMES']
        lines += [', '.join(f'{t:.12g}' for t in times[i:i+8]) for i in range(0,len(times),8)]

    damage = config.get('ductile_damage')
    damage_info = None
    material_data = None
    if damage:
        cards, material_data, damage_info = material_cards(nodes, elements, config, plastic)
        lines += cards
    else:
        lines += ["*MATERIAL,NAME=ALSI10MG", "*ELASTIC",
                  f"{config['young']:.12g}, {config['poisson']:.12g}", "*PLASTIC"]
        lines.extend(f"{s:.12g}, {p:.12g}" for s,p in plastic)
        lines += ["*SOLID SECTION,ELSET=CORE,MATERIAL=ALSI10MG"]
    lines += ["*MATERIAL,NAME=STEEL","*ELASTIC",
              f"{config.get('steel_young',210000.0):.12g}, {config.get('steel_poisson',0.3):.12g}",
              "*SOLID SECTION,ELSET=PLATENS,MATERIAL=STEEL",
              "*SURFACE INTERACTION,NAME=SLIDING",
              "*SURFACE BEHAVIOR,PRESSURE-OVERCLOSURE=LINEAR",
              str(config.get("normal_penalty", 1.0e5)),
              ]
    if friction > 0:
        lines += ["*FRICTION", f"{friction:.12g}, {config.get('stick_slope', 100000.0):.12g}"]
    lines += [f"*CONTACT PAIR,INTERACTION=SLIDING,TYPE={contact_type}",
              "CORE_BOTTOM,BOTTOM_PLATEN_CONTACT",
              f"*CONTACT PAIR,INTERACTION=SLIDING,TYPE={contact_type}",
              "CORE_TOP,TOP_PLATEN_CONTACT"]

    for i in range(len(wall_pairs)):
        lines += [f'*CONTACT PAIR,INTERACTION=SLIDING,TYPE={contact_type}',
                  f'WALL_SLAVE_{i},WALL_MASTER_{i}']

    lines += ["*BOUNDARY","BOTTOM_PLATEN_NODES,1,3,0.",
              "TOP_PLATEN_NODES,1,1,0.","TOP_PLATEN_NODES,3,3,0.",
              "CORE_PLANE_STRAIN,3,3,0.", "CORE_X_ANCHOR,1,1,0.",
              "*STEP,NLGEOM=YES,INC=100000",
              "*STATIC",
              f"{config.get('initial_increment',1e-5):.12g}, 1., {config.get('minimum_increment',1e-6):.12g}, {config.get('maximum_increment',.002):.12g}",
              "*BOUNDARY",
              f"TOP_PLATEN_NODES,2,2,{-float(config['compression_strain'])*height:.12g}",
              ("*NODE FILE,GLOBAL=YES,TIME POINTS=TARGET_FRAMES" if times and target_output_only else "*NODE FILE,GLOBAL=YES"),"U,RF",
              "*EL FILE",("S,E,SDV" if damage else "S,E,PEEQ"),
              "*NODE PRINT,NSET=TOP_PLATEN_NODES,TOTALS=ONLY","RF",
              "*EL PRINT,ELSET=CORE", ("S,SDV" if damage else "S,PEEQ"),
              "*CONTACT FILE", "CDIS,CSTR",
              "*END STEP"]
    path = Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text("\n".join(lines)+"\n",encoding="utf-8")
    if material_data:
        path.with_suffix('.material.dat').write_text(material_data,encoding='utf8')
    return {"path":str(path),"core_nodes":ncore,"core_elements":len(elements), "element_type":element_type,
            "top_faces":len(core_top),"bottom_faces":len(core_bottom),
            "initial_height":height,"nominal_area":float(config.get("reference_area",76.0*25.0)),
            "self_contact_experimental":bool(self_contact), "wall_contact_pairs":len(wall_pairs), "rigid_platen_kinematics":rigid_platens,
            "plane_strain":bool(config.get("plane_strain",True)), "x_anchor_node":int(anchor+1),
            "target_strains":[.12,.35,.65], "damage_model":("compiled finite-strain J2 crack-band damage; retained residual elements" if damage else "J2 plasticity only"), "damage_parameters":damage_info}
