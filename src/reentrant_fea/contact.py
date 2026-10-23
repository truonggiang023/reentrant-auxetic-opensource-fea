"""Topology-based exterior wedge faces; supports contact-deck experiments.

CalculiX self-contact on a single deformable surface is version-dependent.
This module DOES NOT claim validated self-contact.
"""
from collections import defaultdict

# C3D6 node order: bottom triangle 0,1,2; top triangle 3,4,5.
WEDGE_FACES = {"S1": (0, 1, 2), "S2": (3, 5, 4),
               "S3": (0, 3, 4, 1), "S4": (1, 4, 5, 2),
               "S5": (2, 5, 3, 0)}

def exterior_faces(elements):
    """Return [(1-based element id, face label, vertex tuple)] of free faces."""
    counts = defaultdict(list)
    for eid, el in enumerate(elements, 1):
        for label, loc in WEDGE_FACES.items():
            verts = tuple(int(el[i]) for i in loc)
            counts[tuple(sorted(verts))].append((eid, label, verts))
    if any(len(v) > 2 for v in counts.values()):
        raise ValueError('Non-manifold mesh: more than two elements share a face')
    return [v[0] for v in counts.values() if len(v) == 1]

def lateral_exterior_faces(elements, nodes, y_min, y_max, tol=1e-7):
    """Exclude end-plane faces from potential self-contact surface."""
    result=[]
    for eid,label,verts in exterior_faces(elements):
        ys=[nodes[i,1] for i in verts]
        if all(abs(y-y_min)<tol for y in ys) or all(abs(y-y_max)<tol for y in ys):
            continue
        result.append((eid,label,verts))
    return result
