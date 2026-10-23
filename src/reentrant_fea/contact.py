"""Topology-based exterior wedge faces; supports contact-deck experiments.

CalculiX self-contact on a single deformable surface is version-dependent.
This module DOES NOT claim validated self-contact.
"""
from collections import defaultdict

# C3D6 node order: bottom triangle 0,1,2; top triangle 3,4,5.
WEDGE_FACES = {"S1": (0, 1, 2), "S2": (3, 5, 4),
               "S3": (0, 3, 4, 1), "S4": (1, 4, 5, 2),
               "S5": (2, 5, 3, 0)}

HEX_FACES = {'S1':(0,1,2,3),'S2':(4,7,6,5),
             'S3':(0,4,5,1),'S4':(1,5,6,2),
             'S5':(2,6,7,3),'S6':(3,7,4,0)}

def exterior_faces(elements):
    """Return [(1-based element id, face label, vertex tuple)] of free faces."""
    counts = defaultdict(list)
    for eid, el in enumerate(elements, 1):
        faces=WEDGE_FACES if len(el)==6 else HEX_FACES if len(el)==8 else None
        if faces is None:raise ValueError('Only C3D6 and C3D8(R) faces are supported')
        for label, loc in faces.items():
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


def wall_contact_patches(elements, nodes, tol=1e-6):
    """Coplanar exterior wall patches for an extruded XY wedge mesh.

    End faces and out-of-plane caps are excluded. Initially back-facing pairs
    and segments passing entirely through solid are excluded. This fixed
    candidate set can miss new contacts after large rotations. Master collections omit the
    slave patch and every patch sharing any node with it. Each unordered patch
    pair occurs once, so there is no same-surface or duplicated contact pair.
    Mesh-specific topology is preserved; this is not Abaqus general contact.
    """
    import numpy as np
    ymin, ymax = nodes[:,1].min(), nodes[:,1].max()
    planes = defaultdict(list)
    for eid, face, verts in lateral_exterior_faces(elements,nodes,ymin,ymax,tol):
        xyz=nodes[list(verts)]
        if len(verts)!=4 or np.ptp(xyz[:,2])<=tol:
            continue
        normal=np.cross(xyz[1]-xyz[0],xyz[2]-xyz[0])
        normal/=np.linalg.norm(normal)
        if np.dot(normal,xyz.mean(axis=0)-nodes[elements[eid-1]].mean(axis=0))<0:
            normal=-normal
        key=tuple(np.round(np.r_[normal[:2],np.dot(normal,xyz.mean(axis=0))]/tol).astype(int))
        planes[key].append((eid,face,verts))
    patches=[]
    # Split disconnected parts of each plane; a union by plane alone could
    # incorrectly suppress contact between distant but initially collinear walls.
    for entries in planes.values():
        remaining=set(range(len(entries)))
        by_node=defaultdict(set)
        for i,entry in enumerate(entries):
            for node in entry[2]:by_node[node].add(i)
        while remaining:
            todo=[remaining.pop()];group=[]
            while todo:
                i=todo.pop();group.append(entries[i])
                for node in entries[i][2]:
                    new=by_node[node]&remaining
                    remaining-=new;todo.extend(new)
            patches.append(group)
    node_sets=[set(v for _,_,verts in p for v in verts) for p in patches]
    # Opposite faces of the SAME solid wall are not a contact gap. Otherwise
    # a one-sided master normal can see the other face through the material
    # as initial penetration. Exclude pairs whose connecting segment lies
    # entirely within the initial solid, in addition to shared-node neighbors.
    from shapely.geometry import Polygon, LineString
    from shapely.ops import unary_union
    zmin=nodes[:,2].min()
    triangles=[Polygon(nodes[el[:len(el)//2],:2]) for el in elements
               if np.all(np.abs(nodes[el[:len(el)//2],2]-zmin)<=tol)]
    solid=unary_union(triangles)
    solid=solid.buffer(tol)
    centres=[np.mean([nodes[list(entry[2])].mean(axis=0) for entry in p],axis=0)[:2]
             for p in patches]
    normals=[]
    for patch in patches:
        eid,_,verts=patch[0];xyz=nodes[list(verts)]
        normal=np.cross(xyz[1]-xyz[0],xyz[2]-xyz[0]);normal/=np.linalg.norm(normal)
        if np.dot(normal,xyz.mean(axis=0)-nodes[elements[eid-1]].mean(axis=0))<0:normal=-normal
        normals.append(normal[:2])
    pairs=[]
    for i,patch in enumerate(patches):
        master=[entry for j in range(i+1,len(patches))
                if not (node_sets[i]&node_sets[j])
                and not solid.covers(LineString([centres[i],centres[j]]))
                and np.dot(normals[i],centres[j]-centres[i])>=-tol
                and np.dot(normals[j],centres[i]-centres[j])>=-tol
                for entry in patches[j]]
        if master:pairs.append((patch,master))
    return pairs
