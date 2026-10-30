"""Gmsh triangulation of planar solid, extruded in Python to 6-node wedge elements."""
import numpy as np

def _loop(gmsh, coords, size):
    coords=list(coords)[:-1]
    ids=[gmsh.model.geo.addPoint(float(x),float(y),0,size) for x,y in coords]
    lines=[gmsh.model.geo.addLine(ids[i],ids[(i+1)%len(ids)]) for i in range(len(ids))]
    return gmsh.model.geo.addCurveLoop(lines)

def mesh_wedges(poly, depth, mesh_size, layers):
    import gmsh
    if layers < 1: raise ValueError('depth_layers must be >=1')
    gmsh.initialize()
    try:
        gmsh.model.add('ReEntrant')
        polys=list(poly.geoms) if poly.geom_type=='MultiPolygon' else [poly]
        for p in polys:
            outer=_loop(gmsh,p.exterior.coords,mesh_size)
            holes=[_loop(gmsh,ring.coords,mesh_size) for ring in p.interiors]
            gmsh.model.geo.addPlaneSurface([outer]+holes)
        gmsh.model.geo.synchronize()
        gmsh.option.setNumber('Mesh.MeshSizeMin',mesh_size)
        gmsh.option.setNumber('Mesh.MeshSizeMax',mesh_size)
        gmsh.option.setNumber('Mesh.Algorithm',6)
        gmsh.option.setNumber('Mesh.RecombineAll',0)
        gmsh.option.setNumber('Mesh.SubdivisionAlgorithm',0)
        gmsh.model.mesh.generate(2)
        tags, xyz, _=gmsh.model.mesh.getNodes()
        pts={int(t):np.array(xyz[3*i:3*i+3]) for i,t in enumerate(tags)}
        types,_,node_groups=gmsh.model.mesh.getElements(2)
        tris=[]
        for typ, flat in zip(types,node_groups):
            if gmsh.model.mesh.getElementProperties(typ)[0].startswith('Triangle'):
                n_per=gmsh.model.mesh.getElementProperties(typ)[3]
                if n_per!=3: raise RuntimeError('Only linear triangles supported')
                tris.extend(np.asarray(flat,dtype=int).reshape(-1,3).tolist())
        if not tris: raise RuntimeError('No triangular elements generated')
        used=sorted({n for tri in tris for n in tri})
        idmap={tag:i for i,tag in enumerate(used)}
        xy=np.array([pts[tag][:2] for tag in used])
        nn=len(xy)
        nodes=np.array([[x,y,depth*k/layers] for k in range(layers+1) for x,y in xy])
        elements=[]
        for k in range(layers):
            for tri in tris:
                a,b,c=[idmap[n] for n in tri]
                ab = xy[b] - xy[a]
                ac = xy[c] - xy[a]

                cross_z = ab[0] * ac[1] - ab[1] * ac[0]

                if cross_z < 0:
                    b, c = c, b
                elements.append([a+k*nn,b+k*nn,c+k*nn,a+(k+1)*nn,b+(k+1)*nn,c+(k+1)*nn])
        elements=np.asarray(elements,dtype=int)
        y0,y1=nodes[:,1].min(),nodes[:,1].max()
        tol=max(1e-5,mesh_size*0.03)
        bottom=np.flatnonzero(abs(nodes[:,1]-y0)<tol)
        top=np.flatnonzero(abs(nodes[:,1]-y1)<tol)
        if not len(top) or not len(bottom): raise RuntimeError('Missing top/bottom boundary nodes')
        return nodes,elements,top,bottom
    finally:
        gmsh.finalize()


def mesh_bricks_subdivision(poly, depth, mesh_size, layers=1):
    """Conforming all-quad subdivision of a linear triangular Gmsh mesh.

    Each triangle becomes three quads using shared edge midpoints. Extrude to
    C3D8(R) connectivity. This changes discretization, not the input geometry.
    """
    if layers < 1 or depth <= 0:
        raise ValueError('Positive depth and depth_layers required')
    base, wedges, _, _ = mesh_wedges(poly, depth, mesh_size, 1)
    nn=len(base)//2
    xy=base[:nn,:2].tolist();midpoints={};quads=[]
    def midpoint(a,b):
        key=tuple(sorted((int(a),int(b))))
        if key not in midpoints:
            midpoints[key]=len(xy)
            xy.append(((np.asarray(xy[a])+np.asarray(xy[b]))/2).tolist())
        return midpoints[key]
    for wedge in wedges:
        a,b,c=map(int,wedge[:3]);ab=midpoint(a,b);bc=midpoint(b,c);ca=midpoint(c,a)
        centre=len(xy);xy.append(np.mean([xy[a],xy[b],xy[c]],axis=0).tolist())
        quads.extend(((a,ab,centre,ca),(b,bc,centre,ab),(c,ca,centre,bc)))
    nn=len(xy)
    nodes=np.array([[x,y,depth*k/layers] for k in range(layers+1) for x,y in xy])
    elements=np.array([[i+k*nn for i in q]+[i+(k+1)*nn for i in q]
                       for k in range(layers) for q in quads],dtype=int)
    y0,y1=nodes[:,1].min(),nodes[:,1].max()
    bottom=np.flatnonzero(abs(nodes[:,1]-y0)<1e-7)
    top=np.flatnonzero(abs(nodes[:,1]-y1)<1e-7)
    return nodes,elements,top,bottom


def mesh_bricks(poly, depth, mesh_size, layers=1, strategy='recombine', config=None):
    """All-quad blossom recombination followed by one conforming refinement.

    Refinement splits residual triangles into quads and gives two elements
    through initially single-quad wall widths. It avoids subdividing every
    skinny triangle into three skew quads, which can overexcite hourglass modes.
    """
    if strategy=='structured_walls':
        from .structured_mesh import mesh_structured_walls
        return mesh_structured_walls(poly,depth,mesh_size,layers,config)
    if strategy=='triangle_subdivision':
        return mesh_bricks_subdivision(poly,depth,mesh_size,layers)
    if strategy!='recombine':raise ValueError('Unknown brick mesh strategy')
    if layers<1 or depth<=0:raise ValueError('Positive depth and layers required')
    import gmsh
    gmsh.initialize()
    try:
        gmsh.model.add('ReEntrantQuads')
        polys=list(poly.geoms) if poly.geom_type=='MultiPolygon' else [poly]
        for p in polys:
            outer=_loop(gmsh,p.exterior.coords,mesh_size)
            holes=[_loop(gmsh,ring.coords,mesh_size) for ring in p.interiors]
            gmsh.model.geo.addPlaneSurface([outer]+holes)
        gmsh.model.geo.synchronize()
        for key,value in [('Mesh.MeshSizeMin',mesh_size),('Mesh.MeshSizeMax',mesh_size),
                          ('Mesh.Algorithm',8),('Mesh.RecombineAll',1),
                          ('Mesh.RecombinationAlgorithm',1),('Mesh.SubdivisionAlgorithm',0)]:
            gmsh.option.setNumber(key,value)
        gmsh.model.mesh.generate(2)
        gmsh.option.setNumber('Mesh.SubdivisionAlgorithm',1)
        gmsh.model.mesh.refine()
        tags,xyz,_=gmsh.model.mesh.getNodes()
        coords={int(t):np.array(xyz[3*i:3*i+2]) for i,t in enumerate(tags)}
        types,_,groups=gmsh.model.mesh.getElements(2)
        quads=[]
        for typ,flat in zip(types,groups):
            props=gmsh.model.mesh.getElementProperties(typ)
            if props[0]!='Quadrilateral 4':raise RuntimeError('All-quad generation failed: '+props[0])
            quads.extend(np.asarray(flat,dtype=int).reshape(-1,4))
        used=sorted({int(v) for q in quads for v in q})
        mapping={tag:i for i,tag in enumerate(used)}
        xy=np.array([coords[tag] for tag in used]);quad_ids=[]
        for q in quads:
            ids=[mapping[int(tag)] for tag in q];v=xy[ids]
            if np.sum(v[:,0]*np.roll(v[:,1],-1)-v[:,1]*np.roll(v[:,0],-1))<0:
                ids=[ids[0],ids[3],ids[2],ids[1]]
            quad_ids.append(ids)
        nn=len(xy)
        nodes=np.array([[x,y,depth*k/layers] for k in range(layers+1) for x,y in xy])
        elements=np.array([[v+k*nn for v in q]+[v+(k+1)*nn for v in q]
                           for k in range(layers) for q in quad_ids],dtype=int)
        ymin,ymax=nodes[:,1].min(),nodes[:,1].max()
        return nodes,elements,np.flatnonzero(abs(nodes[:,1]-ymax)<1e-7),np.flatnonzero(abs(nodes[:,1]-ymin)<1e-7)
    finally:
        gmsh.finalize()
