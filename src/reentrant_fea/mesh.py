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
