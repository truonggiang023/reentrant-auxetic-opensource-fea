"""Geometry-preserving wall-aligned blocks for the reentrant centerline model.

Straight member interiors are mapped rectangles. Junctions retain the exact
union polygon and use unstructured recombination. One conforming refinement
makes both regions all-quad. This is a discretization choice, not a material
or geometry fit. Elastic qualification is separate for each specimen/mesh.
config['structured_refinement_levels'] (default 1) adds uniform refinements;
mapped straight walls have 2**levels elements across their thickness.
"""
import math
from collections import Counter
import numpy as np
from shapely.geometry import LineString
from shapely.ops import unary_union
from .geometry import core_polygon, parameters


def _segments(c):
    h, shift, pitch = parameters(c)
    segments = []
    def append(points):
        segments.extend(LineString([a,b]) for a,b in zip(points,points[1:]))
    for j in range(c['cols']):
        x = j*pitch
        pts = [(x,0.)]
        for i in range(c['rows']):
            pts.extend([(x+shift,(i+.5)*h),(x,(i+1)*h)])
        pts.append((x+c['l1'],c['height']))
        for i in reversed(range(c['rows'])):
            pts.extend([(x+c['l1']-shift,(i+.5)*h),(x+c['l1'],i*h)])
        pts.append(pts[0]); append(pts)
        for i in range(1,c['rows']):
            append([(x,i*h),(x+c['l1'],i*h)])
        if j < c['cols']-1:
            for i in range(c['rows']):
                y = (i+.5)*h
                append([(x+c['l1']-shift,y),((j+1)*pitch+shift,y)])
    return unary_union(segments)


def mesh_structured_walls(poly, depth, mesh_size, layers, config):
    if config is None:
        raise ValueError('structured_walls requires the reentrant geometry config')
    if (not np.isfinite([depth,mesh_size]).all() or min(depth,mesh_size)<=0
            or isinstance(layers,bool) or not isinstance(layers,int) or layers<1):
        raise ValueError('Positive depth/mesh size and integer depth layers required')
    refinement_levels = config.get('structured_refinement_levels', 1)
    if (isinstance(refinement_levels, bool) or not isinstance(refinement_levels, int)
            or refinement_levels < 1):
        raise ValueError('structured_refinement_levels must be a positive integer')
    expected = core_polygon(config)
    if poly.symmetric_difference(expected).area > 1e-8*expected.area:
        raise ValueError('structured_walls polygon does not match its geometry config')
    width = float(config['wall'])
    strips = []
    for line in _segments(config).geoms:
        coords = list(line.coords)
        for first,second in zip(coords,coords[1:]):
            a,b = np.asarray(first),np.asarray(second)
            length = np.linalg.norm(b-a); direction = (b-a)/length
            trim = 1.5*width
            if length <= 2*trim:
                continue
            rectangle = LineString([a+trim*direction,b-trim*direction]).buffer(width/2,cap_style=2)
            if rectangle.difference(poly).area > 1e-9:
                raise ValueError('Mapped wall block extends outside core geometry')
            strips.append(rectangle)
    merged = unary_union(strips)
    if abs(sum(s.area for s in strips)-merged.area)>1e-8:
        raise ValueError('Mapped wall blocks overlap')
    rest = poly.difference(merged)
    patches = strips + list(rest.geoms if rest.geom_type=='MultiPolygon' else [rest])
    import gmsh
    gmsh.initialize()
    try:
        gmsh.model.add('ReEntrantWallBlocks')
        points,curves,structured = {},{},[]
        def point(xy):
            key = tuple(round(float(v),9) for v in xy)
            if key not in points:
                points[key] = gmsh.model.geo.addPoint(*key,0,mesh_size)
            return points[key]
        def loop(coords):
            vertices = list(coords)
            ids = [point(v) for v in vertices[:-1]]
            edges,lengths = [],[]
            for a,b,xy0,xy1 in zip(ids,ids[1:]+ids[:1],vertices[:-1],vertices[1:]):
                if a==b: continue
                key = tuple(sorted((a,b)))
                if key not in curves:
                    curves[key] = (gmsh.model.geo.addLine(a,b),a,b)
                tag,origin,_ = curves[key]
                edges.append(tag if a==origin else -tag)
                lengths.append(float(np.linalg.norm(np.asarray(xy0)-xy1)))
            return gmsh.model.geo.addCurveLoop(edges),edges,lengths
        for i,patch in enumerate(patches):
            outer,edges,lengths = loop(patch.exterior.coords)
            holes = [loop(r.coords)[0] for r in patch.interiors]
            surface = gmsh.model.geo.addPlaneSurface([outer]+holes)
            if i < len(strips):
                if len(edges)!=4: raise ValueError('Expected rectangular member block')
                structured.append((surface,edges,lengths))
        gmsh.model.geo.synchronize()
        requests = {}
        for surface,edges,lengths in structured:
            for edge,length in zip(edges,lengths):
                count = 2 if abs(length-width)<1e-6 else max(2,int(math.ceil(length/mesh_size))+1)
                if abs(edge) in requests and requests[abs(edge)]!=count:
                    raise ValueError('Conflicting shared-edge discretization')
                requests[abs(edge)] = count
            gmsh.model.mesh.setTransfiniteSurface(surface)
            gmsh.model.mesh.setRecombine(2,surface)
        for edge,count in requests.items():
            gmsh.model.mesh.setTransfiniteCurve(edge,count)
        for name,value in [('Mesh.MeshSizeMin',width/2),('Mesh.MeshSizeMax',mesh_size),
                           ('Mesh.Algorithm',8),('Mesh.RecombineAll',1),('Mesh.RecombinationAlgorithm',1)]:
            gmsh.option.setNumber(name,value)
        gmsh.model.mesh.generate(2)
        gmsh.option.setNumber('Mesh.SubdivisionAlgorithm',1)
        # The first subdivision makes the junctions all-quad. Additional
        # uniform subdivisions refine both along and across every wall while
        # retaining the same conforming interfaces and exact polygon.
        for _ in range(refinement_levels):
            gmsh.model.mesh.refine()
        tags,xyz,_ = gmsh.model.mesh.getNodes()
        coordinates = {int(tag):np.asarray(xyz[3*i:3*i+2]) for i,tag in enumerate(tags)}
        types,_,groups = gmsh.model.mesh.getElements(2)
        quads = []
        for kind,flat in zip(types,groups):
            if gmsh.model.mesh.getElementProperties(kind)[0]!='Quadrilateral 4':
                raise RuntimeError('All-quad block meshing failed')
            quads.extend(np.asarray(flat,int).reshape(-1,4))
        used = sorted({int(v) for q in quads for v in q})
        mapping = {tag:i for i,tag in enumerate(used)}
        xy = np.asarray([coordinates[tag] for tag in used])
        connectivity = []
        for q in quads:
            ids = [mapping[int(v)] for v in q]; v = xy[ids]
            area = .5*np.sum(v[:,0]*np.roll(v[:,1],-1)-v[:,1]*np.roll(v[:,0],-1))
            if area<0: ids = [ids[0],ids[3],ids[2],ids[1]]
            connectivity.append(ids)
        q = np.asarray(connectivity,int); v = xy[q]
        edge = np.roll(v,-1,axis=1)-v; previous = -np.roll(edge,1,axis=1)
        jacobian = edge[:,:,0]*previous[:,:,1]-edge[:,:,1]*previous[:,:,0]
        area = .5*np.sum(v[:,:,0]*np.roll(v[:,:,1],-1,axis=1)-v[:,:,1]*np.roll(v[:,:,0],-1,axis=1))
        boundary = Counter(tuple(sorted((int(a),int(b)))) for face in q for a,b in zip(face,np.roll(face,-1)))
        perimeter = sum(np.linalg.norm(xy[a]-xy[b]) for (a,b),count in boundary.items() if count==1)
        if (np.min(jacobian)<=0 or max(boundary.values())>2
                or not np.isclose(area.sum(),poly.area,rtol=1e-8)
                or not np.isclose(perimeter,poly.length,rtol=1e-8)):
            raise ValueError('Invalid or nonconforming wall-block mesh geometry')
        nn = len(xy)
        nodes = np.asarray([[x,y,depth*k/layers] for k in range(layers+1) for x,y in xy])
        elements = np.asarray([[v+k*nn for v in face]+[v+(k+1)*nn for v in face]
                               for k in range(layers) for face in connectivity],int)
        bottom = np.flatnonzero(abs(nodes[:,1]-nodes[:,1].min())<1e-7)
        top = np.flatnonzero(abs(nodes[:,1]-nodes[:,1].max())<1e-7)
        return nodes,elements,top,bottom
    finally:
        gmsh.finalize()
