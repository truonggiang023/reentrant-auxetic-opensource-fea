"""Standalone geometric model: XY honeycomb, Z extrusion. mm units."""
import math
from shapely.geometry import LineString
from shapely.ops import unary_union

def parameters(c):
    h = c['height'] / c['rows']
    shift = c['lambda_m'] * math.sqrt(c['bevel']**2 - (h / 2)**2)
    pitch = c['l1'] + c['l2'] - 2*shift
    if h <= 0 or c['bevel'] <= h/2 or pitch <= 0 or c['wall'] <= 0:
        raise ValueError('Invalid geometric parameters')
    return h, shift, pitch

def core_polygon(c):
    h, shift, pitch = parameters(c)
    strips = []
    t = c['wall']
    for j in range(c['cols']):
        x = j*pitch
        # One continuous closed zigzag ring per column.
        pts = [(x, 0.0)]
        for i in range(c['rows']):
            pts.extend([(x+shift, (i+0.5)*h), (x, (i+1)*h)])
        pts.append((x+c['l1'], c['height']))
        for i in reversed(range(c['rows'])):
            pts.extend([(x+c['l1']-shift, (i+0.5)*h), (x+c['l1'], i*h)])
        pts.append(pts[0])
        strips.append(LineString(pts).buffer(t/2, cap_style=2, join_style=2))
        for i in range(1, c['rows']):
            strips.append(LineString([(x,i*h),(x+c['l1'],i*h)]).buffer(t/2, cap_style=2, join_style=2))
        if j < c['cols']-1:
            for i in range(c['rows']):
                y=(i+0.5)*h
                strips.append(LineString([(x+c['l1']-shift,y),((j+1)*pitch+shift,y)]).buffer(t/2,cap_style=2,join_style=2))
    geom = unary_union(strips)
    if not geom.is_valid or geom.is_empty:
        raise ValueError('Invalid unioned geometry')
    return geom
