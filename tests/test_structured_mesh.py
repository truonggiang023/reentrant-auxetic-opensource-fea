from collections import Counter
from copy import deepcopy
import numpy as np
import pytest
from shapely.geometry import box
from reentrant_fea.geometry import core_polygon
from reentrant_fea.mesh import mesh_bricks


def config():
    return dict(l1=12.,l2=12.,bevel=5.66,wall=.3,rows=2,cols=2,
                height=15.7,depth=25.,lambda_m=1.)


def test_wall_blocks_preserve_geometry_and_conforming_boundaries():
    c=config();original=deepcopy(c);p=core_polygon(c)
    n,e,top,bottom=mesh_bricks(p,25.,.8,1,strategy='structured_walls',config=c)
    assert c==original
    assert len(top)>0 and len(bottom)>0 and e.shape[1]==8
    v=n[e[:,:4],:2];a=np.roll(v,-1,axis=1)-v;b=-np.roll(a,1,axis=1)
    jac=a[:,:,0]*b[:,:,1]-a[:,:,1]*b[:,:,0]
    assert jac.min()>0
    area=.5*np.sum(v[:,:,0]*np.roll(v[:,:,1],-1,axis=1)-v[:,:,1]*np.roll(v[:,:,0],-1,axis=1))
    assert area.sum()==pytest.approx(p.area,rel=1e-8)
    edges=Counter(tuple(sorted((int(i),int(j)))) for f in e[:,:4] for i,j in zip(f,np.roll(f,-1)))
    assert max(edges.values())==2
    perimeter=sum(np.linalg.norm(n[i,:2]-n[j,:2]) for (i,j),count in edges.items() if count==1)
    assert perimeter==pytest.approx(p.length,rel=1e-8)
    assert np.ptp(n[:,2])==25.
    assert np.allclose(n[e[:,4:],:2],n[e[:,:4],:2])


def test_structured_mesh_requires_matching_reentrant_config():
    with pytest.raises(ValueError,match='requires'):
        mesh_bricks(box(0,0,1,1),25,.8,1,strategy='structured_walls')
    with pytest.raises(ValueError,match='does not match'):
        mesh_bricks(box(0,0,1,1),25,.8,1,strategy='structured_walls',config=config())
