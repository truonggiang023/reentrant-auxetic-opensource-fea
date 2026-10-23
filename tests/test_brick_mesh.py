import numpy as np
from shapely.geometry import box
from reentrant_fea.mesh import mesh_bricks
from reentrant_fea.contact import exterior_faces
from reentrant_fea.paper_contact_solver import write_paper_contact_deck
from reentrant_fea.native_damage import material_cards


def test_conforming_positive_brick_mesh_and_material_lengths(tmp_path):
    nodes,elements,top,bottom=mesh_bricks(box(0,0,1,1),2,.5,1)
    assert elements.shape[1]==8 and len(top)>0 and len(bottom)>0
    volume=0.
    for el in elements:
        x=nodes[el]
        assert np.linalg.det(np.column_stack([x[1]-x[0],x[3]-x[0],x[4]-x[0]]))>0
        q=x[:4,:2]
        area=.5*abs(np.sum(q[:,0]*np.roll(q[:,1],-1)-q[:,1]*np.roll(q[:,0],-1)))
        volume+=area*2
    assert abs(volume-2)<1e-10
    faces=exterior_faces(elements)
    assert faces and all(len(v)==4 for _,_,v in faces)
    cfg={'young':69000.,'poisson':.3,'compression_strain':.12,'element_type':'C3D8R',
         'ductile_damage':{'fracture_locus':[[0,.08],[.333,.047],[1,.018]],
                           'failure_displacement':.001,'provenance':'Synthetic test only'}}
    deck=tmp_path/'bricks.inp'
    info=write_paper_contact_deck(deck,nodes,elements,cfg,[(225,0),(300,1)])
    assert '*ELEMENT,TYPE=C3D8R,ELSET=CORE' in deck.read_text()
    assert '*PLASTIC' not in deck.read_text()
    assert '*USER MATERIAL' in deck.read_text()
    assert info['damage_parameters']['length_min']>0
