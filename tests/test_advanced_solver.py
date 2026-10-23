import numpy as np
import pytest
from reentrant_fea.contact import exterior_faces,lateral_exterior_faces
from reentrant_fea.advanced_solver import write_advanced_deck

N=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,0,1],[0,1,1]],float)
E=np.array([[0,1,2,3,4,5]])
C={'young':69000,'poisson':.3,'compression_strain':.6}
P=[(200,0),(210,.1)]

def test_wedge_exterior_faces():
    assert len(exterior_faces(E))==5
    assert len(lateral_exterior_faces(E,N,0,1))==4

def test_bonded_rigid_grips(tmp_path):
    p=write_advanced_deck(tmp_path/'a.inp',N,E,[2,5],[0,1,3,4],C,P)
    s=p.read_text()
    assert '*RIGID BODY,NSET=TOP_GRIP' in s
    assert '*RIGID BODY,NSET=BOTTOM_GRIP' in s
    assert 'TOP_RP,2,2,-0.6' in s
    assert '*CONTACT PAIR' not in s

def test_damage_not_misrepresented(tmp_path):
    with pytest.raises(NotImplementedError):
        write_advanced_deck(tmp_path/'a.inp',N,E,[2,5],[0,1,3,4],C,P,ductile_damage=True)

def test_self_contact_opt_in(tmp_path):
    with pytest.raises(NotImplementedError):
        write_advanced_deck(tmp_path/'a.inp',N,E,[2,5],[0,1,3,4],C,P,self_contact=True)
    c={**C,'experimental_self_contact_ack':True}
    with pytest.raises(NotImplementedError, match='Same-surface'):
        write_advanced_deck(tmp_path/'b.inp',N,E,[2,5],[0,1,3,4],c,P,self_contact=True)
