import numpy as np
from reentrant_fea.solver import write_deck

def test_deck(tmp_path):
    nodes=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,0,1],[0,1,1]],float)
    c={'young':69000,'poisson':.3,'compression_strain':.6}
    path=write_deck(tmp_path/'x.inp',nodes,np.array([[0,1,2,3,4,5]]),np.array([2,5]),np.array([0,1,3,4]),c,[(200,0),(210,.1)])
    s=path.read_text()
    assert '*ELEMENT,TYPE=C3D6' in s and '*PLASTIC' in s and 'TOP,2,2,-0.6' in s
