import numpy as np
from scripts.probe_matched_elastic_lattice import make_deck


def test_tiny_negative_coordinates_fit_calculix_20_character_fields(tmp_path):
    n=np.array([[x,y,z] for z in (0.,25.) for y in (-.15,-1.82340254006874e-13)
                for x in (.45,.846428571429636)])
    # Standard C3D8 connectivity for the coordinate enumeration above.
    m=dict(nodes=[[i+1,*v] for i,v in enumerate(n)],elements=[[1,1,2,4,3,5,6,8,7]],
           boundary=dict(all_z_fixed=list(range(1,9)),bottom_y_fixed=[1,2,5,6],
                         top_y_displaced=[3,4,7,8],one_bottom_x_pin=1))
    p=tmp_path/'tiny.inp';make_deck(m,p,'C3D8I')
    node_lines=p.read_text().split('*NODE\n',1)[1].split('*ELEMENT',1)[0].strip().splitlines()
    numeric=[[f.strip() for f in line.split(',')[1:]] for line in node_lines]
    assert max(len(f) for row in numeric for f in row)<=20
    np.testing.assert_allclose(np.array(numeric,float),n,rtol=1e-11,atol=0.)
