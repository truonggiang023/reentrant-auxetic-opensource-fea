import numpy as np
import pytest
from reentrant_fea.native_rp_history import rp_reaction_history, VARIABLES


def fixture():
    h={'time':np.array([0.,1.,2.])}
    for node,offset,dy,impulse in [(25,69,[0.,-.1,-.2],[0.,-20.,-50.]),(26,75,[0.,.02,.04],[0.,19.,48.])]:
        for j in range(6):
            h[f'Loading face nodal histories {node} load-node-{node} var {offset+j}']=np.asarray(dy if j==1 else impulse if j==4 else [0.,0.,0.])
    m=dict(platen_reference_points=True,loading_node_variables=VARIABLES.copy(),
           loading_face_nodes=[25],bottom_support_nodes=[26],height_mm=2.,reference_area_mm2=10.,target_nominal_strain=.12)
    return h,m


def test_signed_rp_impulses_and_observed_relative_displacement():
    h,m=fixture();r=rp_reaction_history(h,m)
    assert r['prescribed_motion_reached']
    assert r['last_observed_nominal_strain']==pytest.approx(.12)
    rows=r['reaction_history']
    assert [x['strain'] for x in rows]==pytest.approx([.03,.09])
    assert [x['nominal_stress_MPa'] for x in rows]==[2.,3.]
    assert [x['reaction_y_N'] for x in rows]==[-20.,-30.]
    assert [x['top_bottom_reaction_sum_y_N'] for x in rows]==[-1.,-1.]


def test_no_amplitude_fallback_when_motion_or_channels_are_missing():
    h,m=fixture();m['target_nominal_strain']=.2
    assert not rp_reaction_history(h,m)['prescribed_motion_reached']
    del h[next(k for k in h if 'var 73' in k)]
    with pytest.raises(ValueError,match='history'):
        rp_reaction_history(h,m)


@pytest.mark.parametrize('field,value',[('platen_reference_points',False),('loading_face_nodes',[25,27]),('height_mm',0),('loading_node_variables',['DY','REACY'])])
def test_unsupported_history_conventions_rejected(field,value):
    h,m=fixture();m[field]=value
    with pytest.raises(ValueError):rp_reaction_history(h,m)
