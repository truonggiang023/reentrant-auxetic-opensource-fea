"""Compare native prescribed-face and RBE2 RP histories (interval mean forces)."""
import argparse,json,re
from pathlib import Path
import numpy as np
from reentrant_fea.radioss_results import read_history_csv,history_energy_screen


def extract(prefix):
    prefix=Path(prefix)
    meta=json.loads(prefix.with_suffix('.radioss.json').read_text())
    h=read_history_csv(str(prefix)+'T01.csv');t=h['time']
    def channels(node):
        keys=[k for k in h if re.search(r'load-node-'+str(node)+r'\s',k)]
        if len(keys)!=6:raise ValueError('Require six declared nodal channels')
        return [h[k] for k in keys]
    # Converter orders nodes numerically; match explicit labels, never group order.
    top=[channels(n) for n in meta['loading_face_nodes']]
    bottom=[channels(n) for n in meta['bottom_support_nodes']]
    strain=-(np.mean([x[1] for x in top],axis=0)-np.mean([x[1] for x in bottom],axis=0))/meta['height_mm']
    force=np.diff(np.sum([x[4] for x in top],axis=0))/np.diff(t)
    bottom_force=np.diff(np.sum([x[4] for x in bottom],axis=0))/np.diff(t)
    rows=np.column_stack(((t[1:]+t[:-1])/2,(strain[1:]+strain[:-1])/2,force,bottom_force,-force/meta['reference_area_mm2']))
    np.savetxt(str(prefix)+'.rp_curve.csv',rows,delimiter=',',header='time_s,strain,top_force_N,bottom_force_N,stress_MPa',comments='')
    return rows,dict(final_sample_strain=float(strain[-1]),energy=history_energy_screen(h,meta['duration_seconds']))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('direct');p.add_argument('rp');p.add_argument('--output',required=True);a=p.parse_args()
    d,dm=extract(a.direct);r,rm=extract(a.rp)
    if d.shape!=r.shape or not np.allclose(d[:,0],r[:,0]):raise ValueError('Comparison requires matching time grids')
    result=dict(direct=dm,rp=rm,max_force_difference_N=float(np.max(abs(d[:,2]-r[:,2]))),relative_force_rms=float(np.linalg.norm(d[:,2]-r[:,2])/np.linalg.norm(d[:,2])),scope='native coupon numerical coupling check; no Abaqus or lattice equivalence established')
    Path(a.output).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
