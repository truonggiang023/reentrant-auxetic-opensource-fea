"""Native coupling coupon diagnostic. Run from repo root with PYTHONPATH=.:src.
Source /workspace/opencourant/env.sh first. Outputs require a fresh directory.
"""
import numpy as np
from reentrant_fea.radioss import write_radioss_deck
from run_radioss import run_generated
n=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)
for rp in (False,True):
 c=dict(young=69000.,poisson=.3,compression_strain=.12,reference_area=1.,platen_thickness=.2,platen_margin=.1,friction=.35,platen_reference_points=rp,platen_contact_surface='inner_y')
 m=write_radioss_deck('results/native_rp/'+('rp' if rp else 'direct'),n,np.array([[0,1,2,3,4,5,6,7]]),c,[(225.,0.),(300.,.2)],duration=.001)
 print(run_generated(m,timeout=120),flush=True)
