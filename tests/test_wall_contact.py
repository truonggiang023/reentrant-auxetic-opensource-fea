import os
import subprocess
import numpy as np
import pytest
from reentrant_fea.contact import wall_contact_patches
from reentrant_fea.paper_contact_solver import write_paper_contact_deck
from reentrant_fea.ccx_results import validate_run
from reentrant_fea.frd import fields, contact_envelope


def mesh():
    a=np.array([[0,0,0],[.3,0,0],[.3,1,0],[0,1,0],
                [0,0,1],[.3,0,1],[.3,1,1],[0,1,1]],float)
    n=np.vstack([a,a+[.301,0,0]])
    e=np.array([[0,1,2,4,5,6],[0,2,3,4,6,7],
                [8,9,10,12,13,14],[8,10,11,12,14,15]])
    return n,e


def test_wall_candidates_exclude_adjacency_and_backfaces():
    n,e=mesh();pairs=wall_contact_patches(e,n)
    assert len(pairs)==1  # Only the two faces bordering the initial void.
    for slave,master in pairs:
        sn=set(x for _,_,v in slave for x in v)
        mn=set(x for _,_,v in master for x in v)
        assert not sn&mn
        assert {eid for eid,_,_ in slave}.isdisjoint({eid for eid,_,_ in master})


def test_closed_wall_gap_transmits_contact(tmp_path):
    ccx=os.environ.get('CCX_EXECUTABLE')
    if not ccx:pytest.skip('CCX_EXECUTABLE not set')
    n,e=mesh();results=[]
    for contact in (False,True):
        c={'young':69000,'poisson':.3,'compression_strain':.02,
           'normal_penalty':1e6,'stick_slope':1e3,'friction':.35,
           'initial_increment':.001,'maximum_increment':.005}
        deck=tmp_path/f'wall{int(contact)}.inp'
        write_paper_contact_deck(deck,n,e,c,[(225,0),(300,1)],self_contact=contact)
        # Two disconnected bodies need one lateral anchor apiece. Anchoring
        # their outer corners forces expansion to close the intervening gap.
        deck.write_text(deck.read_text().replace('CORE_X_ANCHOR,1,1,0.',
                                                  '1,1,1,0.\n10,1,1,0.'))
        log=deck.with_suffix('.log')
        with log.open('w') as stream:
            subprocess.run([ccx,'-i',deck.stem],cwd=tmp_path,stdout=stream,
                           stderr=subprocess.STDOUT,check=True,timeout=90)
        validate_run(deck,log)
        for t,name,data in fields(deck.with_suffix('.frd')):
            if name=='DISP':u=data
        gap=.001+u[9][0]-u[2][0]
        from reentrant_fea.rp_curve import ccx_rp_history
        rows, _ = ccx_rp_history(deck)
        forces = [row['reaction_y_N'] for row in rows]
        results.append((gap,forces,contact_envelope(deck.with_suffix('.frd'))))
    assert results[0][0]<-.001  # Closure exceeds the initial 0.001 mm void.
    assert results[1][0]>results[0][0]/2
    # Pre-contact reactions agree: no artificial spring through the solid.
    assert results[0][1][0]==pytest.approx(results[1][1][0],rel=1e-5)
    assert results[1][2]['max_reported_penetration_mm']<.002
    assert results[1][2]['max_reported_pressure_MPa']>0


def test_smooth_imperfection_preserves_reference_patch_partition(tmp_path):
    base=np.array([[x,y,z] for z in [0,1] for y in [0,.5,1] for x in [0,.3]],float)
    def idx(x,y,z):return z*6+y*2+x
    wedge=[]
    for y in [0,1]:
        a,b,c,d=idx(0,y,0),idx(1,y,0),idx(1,y+1,0),idx(0,y+1,0)
        wedge.extend([[a,b,c,a+6,b+6,c+6],[a,c,d,a+6,c+6,d+6]])
    n=np.vstack([base,base+[.301,0,0]])
    e=np.vstack([wedge,np.array(wedge)+12])
    perturb=np.zeros((len(n),2));perturb[n[:,1]==.5,0]=.003
    cfg={'young':69000,'poisson':.3,'compression_strain':.02}
    normal=write_paper_contact_deck(tmp_path/'normal.inp',n,e,cfg,[(225,0),(300,1)],self_contact=True)
    perturbed=write_paper_contact_deck(tmp_path/'perturbed.inp',n,e,cfg,[(225,0),(300,1)],self_contact=True,imperfections=perturb)
    assert normal['wall_contact_pairs']==perturbed['wall_contact_pairs']==1
    coordinate=next(s for s in (tmp_path/'perturbed.inp').read_text().splitlines() if s.startswith('3,'))
    assert float(coordinate.split(',')[1])==.003
