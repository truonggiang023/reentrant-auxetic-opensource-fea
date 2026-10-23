"""Compiled constitutive checks: Python only drives the C material tests."""
import ctypes
import os
from pathlib import Path
import subprocess
import numpy as np
import pytest
from reentrant_fea.native_damage import material_cards
from reentrant_fea.paper_contact_solver import write_paper_contact_deck
from reentrant_fea.ccx_results import validate_run

PROPS=np.array([69000,.3,.001,1e-4,1,2,3,-1/3,
                225,0,300,1,0,.08,.333,.047,1,.018,0.],dtype='d')

@pytest.fixture(scope='module')
def kernel(tmp_path_factory):
    so=tmp_path_factory.mktemp('native')/'material.so'
    source=Path(__file__).resolve().parents[1]/'native/ductile.c'
    subprocess.run(['gcc','-O2','-shared','-fPIC',str(source),
                    '/lib/x86_64-linux-gnu/liblapack.so.3','-lm','-o',str(so)],check=True)
    lib=ctypes.CDLL(str(so));f=lib.reentrant_response_dt
    ptr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    f.argtypes=[ptr,ptr,ptr,ctypes.c_int,ctypes.c_double,ptr,ptr];f.restype=ctypes.c_int
    def response(F,old=None,props=PROPS,elastic=0,dt=1.):
        C=np.ascontiguousarray(F.T@F);old=np.zeros(10) if old is None else old.copy()
        stress=np.zeros(6);state=np.zeros(10)
        assert f(C,old,props,elastic,dt,stress,state)==0
        return stress,state
    return response


def test_objectivity_and_plastic_volume(kernel):
    F=np.diag([1.11,1/np.sqrt(1.11),1/np.sqrt(1.11)])
    a=.73;R=np.array([[np.cos(a),-np.sin(a),0],[np.sin(a),np.cos(a),0],[0,0,1]])
    s,x=kernel(F);sr,xr=kernel(R@F)
    np.testing.assert_allclose(s,sr,atol=1e-8)
    np.testing.assert_allclose(x,xr,atol=1e-10)
    cp=np.array([[x[0],x[3],x[4]],[x[3],x[1],x[5]],[x[4],x[5],x[2]]])
    assert np.linalg.det(cp)==pytest.approx(1.,abs=1e-10)
    assert x[6]>0 and x[9]>0


def test_commit_unloading_and_compression_cutoff(kernel):
    old=np.zeros(10)
    F=np.diag([1.08,1/np.sqrt(1.08),1/np.sqrt(1.08)])
    s,x=kernel(F,old)
    np.testing.assert_array_equal(old,np.zeros(10))
    # Re-evaluation from the same committed state must not accumulate history.
    np.testing.assert_allclose(kernel(F,old)[1],x)
    np.testing.assert_allclose(kernel(F,x)[1],x,atol=1e-10)
    compressed=np.diag([.85,.98,.98])
    _,xc=kernel(compressed)
    assert xc[6]>0 and xc[7]==0 and xc[9]==0


def test_crack_band_softening_displacement(kernel):
    for length in [.5,1.,2.]:
        p=PROPS.copy();p[4]=length
        old=np.zeros(10);init_alpha=None
        for alpha in np.linspace(0,.12,2401):
            # Isochoric uniaxial stretch, positive shear-dominated triaxiality=0.
            F=np.diag([np.exp(alpha),np.exp(-alpha/2),np.exp(-alpha/2)])
            _,old=kernel(F,old,p)
            if old[7]>=1 and init_alpha is None:init_alpha=old[6]
            if old[9]>=.999:
                assert old[8]==pytest.approx(.001,abs=length*.00006)
                break
        else:
            pytest.fail('Damage did not reach residual stiffness')
        assert init_alpha is not None


@pytest.mark.parametrize('full_table',[False,True])
def test_solver_side_softening(tmp_path,full_table):
    ccx=os.environ.get('CCX_DAMAGE_EXECUTABLE')
    if not ccx:pytest.skip('CCX_DAMAGE_EXECUTABLE not set')
    source=Path(__file__).resolve().parents[1]/'tests/fixtures/damage_coupon.inp'
    deck=tmp_path/'coupon.inp';text=source.read_text()
    env=os.environ.copy()
    if full_table:
        from reentrant_fea.plastic import from_reference
        n=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],
                    [0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)
        cfg={'young':69000,'poisson':.3,'compression_strain':.1,
             'ductile_damage':{'fracture_locus':[[0,.08],[.333,.047],[1,.018]],
                               'failure_displacement':.001,'provenance':'Test using full repository hardening; assumed locus'}}
        root=Path(__file__).resolve().parents[1]
        cards,data,_=material_cards(n,np.array([[0,1,2,3,4,5,6,7]]),cfg,
                                    from_reference(root/'reference_setup_numeric_source.py'))
        assert all(int(line.split('=')[-1])<=21 for line in cards if line.startswith('*USER MATERIAL'))
        text=text[:text.index('*MATERIAL')]+'\n'.join(cards)+'\n'+text[text.index('*BOUNDARY'):]
        sidecar=tmp_path/'material.txt';sidecar.write_text(data)
        env['REENTRANT_MATERIAL_DATA']=str(sidecar)
    deck.write_text(text)
    log=tmp_path/'coupon.log'
    with log.open('w') as stream:
        subprocess.run([ccx,'-i','coupon'],cwd=tmp_path,stdout=stream,
                       stderr=subprocess.STDOUT,check=True,timeout=90,env=env)
    validate_run(deck,log)
    assert 'REENTRANT_MATERIAL_ABI=3' in log.read_text()
    lines=deck.with_suffix('.dat').read_text().splitlines()
    force=[float(lines[i+2].split()[1]) for i,s in enumerate(lines) if 'total force' in s]
    assert max(force)>200 and 0<force[-1]<.1
    sdv=[list(map(float,s.split()[2:])) for s in lines if len(s.split())==12 and s.split()[0]=='1']
    assert sdv[-1][6]>.08 and sdv[-1][7]>1 and sdv[-1][9]==pytest.approx(.9999)


def test_viscous_damage_has_explicit_time_scale(kernel):
    p=PROPS.copy();p[-1]=.005
    F=np.diag([1.11,1/np.sqrt(1.11),1/np.sqrt(1.11)])
    _,x=kernel(F,props=p,dt=.001)
    assert x[9]==pytest.approx(.9999/6)
    _,frozen=kernel(F,x,p,elastic=1,dt=.001)
    assert frozen[9]==x[9]
    _,relaxed=kernel(F,x,p,dt=.001)
    assert relaxed[9]==pytest.approx((.005*x[9]+.001*.9999)/.006)


def test_softening_energy_per_band_area(kernel):
    # Perfect plasticity removes hardening-induced energy differences, so the
    # linear displacement law has Gf = sigma_y * uf / 2 for all band lengths.
    for length in [.5,1.,2.]:
        p=PROPS.copy();p[4]=length;p[10]=225.
        alpha0=.08+225/(3*(69000/(2*1.3)))
        old=None;energy=0.;previous_q=225.;previous_u=0.
        for u in np.linspace(0,.001,501):
            strain=alpha0+u/length
            F=np.diag([np.exp(strain),np.exp(-strain/2),np.exp(-strain/2)])
            stress,old=kernel(F,old,p)
            S=np.array([[stress[0],stress[3],stress[4]],
                        [stress[3],stress[1],stress[5]],
                        [stress[4],stress[5],stress[2]]])
            sigma=F@S@F.T
            dev=sigma-np.eye(3)*np.trace(sigma)/3
            q=np.sqrt(1.5*np.sum(dev*dev))
            energy+=(q+previous_q)/2*(u-previous_u)
            previous_q,previous_u=q,u
        assert energy==pytest.approx(225*.001/2,rel=.005)


def test_reduced_brick_hourglass_matches_elastic_reference(tmp_path):
    ccx=os.environ.get('CCX_DAMAGE_EXECUTABLE')
    if not ccx:pytest.skip('CCX_DAMAGE_EXECUTABLE not set')
    fixture=(Path(__file__).resolve().parents[1]/'tests/fixtures/damage_coupon.inp').read_text()
    head=fixture.split('*NSET')[0].replace('TYPE=C3D8,','TYPE=C3D8R,')
    native=fixture[fixture.index('*MATERIAL'):fixture.index('*BOUNDARY')]
    linear='*MATERIAL,NAME=REFERENCE\n*ELASTIC\n69000,.3\n*SOLID SECTION,ELSET=CORE,MATERIAL=REFERENCE\n'
    signs=[1,-1,1,-1,-1,1,-1,1];projections=[]
    for name,material in [('reference',linear),('native',native)]:
        deck=tmp_path/f'{name}.inp'
        text=head+'*NSET,NSET=ALL,GENERATE\n1,8,1\n'+material
        text+='*BOUNDARY\nALL,3,3,0\n1,2,2,0\n2,2,2,0\n3,2,2,0\n4,2,2,0\n5,2,2,0\n6,2,2,0\n7,2,2,0\n'
        text+='*STEP,NLGEOM\n*STATIC\n1,1\n*BOUNDARY\n'
        text+=''.join(f'{i+1},1,1,{s*.001}\n' for i,s in enumerate(signs))
        text+='*NODE PRINT,NSET=ALL\nRF\n*NODE FILE\nU,RF\n*END STEP\n'
        deck.write_text(text);log=deck.with_suffix('.log')
        with log.open('w') as out:
            subprocess.run([ccx,'-i',name],cwd=tmp_path,stdout=out,
                           stderr=subprocess.STDOUT,check=True,timeout=90)
        validate_run(deck,log)
        rows=[s.split() for s in deck.with_suffix('.dat').read_text().splitlines()]
        force={int(p[0]):float(p[1]) for p in rows if len(p)==4 and p[0].isdigit()}
        projections.append(sum(signs[i-1]*force[i] for i in range(1,9)))
    assert projections[0]>1.
    assert projections[1]==pytest.approx(projections[0],rel=1e-5)


def test_unloading_keeps_plastic_state_and_damage(kernel):
    strain=.0835;delta=.0005
    def F(e):return np.diag([np.exp(e),np.exp(-e/2),np.exp(-e/2)])
    def q(stress,f):
        S=np.array([[stress[0],stress[3],stress[4]],
                    [stress[3],stress[1],stress[5]],
                    [stress[4],stress[5],stress[2]]])
        sigma=f@S@f.T
        dev=sigma-np.eye(3)*np.trace(sigma)/3
        return np.sqrt(1.5*np.sum(dev*dev))
    peak,old=kernel(F(strain))
    assert 0<old[9]<.9999
    unloaded,state=kernel(F(strain-delta),old)
    np.testing.assert_allclose(state,old,atol=1e-10)
    G=69000/(2*1.3)
    assert q(peak,F(strain))-q(unloaded,F(strain-delta))==pytest.approx(3*G*delta*(1-old[9]),rel=1e-5)


def test_two_brick_bending_matches_elastic_solver(tmp_path):
    """Exercise free bending DOFs, full table transport and hourglass tangent."""
    ccx=os.environ.get('CCX_DAMAGE_EXECUTABLE')
    if not ccx:pytest.skip('CCX_DAMAGE_EXECUTABLE not set')
    from reentrant_fea.plastic import from_reference
    n=np.array([[x,y,z] for z in [0,1] for y in [0,1] for x in [0,1,2]],float)
    def idx(x,y,z):return z*6+y*3+x
    e=np.array([[idx(x,0,0),idx(x+1,0,0),idx(x+1,1,0),idx(x,1,0),
                 idx(x,0,1),idx(x+1,0,1),idx(x+1,1,1),idx(x,1,1)] for x in [0,1]])
    cfg={'young':69000,'poisson':.3,'compression_strain':.1,
         'ductile_damage':{'fracture_locus':[[0,.08],[.333,.047],[1,.018]],
                           'failure_displacement':.001,'provenance':'Synthetic bending test'}}
    root=Path(__file__).resolve().parents[1]
    cards,data,_=material_cards(n,e,cfg,from_reference(root/'reference_setup_numeric_source.py'))
    sidecar=tmp_path/'tables.txt';sidecar.write_text(data)
    forces=[]
    for name,mat in [('native',cards),('reference',['*MATERIAL,NAME=REF','*ELASTIC','69000,.3',
                                                  '*SOLID SECTION,ELSET=CORE,MATERIAL=REF'])]:
        lines=['*NODE']+[f'{i+1},'+','.join(map(str,p)) for i,p in enumerate(n)]
        lines+=['*ELEMENT,TYPE=C3D8R,ELSET=CORE']+[f'{i+1},'+','.join(str(v+1) for v in p) for i,p in enumerate(e)]
        lines+=['*NSET,NSET=ALL,GENERATE','1,12,1','*NSET,NSET=LEFT',
                ','.join(str(i+1) for i,p in enumerate(n) if p[0]==0),'*NSET,NSET=RIGHT',
                ','.join(str(i+1) for i,p in enumerate(n) if p[0]==2)]
        lines+=mat+['*BOUNDARY','ALL,3,3,0','LEFT,1,2,0','*STEP,NLGEOM','*STATIC','.1,1',
                    '*BOUNDARY','RIGHT,2,2,.00001','*NODE PRINT,NSET=RIGHT,TOTALS=ONLY','RF',
                    '*NODE FILE','U,RF','*END STEP']
        deck=tmp_path/f'{name}.inp';deck.write_text('\n'.join(lines)+'\n')
        log=deck.with_suffix('.log')
        with log.open('w') as stream:
            subprocess.run([ccx,'-i',name],cwd=tmp_path,stdout=stream,stderr=subprocess.STDOUT,
                           env={**os.environ,'REENTRANT_MATERIAL_DATA':str(sidecar)},check=True,timeout=90)
        validate_run(deck,log)
        rows=deck.with_suffix('.dat').read_text().splitlines()
        forces.append([float(rows[i+2].split()[1]) for i,s in enumerate(rows) if 'total force' in s][-1])
    assert forces[0]>0
    assert forces[0]==pytest.approx(forces[1],rel=1e-5)
