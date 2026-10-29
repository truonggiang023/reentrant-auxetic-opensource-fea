"""Actual solver regressions; enable with CCX_EXECUTABLE (stock CCX 2.20)."""
import os
from pathlib import Path
import subprocess
import numpy as np
import pytest
from reentrant_fea.paper_contact_solver import write_paper_contact_deck
from reentrant_fea.ccx_results import validate_run, integration_history, initiation_indices
from reentrant_fea.plastic import validate_material

N=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)
E=np.array([[0,1,2,4,5,6],[0,2,3,4,6,7]])

@pytest.mark.parametrize('strain',[.12,.35,.65])
@pytest.mark.parametrize('contact_type',['SURFACE TO SURFACE','MORTAR'])
def test_contact_compression_solver(tmp_path, strain, contact_type):
    ccx=os.environ.get('CCX_EXECUTABLE')
    if not ccx:
        pytest.skip('CCX_EXECUTABLE not set: solver regression not executed')
    cfg={'young':69000.,'poisson':.3,'compression_strain':strain,
         'initial_increment':.001,'maximum_increment':.02, 'contact_type':contact_type}
    deck=tmp_path/'slab.inp'
    write_paper_contact_deck(deck,N,E,cfg,[(225,0),(300,2.)])
    log=tmp_path/'slab.log'
    with log.open('w') as stream:
        run = subprocess.run([ccx,'-i','slab'],cwd=tmp_path,stdout=stream,
                             stderr=subprocess.STDOUT,timeout=90)
    if contact_type == 'MORTAR' and run.returncode:
        assert any(message in log.read_text() for message in
                   ('*ERROR: too many cutbacks', '*ERROR: increment size smaller than minimum'))
        pytest.xfail('CCX 2.20 mortar contact loses convergence on this coarse slab; penalty regression passes')
    run.check_returncode()
    inc=validate_run(deck,log)
    assert inc[-1][2] == pytest.approx(1.)
    history=list(integration_history(deck.with_suffix('.dat')))
    assert len(history)>=4 and max(row[3] for row in history)>0
    # Synthetic locus for plumbing only, not an AlSi10Mg calibration.
    indices=initiation_indices(deck.with_suffix('.dat'),[(-100.,.5),(100.,.5)])
    assert max(row[3] for row in indices)>0
    from reentrant_fea.rp_curve import ccx_rp_history
    rows, area = ccx_rp_history(deck)
    assert rows[-1]['strain'] == pytest.approx(strain)
    assert rows[-1]['reaction_y_N'] < -1.
    assert rows[-1]['nominal_stress_MPa'] == pytest.approx(-rows[-1]['reaction_y_N']/area)
    assert abs(rows[-1]['force_balance_y_N']) < abs(rows[-1]['reaction_y_N'])*1e-4

@pytest.mark.parametrize('table', [[(200,.01),(210,.1)],[(200,0),(190,.1)],[(200,0),(float('nan'),.1)]])
def test_invalid_hardening(table):
    with pytest.raises(ValueError):
        validate_material({'young':69000,'poisson':.3,'compression_strain':.65},table)


def test_zero_exit_is_not_convergence(tmp_path):
    log=tmp_path/'x.log'; log.write_text('Job finished\n')
    (tmp_path/'x.sta').write_text('1 1 1U 20 0.0 0.0 0.1\n')
    with pytest.raises(RuntimeError,match='final step time'):
        validate_run(tmp_path/'x.inp',log)


def test_damage_requires_calibrated_locus(tmp_path):
    with pytest.raises(ValueError, match='calibrated'):
        initiation_indices(tmp_path/'missing.dat',[(.333,.047)])


def test_damage_preserves_integration_point_identity(tmp_path):
    dat=tmp_path/'history.dat'
    dat.write_text('''stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sxz,syz) for set CORE and time 1.0
1 1 100 0 0 0 0 0
1 2 -100 0 0 0 0 0
equivalent plastic strain (elem, integ.pnt.,pe)for set CORE and time 1.0
1 1 .1
1 2 .2
''')
    rows=initiation_indices(dat,[(-1.,.5),(1.,.5)])
    assert rows[0][3] == pytest.approx(.2)
    assert rows[1][3] == pytest.approx(.4)
    with pytest.raises(ValueError, match='outside calibrated'):
        initiation_indices(dat,[(0.,.047),(.5,.05)])


def test_rigid_grip_rotation_nodes_with_solver(tmp_path):
    from reentrant_fea.advanced_solver import write_advanced_deck
    ccx=os.environ.get('CCX_EXECUTABLE')
    if not ccx:
        pytest.skip('CCX_EXECUTABLE not set')
    deck=tmp_path/'grip.inp'
    write_advanced_deck(deck,N,E,[2,3,6,7],[0,1,4,5],
                        {'young':69000.,'poisson':.3,'compression_strain':.01},
                        [(225.,0.),(300.,.2)])
    log=tmp_path/'grip.log'
    with log.open('w') as stream:
        subprocess.run([ccx,'-i','grip'],cwd=tmp_path,stdout=stream,
                       stderr=subprocess.STDOUT,check=True,timeout=90)
    validate_run(deck,log)
