"""Native free-lateral coupons: true-flow input and the located fracture point.

The existing bonded-end erosion coupon deliberately has constrained strain.
These tests instead leave lateral faces free and remove only rigid modes.
They qualify the material input convention, not the missing fracture locus,
mesh objectivity, or experimental honeycomb agreement.
"""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from reentrant_fea.plastic import from_reference
from reentrant_fea.radioss import write_radioss_deck


NODES = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
ELEMENTS = np.array([list(range(8))])
VARIABLES = ['OFF', 'SX', 'SY', 'SZ', 'SXY', 'SYZ', 'SXZ', 'PLAS']


def _uniaxial_run(directory, name, plastic, damage=None, *,
                  lateral_free=True, compression=False, end_strain=.08):
    executables = [shutil.which(executable) for executable in
                   ('starter_linux64_gf', 'engine_linux64_gf', 'th_to_csv_linux64_gf')]
    if not all(executables):
        pytest.skip('Source the installed OpenCourant env.sh to run native uniaxial qualification')
    config = {'young': 69000., 'poisson': .3, 'compression_strain': end_strain}
    if damage:
        config['ductile_damage'] = damage
    meta = write_radioss_deck(directory/name, NODES, ELEMENTS, config, plastic,
                              duration=.001, tensile_coupon=True)
    # The writer's bonded-end coupon is useful for erosion, but is not a
    # uniaxial-stress material test. Here bottom Y is fixed, top Y prescribed,
    # node 1 fixes XZ, and node 2 fixes Z. All other lateral DOFs are free.
    deck = Path(meta['starter'])
    original = deck.read_text().splitlines()
    lines, index = [], 0
    while index < len(original):
        if lateral_free and original[index] in ('/BCS/1', '/BCS/2'):
            if original[index] == '/BCS/1':
                lines += ['/BCS/1', 'Bottom Y restraint only',
                          f'{"010":>6} 000'+f'{0:10d}{1:10d}']
            index += 3
        else:
            lines.append(original[index])
            index += 1
    assert lines.pop() == '/END'
    if lateral_free:
        lines += ['/GRNOD/NODE/4', 'One XZ anchor', f'{1:10d}',
                  '/BCS/4', 'Node 1 XZ only', f'{"101":>6} 000'+f'{0:10d}{4:10d}',
                  '/GRNOD/NODE/5', 'One Z anchor', f'{2:10d}',
                  '/BCS/5', 'Node 2 Z only', f'{"001":>6} 000'+f'{0:10d}{5:10d}']
    if compression:
        motion = lines.index('/IMPDISP/1')
        lines[motion+3] = ''.join(f'{value:20.12g}' for value in (1., -end_strain, 0., 0.))
    lines += ['/TH/BRIC/4', 'Uniaxial material point',
              ''.join(f'{variable:<10}' for variable in VARIABLES),
              f'{1:10d}          solid1', '/END']
    deck.write_text('\n'.join(lines)+'\n')
    # The uf=.001 softening interval is short. Dense histories are necessary
    # to resolve its stress loss; coarse lattice output can miss it entirely.
    engine_deck = Path(meta['engine'])
    engine_lines = engine_deck.read_text().splitlines()
    interval = meta['duration_seconds']/2000
    engine_lines[engine_lines.index('/TFILE/4')+1] = f'{interval:20.12g}'
    engine_deck.write_text('\n'.join(engine_lines)+'\n')
    meta['time_history_interval_seconds'] = interval
    meta['coupon_boundary_conditions'] = ('free lateral faces with minimal rigid-body restraints' if lateral_free
                                          else 'constrained lateral strain; diagnostic of high stress triaxiality only')
    meta['coupon_loading_sign'] = -1 if compression else 1
    meta['sha256'][deck.name] = hashlib.sha256(deck.read_bytes()).hexdigest()
    meta['sha256'][engine_deck.name] = hashlib.sha256(engine_deck.read_bytes()).hexdigest()
    (directory/(name+'.radioss.json')).write_text(json.dumps(meta, indent=2)+'\n')
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '2'
    for executable, args in [(executables[0], ['-i', deck.name, '-np', '1']),
                             (executables[1], ['-i', Path(meta['engine']).name]),
                             (executables[2], [name+'T01'])]:
        with (directory/(name+'-'+Path(executable).name+'.log')).open('w') as stream:
            subprocess.run([executable, *args], cwd=directory, env=env, stdout=stream,
                           stderr=subprocess.STDOUT, check=True, timeout=40)
    starter = deck.with_suffix('.out').read_text()
    engine = (directory/(name+'_0001.out')).read_text()
    assert '0 ERROR(S)' in starter and '0 WARNING(S)' in starter
    assert 'NORMAL TERMINATION' in starter and 'NORMAL TERMINATION' in engine
    with (directory/(name+'T01.csv')).open() as stream:
        reader = csv.reader(stream)
        headers = next(reader)
        values = np.array([[float(value) for value in row] for row in reader])
    columns = [i for i, label in enumerate(headers) if label.startswith('Uniaxial material point')]
    assert len(columns) == len(VARIABLES) and len(values) > 200
    point = values[:, columns]
    assert np.isfinite(point).all()
    stress, peeq = point[:, 1:7], point[:, 7]
    mean = stress[:, :3].mean(axis=1)
    mises = np.sqrt(1.5*np.sum((stress[:, :3]-mean[:, None])**2, axis=1)
                    + 3*np.sum(stress[:, 3:]**2, axis=1))
    eta = np.divide(mean, mises, out=np.zeros_like(mises), where=mises > 1e-9)
    return dict(active=point[:, 0], stress=stress, peeq=peeq, mises=mises,
                eta=eta, values=values, headers=headers, engine=engine)


def test_full_table_uses_true_flow_stress_and_clamps_beyond_final_data(tmp_path):
    root = Path(__file__).resolve().parents[1]
    plastic = np.asarray(from_reference(root/'reference_setup_numeric_source.py'))
    result = _uniaxial_run(tmp_path, 'hardening', plastic)
    selected = result['peeq'] > .005
    assert selected.sum() > 50
    np.testing.assert_allclose(result['eta'][selected], 1/3, atol=5e-4)
    expected = np.interp(result['peeq'], plastic[:, 1], plastic[:, 0])
    np.testing.assert_allclose(result['mises'][selected], expected[selected], rtol=5e-4)
    # This specifically catches the native FUNCT linear-extrapolation defect:
    # the source table ends at .03814; the coupon must actually exceed it.
    beyond = result['peeq'] > plastic[-1, 1]+.01
    assert beyond.sum() > 20
    np.testing.assert_allclose(result['mises'][beyond], plastic[-1, 0], rtol=5e-4)


def test_uniaxial_anchor_initiates_then_softens_and_erodes(tmp_path):
    plastic = [(225., 0.), (300., 1.)]
    damage = {'fracture_locus': [[0., .08], [.333, .047], [1., .018]],
              'failure_displacement': .001}
    result = _uniaxial_run(tmp_path, 'anchor', plastic, damage)
    peeq, active = result['peeq'], result['active'] > .5
    preinit = (peeq > .005) & (peeq < .04) & active
    np.testing.assert_allclose(result['eta'][preinit], 1/3, atol=5e-4)
    # Here linear hardening stays within its supplied interval; q is Cauchy
    # flow stress and PLAS is equivalent plastic strain, not engineering strain.
    flow = 225.+75.*peeq
    ratio = result['mises']/flow
    np.testing.assert_allclose(ratio[preinit], 1., rtol=5e-4)
    softening = active & (peeq > .04) & (ratio > .05) & (ratio < .95)
    assert softening.sum() >= 2
    # For this unit cube, native initial characteristic length is 1 mm.
    # Linear uf=.001 gives D=(PEEQ-PEEQ_start)/.001. Infer initiation from
    # actual stress loss; the sample interval need not land exactly on onset.
    observed_onset = np.median(peeq[softening]-.001*(1-ratio[softening]))
    expected_onset = np.interp(1/3, [0., .333, 1.], [.08, .047, .018])
    assert observed_onset == pytest.approx(expected_onset, abs=5e-5)
    assert peeq[-1] == pytest.approx(expected_onset+.001, abs=5e-5)
    assert result['active'][0] == 1 and result['active'][-1] == 0
    assert 'RUPTURE OF SOLID ELEMENT' in result['engine']


@pytest.mark.parametrize('case', ['high_eta', 'negative_eta'])
def test_fracture_endpoints_are_clamped_outside_supplied_eta_range(tmp_path, case):
    """No unmeasured slope extrapolation may create negative fracture strains.

    The negative-eta case also makes the native difference explicit: INIEVO
    accumulates compression damage here, whereas the CCX example has a cutoff.
    This verifies a declared boundary convention, not missing fracture data.
    """
    damage = {'fracture_locus': [[0., .08], [.333, .047], [1., .018]],
              'failure_displacement': .001}
    high = case == 'high_eta'
    result = _uniaxial_run(tmp_path, case, [(225., 0.), (300., 1.)], damage,
                           lateral_free=not high, compression=not high,
                           end_strain=.08 if high else .15)
    peeq, active = result['peeq'], result['active'] > .5
    ratio = result['mises']/(225.+75.*peeq)
    before = active & (peeq > .005) & (peeq < (.015 if high else .06))
    assert before.sum() > 10
    if high:
        # Fully restrained lateral strain starts above eta=1 and climbs; its
        # old linear extension would make epsilon_f negative by eta~1.414.
        assert np.min(result['eta'][before]) > 1.4
    else:
        np.testing.assert_allclose(result['eta'][before], -1/3, atol=1e-3)
    expected_onset = .018 if high else .08
    soft = active & (peeq > expected_onset-.003) & (ratio > .05) & (ratio < .95)
    assert soft.sum() >= 2
    observed = np.median(peeq[soft]-.001*(1-ratio[soft]))
    assert observed == pytest.approx(expected_onset, abs=1e-4)
    assert peeq[-1] == pytest.approx(expected_onset+.001, abs=1e-4)
    assert result['active'][-1] == 0
