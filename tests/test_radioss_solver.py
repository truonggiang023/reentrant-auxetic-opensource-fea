"""Real native-solver regressions; source /workspace/opencourant/env.sh first.

These are component qualification, not experimental calibration or a proof of
complete fracture-contact coverage.
"""
import csv
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
import pytest

from reentrant_fea.radioss import write_radioss_deck


NODES = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
ELEMENTS = np.array([list(range(8))])


def _run(tmp_path, name, config, *, tensile_coupon, ams_time_step=None):
    executables = [shutil.which(name) for name in
                   ('starter_linux64_gf', 'engine_linux64_gf', 'th_to_csv_linux64_gf')]
    if not all(executables):
        pytest.skip('OpenCourant runtime not on PATH; source its env.sh for native qualification')
    meta = write_radioss_deck(tmp_path/name, NODES, ELEMENTS, config,
                              [(225., 0.), (300., 1.)], duration=.0001,
                              tensile_coupon=tensile_coupon, ams_time_step=ams_time_step)
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '2'
    for executable, arguments in [(executables[0], ['-i', name+'_0000.rad', '-np', '1']),
                                  (executables[1], ['-i', name+'_0001.rad']),
                                  (executables[2], [name+'T01'])]:
        log = tmp_path/(Path(executable).name+'.log')
        with log.open('w') as stream:
            subprocess.run([executable, *arguments], cwd=tmp_path, env=env,
                           stdout=stream, stderr=subprocess.STDOUT,
                           check=True, timeout=45)
    starter = (tmp_path/(name+'_0000.out')).read_text()
    engine = (tmp_path/(name+'_0001.out')).read_text()
    assert re.search(r'0 ERROR\(S\)', starter)
    assert re.search(r'0 WARNING\(S\)', starter)
    assert 'NORMAL TERMINATION' in engine
    assert re.search(r'TOTAL NUMBER OF CYCLES\s*:\s*[1-9]', engine)
    with (tmp_path/(name+'T01.csv')).open() as stream:
        reader = csv.reader(stream)
        headers = next(reader)
        values = np.array([[float(x) for x in row] for row in reader])
    assert len(values) > 100 and np.isfinite(values).all()
    assert values[-1, 0] >= .95*meta['duration_seconds']
    assert np.all(np.diff(values[:, 0]) > 0)
    return meta, headers, values, engine


def test_native_inievo_erodes_element_and_releases_bonded_end_force(tmp_path):
    config = {'young': 69000., 'poisson': .3, 'compression_strain': .15,
              'ductile_damage': {'fracture_locus': [[-1., .04], [1., .04]],
                                 'failure_displacement': .001}}
    meta, headers, values, engine = _run(tmp_path, 'tension', config, tensile_coupon=True)
    # The solver reports volume-element removal, not only stress degradation.
    match = re.search(r'RUPTURE OF SOLID ELEMENT\s*:\s*1 AT TIME\s*:\s*([\d.E+-]+)', engine)
    assert match, engine[-3000:]
    failed_time = float(match.group(1))
    assert 0 < failed_time < meta['duration_seconds']
    eroded = next(i for i, h in enumerate(headers) if h.strip().endswith('ERODED'))
    assert values[0, eroded] == 0 and values[-1, eroded] == 1
    # TH/REACY channels are accumulated impulses. Compare interval forces after
    # differentiation; interpreting their raw plateau as residual force is wrong.
    node_columns = [i for i, h in enumerate(headers) if h.startswith('Loading face nodal histories')]
    assert len(node_columns) == 6*len(meta['loading_face_nodes'])
    impulse = values[:, node_columns[4::6]].sum(axis=1)
    time = values[:, 0]
    force = np.diff(impulse)/np.diff(time)
    midtime = .5*(time[:-1]+time[1:])
    before = np.max(np.abs(force[midtime < .9*failed_time]))
    after = np.max(np.abs(force[midtime > 1.2*failed_time]))
    assert before > 100.  # Nontrivial tension load actually developed.
    assert after < .01*before


def test_type25_transmits_compression_between_separate_elastic_platens(tmp_path):
    config = {'young': 69000., 'poisson': .3, 'compression_strain': .12,
              'platen_thickness': 1., 'platen_margin': .1, 'friction': .35}
    meta, headers, values, _ = _run(tmp_path, 'compression', config, tensile_coupon=False)
    contact_columns = [i for i, h in enumerate(headers) if h.startswith('Contact force and energy')]
    assert len(contact_columns) == 18  # Three interfaces, six components each.
    top_impulse = values[:, contact_columns[13]]
    bottom_impulse = values[:, contact_columns[7]]
    assert top_impulse[-1] < -.001 and bottom_impulse[-1] > .001
    assert abs(top_impulse[-1]+bottom_impulse[-1]) < .05*abs(top_impulse[-1])
    reaction_force = np.diff(top_impulse)/np.diff(values[:, 0])
    assert -np.min(reaction_force) > 100
    loading_columns = [i for i, h in enumerate(headers) if h.startswith('Loading face nodal histories')]
    top_displacements = values[-1, loading_columns[1::6]]
    expected = -meta['height_mm']*meta['target_nominal_strain']
    assert np.allclose(top_displacements, expected, rtol=.005)


def test_native_ams_keeps_inievo_erosion_and_reports_actual_operator(tmp_path):
    config = {'young': 69000., 'poisson': .3, 'compression_strain': .15,
              'ductile_damage': {'fracture_locus': [[-1., .04], [1., .04]],
                                 'failure_displacement': .001}}
    meta, headers, values, engine = _run(tmp_path, 'tension_ams', config,
                                        tensile_coupon=True, ams_time_step=5e-7)
    assert 'ADVANCED MASS SCALING IS ON' in engine
    assert 'RUPTURE OF SOLID ELEMENT' in engine
    eroded = next(i for i, h in enumerate(headers) if h.strip().endswith('ERODED'))
    assert values[-1, eroded] == 1
    # No net added mass does not mean no mass-operator modification: retain
    # explicit provenance and the independent sensitivity requirement.
    assert meta['ams']['modifies_inertial_operator']
    assert meta['ams']['timestep_sensitivity_required']
