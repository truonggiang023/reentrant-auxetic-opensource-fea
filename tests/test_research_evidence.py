"""Reject stale provenance and misleading explicit-energy readiness claims."""
import json

import pytest

from reentrant_fea.provenance import record_generation, read_generation
from reentrant_fea.radioss_results import (read_engine_out, energy_screen,
                                         read_history_csv, history_energy_screen)


def test_generation_rejects_replaced_deck_or_material(tmp_path):
    deck = tmp_path/'core.inp'
    deck.write_text('original mesh')
    material = deck.with_suffix('.material.dat')
    material.write_text('original law')
    config, hardening = tmp_path/'config.json', tmp_path/'plastic.py'
    config.write_text(json.dumps({'wall': .3}))
    hardening.write_text('original hardening')
    record_generation(deck, config, hardening, config={'wall': .3}, metadata={})
    assert read_generation(deck)['config_wall_mm'] == .3
    material.write_text('different law')
    with pytest.raises(ValueError, match='sidecar'):
        read_generation(deck)
    material.write_text('original law')
    deck.write_text('wrong mesh')
    with pytest.raises(ValueError, match='input deck'):
        read_generation(deck)


def engine_output(tmp_path, kinetic=.01, *, normal=True, last=1., mass=True):
    path = tmp_path/'engine.out'
    header = 'CYCLE TIME TIME-STEP ELEMENT ERROR I-ENERGY K-ENERGY T K-ENERGY R EXT-WORK MAS.ERR'
    lines = ['FINAL TIME ... 1.0', header]
    for i in range(11):
        suffix = ' 1.000 0.000' if mass else ''
        lines.append(f'{i*100} {i/10*last:.4f} 1.0D-04 SOLID 1 -1.0% 1.0 {kinetic} 0.0 1.01 0.0{suffix}')
    lines += ['DELETE SOLID ELEMENT NUMBER 12']
    if normal:
        lines += ['NORMAL TERMINATION']
    path.write_text('\n'.join(lines))
    return read_engine_out(path)


def test_dynamic_energy_cannot_be_passed_by_normal_termination(tmp_path):
    evidence = engine_output(tmp_path, kinetic=.5)
    assert evidence['normal_termination'] and evidence['deleted_solid_elements'] == [12]
    result = energy_screen(evidence, 1.)
    assert not result['sampled_energy_screen_passed']
    assert any('kinetic/internal' in reason for reason in result['reasons'])


def test_header_time_does_not_supply_missing_final_evidence(tmp_path):
    evidence = engine_output(tmp_path, last=.5)
    assert evidence['last_printed_time'] == .5
    result = energy_screen(evidence, 1.)
    assert not result['sampled_energy_screen_passed']
    assert any('end of loading' in reason for reason in result['reasons'])


def test_sampled_energy_is_only_a_screen_and_missing_mass_blocks_it(tmp_path):
    result = energy_screen(engine_output(tmp_path), 1.)
    assert result['sampled_energy_screen_passed']
    assert result['quasi_static_validated'] is False
    result = energy_screen(engine_output(tmp_path, mass=False), 1.)
    assert not result['sampled_energy_screen_passed']
    assert 'missing mass-scaling evidence' in result['reasons']


def test_csv_corruption_and_restart_are_not_silently_read(tmp_path):
    path = tmp_path/'history.csv'
    path.write_text('time,INTERNAL ENERGY\n0,1\n.1,2\n.05,3\n')
    with pytest.raises(ValueError, match='increasing'):
        read_history_csv(path)
    path.write_text('time,INTERNAL ENERGY\n0,1\n.1\n')
    with pytest.raises(ValueError, match='Incomplete'):
        read_history_csv(path)


def test_dense_history_catches_energy_peak_hidden_between_prints(tmp_path):
    import numpy as np
    times = np.linspace(0, 1, 101)
    history = {name: np.ones(101) for name in ('INTERNAL ENERGY', 'MASS')}
    history.update({name: np.zeros(101) for name in (
        'KINETIC ENERGY', 'ROTATION ENERGY', 'HOURGLASS ENERGY',
        'ADDED MASS', 'EXTERNAL WORK', 'CONTACT ENERGY')})
    history['time'] = times
    history['ADDED MASS'][:] = -1e-16  # ordinary converter roundoff
    assert history_energy_screen(history, 1)['sampled_energy_screen_passed']
    history['KINETIC ENERGY'][51] = .5
    result = history_energy_screen(history, 1)
    assert not result['sampled_energy_screen_passed']
    assert result['max_kinetic_internal_ratio'] == .5


def test_solver_runner_preserves_exit_status_and_bounds_execution(tmp_path):
    import subprocess
    import sys
    from run_radioss import run_logged
    log = tmp_path/'solver.log'
    with pytest.raises(subprocess.CalledProcessError) as error:
        run_logged([sys.executable, '-c', 'print("solver failure"); raise SystemExit(7)'],
                   tmp_path, log, 2)
    assert error.value.returncode == 7
    assert 'solver failure' in log.read_text()
    with pytest.raises(RuntimeError, match='timeout'):
        run_logged([sys.executable, '-c', 'import time; time.sleep(20)'],
                   tmp_path, log, .1)
