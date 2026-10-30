"""Native contact on a newly exposed solid face, with an exterior-only control.

The retained and sacrificial bricks share four nodes and one initially internal
face. After the upper brick erodes, a smaller prescribed impactor crosses its
former volume and loads the retained brick. Keeping the impactor away from all
sidewalls rules out side/corner contact as an explanation. This is a component
test of erosion-face activation, not complete general-contact qualification.
"""
import csv
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
import pytest

from reentrant_fea.radioss import _ints, _reals, _group, _bcs, _contact, _material


def _fresh_face_run(directory, mode):
    executables = [shutil.which(name) for name in
                   ('starter_linux64_gf', 'engine_linux64_gf', 'th_to_csv_linux64_gf')]
    if not all(executables):
        pytest.skip('Source the installed OpenCourant env.sh for native contact qualification')
    directory.mkdir(exist_ok=True)
    name = 'fresh_'+mode.lower()
    nodes = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
        [0, 2, 0], [1, 2, 0], [0, 2, 1], [1, 2, 1],
        [0, 3, 0], [1, 3, 0], [1, 4, 0], [0, 4, 0],
        [0, 3, 1], [1, 3, 1], [1, 4, 1], [0, 4, 1]], float)
    nodes[12:, 0] = .25+.5*nodes[12:, 0]
    nodes[12:, 2] = .25+.5*nodes[12:, 2]
    config = {'young': 69000., 'poisson': .3, 'density': 2.67e-9}
    law = {'fracture_locus': [[-1., .04], [1., .04]], 'failure_displacement': .001}
    lines = ['#RADIOSS STARTER', '/BEGIN', name, _ints(2022, 0),
             f'{"Mg":>20}{"mm":>20}{"s":>20}', f'{"Mg":>20}{"mm":>20}{"s":>20}',
             '/TITLE', 'Native internal-face erosion/contact component test']
    _material(lines, config, np.array([[225., 0.], [300., 1.]]), law)
    lines += ['/MAT/ELAST/2', 'Elastic impactor', _reals(2.67e-9), _reals(69000., .3), '/NODE']
    lines += [_ints(i)+_reals(*xyz) for i, xyz in enumerate(nodes, 1)]
    lines += ['/PROP/SOLID/1', 'HEPH solid',
              _ints(24, 0, None, 0, None, 0, 0, 0)+_reals(0),
              _reals(0, 0, 0), _reals(0)+_ints(1),
              '/PART/1', 'Connected lower and sacrificial solids', _ints(1, 1),
              '/BRICK/1', _ints(1, *range(1, 9)), _ints(2, 4, 3, 10, 9, 8, 7, 12, 11),
              '/PART/2', 'Impactor', _ints(1, 2), '/BRICK/2', _ints(3, *range(13, 21))]
    _group(lines, 1, 'Fixed retained lower solid', range(1, 9))
    _bcs(lines, 1, 'Fixed lower solid', 1, '111')
    _group(lines, 2, 'Sacrificial solid top nodes', range(9, 13))
    _bcs(lines, 2, 'Upper top XZ', 2, '101')
    _group(lines, 3, 'Prescribed impactor', range(13, 21))
    _bcs(lines, 3, 'Impactor XZ', 3, '101')
    for ident, points in [(20, [(0, 0), (.0003, 1), (.001, 1)]),
                          (21, [(0, 0), (.0005, 0), (.001, 1)])]:
        lines += [f'/FUNCT/{ident}', f'Motion fraction {ident}']
        lines += [_reals(time, fraction) for time, fraction in points]
    for ident, group, function, amplitude in [(1, 2, 20, .2), (2, 3, 21, -2.2)]:
        lines += [f'/IMPDISP/{ident}', f'Motion {ident}',
                  _ints(function)+'         Y'+_ints(0, 0, group, 0, 0),
                  _reals(1, amplitude, 0, 0)]
    # ALL retains the internal face; EXT is a deliberate missing-face control.
    # Both connected bricks are in the SAME part to make that distinction real.
    lines += [f'/SURF/PART/{mode}/1', 'Connected eroding part', _ints(1),
              '/SURF/PART/ALL/2', 'Impactor', _ints(2)]
    _contact(lines, 1, 2, 1, 0.)
    contact = lines.index('/INTER/TYPE25/1')
    lines[contact+3] = _ints(0, None)+_reals(0, 0, .005, .005)
    lines += ['/TH/INTER/1', 'Fresh-face impact contact', 'DEF', _ints(1),
              '/TH/BRIC/2', 'Active element states',
              ''.join(f'{v:<10}' for v in ('OFF', 'PLAS')),
              _ints(1)+' '*10+'lower', _ints(2)+' '*10+'sacrificial',
              _ints(3)+' '*10+'impactor', '/TH/NODE/3', 'Impactor geometry',
              ''.join(f'{v:<10}' for v in ('DX', 'DY', 'DZ')),
              _ints(13, 0)+'impact-base', '/END']
    starter = directory/(name+'_0000.rad')
    engine = directory/(name+'_0001.rad')
    starter.write_text('\n'.join(lines)+'\n')
    engine.write_text('\n'.join([
        '#RADIOSS ENGINE', f'/RUN/{name}/1', _reals(.001),
        '/PRINT/-100', '/TFILE/4', _reals(.001/500),
        '/ANIM/DT', _reals(0, .0001), '/ANIM/BRICK/DAMA', '/ANIM/VECT/DISP'])+'\n')
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '2'
    for executable, args in [(executables[0], ['-i', starter.name, '-np', '1']),
                             (executables[1], ['-i', engine.name]),
                             (executables[2], [name+'T01'])]:
        with (directory/(Path(executable).name+'.log')).open('w') as stream:
            subprocess.run([executable, *args], cwd=directory, env=env,
                           stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=40)
    starter_text, engine_text = starter.with_suffix('.out').read_text(), engine.with_suffix('.out').read_text()
    assert '0 ERROR(S)' in starter_text and '0 WARNING(S)' in starter_text
    assert 'NORMAL TERMINATION' in starter_text and 'NORMAL TERMINATION' in engine_text
    match = re.search(r'RUPTURE OF SOLID ELEMENT\s*:\s*2 AT TIME\s*:\s*([\d.E+-]+)', engine_text)
    assert match and 0 < float(match.group(1)) < .0005
    with (directory/(name+'T01.csv')).open() as stream:
        reader = csv.reader(stream)
        header = next(reader)
        values = np.array([[float(value) for value in row] for row in reader])
    assert np.isfinite(values).all() and len(values) > 400
    t = values[:, 0]
    assert np.all(np.diff(t) > 0) and t[-1] > .00099
    contact_cols = [i for i, label in enumerate(header) if label.startswith('Fresh-face impact contact')]
    assert len(contact_cols) == 6
    # The Y channel is cumulative impulse; use interval-averaged force.
    impulse = values[:, contact_cols[1]]
    force = np.diff(impulse)/np.diff(t)
    geometry = [i for i, label in enumerate(header) if label.startswith('Impactor geometry')]
    assert len(geometry) == 3
    y = 3.+values[:, geometry[1]]
    mid_y = .5*(y[1:]+y[:-1])
    for eid, title, final in [(1, 'lower', 1), (2, 'sacrificial', 0), (3, 'impactor', 1)]:
        cols = [i for i, label in enumerate(header)
                if label.startswith('Active element states') and re.search(rf'\s{eid}\s+{title}\s', label)]
        assert len(cols) == 2
        assert values[0, cols[0]] == 1 and values[-1, cols[0]] == final
    return dict(time=t, y=y, mid_y=mid_y, force=force, impulse=impulse,
                rupture_time=float(match.group(1)))


def test_new_internal_face_carries_contact_after_neighbor_erodes(tmp_path):
    all_faces = _fresh_face_run(tmp_path/'all', 'ALL')
    exterior = _fresh_face_run(tmp_path/'exterior_control', 'EXT')
    # Both motions pass through the vanished upper volume before reaching y=1.
    # No traction may remain on the deleted upper face at its former y=2.2.
    assert np.max(np.abs(all_faces['force'][all_faces['mid_y'] > 1.3])) < 1e-6
    assert all_faces['y'][-1] < .85 and exterior['y'][-1] < .85
    # Actual traction appears only when penetrating the newly exposed face.
    assert np.max(np.abs(all_faces['force'][all_faces['mid_y'] < .98])) > 1000.
    assert abs(all_faces['impulse'][-1]) > .01
    # With initially internal facets omitted, the same centered impactor has
    # no face to hit. This negative control catches a side/corner explanation.
    assert np.max(np.abs(exterior['force'])) < 1e-6
    assert exterior['rupture_time'] == pytest.approx(all_faces['rupture_time'], rel=.01)


def _rotation_run(directory, contact_enabled):
    executables = [shutil.which(name) for name in
                   ('starter_linux64_gf', 'engine_linux64_gf', 'th_to_csv_linux64_gf')]
    if not all(executables):
        pytest.skip('Source the installed OpenCourant env.sh for native contact qualification')
    directory.mkdir(exist_ok=True)
    name = 'rotation_contact' if contact_enabled else 'rotation_control'
    unit = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                     [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
    moving = .4*(unit-.5)+[.5, 2.5, .5]
    nodes = np.vstack([unit, moving])
    rotation_time, duration = .0002, .00025
    lines = ['#RADIOSS STARTER', '/BEGIN', name, _ints(2022, 0),
             f'{"Mg":>20}{"mm":>20}{"s":>20}', f'{"Mg":>20}{"mm":>20}{"s":>20}',
             '/TITLE', 'Native initially back-facing rotation self-contact component',
             '/MAT/ELAST/1', 'Elastic solids', _reals(2.67e-9), _reals(69000., .3), '/NODE']
    lines += [_ints(i)+_reals(*xyz) for i, xyz in enumerate(nodes, 1)]
    lines += ['/PROP/SOLID/1', 'HEPH solid',
              _ints(24, 0, None, 0, None, 0, 0, 0)+_reals(0),
              _reals(0, 0, 0), _reals(0)+_ints(1),
              '/PART/1', 'Two disconnected solids', _ints(1, 1), '/BRICK/1',
              _ints(1, *range(1, 9)), _ints(2, *range(9, 17))]
    _group(lines, 1, 'Fixed target', range(1, 9))
    _bcs(lines, 1, 'Fixed target', 1, '111')
    _group(lines, 20, 'Moving Z restraint', range(9, 17))
    _bcs(lines, 20, 'Moving Z restraint', 20, '001')
    # Per-node prescribed paths are sampled densely. Actual output geometry
    # below must retain all distances and reach the intended reversed face.
    for index, original in enumerate(moving):
        node, group = index+9, index+2
        _group(lines, group, f'Moving node {node}', [node])
        relative = original-np.array([.5, 2.5, .5])
        for axis, component in [('X', 0), ('Y', 1)]:
            function = 100+2*index+component
            lines += [f'/FUNCT/{function}', f'{axis} displacement node {node}']
            for fraction in np.r_[np.linspace(0, 1, 201), duration/rotation_time]:
                f = min(fraction, 1)
                angle = np.pi*f
                rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                                     [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
                position = rotation@relative+[.5, 2.5-1.4*f, .5]
                lines.append(_reals(rotation_time*fraction,
                                    position[component]-original[component]))
            lines += [f'/IMPDISP/{function}', f'{axis} motion node {node}',
                      _ints(function)+f'{axis:>10}'+_ints(0, 0, group, 0, 0),
                      _reals(1, 1, 0, 0)]
    if contact_enabled:
        # SAME part, two disconnected bricks: native self-contact, no initial
        # orientation filter or manually selected facing surface pair.
        lines += ['/SURF/PART/ALL/1', 'All solid faces', _ints(1)]
        _contact(lines, 1, 1, 0, 0.)
        k = lines.index('/INTER/TYPE25/1')
        lines[k+3] = _ints(0, None)+_reals(0, 0, .005, .005)
        lines += ['/TH/INTER/1', 'Self contact resultant', 'DEF', _ints(1)]
    lines += ['/TH/NODE/2', 'Fixed target reaction', f'{"REACY":<10}']
    lines += [_ints(i, 0)+f'target-{i}' for i in range(1, 9)]
    lines += ['/TH/NODE/3', 'Moving solid geometry',
              ''.join(f' {v:<9}' for v in ('DX', 'DY', 'DZ'))]
    lines += [_ints(i, 0)+f'moving-{i}' for i in range(9, 17)]
    lines += ['/END']
    starter = directory/(name+'_0000.rad')
    engine = directory/(name+'_0001.rad')
    starter.write_text('\n'.join(lines)+'\n')
    engine.write_text('\n'.join([
        '#RADIOSS ENGINE', f'/RUN/{name}/1', _reals(duration),
        '/PRINT/-100', '/TFILE/4', _reals(duration/500),
        '/ANIM/DT', _reals(0, duration/10), '/ANIM/VECT/DISP'])+'\n')
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '2'
    for executable, args in [(executables[0], ['-i', starter.name, '-np', '1']),
                             (executables[1], ['-i', engine.name]),
                             (executables[2], [name+'T01'])]:
        with (directory/(Path(executable).name+'.log')).open('w') as stream:
            subprocess.run([executable, *args], cwd=directory, env=env,
                           stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=40)
    starter_text, engine_text = starter.with_suffix('.out').read_text(), engine.with_suffix('.out').read_text()
    assert '0 ERROR(S)' in starter_text and '0 WARNING(S)' in starter_text
    assert 'NORMAL TERMINATION' in starter_text and 'NORMAL TERMINATION' in engine_text
    with (directory/(name+'T01.csv')).open() as stream:
        reader = csv.reader(stream)
        header = next(reader)
        values = np.array([[float(value) for value in row] for row in reader])
    t = values[:, 0]
    assert len(t) > 400 and np.isfinite(values).all() and np.all(np.diff(t) > 0)
    reaction = [i for i, label in enumerate(header) if label.startswith('Fixed target reaction')]
    assert len(reaction) == 8
    force = np.diff(values[:, reaction].sum(axis=1))/np.diff(t)
    positions = []
    for node in range(9, 17):
        columns = [i for i, label in enumerate(header)
                   if label.startswith('Moving solid geometry')
                   and re.search(rf'\s{node}\s+moving-{node}\s', label)]
        assert len(columns) == 3
        positions.append(values[:, columns]+moving[node-9])
    positions = np.stack(positions, axis=1)
    rigidity_error = 0.
    for a in range(8):
        for b in range(a+1, 8):
            original_distance = np.linalg.norm(moving[a]-moving[b])
            distances = np.linalg.norm(positions[:, a]-positions[:, b], axis=1)
            rigidity_error = max(rigidity_error, np.max(abs(distances-original_distance))/original_distance)
    assert rigidity_error < 1e-4
    # The original TOP face points away from the target initially. Measured
    # final coordinates must put this same face below the target's y=1 face.
    initial_normal = np.cross(moving[3]-moving[2], moving[6]-moving[2])
    final_normal = np.cross(positions[-1, 3]-positions[-1, 2],
                            positions[-1, 6]-positions[-1, 2])
    assert np.dot(initial_normal, np.array([.5, .5, .5])-moving[[2, 3, 6, 7]].mean(axis=0)) < 0
    assert np.dot(final_normal, np.array([.5, .5, .5])-positions[-1, [2, 3, 6, 7]].mean(axis=0)) > 0
    np.testing.assert_allclose(positions[-1, [2, 3, 6, 7], 1], .9, atol=1e-4)
    return dict(force=force, midtime=.5*(t[1:]+t[:-1]),
                rigidity_error=float(rigidity_error), positions=positions)


def test_self_contact_detects_originally_back_facing_face_after_half_turn(tmp_path):
    contact = _rotation_run(tmp_path/'contact', True)
    control = _rotation_run(tmp_path/'no_contact', False)
    # The target is fixed: its summed reaction has no target inertial term.
    # A no-contact control has exactly the same imposed moving-solid motion.
    assert np.max(abs(control['force'])) < 1e-6
    assert np.max(abs(contact['force'][contact['midtime'] < .00016])) < 1e-6
    active = np.flatnonzero(abs(contact['force']) > 1e-4)
    assert len(active) > 20
    first_angle = contact['midtime'][active[0]]/.0002*180
    assert 150 < first_angle < 180
    held_force = abs(contact['force'][contact['midtime'] > .00021])
    assert len(held_force) > 20 and np.median(held_force) > 1000.
    assert np.std(held_force)/np.mean(held_force) < .02
