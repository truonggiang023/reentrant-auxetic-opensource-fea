"""Check contact-channel arithmetic using synthetic, unbound solver evidence."""
import csv
import hashlib
import json
import sys

import pytest

from assess_courant import assess, main


@pytest.fixture
def contact_history(tmp_path):
    """Nonzero FT and unequal intervals distinguish sum/difference mistakes."""
    prefix = 'contact'
    deck = tmp_path / (prefix + '_0000.rad')
    deck.write_text('synthetic input deck\n')
    (tmp_path / (prefix + '_0001.out')).write_text('NORMAL TERMINATION\n')
    history = tmp_path / (prefix + 'T01.csv')
    with history.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(('time', 'fn', 'ft'))
        writer.writerows(((0, 0, 0), (.1, -2, 1), (.4, -14, 4), (1, -50, 10)))
    metadata = tmp_path / (prefix + '.radioss.json')
    metadata.write_text(json.dumps({
        'prefix': prefix,
        'backend': 'OpenCourant',
        'sha256': {deck.name: hashlib.sha256(deck.read_bytes()).hexdigest()},
        'duration_seconds': 1.,
        'target_nominal_strain': .6,
        'reference_area_mm2': 2.,
        'displacement_fraction_samples': [[0., 0.], [1., 1.]],
    }))
    return metadata, history


@pytest.mark.parametrize('force_sign', [-1, 1])
def test_total_impulse_is_differentiated_once_over_actual_intervals(contact_history, force_sign):
    result = assess(*contact_history, force_column='fn', tangential_column='ft',
                    force_sign=force_sign, force_kind='impulse')
    rows = result['reaction_history']
    # -(FN - FT) = [0, 3, 18, 60], giving [30, 50, 70] N on
    # the actual 0.1, 0.3 and 0.6 second intervals. A second derivative
    # would lose another sample; FN + FT gives different forces.
    expected_forces = [-force_sign * value for value in (30., 50., 70.)]
    assert len(rows) == 3
    assert [row['time'] for row in rows] == pytest.approx([.05, .25, .7])
    assert [row['reaction_y_N'] for row in rows] == pytest.approx(expected_forces)
    assert [row['nominal_stress_MPa'] for row in rows] == pytest.approx(
        [force / 2. for force in expected_forces])
    assert [row['strain'] for row in rows] == pytest.approx([.03, .15, .42])
    mapping = result['force_mapping']
    assert mapping['contact_force_quantity'] == 'total_force_on_secondary'
    assert mapping['component_operation'] == 'FN-FT'
    assert mapping['csv_column'] == 'fn'
    assert mapping['tangential_csv_column'] == 'ft'
    assert mapping['compressive_sign'] == force_sign
    assert mapping['channel_kind'] == 'impulse'
    assert mapping['sampling'] == 'interval average'
    assert result['engine_normal_termination']
    assert not result['run_binding_verified']
    assert not result['completed']
    assert not result['quasi_static_validated']


def test_total_force_samples_are_not_differentiated(contact_history):
    result = assess(*contact_history, force_column='fn', tangential_column='ft',
                    force_sign=-1, force_kind='force')
    rows = result['reaction_history']
    assert len(rows) == 4
    assert [row['time'] for row in rows] == pytest.approx([0., .1, .4, 1.])
    assert [row['reaction_y_N'] for row in rows] == pytest.approx([0., 3., 18., 60.])
    assert [row['nominal_stress_MPa'] for row in rows] == pytest.approx([0., 1.5, 9., 30.])
    assert result['force_mapping']['channel_kind'] == 'force'
    assert result['force_mapping']['sampling'] == 'sampled force'


@pytest.mark.parametrize('force_kind, expected', [
    ('force', [0., 2., 14., 50.]),
    ('impulse', [20., 40., 60.]),
])
def test_legacy_single_channel_is_explicitly_normal_only(contact_history, force_kind, expected):
    result = assess(*contact_history, force_column='fn', force_sign=-1, force_kind=force_kind)
    assert [row['reaction_y_N'] for row in result['reaction_history']] == pytest.approx(expected)
    assert result['force_mapping']['contact_force_quantity'] == 'normal_component_only'


@pytest.mark.parametrize('force_column, tangential_column', [
    ('fn', 'missing'),
    ('missing', 'ft'),
    (None, 'ft'),
    ('fn', 'fn'),
], ids=['missing-tangential-channel', 'missing-normal-channel',
        'tangential-without-normal-selection', 'same-channel-twice'])
def test_invalid_contact_channel_pairs_are_rejected(contact_history, force_column, tangential_column):
    with pytest.raises(ValueError):
        assess(*contact_history, force_column=force_column, tangential_column=tangential_column,
               force_sign=-1, force_kind='impulse')


def test_cli_accepts_and_forwards_tangential_column(contact_history, tmp_path, monkeypatch, capsys):
    metadata, history = contact_history
    output = tmp_path / 'assessment.json'
    monkeypatch.setattr(sys, 'argv', [
        'assess_courant.py', str(metadata), str(history),
        '--force-column', 'fn', '--tangential-column', 'ft',
        '--force-sign', '-1', '--force-kind', 'impulse', '--output', str(output),
    ])
    main()
    capsys.readouterr()
    result = json.loads(output.read_text())
    assert [row['reaction_y_N'] for row in result['reaction_history']] == pytest.approx([30., 50., 70.])
    assert result['force_mapping']['tangential_csv_column'] == 'ft'
