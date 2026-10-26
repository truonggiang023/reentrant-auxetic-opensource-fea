import json
import numpy as np
import pytest

from reentrant_fea.radioss import write_radioss_deck


NODES = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
ELEMENTS = np.array([list(range(8))])
CONFIG = {'young': 69000., 'poisson': .3, 'compression_strain': .65,
          'wall': .3, 'ductile_damage': {'fracture_locus': [[0., .08], [.333, .047], [1., .018]],
                                       'failure_displacement': .001, 'viscosity_time': .005}}


def test_requires_explicit_viscosity_omission_and_preserves_provenance(tmp_path):
    with pytest.raises(ValueError, match='pseudo-time viscosity'):
        write_radioss_deck(tmp_path/'case', NODES, ELEMENTS, CONFIG,
                           [(225., 0.), (300., 1.)], duration=.01)
    meta = write_radioss_deck(tmp_path/'case', NODES, ELEMENTS, CONFIG,
                              [(225., 0.), (300., 1.)], duration=.01,
                              drop_ccx_viscosity=True)
    assert not meta['calibrated'] and not meta['experimental_agreement']
    assert not meta['completion_verified']
    assert meta['ccx_viscosity_omitted'] == .005
    assert meta['physical_duration_seconds'] == pytest.approx(19.5)
    saved = json.loads((tmp_path/'case.radioss.json').read_text())
    assert saved['fracture_locus'] == CONFIG['ductile_damage']['fracture_locus']
    assert saved['acceleration_requires_energy_and_speed_sensitivity']
    # A repeated generator must not overwrite inputs/results from a prior run.
    with pytest.raises(FileExistsError):
        write_radioss_deck(tmp_path/'case', NODES, ELEMENTS, CONFIG,
                           [(225., 0.), (300., 1.)], duration=.01,
                           drop_ccx_viscosity=True)


@pytest.mark.parametrize('bad', [np.array([[0, 1, 2, 3, 4, 5]]),
                                np.array([[0, 1, 2, 3, 4, 5, 6, 6]]),
                                np.array([[0, 1, 2, 3, 4, 5, 6, 8]])])
def test_reject_unqualified_or_invalid_connectivity(tmp_path, bad):
    with pytest.raises(ValueError, match='8-node brick'):
        write_radioss_deck(tmp_path/'case', NODES, bad, CONFIG,
                           [(225., 0.), (300., 1.)], duration=.01)


def test_invalid_physical_inputs_cannot_create_a_deck(tmp_path):
    for duration in [0, -1, float('nan')]:
        with pytest.raises(ValueError, match='positive duration'):
            write_radioss_deck(tmp_path/'case', NODES, ELEMENTS, CONFIG,
                               [(225., 0.), (300., 1.)], duration=duration)
    assert not list(tmp_path.iterdir())


def test_advanced_mass_scaling_requires_explicit_valid_target_and_provenance(tmp_path):
    for target in [0, -1, float('nan'), .001]:
        with pytest.raises(ValueError, match='AMS time step'):
            write_radioss_deck(tmp_path/'bad', NODES, ELEMENTS, CONFIG,
                               [(225., 0.), (300., 1.)], duration=.01,
                               drop_ccx_viscosity=True, ams_time_step=target)
    meta = write_radioss_deck(tmp_path/'ams', NODES, ELEMENTS, CONFIG,
                              [(225., 0.), (300., 1.)], duration=.01,
                              drop_ccx_viscosity=True, ams_time_step=1e-6)
    assert meta['ams']['enabled']
    assert meta['ams']['modifies_inertial_operator']
    assert meta['ams']['timestep_sensitivity_required']
    assert meta['mass_scaling'] == 'native Advanced Mass Scaling'
    assert not meta['experimental_agreement']
    assert meta['numerical_damping']['native_solid_defaults_active']
    assert meta['numerical_damping']['starter_reported_contact_critical_damping'] == .05
