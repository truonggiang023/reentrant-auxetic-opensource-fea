import json
from pathlib import Path
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


def test_inner_y_platen_control_preserves_core_erosion_surface(tmp_path):
    config = dict(CONFIG, platen_contact_surface='inner_y')
    meta = write_radioss_deck(tmp_path/'inner', NODES, ELEMENTS, config,
                              [(225., 0.), (300., 1.)], duration=.01,
                              drop_ccx_viscosity=True)
    text = (tmp_path/'inner_0000.rad').read_text()
    assert '/SURF/PART/ALL/1' in text
    assert '/SURF/PART/ALL/2' not in text
    assert '/SURF/PART/ALL/3' not in text
    lines = text.splitlines()
    for part, expected in [(2, [2, 11, 12, 16, 15]), (3, [3, 17, 18, 22, 21])]:
        assert list(map(int, lines[lines.index('/SURF/SEG/'+str(part))+2].split())) == expected
    assert meta['platen_contact_surface'] == 'inner_y'


def test_rp_coupling_moves_only_outer_platen_faces(tmp_path):
    cfg=dict(young=69000.,poisson=.3,compression_strain=.12,
             platen_reference_points=True)
    meta=write_radioss_deck(tmp_path/'rp',NODES,ELEMENTS,cfg,
                           [(225.,0.),(300.,.2)],duration=.001)
    text=(tmp_path/'rp_0000.rad').read_text()
    assert '/RBE2/1' in text and '/RBE2/2' in text
    assert meta['loading_face_nodes']==[26]
    assert meta['bottom_support_nodes']==[25]
    assert '/BCS/11' in text and '/BCS/12' in text
    assert 'Top face XZ restraint' not in text
    assert meta['core_nodes']==8
    # All original solid elements retain their connectivity; RPs add no solids.
    assert meta['core_elements']==1


@pytest.mark.parametrize('coupled,support_requested,expected_support', [
    (False, False, False), (False, True, True),
    (True, False, True), (True, True, True)])
def test_support_history_is_opt_in_except_for_platen_rps(tmp_path, coupled, support_requested, expected_support):
    cfg = dict(young=69000., poisson=.3, compression_strain=.12,
               platen_reference_points=coupled, support_reaction_history=support_requested)
    meta = write_radioss_deck(tmp_path/'case', NODES, ELEMENTS, cfg,
                              [(225., 0.), (300., .2)], duration=.001)
    lines = (tmp_path/'case_0000.rad').read_text().splitlines()
    start = lines.index('/TH/NODE/3')+3
    stop = next(i for i in range(start, len(lines)) if lines[i].startswith('/'))
    history_nodes = [int(line[:10]) for line in lines[start:stop]]
    expected = meta['loading_face_nodes'] + (meta['bottom_support_nodes'] if expected_support else [])
    assert history_nodes == expected
    assert meta['support_reaction_history'] is expected_support
    # Support-node identities remain available even when not requested in T01.
    assert meta['bottom_support_nodes']


@pytest.mark.parametrize('scope,group,parts', [('all',0,[1,2,3]), ('core',50,[1])])
def test_ams_scope_selects_declared_parts_without_changing_platen_coupling(tmp_path,scope,group,parts):
    cfg=dict(young=69000.,poisson=.3,compression_strain=.12,
             platen_reference_points=True,ams_scope=scope)
    meta=write_radioss_deck(tmp_path/'rp',NODES,ELEMENTS,cfg,
                          [(225.,0.),(300.,.2)],duration=.001,ams_time_step=1e-6)
    lines=(tmp_path/'rp_0000.rad').read_text().splitlines()
    assert int(lines[lines.index('/AMS')+1])==group
    assert ('/GRPART/PART/50' in lines)==(scope=='core')
    if scope=='core':
        assert lines[lines.index('/GRPART/PART/50')+2].split()==['1']
    assert meta['ams']['scope']==scope
    assert meta['ams']['selected_part_ids']==parts
    assert meta['ams']['part_group_id']==group
    assert '/RBE2/1' in lines and '/RBE2/2' in lines
    assert meta['loading_face_nodes']==[26]
    # Native AMS modifies inertia, but this selection adds no nodal mass card.
    assert not any(line.startswith('/ADMAS') for line in lines)


def test_default_and_invalid_ams_scope(tmp_path):
    cfg=dict(young=69000.,poisson=.3,compression_strain=.12)
    meta=write_radioss_deck(tmp_path/'default',NODES,ELEMENTS,cfg,
                          [(225.,0.),(300.,.2)],duration=.001,ams_time_step=1e-6)
    assert meta['ams']['scope']=='all'
    assert meta['ams']['selected_part_ids']==[1,2,3]
    assert int((tmp_path/'default_0000.rad').read_text().split('/AMS\n')[1].splitlines()[0])==0
    for value in ('steel','',None,True,1):
        with pytest.raises(ValueError,match='ams_scope'):
            write_radioss_deck(tmp_path/'invalid',NODES,ELEMENTS,dict(cfg,ams_scope=value),
                              [(225.,0.),(300.,.2)],duration=.001,ams_time_step=1e-6)
    assert not (tmp_path/'invalid_0000.rad').exists()


def test_optional_restart_checkpoints_preserve_default_and_validate_input(tmp_path):
    cfg=dict(young=69000.,poisson=.3,compression_strain=.12)
    for name,interval in [('default',None),('checkpoint',1000)]:
        c=dict(cfg)
        if interval is not None:c['restart_cycle_interval']=interval
        m=write_radioss_deck(tmp_path/name,NODES,ELEMENTS,c,[(225.,0.),(300.,.2)],duration=.01)
        engine=Path(m['engine']).read_text()
        assert ('/RFILE/2' in engine)==(interval is not None)
        assert not m['restart_policy']['continuation_verified']
    for invalid in (0,-1,1.5,True):
        with pytest.raises(ValueError,match='restart_cycle_interval'):
            write_radioss_deck(tmp_path/'invalid',NODES,ELEMENTS,dict(cfg,restart_cycle_interval=invalid),[(225.,0.),(300.,.2)],duration=.01)


def test_linear_profile_means_constant_velocity_instead_of_smoothstep(tmp_path):
    cfg=dict(young=69000.,poisson=.3,compression_strain=.12,loading_profile='linear')
    m=write_radioss_deck(tmp_path/'linear',NODES,ELEMENTS,cfg,[(225.,0.),(300.,.2)],duration=3.6)
    assert m['displacement_fraction_samples']==[[0.,0.],[3.6,1.]]
    assert m['peak_imposed_speed_mm_per_s']==pytest.approx(2./60.)
    assert m['mean_imposed_speed_mm_per_s']==pytest.approx(2./60.)
    assert not m['constitutive_rate_dependence']
    assert m['hardening_curve_count']==1
    with pytest.raises(ValueError,match='loading_profile'):
        write_radioss_deck(tmp_path/'bad_profile',NODES,ELEMENTS,dict(cfg,loading_profile='bad'),[(225.,0.),(300.,.2)],duration=1.)
