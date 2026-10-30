import numpy as np
import pytest

import plot_courant_frames as mod
from plot_courant_frames import front_polygons, frame_time, projected_diagnostics, read_frame, strain_at_time


def sample_frame():
    initial = np.array([[0,0,0], [1,0,0], [1,1,0], [0,1,0],
                        [0,0,1], [1,0,1], [1,1,1], [0,1,1]], float)
    return {"initial_points": initial, "points": initial+[2,3,0],
            "cells": np.array([range(8), range(8)]), "active": np.array([True, False]),
            "part_id": np.array([1,1]), "element_id": np.array([101,102])}


def test_eroded_cells_are_not_rendered_and_displacement_is_not_applied_twice():
    frame = sample_frame()
    polygons, indices = front_polygons(frame)
    assert len(polygons) == 1
    assert indices.tolist() == [0]
    np.testing.assert_array_equal(polygons[0], [[2,3],[3,3],[3,4],[2,4]])


def test_frame_strain_uses_imposed_samples_instead_of_linear_time():
    meta = {"displacement_fraction_samples": [[0,0], [.5,.2], [1,1]], "target_nominal_strain": .65}
    assert strain_at_time(.5, meta) == pytest.approx(.13)
    assert strain_at_time(.75, meta) == pytest.approx(.39)
    with pytest.raises(ValueError, match="outside"):
        strain_at_time(1.1, meta)
    assert frame_time("TIME 1 1 double\n0.0052\nCYCLE 1 1 int\n0\n") == .0052


def test_overlap_and_inversion_diagnostic_preserves_original_element_ids():
    square = np.array([[0,0],[1,0],[1,1],[0,1]], float)
    result = projected_diagnostics([square, square+.5, square[::-1]], [41,77,109])
    assert result["invalid_projected_element_ids"] == [109]
    assert result["overlap_pairs"] == 1
    assert result["max_overlap_mm2"] == pytest.approx(.25)


def test_parser_rejects_non_ascii_or_missing_time():
    with pytest.raises(ValueError, match="ASCII"):
        read_frame("# vtk\nBINARY\n")
    with pytest.raises(ValueError, match="TIME"):
        frame_time("TIME absent")


def test_rp_frame_sampling_uses_verified_raw_relative_motion_and_rejects_extrapolation(tmp_path, monkeypatch):
    metadata = dict(prefix='case', platen_reference_points=True, height_mm=10.,
                    displacement_fraction_samples=[[0., 0.], [.5, .8], [1., 1.]],
                    target_nominal_strain=.12)
    metadata_path = tmp_path/'case.radioss.json'
    metadata_path.write_text('{}')
    (tmp_path/'caseT01.csv').write_text('fixture history')
    (tmp_path/'case.run.json').write_text('{}')
    def assess(mp, hp, **kwargs):
        assert mp == metadata_path and hp == tmp_path/'caseT01.csv'
        assert kwargs == {'force_source': 'rp'}
        return dict(completed=True, run_binding_verified=True, metadata_binding_verified=True,
                    force_mapping=dict(top_displacement_column='actual top DY',
                                       bottom_displacement_column='actual bottom DY'))
    monkeypatch.setattr(mod, 'assess', assess)
    monkeypatch.setattr(mod, 'read_history_csv', lambda path:
                        {'time': np.array([0., .5, .9]),
                         'actual top DY': np.array([0., -.3, -.8]),
                         'actual bottom DY': np.array([0., .1, .2])})
    samples, source = mod.strain_source(metadata_path, metadata)
    assert source['observed_RP_measurement'] and source['completion_verified_by_strain_source']
    assert mod.observed_strain_at_time(.5, samples) == pytest.approx(.04)
    assert strain_at_time(.5, metadata) == pytest.approx(.096)
    assert mod.observed_strain_at_time(.7, samples) == pytest.approx(.07)
    assert mod.observed_strain_at_time(.9, samples) == pytest.approx(.1)
    with pytest.raises(ValueError, match='outside observed'):
        mod.observed_strain_at_time(1., samples)
    with pytest.raises(ValueError, match='outside observed'):
        mod.observed_strain_at_time(-1e-6, samples)


def test_rp_frame_source_rejects_unfinished_run(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, 'assess', lambda *a, **k:
                        dict(completed=False, run_binding_verified=True, metadata_binding_verified=True))
    with pytest.raises(ValueError, match='completed receipt-bound'):
        mod.strain_source(tmp_path/'case.radioss.json',
                          dict(prefix='case', platen_reference_points=True))


def test_thick_native_frames_are_allowed_but_thin_paper_crops_are_guarded():
    assert mod.wall_for_plot({'wall_mm': 1.}) == 1.
    with pytest.raises(ValueError, match='thick-wall paper-pdf comparison is unsupported'):
        mod.wall_for_plot({'wall_mm': 1.}, paper_pdf='paper.pdf')
    assert mod.wall_for_plot({'wall_mm': .3}, paper_pdf='paper.pdf') == .3
    with pytest.raises(ValueError, match='wall_mm'):
        mod.wall_for_plot({'wall_mm': .5})


def test_legacy_non_rp_source_is_explicitly_labelled():
    samples, source = mod.strain_source('not-read.json', {})
    assert samples is None and not source['observed_RP_measurement']
    assert 'Legacy non-RP fallback' in source['source']


def test_rp_connector_lines_are_excluded_without_removing_physical_hexes():
    points = '0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1  2 0 0'
    text = '\n'.join([
        '# vtk DataFile Version 3.0', 'fixture', 'ASCII', 'DATASET UNSTRUCTURED_GRID',
        'FIELD FieldData 1', 'TIME 1 1 double', '.5', 'POINTS 9 float', points,
        'CELLS 2 12', '8 0 1 2 3 4 5 6 7', '2 8 1', 'CELL_TYPES 2', '12 3',
        'POINT_DATA 9', 'VECTORS Displacement float', ' '.join(['0']*27),
        'CELL_DATA 2', 'SCALARS ELEMENT_ID int 1', 'LOOKUP_TABLE default', '41 0',
        'SCALARS PART_ID int 1', 'LOOKUP_TABLE default', '1 1',
        'SCALARS EROSION_STATUS int 1', 'LOOKUP_TABLE default', '1 1',
        'SCALARS 3DELEM_MAX_DAMAGE_ELEMENT float 1', 'LOOKUP_TABLE default', '.2 0',
        ''])
    frame = read_frame(text)
    assert frame['cells'].tolist() == [list(range(8))]
    assert frame['element_id'].tolist() == [41]
    assert frame['damage'].tolist() == [.2]
    assert frame['excluded_rp_connector_cells'] == 1
    assert len(front_polygons(frame)[0]) == 1
    with pytest.raises(ValueError, match='zero-ELEMENT_ID RP connector'):
        read_frame(text.replace('41 0', '41 99'))
