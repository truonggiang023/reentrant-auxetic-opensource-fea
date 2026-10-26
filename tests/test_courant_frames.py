import numpy as np
import pytest

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
