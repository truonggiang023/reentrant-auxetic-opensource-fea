"""Reference extraction must not turn axes, legends, or bad histories into data."""
import json
from pathlib import Path

import numpy as np
import pytest

from scripts.digitize_dong2019 import PANELS, digitize_panel
from plot_paper_comparison import read_history, read_reference, resolve_wall_mm, solver_label


def test_digitizer_maps_pixels_and_excludes_legend_and_right_ticks():
    image = np.full((772, 1884, 3), 255, dtype=np.uint8)
    # Thin-wall origin=(1094,591); 384 columns is strain .4; y=428 is .4 MPa.
    image[427:430, 1478] = [240, 5, 5]
    image[74:78, 1145:1200] = [240, 5, 5]  # legend, never observations
    image[345:348, 1849:1862] = 0  # inward right-axis tick
    rows = digitize_panel(image, PANELS["thin"])
    assert len(rows) == 1
    assert rows[0]["curve"] == "experiment_1"
    assert rows[0]["strain"] == pytest.approx(.4)
    assert rows[0]["stress_MPa"] == pytest.approx((591-428)*1.4/570)
    assert rows[0]["y_span_MPa"] == pytest.approx(2*1.4/570)


def test_digitizer_rejects_resized_input():
    with pytest.raises(ValueError, match="original"):
        digitize_panel(np.zeros((200, 200, 3), dtype=np.uint8), PANELS["thin"])


def test_published_reference_preserves_panel_labels_and_range():
    root = Path(__file__).resolve().parents[1]
    for panel in ("thin", "thick"):
        curves = read_reference(root/f"validation/paper_reference/dong2019_fig12_{panel}.csv")
        assert set(curves) == {"experiment_1", "experiment_2", "paper_fe"}
        assert all(len(data) > 500 for data in curves.values())
        for data in curves.values():
            assert np.all(data[:, 0] > 0)
            assert np.all(data[:, 0] < PANELS[panel]["max_strain"])
            assert np.all(data[:, 1] < PANELS[panel]["max_stress_MPa"])


@pytest.mark.parametrize("rows", [[], [{"strain": 0.1, "nominal_stress_MPa": 1}],
    [{"strain": .2, "nominal_stress_MPa": 1}, {"strain": .1, "nominal_stress_MPa": 2}],
    [{"strain": .1, "nominal_stress_MPa": 1}, {"strain": .2, "nominal_stress_MPa": float("nan")}]] )
def test_invalid_assessment_cannot_be_plotted_as_history(tmp_path, rows):
    source = tmp_path/"bad.json"
    source.write_text(json.dumps({"reaction_history": rows}))
    with pytest.raises(ValueError):
        read_history(source)


def test_unknown_assessment_requires_explicit_wall_assertion(tmp_path):
    with pytest.raises(ValueError, match="no verified wall thickness"):
        resolve_wall_mm({}, root=tmp_path)
    wall, evidence = resolve_wall_mm({}, .3, root=tmp_path)
    assert wall == .3
    assert "caller assertion" in evidence


def test_thick_wall_or_conflicting_metadata_cannot_be_labelled_thin(tmp_path):
    with pytest.raises(ValueError, match="requires t=0.3"):
        resolve_wall_mm({"config_wall_mm": 1}, root=tmp_path)
    with pytest.raises(ValueError, match="conflicts"):
        resolve_wall_mm({"generation": {"config_wall_mm": 1}}, .3, root=tmp_path)
    with pytest.raises(ValueError, match="Conflicting"):
        resolve_wall_mm({"config_wall_mm": .3, "generation": {"config_wall_mm": 1}}, root=tmp_path)
    assert resolve_wall_mm({"generation": {"config_wall_mm": .3}}, root=tmp_path)[0] == .3


def test_known_legacy_assessment_matches_qualification_and_config_hash():
    root = Path(__file__).resolve().parents[1]
    assessment = json.loads((root/"validation/native-final.json").read_text())
    wall, evidence = resolve_wall_mm(assessment, root=root)
    assert wall == .3
    assert "SHA256-verified" in evidence
    assessment["input_sha256"] = "different-input"
    with pytest.raises(ValueError, match="no verified wall thickness"):
        resolve_wall_mm(assessment, root=root)


def test_native_plot_distinguishes_total_contact_and_legacy_normal_component():
    common = dict(backend='OpenCourant / Radioss format experimental', completed=True,
                  quasi_static_validated=False)
    total = dict(common, force_mapping=dict(contact_force_quantity='total_force_on_secondary',
                                           component_operation='FN-FT'))
    assert 'total FN-FT' in solver_label(total)[1]
    assert 'quasi-static unverified' in solver_label(total)[1]
    legacy = dict(common, force_mapping=dict(csv_column='Contact force and energy 3 var 64'))
    assert 'normal component only' in solver_label(legacy)[1]
    explicit = dict(common, force_mapping=dict(contact_force_quantity='normal_component_only'))
    assert 'normal component only' in solver_label(explicit)[1]
