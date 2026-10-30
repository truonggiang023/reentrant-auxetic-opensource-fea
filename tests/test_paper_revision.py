import json
from pathlib import Path
import pytest
from revise_paper_calibration import build_report

ROOT = Path(__file__).resolve().parents[1]


def assessment(tmp_path, end=.12, **extra):
    path=tmp_path/'thin.json'
    result=dict(completed=True, config_wall_mm=.3, stress_area_convention_mm2=1900.,
                reaction_history=[dict(strain=.001,nominal_stress_MPa=.01),
                                  dict(strain=end,nominal_stress_MPa=.3)])
    result.update(extra);path.write_text(json.dumps(result));return path


def test_completed_early_run_cannot_supply_full_curve_or_calibration(tmp_path):
    report=build_report([assessment(tmp_path)],ROOT/'validation/paper_reference')
    run=report['assessments'][0]
    assert run['solver_completed'] and not run['calibration_eligible']
    assert not run['metrics']['densification']['available']
    assert not run['discrepancy_to_both_experiments']['full_visible']['full_interval_covered']
    assert report['selected_material'] is None and not report['calibrated']
    assert set(report['reference_diagnostics']['thin'])=={'experiment_1','experiment_2'}
    assert report['paper']['paper_fe_is_experimental_target'] is False


@pytest.mark.parametrize('extra',[{'stress_area_convention_mm2':1.}, {'config_wall_mm':1.}])
def test_incompatible_reference_conventions_rejected(tmp_path,extra):
    with pytest.raises(ValueError):
        build_report([assessment(tmp_path,**extra)],ROOT/'validation/paper_reference')
