import numpy as np
import pytest
from run_radioss import physical_loading_duration


def test_paper_speed_uses_actual_height_not_centerline_height():
    n=np.array([[0,-.15,0],[76,47.85,25]])
    assert physical_loading_duration(n,.12,2.)==pytest.approx(172.8)
    assert physical_loading_duration(n,.65,2.)==pytest.approx(936.)
    assert physical_loading_duration(n,.78,2.)==pytest.approx(1123.2)


@pytest.mark.parametrize('rate',[0.,-2.,float('nan')])
def test_invalid_speed_rejected(rate):
    with pytest.raises(ValueError):physical_loading_duration([[0,0,0],[1,48,25]],.12,rate)
