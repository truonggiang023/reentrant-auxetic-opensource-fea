import pytest
from reentrant_fea.material import DuctilePoint

def test_initiation_and_failure():
    p=DuctilePoint()
    p.update(.047,.333,.2,[(.333,.047)],.001)
    assert p.initiation_index==pytest.approx(1)
    assert p.damage==0
    p.update(.005,.333,.2,[(.333,.047)],.001)
    assert p.failed and p.damage==1

def test_partial_crossing():
    p=DuctilePoint()
    p.update(.05,.333,.1,[(.333,.047)],.001)
    assert p.damage==pytest.approx(.3)
