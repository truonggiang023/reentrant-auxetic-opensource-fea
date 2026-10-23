from reentrant_fea.plastic import from_reference

def test_reference_table():
    p=from_reference('reference_setup_numeric_source.py')
    assert len(p)==178
    assert p[0][1]==0.0
    assert p[-1][1]>p[0][1]
