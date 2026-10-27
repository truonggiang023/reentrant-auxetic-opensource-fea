from pathlib import Path

import pytest

from reentrant_fea.plastic import from_reference, from_abaqus_plastic, load_plastic

def test_reference_table():
    p=from_reference('reference_setup_numeric_source.py')
    assert len(p)==178
    assert p[0][1]==0.0
    assert p[-1][1]>p[0][1]


def test_supplied_abaqus_table_and_legacy_dispatch_without_execution(tmp_path):
    root = Path(__file__).resolve().parents[1]
    supplied = root/'validation/materials/user_plastic_table_2026-10-26.inp'
    values = load_plastic(supplied)
    assert len(values) == 178
    assert values[0] == (224.668, 0.)
    assert values[-1] == (351.151, .0381444)
    sentinel = tmp_path/'executed'
    legacy = tmp_path/'source.py'
    legacy.write_text(f"from pathlib import Path\nPath({str(sentinel)!r}).write_text('bad')\n"
                      'PLASTIC_DATA=((100,0),(120,.1))\n')
    assert load_plastic(legacy) == ((100, 0), (120, .1))
    assert not sentinel.exists()


@pytest.mark.parametrize('text, reason', [
    ('*Elastic\n69000,.3\n', 'exactly one'),
    ('*Plastic\n100,0\n120,.1\n*Plastic\n100,0\n120,.1\n', 'exactly one'),
    ('*Plastic, DEPENDENCIES=1\n100,0\n120,.1\n', 'rate-independent isotropic'),
    ('*Plastic\n100,0,293\n120,.1,293\n', 'true stress and plastic strain only'),
    ('*Plastic\n100,0\nNaN,.1\n', 'positive nondecreasing'),
    ('*Plastic\n100,.01\n120,.1\n', 'starting at zero'),
    ('*Plastic\n100,0\n90,.1\n', 'nondecreasing'),
    ('*Plastic\n100,0\n120,0\n', 'increasing plastic strain'),
    ("*Plastic\n100,0\n__import__('os').system('bad'),.1\n", 'Nonnumeric'),
])
def test_abaqus_reader_rejects_bad_or_ambiguous_curves(tmp_path, text, reason):
    path = tmp_path/'material.inp'
    path.write_text(text)
    with pytest.raises(ValueError, match=reason):
        from_abaqus_plastic(path)


def test_abaqus_reader_comments_and_fortran_exponents(tmp_path):
    path = tmp_path/'material.inp'
    path.write_text('*Plastic, HARDENING=ISOTROPIC\n** stress, plastic strain\n'
                    '1.0D2,0.\n1.2D2,1.0D-1\n*Density\n2.67e-9\n')
    assert from_abaqus_plastic(path) == ((100., 0.), (120., .1))
