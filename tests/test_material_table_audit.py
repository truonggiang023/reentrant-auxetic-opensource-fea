from pathlib import Path

import pytest

from scripts.audit_plastic_table import audit, canonical, supplied_rows


def test_audit_compares_interior_rows_not_only_endpoints(tmp_path):
    supplied = tmp_path/'user.inp'
    source = tmp_path/'source.py'
    supplied.write_text('*Plastic\n100,0\n121,0.1\n140,0.2\n')
    source.write_text('PLASTIC_DATA=((100,0.0),(120,.1),(140,.2))\n')
    result = audit(supplied, source)
    assert result['mismatching_row_count'] == 1
    assert result['mismatching_rows'][0]['row_1based'] == 2
    assert not result['numerically_identical']
    assert not result['six_significant_digit_diagnostic']['all_supplied_values_match']
    assert not result['audit_changes_solver_inputs']
    assert canonical(supplied_rows('*Plastic\n100.00,0.00\n140,.200\n')) == b'100,0\n140,0.2\n'


def test_actual_user_table_is_rounded_source_not_exactly_equal():
    root = Path(__file__).resolve().parents[1]
    result = audit(root/'validation/materials/user_plastic_table_2026-10-26.inp',
                   root/'reference_setup_numeric_source.py')
    assert result['supplied_row_count'] == result['existing_row_count'] == 178
    assert result['compared_numeric_value_count'] == 356
    assert result['mismatching_row_count'] == 178
    assert not result['numerically_identical']
    assert result['six_significant_digit_diagnostic']['all_supplied_values_match']
    assert [row['differing_value_count'] for row in result['column_comparison']] == [178, 170]
    assert result['supplied_numeric_sha256'] != result['existing_numeric_sha256']


def test_audit_rejects_ambiguous_sections_and_nonnumeric_data():
    with pytest.raises(ValueError, match='exactly one'):
        supplied_rows('*Plastic\n100,0\n120,.1\n*Plastic\n100,0\n120,.1\n')
    with pytest.raises(ValueError, match='finite decimal'):
        supplied_rows('*Plastic\n100,0\nNaN,.1\n')
