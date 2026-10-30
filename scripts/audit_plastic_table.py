"""Audit a supplied Abaqus *Plastic block against the source, using exact decimals.

This records equality and provenance; it does not certify experimental origin,
convert engineering data or change any material/solver input.
"""
import argparse
import ast
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re


NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'


def decimal_value(token):
    if not re.fullmatch(NUMBER, token):
        raise ValueError('Expected a finite decimal number: '+token)
    return Decimal(token.replace('D', 'E').replace('d', 'e'))


def supplied_rows(text):
    lines = text.splitlines()
    headers = [i for i, line in enumerate(lines) if re.match(r'^\s*\*Plastic(?:\s|,|$)', line, re.I)]
    if len(headers) != 1:
        raise ValueError('Require exactly one explicitly labelled *Plastic block')
    rows = []
    for line in lines[headers[0]+1:]:
        line = line.strip()
        if not line or line.startswith('**'):
            continue
        if line.startswith('*'):
            break
        fields = [value.strip() for value in line.split(',')]
        if len(fields) != 2:
            raise ValueError('Every plastic row requires exactly two comma-separated values')
        rows.append(tuple(decimal_value(value) for value in fields))
    if len(rows) < 2:
        raise ValueError('Require at least two material points')
    return rows


def source_rows(text):
    match = re.search(r'\bPLASTIC_DATA\s*=\s*\(', text)
    if match is None:
        raise ValueError('PLASTIC_DATA assignment missing')
    start = text.index('(', match.start())
    depth = 0
    for index in range(start, len(text)):
        depth += (text[index] == '(') - (text[index] == ')')
        if depth == 0:
            literal = text[start:index+1]
            break
    else:
        raise ValueError('Unclosed PLASTIC_DATA tuple')
    tree = ast.parse(literal, mode='eval').body
    if not isinstance(tree, ast.Tuple):
        raise ValueError('Expected a literal material tuple')
    rows = []
    for row in tree.elts:
        if not isinstance(row, ast.Tuple) or len(row.elts) != 2:
            raise ValueError('Every source row requires two numeric literals')
        rows.append(tuple(decimal_value(ast.get_source_segment(literal, value)) for value in row.elts))
    return rows


def canonical(rows):
    def number(value):
        if value == 0:
            return '0'
        result = format(value, 'f')
        return result.rstrip('0').rstrip('.') if '.' in result else result
    return ''.join(','.join(number(value) for value in row)+'\n' for row in rows).encode('ascii')


def audit(supplied_path, source_path):
    supplied_path, source_path = Path(supplied_path), Path(source_path)
    user = supplied_rows(supplied_path.read_text(encoding='utf-8-sig'))
    existing = source_rows(source_path.read_text(encoding='utf-8-sig'))
    mismatches = [dict(row_1based=i+1, user=[str(value) for value in a],
                       existing=[str(value) for value in b])
                  for i, (a, b) in enumerate(zip(user, existing)) if a != b]
    identical = len(user) == len(existing) and not mismatches
    rounded_matches = len(user) == len(existing) and all(
        a == Decimal(format(float(b), '.6g'))
        for user_row, source_row in zip(user, existing) for a, b in zip(user_row, source_row))
    columns = []
    for column, name in enumerate(('true_stress_MPa', 'true_plastic_strain')):
        pairs = [(a[column], b[column]) for a, b in zip(user, existing)]
        deltas = [abs(a-b) for a, b in pairs]
        maximum = max(deltas, default=Decimal(0))
        columns.append(dict(column=name, exact_matching_value_count=sum(a == b for a, b in pairs),
                            differing_value_count=sum(a != b for a, b in pairs),
                            maximum_absolute_difference=str(maximum),
                            maximum_difference_rows_1based=[i+1 for i, value in enumerate(deltas) if value == maximum],
                            matches_existing_binary64_six_significant_digits=all(
                                a == Decimal(format(float(b), '.6g')) for a, b in pairs)))
    digest = lambda data: hashlib.sha256(data).hexdigest()
    return dict(
        attribution='User supplied the complete *Plastic block in this conversation (2026-10-26).',
        data_interpretation='Abaqus *Plastic ordering: true stress (MPa), equivalent true plastic strain (dimensionless), as supplied by user and already labelled in the repository.',
        independent_experimental_provenance='Not independently established; this audit does not claim a paper page, raw-test archive or digitization source for the table.',
        comparison_method='Exact decimal numeric comparison of every row and both columns; no float tolerance or endpoint-only check.',
        numeric_sha256_encoding='ASCII stress,plastic_strain with one LF per row; fixed-point decimals, insignificant trailing fractional zeros removed, zero normalized to 0; no header.',
        supplied_file=str(supplied_path), existing_source=str(source_path),
        supplied_file_sha256=digest(supplied_path.read_bytes()), existing_source_file_sha256=digest(source_path.read_bytes()),
        supplied_numeric_sha256=digest(canonical(user)), existing_numeric_sha256=digest(canonical(existing)),
        supplied_row_count=len(user), existing_row_count=len(existing),
        compared_numeric_value_count=2*min(len(user), len(existing)),
        row_count_matches=len(user) == len(existing), mismatching_rows=mismatches,
        mismatching_row_count=len(mismatches), numerically_identical=identical,
        column_comparison=columns,
        six_significant_digit_diagnostic=dict(all_supplied_values_match=rounded_matches,
            method="For every source decimal value b compare supplied decimal to Decimal(format(float(b), '.6g')); float is the source's binary64 representation.",
            interpretation='A formatting-consistency check, not exact numeric equality or independent experimental verification.'),
        audit_changes_solver_inputs=False,
        result_interpretation=('All constitutive input values match; this is a provenance audit and does not require a physics rerun.'
                               if identical else ('The supplied table is not exactly identical, but every value matches six-significant-digit binary64 formatting of the existing table. This audit preserves the existing higher-precision solver input.'
                                                  if rounded_matches else 'Material data differs; review before replacing any existing input.')),
        unresolved=['Fracture-locus coordinates outside the located point are not supplied by a hardening table.',
                    'Hardening-data provenance does not establish full-lattice calibration or experimental agreement.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('supplied_table', type=Path)
    parser.add_argument('--source', type=Path, default=Path('reference_setup_numeric_source.py'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.supplied_table, args.source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ('supplied_row_count', 'existing_row_count',
                                                  'mismatching_row_count', 'numerically_identical',
                                                  'supplied_numeric_sha256')}, indent=2))


if __name__ == '__main__':
    main()
