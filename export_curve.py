"""Export accepted CalculiX RP history to an engineering stress-strain CSV."""
import argparse
import csv
from pathlib import Path
from reentrant_fea.rp_curve import ccx_rp_history


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('deck', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rows, area = ccx_rp_history(args.deck)
    output = args.output or args.deck.with_suffix('.curve.csv')
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f'{output}: {len(rows)} accepted RP samples; A0={area:g} mm^2')


if __name__ == '__main__':
    main()
