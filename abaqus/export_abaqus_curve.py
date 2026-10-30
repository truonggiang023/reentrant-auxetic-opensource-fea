"""Run with Abaqus Python: export RF2/U2 histories at coupled platen RPs.

Usage: abaqus python export_abaqus_curve.py job.odb --metadata rp_curve_metadata.json
Alternatively supply --height H0 --area A0 --step STEP (mm, mm^2).
Uses standard library only; no editable package installation inside Abaqus needed.
"""
from __future__ import print_function
import argparse
import csv
import json
import math
import sys


def find_region(step, node_set):
    labels = set()
    def collect(items):
        for item in items:
            if hasattr(item, 'label'):
                labels.add(item.label)
            else:
                collect(item)
    collect(node_set.nodes)
    if len(labels) != 1:
        raise ValueError('Expected exactly one node in each reference-point set')
    label = next(iter(labels))
    matches = [region for key, region in step.historyRegions.items()
               if key.startswith('Node ASSEMBLY.') and key.split('.')[-1] == str(label)]
    if len(matches) != 1:
        raise ValueError('Missing or ambiguous assembly RP history region')
    return matches[0]


def curve_rows(top, bottom, height, area):
    if not all((not math.isnan(v) and not math.isinf(v)) for v in (height, area)) or min(height, area) <= 0:
        raise ValueError('Positive finite initial height and reference area required')
    histories = [dict(region.historyOutputs[name].data)
                 for region, name in ((top, 'U2'), (bottom, 'U2'), (top, 'RF2'), (bottom, 'RF2'))]
    times = set(histories[0])
    if not times or any(set(h) != times for h in histories):
        raise ValueError('RF2 and U2 must share the same nonempty sampling times at both RPs')
    rows = []
    for time in sorted(times):
        ut, ub, ft, fb = [h[time] for h in histories]
        if not all((not math.isnan(v) and not math.isinf(v)) for v in (time, ut, ub, ft, fb)):
            raise ValueError('Nonfinite RP history')
        rows.append(dict(time=time, strain=-(ut-ub)/height, reaction_y_N=ft,
                         bottom_reaction_y_N=fb, top_u2_mm=ut, bottom_u2_mm=ub,
                         nominal_stress_MPa=-ft/area, force_balance_y_N=ft+fb))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('odb')
    parser.add_argument('--metadata')
    parser.add_argument('--height', type=float)
    parser.add_argument('--area', type=float)
    parser.add_argument('--step')
    parser.add_argument('--output')
    args = parser.parse_args()
    metadata = {}
    if args.metadata:
        with open(args.metadata) as stream:
            metadata = json.load(stream)
    height = args.height if args.height is not None else metadata.get('initial_height_mm')
    area = args.area if args.area is not None else metadata.get('reference_area_mm2')
    step_name = args.step or metadata.get('step')
    if height is None or area is None or not step_name:
        parser.error('Supply metadata or explicit --height, --area and --step')
    from odbAccess import openOdb
    odb = openOdb(path=args.odb, readOnly=True)
    try:
        step = odb.steps[step_name]
        top = find_region(step, odb.rootAssembly.nodeSets['RP_TOP_PLATEN'])
        bottom = find_region(step, odb.rootAssembly.nodeSets['RP_BOTTOM_PLATEN'])
        rows = curve_rows(top, bottom, height, area)
        output = args.output or args.odb.rsplit('.', 1)[0] + '.curve.csv'
        with (open(output, 'wb') if sys.version_info[0] < 3 else open(output, 'w', newline='')) as stream:
            writer = csv.DictWriter(stream, fieldnames=['time', 'strain', 'reaction_y_N',
                'bottom_reaction_y_N', 'top_u2_mm', 'bottom_u2_mm', 'nominal_stress_MPa', 'force_balance_y_N'])
            writer.writeheader()
            writer.writerows(rows)
        print('%s: %d RP samples; A0=%g mm^2' % (output, len(rows), area))
    finally:
        odb.close()


if __name__ == '__main__':
    main()
