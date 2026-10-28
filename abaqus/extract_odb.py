#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Export signed loading-face reactions from an Abaqus/Explicit ODB.

Run with Abaqus Python (Python 2.7 or 3), not the repository's Python solver:
  abaqus python extract_odb.py model.odb --output curve.csv

Compression is positive: stress=-sum(top RF2)/area, strain=-mean(top U2)/height.
RF2 is an instantaneous nodal constraint reaction, not a cumulative impulse.
An elastic moving platen contributes inertia to its loading-face reaction;
this extraction does not turn a dynamic run into a quasi-static measurement.
"""
from __future__ import print_function

import argparse
import csv
import io
import json
import math
import os
import re
import sys


ENERGIES = ('ALLIE', 'ALLKE', 'ALLAE', 'ALLVD', 'ALLFD', 'ALLWK', 'ETOTAL')
CSV_FIELDS = ('strain', 'stress_mpa', 'time_s', 'compression_mm',
              'top_reaction_y_n', 'bottom_reaction_y_n') + ENERGIES


def finite(value):
    return not (math.isnan(value) or math.isinf(value))


def positive(value, label):
    value = float(value)
    if not finite(value) or value <= 0:
        raise ValueError(label + ' must be finite and positive')
    return value


def identity(instance, label):
    return (str(instance).upper(), int(label))


def node_set_members(odb, name):
    """Return (instance, label) keys, keeping repeated labels in separate parts."""
    sets = odb.rootAssembly.nodeSets
    lookup = dict((str(key).upper(), key) for key in sets.keys())
    if name.upper() not in lookup:
        raise ValueError('ODB has no assembly node set ' + name)
    node_set = sets[lookup[name.upper()]]
    raw = node_set.nodes
    if not len(raw):
        raise ValueError('Empty assembly node set ' + name)
    groups = [raw] if hasattr(raw[0], 'label') else list(raw)
    instances = list(getattr(node_set, 'instanceNames', ()) or ())
    root_instances = list(odb.rootAssembly.instances.keys())
    result = set()
    for index, nodes in enumerate(groups):
        for node in nodes:
            instance = getattr(node, 'instanceName', None)
            if not instance and len(instances) == len(groups):
                instance = instances[index]
            if not instance and len(root_instances) == 1:
                instance = root_instances[0]
            if not instance and not root_instances:
                instance = 'ASSEMBLY'
            if not instance:
                raise ValueError('Ambiguous node instance in set ' + name)
            key = identity(instance, node.label)
            if key in result:
                raise ValueError('Duplicate node in set ' + name + ': ' + str(key))
            result.add(key)
    if not result:
        raise ValueError('Empty assembly node set ' + name)
    return result


def history_identity(name, region):
    match = re.match(r'^Node\s+(.+)\.(\d+)\s*$', str(name), re.I)
    from_name = identity(match.group(1), match.group(2)) if match else None
    point = getattr(region, 'point', None)
    node = getattr(point, 'node', None)
    if node is not None and hasattr(node, 'label'):
        instance = getattr(getattr(point, 'instance', None), 'name', None)
        instance = instance or getattr(node, 'instanceName', None)
        if instance:
            from_point = identity(instance, node.label)
            if from_name is not None and from_name != from_point:
                raise ValueError('Conflicting ODB history node identity: ' + str(name))
            return from_point
    return from_name


def series(output, label):
    data = [(float(t), float(value)) for t, value in output.data]
    if len(data) < 2 or any(not finite(t) or not finite(v) for t, v in data):
        raise ValueError('Need at least two finite samples for ' + label)
    if any(data[i + 1][0] <= data[i][0] for i in range(len(data) - 1)):
        raise ValueError('History times must increase strictly for ' + label)
    return data


def aligned(data, times, label):
    if len(data) != len(times) or any(
            abs(pair[0] - time) > 1.e-10 * max(1., abs(time))
            for pair, time in zip(data, times)):
        raise ValueError('Mismatched history time grids for ' + label +
                         '; request all histories at one common output interval')
    return [pair[1] for pair in data]


def node_histories(step, members, required):
    records = {}
    for name, region in step.historyRegions.items():
        key = history_identity(name, region)
        if key not in members:
            continue
        if key in records:
            raise ValueError('Multiple history regions for node ' + str(key))
        missing = [variable for variable in required if variable not in region.historyOutputs]
        if missing:
            raise ValueError('Missing ' + ', '.join(missing) + ' for node ' + str(key))
        records[key] = dict((variable, series(region.historyOutputs[variable],
                                            str(key) + ' ' + variable))
                            for variable in required)
    missing = members - set(records)
    if missing:
        raise ValueError('Missing nodal history for ' + str(sorted(missing)))
    return records


def extract(odb, step_name='COMPRESS', height_mm=48., reference_area_mm2=1900.,
            wall_mm=.3):
    """Return CSV rows and an assessment; does not import proprietary modules."""
    height_mm = positive(height_mm, 'height_mm')
    reference_area_mm2 = positive(reference_area_mm2, 'reference_area_mm2')
    wall_mm = positive(wall_mm, 'wall_mm')
    if step_name not in odb.steps:
        raise ValueError('ODB has no step ' + step_name)
    step = odb.steps[step_name]
    top = node_set_members(odb, 'TOP_LOAD')
    bottom = node_set_members(odb, 'BOTTOM_LOAD')
    if top & bottom:
        raise ValueError('TOP_LOAD and BOTTOM_LOAD overlap')
    top_history = node_histories(step, top, ('U2', 'RF2'))
    bottom_history = node_histories(step, bottom, ('RF2',))
    first = top_history[sorted(top)[0]]['U2']
    times = [pair[0] for pair in first]
    top_force, bottom_force, displacement = ([0.] * len(times) for _ in range(3))
    for key in sorted(top):
        force = aligned(top_history[key]['RF2'], times, str(key) + ' RF2')
        motion = aligned(top_history[key]['U2'], times, str(key) + ' U2')
        for i in range(len(times)):
            top_force[i] += force[i]
            displacement[i] -= motion[i] / len(top)
    for key in sorted(bottom):
        force = aligned(bottom_history[key]['RF2'], times, str(key) + ' RF2')
        for i in range(len(times)):
            bottom_force[i] += force[i]
    energy_candidates = []
    for name, region in step.historyRegions.items():
        if history_identity(name, region) is not None:
            continue
        outputs = region.historyOutputs
        if all(variable in outputs for variable in ENERGIES):
            energy_candidates.append((name, outputs))
    if len(energy_candidates) != 1:
        raise ValueError('Require exactly one whole-model region containing ' +
                         ', '.join(ENERGIES) + '; found ' + str(len(energy_candidates)))
    energy_name, outputs = energy_candidates[0]
    energies = dict((variable, aligned(series(outputs[variable], variable), times, variable))
                    for variable in ENERGIES)
    rows = []
    for i, time in enumerate(times):
        row = dict(strain=displacement[i] / height_mm,
                   stress_mpa=-top_force[i] / reference_area_mm2,
                   time_s=time, compression_mm=displacement[i],
                   top_reaction_y_n=top_force[i], bottom_reaction_y_n=bottom_force[i])
        row.update((variable, energies[variable][i]) for variable in ENERGIES)
        rows.append(row)
    requested_time = positive(step.timePeriod, 'step.timePeriod')
    last_frame = float(step.frames[-1].frameValue) if step.frames else None
    tolerance = max(1.e-10, requested_time * 1.e-6)
    reached = (last_frame is not None and finite(last_frame) and
               abs(last_frame - requested_time) <= tolerance and
               abs(times[-1] - requested_time) <= tolerance)
    job = getattr(odb, 'jobData', None)
    job_status = str(getattr(job, 'status', 'unavailable'))
    assessment = dict(
        backend='Abaqus/Explicit',
        completed=reached,
        completion_basis='Last ODB frame and history reach the requested step duration; inspect Abaqus status/log separately',
        job_status=job_status, step=step_name, requested_step_time_s=requested_time,
        last_frame_time_s=last_frame, last_history_time_s=times[-1],
        last_nominal_strain=rows[-1]['strain'],
        height_mm=height_mm, reference_area_mm2=reference_area_mm2,
        config_wall_mm=wall_mm,
        dimension_provenance='Caller-supplied dimensions; defaults correspond to the exported t=0.3 mm deck',
        sample_count=len(rows), top_loading_nodes=len(top), bottom_loading_nodes=len(bottom),
        whole_model_energy_history_region=str(energy_name),
        force_quantity='Sum of instantaneous top loading-face nodal RF2; N',
        force_sign_convention='Compression positive: nominal_stress_MPa = -top_reaction_y_N / reference_area_mm2',
        calibrated=False, experimental_agreement=False, quasi_static_validated=False,
        limitations=[
            'Loading-face reaction includes moving elastic-platen inertia; assess force balance and energy histories.',
            'Reaching step duration does not establish numerical convergence or quasi-static behavior.',
            'RF2 is force, not impulse. No time differentiation and no absolute-value sign correction are applied.',
            'Dimensions are caller assertions; this exporter does not verify the ODB geometry against a generated deck.'],
        reaction_history=[dict(strain=row['strain'], nominal_stress_MPa=row['stress_mpa'],
                               time_s=row['time_s'], top_reaction_y_N=row['top_reaction_y_n'],
                               bottom_reaction_y_N=row['bottom_reaction_y_n']) for row in rows])
    return rows, assessment


def write_csv(path, rows):
    # Python 2's csv module needs bytes; Python 3 requires newline=''.
    if sys.version_info[0] < 3:
        stream = open(path, 'wb')
    else:
        stream = io.open(path, 'w', newline='', encoding='utf-8')
    with stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('odb')
    parser.add_argument('--output', default='curve.csv')
    parser.add_argument('--height-mm', type=float, default=48.)
    parser.add_argument('--reference-area-mm2', type=float, default=1900.)
    parser.add_argument('--wall-mm', type=float, default=.3)
    parser.add_argument('--step', default='COMPRESS')
    args = parser.parse_args(argv)
    summary_path = os.path.splitext(args.output)[0] + '.assessment.json'
    if os.path.exists(args.output) or os.path.exists(summary_path):
        parser.error('Output already exists; choose a fresh --output to preserve prior results')
    from odbAccess import openOdb
    odb = openOdb(path=args.odb, readOnly=True)
    try:
        rows, assessment = extract(odb, args.step, args.height_mm,
                                   args.reference_area_mm2, args.wall_mm)
    finally:
        odb.close()
    assessment['source_odb'] = os.path.abspath(args.odb)
    write_csv(args.output, rows)
    with io.open(summary_path, 'w', encoding='utf-8') as stream:
        stream.write(json.dumps(assessment, indent=2, ensure_ascii=True) + u'\n')
    print('Wrote ' + args.output + ' and ' + summary_path)
    print('Step duration reached: ' + str(assessment['completed']) +
          '; quasi-static behavior and experimental agreement are unverified.')


if __name__ == '__main__':
    main()
