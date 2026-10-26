"""Read explicit-solver evidence without equating termination with validation.

Engine print tables are sampled diagnostics. They cannot prove quasi-static
response, contact completeness, or experimental agreement. In particular the
configured FINAL TIME in a header is not an observed result timestamp.
"""
import re
import csv
from pathlib import Path

import numpy as np


def read_history_csv(path):
    """Read the official OpenCourant th_to_csv output without dropping columns."""
    with Path(path).open(newline='') as stream:
        reader = csv.reader(stream)
        names = [name.strip() for name in next(reader)]
        if not names or names[0].lower() != 'time' or len(set(names)) != len(names):
            raise ValueError('Expected uniquely named time-history columns beginning with time')
        rows = []
        for row in reader:
            if not row:
                continue
            if len(row) != len(names):
                raise ValueError('Incomplete time-history CSV row')
            rows.append([float(value.replace('D', 'E')) for value in row])
    values = np.asarray(rows, dtype=float)
    if (values.ndim != 2 or len(values) < 2 or not np.isfinite(values).all()
            or values[0, 0] < 0 or np.any(np.diff(values[:, 0]) <= 0)):
        raise ValueError('Require finite increasing time-history samples')
    return {name: values[:, i] for i, name in enumerate(names)}


def history_energy_screen(history, end_time, *, kinetic_ratio_limit=.05,
                          hourglass_ratio_limit=.05, added_mass_ratio_limit=.01,
                          startup_fraction=.05):
    """Dense global-history screen, independent of the engine print frequency.

    Completion must be checked against engine termination separately. The
    core/platen partition, energy balance and speed sensitivity still need review.
    """
    controls = [end_time, kinetic_ratio_limit, hourglass_ratio_limit,
                added_mass_ratio_limit, startup_fraction]
    if (not np.isfinite(controls).all() or end_time <= 0 or kinetic_ratio_limit <= 0
            or hourglass_ratio_limit <= 0 or added_mass_ratio_limit < 0
            or not 0 <= startup_fraction < .9):
        raise ValueError('Invalid history-screen controls')
    required = ['time', 'INTERNAL ENERGY', 'KINETIC ENERGY', 'ROTATION ENERGY',
                'HOURGLASS ENERGY', 'MASS', 'ADDED MASS', 'EXTERNAL WORK', 'CONTACT ENERGY']
    missing = [name for name in required if name not in history]
    if missing:
        return {'sampled_energy_screen_passed': False,
                'reasons': ['missing global histories: '+', '.join(missing)],
                'quasi_static_validated': False}
    t = history['time']
    selected = t >= startup_fraction*end_time
    reasons = []
    if selected.sum() < 20:
        reasons.append('insufficient dense energy history')
    if t[-1] < .995*end_time:
        reasons.append('time history does not cover the end of loading')
    if np.any(t > end_time*1.01):
        reasons.append('time history exceeds requested run time')
    internal = history['INTERNAL ENERGY'][selected]
    kinetic = (history['KINETIC ENERGY']+history['ROTATION ENERGY'])[selected]
    hg = history['HOURGLASS ENERGY'][selected]
    mass = history['MASS'][selected]
    added = history['ADDED MASS'][selected]
    # Converters may report roundoff-size negative added mass (about 1e-16).
    if np.any(added < -np.maximum(np.abs(mass)*1e-10, 1e-30)):
        reasons.append('negative added mass exceeds roundoff')
    added = np.maximum(added, 0)
    valid = (internal > 0) & (mass > added) & (kinetic >= 0)
    if not np.all(valid) or not len(internal):
        reasons.append('undefined energy or added-mass ratio')
    peak = lambda values: float(np.max(values)) if len(values) else None
    kr = kinetic[valid]/internal[valid]
    hr = np.abs(hg[valid])/internal[valid]
    mr = added[valid]/(mass[valid]-added[valid])
    if len(kr) and max(kr) > kinetic_ratio_limit:
        reasons.append('global kinetic/internal energy ratio exceeds limit')
    if len(hr) and max(hr) > hourglass_ratio_limit:
        reasons.append('global hourglass/internal energy ratio exceeds limit')
    if len(mr) and max(mr) > added_mass_ratio_limit:
        reasons.append('added/original mass ratio exceeds limit')
    return {'sampled_energy_screen_passed': not reasons, 'reasons': reasons,
            'samples_after_startup': int(selected.sum()), 'last_observed_time': float(t[-1]),
            'max_kinetic_internal_ratio': peak(kr),
            'max_hourglass_internal_ratio': peak(hr), 'max_added_original_mass_ratio': peak(mr),
            'thresholds': dict(kinetic_internal=kinetic_ratio_limit,
                               hourglass_internal=hourglass_ratio_limit,
                               added_original_mass=added_mass_ratio_limit,
                               ignored_startup_fraction=startup_fraction),
            'quasi_static_validated': False,
            'qualification': 'global screen only; review core histories, contact/energy balance and loading-speed sensitivity'}


def part_energy_screen(history, part_id, end_time, *, kinetic_ratio_limit=.05,
                       hourglass_ratio_limit=.05, startup_fraction=.05):
    """Screen core energies separately from prescribed platen kinetic energy."""
    if (not np.isfinite([end_time, kinetic_ratio_limit, hourglass_ratio_limit, startup_fraction]).all()
            or end_time <= 0 or min(kinetic_ratio_limit, hourglass_ratio_limit) <= 0
            or not 0 <= startup_fraction < .9):
        raise ValueError('Invalid part-energy controls')
    columns = {}
    for name in history:
        fields = name.split()
        if len(fields) == 3 and fields[:2] == ['Part', str(part_id)]:
            columns[fields[2]] = name
    required = ('IE', 'KE', 'HE', 'ERODED')
    missing = [label for label in required if label not in columns]
    if missing:
        return {'part_id': part_id, 'sampled_energy_screen_passed': False,
                'reasons': ['missing part histories: '+', '.join(missing)],
                'quasi_static_validated': False}
    selected = history['time'] >= startup_fraction*end_time
    ie, ke, he = (history[columns[label]][selected] for label in ('IE', 'KE', 'HE'))
    reasons = []
    if selected.sum() < 20:
        reasons.append('insufficient dense part-energy history')
    if history['time'][-1] < .995*end_time:
        reasons.append('part history does not cover the end of loading')
    valid = (ie > 0) & (ke >= 0)
    if not valid.all() or not len(ie):
        reasons.append('undefined part kinetic/internal ratio')
    kr, hr = ke[valid]/ie[valid], np.abs(he[valid])/ie[valid]
    if len(kr) and max(kr) > kinetic_ratio_limit:
        reasons.append('core kinetic/internal energy ratio exceeds limit')
    if len(hr) and max(hr) > hourglass_ratio_limit:
        reasons.append('core hourglass/internal energy ratio exceeds limit')
    return {'part_id': part_id, 'sampled_energy_screen_passed': not reasons, 'reasons': reasons,
            'max_kinetic_internal_ratio': float(max(kr)) if len(kr) else None,
            'max_hourglass_internal_ratio': float(max(hr)) if len(hr) else None,
            'last_recorded_eroded_elements': float(history[columns['ERODED']][-1]),
            'quasi_static_validated': False,
            'qualification': 'core screen only; balance, contact/erosion work and speed sensitivity remain required'}


NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EDed][+-]?\d+)?"
ROW = re.compile(
    rf"^\s*(\d+)\s+({NUMBER})\s+({NUMBER})\s+\S+\s+\d+\s+"
    rf"({NUMBER})%\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})\s+"
    rf"({NUMBER})\s+({NUMBER})(?:\s+({NUMBER})\s+({NUMBER}))?\s*$")


def normal_termination(text):
    """Require the solver's standalone success line, never ABNORMAL TERMINATION."""
    return bool(re.search(r"^\s*NORMAL TERMINATION\s*$", text, re.MULTILINE))


def solver_errors(text):
    """Some native failures retain OS exit zero and even a normal trailer."""
    pattern = r'^\s*(?:\*+\s*ERROR\b|ERROR\s+ID\b|(?:ERROR|ABNORMAL)\s+TERMINATION\b)'
    return [line.strip() for line in text.splitlines() if re.search(pattern, line, re.I)]


def read_engine_out(path):
    """Parse Radioss I/K-energy print tables and explicit deletion messages.

    Energy/error column positions follow the engine's CYCLE/TIME table, not
    arbitrary numeric lines. The optional mass columns are retained when present.
    """
    text = Path(path).read_text(errors='replace')
    samples, deleted = [], []
    active = False
    for line in text.splitlines():
        if all(label in line for label in ('CYCLE', 'TIME-STEP', 'I-ENERGY', 'K-ENERGY T')):
            active = True
            continue
        event = re.search(r'DELETE SOLID ELEMENT NUMBER\s+(\d+)', line)
        if not event:
            event = re.search(r'RUPTURE OF SOLID ELEMENT\s*:\s*(\d+)', line)
        if event:
            deleted.append(int(event.group(1)))
        if not active:
            continue
        match = ROW.match(line)
        if not match:
            continue
        values = [float(value.replace('D', 'E').replace('d', 'e'))
                  if value is not None else None for value in match.groups()]
        cycle, time, step, error, internal, kt, kr, work, mass_error, mass, added = values
        numeric = [v for v in values if v is not None]
        if not np.isfinite(numeric).all() or time < 0 or step <= 0:
            raise ValueError('Nonfinite or invalid engine energy row')
        # Printed times have only a few significant figures; successive cycles
        # may legitimately have the same rounded timestamp during small steps.
        if samples and (cycle <= samples[-1]['cycle'] or time < samples[-1]['time']):
            raise ValueError('Engine rows are not one increasing run; assess restarts separately')
        samples.append(dict(cycle=int(cycle), time=time, time_step=step,
                            energy_error_percent=error, internal_energy=internal,
                            kinetic_energy_translation=kt, kinetic_energy_rotation=kr,
                            external_work=work, mass_error=mass_error,
                            total_mass=mass, added_mass=added))
    return {'source': str(Path(path).resolve()), 'normal_termination': normal_termination(text),
            'errors': solver_errors(text),
            'samples': samples, 'deleted_solid_elements': sorted(set(deleted)),
            'last_printed_time': samples[-1]['time'] if samples else None,
            'limitations': ['Print-table sampling can miss energy peaks and the final state.',
                            'Deletion messages alone do not verify fresh-surface contact.']}


def energy_screen(evidence, end_time, *, kinetic_ratio_limit=.05,
                  energy_error_limit_percent=5., mass_added_ratio_limit=.01,
                  startup_fraction=.05, minimum_samples=5):
    """Conservative screen of sampled global energies; never a validation claim.

    Thresholds are investigator choices, not values measured in Dong's paper.
    Whole-model energies include platen motion: a part-level core audit, dense
    histories, hourglass/contact energies and loading-speed sensitivity remain
    necessary even when this screen passes.
    """
    controls = [end_time, kinetic_ratio_limit, energy_error_limit_percent,
                mass_added_ratio_limit, startup_fraction]
    if (not np.isfinite(controls).all() or end_time <= 0 or kinetic_ratio_limit <= 0
            or energy_error_limit_percent <= 0 or mass_added_ratio_limit < 0
            or not 0 <= startup_fraction < .9 or minimum_samples < 2):
        raise ValueError('Invalid explicit energy-screen controls')
    rows = [row for row in evidence['samples'] if row['time'] >= startup_fraction*end_time]
    reasons = []
    if not evidence['normal_termination']:
        reasons.append('engine has no normal termination')
    if evidence.get('errors'):
        reasons.append('engine reports errors')
    if len(rows) < minimum_samples:
        reasons.append('insufficient sampled energy history')
    if not rows or rows[-1]['time'] < .95*end_time:
        reasons.append('printed energy history does not cover the end of loading')
    if rows and any(row['time'] > end_time*(1+1e-3) for row in rows):
        reasons.append('energy history exceeds the requested run time')
    ratios, errors, added = [], [], []
    for row in rows:
        internal = row['internal_energy']
        kinetic = row['kinetic_energy_translation'] + row['kinetic_energy_rotation']
        if internal <= 0 or kinetic < 0:
            reasons.append('undefined kinetic/internal ratio or negative kinetic energy')
        else:
            ratios.append(kinetic/internal)
        errors.append(abs(row['energy_error_percent']))
        if row['total_mass'] is None or row['added_mass'] is None:
            reasons.append('missing mass-scaling evidence')
        elif row['total_mass'] <= 0 or not 0 <= row['added_mass'] < row['total_mass']:
            reasons.append('invalid mass-scaling evidence')
        else:
            added.append(row['added_mass']/(row['total_mass']-row['added_mass']))
    peak = lambda values: max(values) if values else None
    if ratios and max(ratios) > kinetic_ratio_limit:
        reasons.append('sampled kinetic/internal energy ratio exceeds limit')
    if errors and max(errors) > energy_error_limit_percent:
        reasons.append('sampled global energy error exceeds limit')
    if added and max(added) > mass_added_ratio_limit:
        reasons.append('added/original mass ratio exceeds limit')
    return {'sampled_energy_screen_passed': not reasons, 'reasons': sorted(set(reasons)),
            'sample_count_after_startup': len(rows),
            'max_kinetic_internal_ratio': peak(ratios),
            'max_abs_energy_error_percent': peak(errors),
            'max_added_original_mass_ratio': peak(added),
            'thresholds': dict(kinetic_internal_ratio=kinetic_ratio_limit,
                               abs_energy_error_percent=energy_error_limit_percent,
                               added_original_mass_ratio=mass_added_ratio_limit,
                               ignored_startup_fraction=startup_fraction),
            'quasi_static_validated': False,
            'qualification': 'sampled screen only; core energies, hourglass/contact work and speed sensitivity still required'}
