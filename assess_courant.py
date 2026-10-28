"""Assess OpenCourant termination, dense energies and optional force history.

Convert T01 with the bundled th_to_csv executable first. The force column and
sign must be selected from verified interface output; opaque column numbers
are not guessed. No curve or quasi-static agreement is manufactured.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from reentrant_fea.radioss_results import (read_engine_out, read_history_csv,
                                          history_energy_screen, part_energy_screen)


def assess(metadata_path, history_path, *, force_column=None, tangential_column=None,
           force_sign=None, force_kind=None):
    """Assess bound output and optionally reconstruct a verified force quantity.

    For this writer's TYPE25 IMS=2 convention, force on the secondary is
    ``FN-FT``. The optional tangential channel is subtracted before sign
    selection and before one time difference of cumulative impulses. A legacy
    single normal channel remains a normal-component diagnostic, not a total
    contact force. No column numbers or force/impulse units are inferred.
    """
    if tangential_column is not None and force_column is None:
        raise ValueError('A tangential column requires a normal force column')
    if tangential_column is not None and tangential_column == force_column:
        raise ValueError('Normal and tangential columns must be different')
    metadata_path = Path(metadata_path)
    metadata = json.loads(metadata_path.read_text())
    directory = metadata_path.parent
    canonical_history = directory/(metadata['prefix']+'T01.csv')
    if Path(history_path).resolve() != canonical_history.resolve():
        raise ValueError('History CSV must belong to the same model prefix')
    for name, expected in metadata['sha256'].items():
        path = directory/name
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError('Input deck changed after generation: '+str(path))
    engine = directory/(metadata['prefix']+'_0001.out')
    evidence = read_engine_out(engine)
    history = read_history_csv(history_path)
    binding_path = directory/(metadata['prefix']+'.run.json')
    bound = binding_path.is_file()
    metadata_bound = False
    if bound:
        binding = json.loads(binding_path.read_text())
        expected_metadata = binding.get('metadata_sha256')
        if expected_metadata:
            if hashlib.sha256(metadata_path.read_bytes()).hexdigest() != expected_metadata:
                raise ValueError('Model metadata changed since execution')
            metadata_bound = True
        if binding.get('input_sha256') != metadata['sha256']:
            raise ValueError('Execution record belongs to different inputs')
        expected_outputs = {metadata['prefix']+suffix for suffix in ('_0000.out', '_0001.out', 'T01', 'T01.csv')}
        if set(binding.get('output_sha256', {})) != expected_outputs:
            raise ValueError('Execution record has missing solver evidence')
        for name, expected in binding['output_sha256'].items():
            if hashlib.sha256((directory/name).read_bytes()).hexdigest() != expected:
                raise ValueError('Solver output changed since execution: '+name)
    duration = metadata['duration_seconds']
    dt = metadata.get('time_history_interval_seconds', duration/400)
    final = float(history['time'][-1])
    reached = abs(final-duration) <= max(1.5*dt, duration*1e-5)
    complete = bound and evidence['normal_termination'] and reached and not evidence['errors']
    result = {
        'backend': metadata['backend'], 'completed': bool(complete),
        'engine_normal_termination': evidence['normal_termination'],
        'errors': evidence['errors'],
        'run_binding_verified': bound,
        'metadata_binding_verified': metadata_bound,
        'duration_seconds': duration, 'last_observed_history_time': final,
        'target_nominal_strain': metadata.get('target_nominal_strain'),
        'last_printed_engine_time': evidence['last_printed_time'],
        'deleted_solid_elements': evidence['deleted_solid_elements'],
        'config_wall_mm': metadata.get('wall_mm'),
        'stress_area_convention_mm2': metadata.get('reference_area_mm2'),
        'history_samples': len(history['time']),
        'global_energy_screen': history_energy_screen(history, duration),
        'core_energy_screen': part_energy_screen(history, 1, duration),
        'input_sha256': metadata['sha256'],
        'history_sha256': hashlib.sha256(Path(history_path).read_bytes()).hexdigest(),
        'metadata_sha256': hashlib.sha256(metadata_path.read_bytes()).hexdigest(),
        'calibrated': False, 'experimental_agreement': 'not established',
        'quasi_static_validated': False,
        'reaction_history': [],
        'qualification': 'solver evidence only; independent material, contact and quasi-static qualification required',
    }
    if force_column is not None:
        if force_column not in history or force_sign not in (-1, 1) or force_kind not in ('force', 'impulse'):
            raise ValueError('Select an existing channel, explicit sign and verified force/impulse kind')
        if tangential_column is not None and tangential_column not in history:
            raise ValueError('Selected tangential column is missing from the history CSV')
        amplitude = np.asarray(metadata.get('displacement_fraction_samples', []), dtype=float)
        if (amplitude.ndim != 2 or amplitude.shape[1] != 2 or len(amplitude) < 2
                or not np.isfinite(amplitude).all() or np.any(np.diff(amplitude[:, 0]) <= 0)):
            raise ValueError('Actual prescribed displacement curve is required for nominal strain')
        area = float(metadata['reference_area_mm2'])
        if not np.isfinite(area) or area <= 0:
            raise ValueError('Require positive nominal reference area')
        # Never extrapolate prescribed motion outside the generated curve.
        times = history['time']
        forces = history[force_column]
        if tangential_column is not None:
            forces = forces-history[tangential_column]
        forces = forces*force_sign
        if force_kind == 'impulse':
            # REACY histories can be cumulative impulse, not force. These
            # differences are interval-averaged forces, not pointwise peaks.
            forces = np.diff(forces)/np.diff(times)
            times = (times[1:]+times[:-1])/2
        selected = (times >= amplitude[0, 0]) & (times <= amplitude[-1, 0])
        times, forces = times[selected], forces[selected]
        strains = np.interp(times, amplitude[:, 0], amplitude[:, 1])*metadata['target_nominal_strain']
        result['reaction_history'] = [dict(time=float(t), strain=float(e),
                                           reaction_y_N=float(f), nominal_stress_MPa=float(f/area))
                                      for t, e, f in zip(times, strains, forces)]
        result['force_mapping'] = dict(csv_column=force_column, compressive_sign=force_sign,
                                       area_mm2=area, channel_kind=force_kind,
                                       sampling='interval average' if force_kind == 'impulse' else 'sampled force',
                                       verification='caller-selected verified channel and physical quantity')
        if tangential_column is not None:
            result['force_mapping'].update(
                tangential_csv_column=tangential_column, component_operation='FN-FT',
                contact_force_quantity='total_force_on_secondary',
                convention='TYPE25 IMS=2 secondary force; subtract FT from FN before sign and impulse differentiation')
        else:
            normal = (force_column.startswith('Contact force and energy')
                      or force_column.lower() in ('fn', 'fny'))
            result['force_mapping'].update(
                component_operation='single_channel',
                contact_force_quantity='normal_component_only' if normal else 'single_selected_channel',
                total_contact_force_verified=False)
        result['last_recorded_force_nominal_strain'] = float(strains[-1]) if len(strains) else None
        sample_times = history['time'][(history['time'] >= amplitude[0, 0])
                                       & (history['time'] <= amplitude[-1, 0])]
        result['last_observed_nominal_strain'] = (
            float(np.interp(sample_times[-1], amplitude[:, 0], amplitude[:, 1])
                  *metadata['target_nominal_strain']) if len(sample_times) else None)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('metadata', type=Path)
    parser.add_argument('history_csv', type=Path)
    parser.add_argument('--force-column')
    parser.add_argument('--tangential-column',
                        help='Verified TYPE25 IMS=2 FT channel; total secondary force is FN-FT')
    parser.add_argument('--force-sign', type=int, choices=(-1, 1))
    parser.add_argument('--force-kind', choices=('force', 'impulse'),
                        help='Required with a channel; nodal REACY may be cumulative impulse')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = assess(args.metadata, args.history_csv,
                    force_column=args.force_column, tangential_column=args.tangential_column,
                    force_sign=args.force_sign, force_kind=args.force_kind)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ('completed', 'history_samples',
                                                 'quasi_static_validated', 'global_energy_screen')}, indent=2))


if __name__ == '__main__':
    main()
