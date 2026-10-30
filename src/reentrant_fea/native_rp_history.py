"""Observed RP strain and signed interval reactions for the native T01 writer."""
import re
import numpy as np
from .force_balance import interval_force

VARIABLES = ['DX', 'DY', 'DZ', 'REACX', 'REACY', 'REACZ']


def rp_reaction_history(history, metadata):
    if metadata.get('platen_reference_points') is not True:
        raise ValueError('RP assessment requires generated platen_reference_points=true')
    if metadata.get('loading_node_variables') != VARIABLES:
        raise ValueError('RP assessment requires the verified six-component nodal request')
    top = metadata.get('loading_face_nodes', [])
    bottom = metadata.get('bottom_support_nodes', [])
    if len(top) != 1 or len(bottom) != 1 or top == bottom:
        raise ValueError('Require distinct single top and bottom reference nodes')
    height, area = float(metadata['height_mm']), float(metadata['reference_area_mm2'])
    if not np.isfinite([height, area]).all() or min(height, area) <= 0:
        raise ValueError('Require positive initial core height and reference area')
    time = np.asarray(history['time'], float)
    def channels(node):
        pattern = re.compile(r'^Loading face nodal histories\s+'+str(node)+r'\s+load-node-'+str(node)+r'\s+var (\d+)$')
        selected = [(int(m.group(1)), key) for key in history if (m := pattern.match(key))]
        selected.sort()
        if len(selected) != 6 or [n for n, _ in selected] != list(range(selected[0][0], selected[0][0]+6)):
            raise ValueError('Missing or ambiguous six-channel history for reference node '+str(node))
        values = [np.asarray(history[key], float) for _, key in selected]
        if any(v.shape != time.shape or not np.isfinite(v).all() for v in values):
            raise ValueError('Nonfinite or mismatched RP history')
        return values, [key for _, key in selected]
    ut, kt = channels(top[0]); ub, kb = channels(bottom[0])
    mid, ft = interval_force(time, ut[4]); _, fb = interval_force(time, ub[4])
    strain = -(ut[1]-ub[1])/height
    if not np.isfinite(strain).all() or len(strain)<2:
        raise ValueError('Invalid observed RP strain')
    emid = (strain[:-1]+strain[1:])/2
    rows = [dict(time=float(t), strain=float(e), reaction_y_N=float(f),
                 bottom_reaction_y_N=float(b), nominal_stress_MPa=float(-f/area),
                 top_bottom_reaction_sum_y_N=float(f+b))
            for t,e,f,b in zip(mid,emid,ft,fb)]
    return dict(reaction_history=rows,
                last_observed_nominal_strain=float(strain[-1]),
                last_recorded_force_nominal_strain=float(emid[-1]),
                prescribed_motion_reached=bool(abs(strain[-1]-metadata['target_nominal_strain'])<=max(1e-6,metadata['target_nominal_strain']*1e-4)),
                force_mapping=dict(quantity='platen_reference_point_reaction', channel_kind='impulse',
                    sampling='interval average; differentiate cumulative REACY once',
                    top_node=top[0],bottom_node=bottom[0],
                    top_reaction_column=kt[4],bottom_reaction_column=kb[4],
                    top_displacement_column=kt[1],bottom_displacement_column=kb[1],
                    strain_convention='-(observed U2_top-observed U2_bottom)/initial_core_height',
                    stress_convention='-signed top reaction / initial_reference_area',area_mm2=area,
                    balance_note='Signed top+bottom reactions include platen inertia; not a momentum-balance residual'))
