"""Cards for the compiled finite-strain CCX 2.20 user material.

The paper does not publish the complete fracture locus. Callers must explicitly
supply their assumptions; no calibrated fracture table is silently invented.
"""
import numpy as np
from .plastic import validate_material


def material_cards(nodes, elements, config, plastic):
    table = validate_material(config, plastic)
    law = config['ductile_damage']
    locus = np.asarray(law['fracture_locus'], dtype=float)
    if (locus.ndim != 2 or locus.shape[1] != 2 or len(locus) < 2
            or not np.isfinite(locus).all() or np.any(np.diff(locus[:, 0]) <= 0)
            or np.any(locus[:, 1] <= 0)):
        raise ValueError('Invalid explicitly supplied fracture locus')
    uf = float(law['failure_displacement'])
    residual = float(law.get('residual_stiffness', 1e-4))
    cutoff = float(law.get('compression_cutoff', -1/3))
    viscosity = float(law.get('viscosity_time', 0.))
    if not np.isfinite(viscosity) or viscosity < 0:
        raise ValueError('viscosity_time must be nonnegative')
    if not np.isfinite([uf, residual, cutoff]).all() or uf <= 0 or not 0 < residual < 1:
        raise ValueError('Invalid ductile softening parameters')
    if not law.get('provenance'):
        raise ValueError('Damage parameters require provenance, including assumed triaxialities')
    groups = {}
    lengths = []
    for eid, el in enumerate(elements, 1):
        xyz = nodes[el]
        volume = 0.
        tets=(((0,1,2,3),(1,2,4,3),(2,4,5,3)) if len(el)==6 else
              ((0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6),(0,5,1,6)))
        for ids in tets:
            tet = xyz[list(ids)]
            volume += abs(np.linalg.det((tet[1:]-tet[0]).T))/6
        if volume <= 0:
            raise ValueError('Zero-volume element in damage mesh')
        length = volume**(1/3)
        # Three significant figures bounds length binning error to 0.5%.
        length = float(f'{length:.3g}')
        lengths.append(length)
        groups.setdefault(length, []).append(eid)
    lines = []
    for index, (length, ids) in enumerate(sorted(groups.items())):
        name = '@REENTRANT@' + str(index)
        # Suffix after @ is a material-instance identifier; CCX loads the same
        # libREENTRANT.so for all characteristic-length groups.
        props = [config['young'], config['poisson'], uf, residual, length,
                 -len(table), len(locus), cutoff, viscosity]
        lines += [f'*MATERIAL,NAME={name}', f'*USER MATERIAL,CONSTANTS={len(props)}']
        lines += [', '.join(f'{v:.12g}' for v in props[i:i+8]) for i in range(0,len(props),8)]
        lines += ['*DEPVAR', '10', f'*ELSET,ELSET=DAMAGE_{index}']
        lines += [', '.join(str(v) for v in ids[i:i+16]) for i in range(0,len(ids),16)]
        lines += [f'*SOLID SECTION,ELSET=DAMAGE_{index},MATERIAL={name}']
    sidecar='REENTRANT_TABLES_V1 '+str(len(table))+' '+str(len(locus))+'\n'
    sidecar+=''.join(f'{a:.12g} {b:.12g}\n' for a,b in [*table.tolist(),*locus.tolist()])
    return lines, sidecar, {'length_min':min(lengths), 'length_max':max(lengths),
                   'material_groups':len(groups), 'provenance':law['provenance']}
