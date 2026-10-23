"""Read solver evidence without equating process exit status with convergence."""
from pathlib import Path
import numpy as np


def completed_increments(path):
    result = []
    for line in Path(path).read_text().splitlines():
        parts = line.split()
        if len(parts) != 7 or not parts[0].isdigit() or not parts[2].isdigit():
            continue  # Failed attempts carry U in the attempt column.
        result.append((int(parts[0]), int(parts[1]), float(parts[5]), float(parts[6])))
    return result


def validate_run(deck, log, end_time=1.):
    deck = Path(deck)
    text = Path(log).read_text(errors='replace')
    if '*ERROR' in text.upper() or 'Job finished' not in text:
        raise RuntimeError('CalculiX did not finish cleanly; inspect ' + str(log))
    increments = completed_increments(deck.with_suffix('.sta'))
    if not increments or abs(increments[-1][2] - end_time) > 1e-6:
        raise RuntimeError('CalculiX did not reach the prescribed final step time')
    for suffix in ('.dat', '.frd'):
        if not deck.with_suffix(suffix).is_file() or deck.with_suffix(suffix).stat().st_size == 0:
            raise RuntimeError('Missing solver output ' + suffix)
    return increments


def integration_history(dat):
    """Yield time, element, integration point, PEEQ, triaxiality from *EL PRINT.

    Uses unaveraged integration-point stresses, never FRD nodal extrapolations.
    Missing/mismatched blocks fail rather than fabricating a stress history.
    """
    blocks = {}
    mode = None
    for line in Path(dat).read_text().splitlines():
        lower = line.lower().strip()
        if 'for set core and time' in lower:
            time = float(lower.split('time')[-1])
            mode = 'stress' if lower.startswith('stresses') else ('peeq' if lower.startswith('equivalent plastic strain') else None)
            continue
        fields = line.split()
        if mode and len(fields) == (8 if mode == 'stress' else 3) and fields[0].isdigit():
            key = (time, int(fields[0]), int(fields[1]))
            blocks.setdefault(key, {})[mode] = np.array([float(x) for x in fields[2:]])
        elif lower and not fields[0].isdigit():
            mode = None
    if not blocks:
        raise ValueError('No CORE integration-point history in DAT')
    for key, values in sorted(blocks.items()):
        if set(values) != {'stress', 'peeq'}:
            raise ValueError('Unpaired stress/PEEQ at ' + str(key))
        if not all(np.isfinite(v).all() for v in values.values()) or values['peeq'][0] < 0:
            raise ValueError('Nonfinite stress or invalid PEEQ')
        s = values['stress']
        mean = s[:3].mean()
        q = np.sqrt(1.5 * (np.sum((s[:3]-mean)**2) + 2*np.sum(s[3:]**2)))
        if q < 1e-12 and values['peeq'][0] > 0:
            raise ValueError('Undefined triaxiality at zero deviatoric stress')
        yield (*key, float(values['peeq'][0]), float(mean/q) if q > 1e-12 else 0.)


def initiation_indices(dat, locus):
    """Uncoupled ductile initiation index ∫dPEEQ/epsilon_f(eta).

    This diagnostic is compatible with stock CCX 2.20 output. It supplies no
    stress degradation or deletion. Requires a calibrated fracture locus;
    extrapolating the single tensile datum to compression is prohibited.
    Output sampling controls integration accuracy; request every increment.
    """
    curve = np.asarray(locus, dtype=float)
    if (curve.ndim != 2 or curve.shape[1] != 2 or len(curve) < 2
            or not np.isfinite(curve).all() or np.any(np.diff(curve[:, 0]) <= 0)
            or np.any(curve[:, 1] <= 0)):
        raise ValueError('Supply a calibrated, increasing triaxiality locus with at least two positive fracture strains')
    states = {}
    output = []
    for time, eid, ip, peeq, eta in integration_history(dat):
        old, omega = states.get((eid, ip), (0., 0.))
        if peeq < old - 1e-9:
            raise ValueError('Nonmonotonic integration-point PEEQ')
        if peeq > old:
            if not curve[0, 0] <= eta <= curve[-1, 0]:
                raise ValueError('Stress triaxiality outside calibrated fracture locus')
            omega += (peeq-old) / np.interp(eta, curve[:, 0], curve[:, 1])
        states[eid, ip] = peeq, omega
        output.append((time, eid, ip, omega))
    return output
