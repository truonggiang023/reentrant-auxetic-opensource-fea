"""Compression-positive engineering curves from coupled platen reference points.

No material calibration or local plastic-strain estimate is performed here.
"""
import math
import re
from pathlib import Path


def engineering_row(time, top_u2, bottom_u2, top_rf2, bottom_rf2, height, area):
    values = (time, top_u2, bottom_u2, top_rf2, bottom_rf2, height, area)
    if not all(math.isfinite(v) for v in values) or height <= 0 or area <= 0:
        raise ValueError('Finite RP histories and positive initial height/area required')
    return {'time': time, 'strain': -(top_u2-bottom_u2)/height,
            'reaction_y_N': top_rf2, 'bottom_reaction_y_N': bottom_rf2,
            'top_u2_mm': top_u2, 'bottom_u2_mm': bottom_u2,
            'nominal_stress_MPa': -top_rf2/area,
            'force_balance_y_N': top_rf2+bottom_rf2}


def ccx_rp_history(deck):
    """Read paired DAT U/RF at both RPs, keeping only accepted increments.

    DAT is used because FRD does not reliably include independent rigid-body RPs.
    Legacy direct-node decks must be rerun; no pseudo-time fallback is used.
    """
    from .ccx_results import completed_increments
    from .frd import read_core_deck
    deck = Path(deck)
    _, _, meta = read_core_deck(deck)
    text = deck.read_text()
    match = re.search(r'^\*\* REFERENCE_AREA_MM2=(\S+)', text, re.M)
    if not match:
        raise ValueError('No RP curve metadata; regenerate and rerun the deck')
    area = float(match.group(1))
    accepted = [r[2] for r in completed_increments(deck.with_suffix('.sta'))]
    blocks = {}
    mode = None
    for line in deck.with_suffix('.dat').read_text().splitlines():
        header = re.match(r'\s*(displacements|forces) .* for set (RP_TOP_PLATEN|RP_BOTTOM_PLATEN) and time\s+(\S+)', line)
        if header:
            kind, region, time = header.groups()
            mode = (float(time), region, kind)
            continue
        fields = line.split()
        if mode and fields:
            if len(fields) != 4 or not fields[0].isdigit():
                raise ValueError('Malformed reference-point DAT block')
            time, region, kind = mode
            key = (region, kind)
            frame = blocks.setdefault(time, {})
            if key in frame:
                raise ValueError('Expected a single node per reference-point set')
            frame[key] = float(fields[2])
            mode = None
    rows = []
    for time, frame in sorted(blocks.items()):
        if not any(abs(time-t) <= max(1e-10, abs(t)*6e-6) for t in accepted):
            continue
        try:
            rows.append(engineering_row(time,
                frame['RP_TOP_PLATEN', 'displacements'], frame['RP_BOTTOM_PLATEN', 'displacements'],
                frame['RP_TOP_PLATEN', 'forces'], frame['RP_BOTTOM_PLATEN', 'forces'],
                meta['height_mm'], area))
        except KeyError as exc:
            raise ValueError('Missing paired RP displacement/reaction at time '+str(time)) from exc
    if not rows:
        raise ValueError('No accepted paired reference-point history; regenerate and rerun')
    return rows, area
