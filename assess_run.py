"""Report accepted solver evidence and contact geometry, without mode claims."""
import argparse
import hashlib
import json
from pathlib import Path
from reentrant_fea.ccx_results import completed_increments
from reentrant_fea.rp_curve import ccx_rp_history
from reentrant_fea.frd import committed_fields,read_core_deck,projected_overlap,contact_envelope


def assess(deck, targets=(.12,.35,.65), tolerance=.002):
    deck=Path(deck);nodes,elements,meta=read_core_deck(deck)
    inc=completed_increments(deck.with_suffix('.sta'))
    log=deck.with_suffix('.log').read_text(errors='replace')
    errors=[s.strip() for s in log.splitlines() if '*ERROR' in s]
    complete=bool(inc and abs(inc[-1][2]-meta['period'])<1e-6 and 'Job finished' in log and not errors)
    if meta['native_material'] and 'REENTRANT_MATERIAL_ABI=3' not in log:
        raise ValueError('Native-material run has no qualified interface signature')
    best={target:(float('inf'),None,None,None) for target in targets}
    latest=None
    for time,name,data in committed_fields(deck.with_suffix('.frd')):
        if name!='DISP':continue
        strain=time/meta['period']*meta['end_nominal_strain']
        if latest is None or time>latest[0]:latest=(time,strain,data)
        for target in targets:
            err=abs(strain-target)
            if err<best[target][0]:best[target]=err,strain,time,data
    frames=[]
    for target,(err,strain,time,data) in best.items():
        row={'target_strain':target,'available':data is not None and err<=tolerance}
        if row['available']:
            row.update(actual_strain=strain,time=time,
                       projected_contact=projected_overlap(nodes,elements,data),
                       experimental_mode_agreement='not established')
        frames.append(row)
    forces, reference_area = ccx_rp_history(deck)
    max_damage=None;max_peeq=None;mode=None
    accepted=[row[2] for row in inc]
    def committed(time):return any(abs(time-t)<=max(1e-10,abs(t)*6e-6) for t in accepted)
    lines=deck.with_suffix('.dat').read_text().splitlines()
    for i,line in enumerate(lines):
        if line.strip().startswith('internal state variables'):
            mode='sdv' if committed(float(line.split()[-1])) else None;continue
        if mode=='sdv':
            fields=line.split()
            if len(fields)==12 and fields[0].isdigit():
                max_damage=max(float(fields[11]),max_damage or 0.)
                max_peeq=max(float(fields[8]),max_peeq or 0.)
            elif fields and not fields[0].isdigit():mode=None
    try:contact=contact_envelope(deck.with_suffix('.frd'))
    except ValueError:contact=None
    result={'input_sha256':hashlib.sha256(deck.read_bytes()).hexdigest(),
            'mesh':meta,'completed':complete,'errors':errors,'converged_increments':len(inc),
            'last_converged_nominal_strain':(inc[-1][2]/meta['period']*meta['end_nominal_strain'] if inc else 0.),
            'frames':frames,'contact':contact,'reaction_history':forces,
            'max_recorded_integration_point_damage':max_damage,
            'max_recorded_integration_point_peeq':max_peeq,
            'stress_area_convention_mm2':reference_area,
            'curve_source':'coupled platen reference-point RF2 and relative U2',
            'qualification':'solver evidence only; calibration, full-contact validity and experimental agreement are not established'}
    if latest is not None:
        time,strain,data=latest
        result['latest_recorded_frame']={'time':time,'nominal_strain':strain,
                                         'projected_contact':projected_overlap(nodes,elements,data)}
    sidecar=deck.with_suffix('.material.dat')
    if sidecar.exists():result['material_sidecar_sha256']=hashlib.sha256(sidecar.read_bytes()).hexdigest()
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('deck',type=Path);p.add_argument('--output',type=Path)
    args=p.parse_args();result=assess(args.deck)
    text=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(text)
    else:print(text,end='')

if __name__=='__main__':main()
