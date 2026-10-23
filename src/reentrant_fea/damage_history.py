"""Postprocess solver-exported PEEQ/triaxiality histories (NOT coupled deletion).

CSV must contain step, element_id, peeq, triaxiality, characteristic_length.
Damage is integrated across successive records of each element.
"""
import csv
from collections import defaultdict
from pathlib import Path
from .material import DuctilePoint

def evaluate_damage_csv(input_csv, output_csv, fracture_curve, failure_displacement):
    rows=list(csv.DictReader(Path(input_csv).open(newline='',encoding='utf8')))
    rows.sort(key=lambda r:(int(r['element_id']),float(r['step'])))
    points=defaultdict(DuctilePoint)
    prev={}
    result=[]
    for row in rows:
        eid=int(row['element_id']);peeq=float(row['peeq'])
        if peeq < prev.get(eid,0)-1e-10:
            raise ValueError(f'Nonmonotonic PEEQ for element {eid}')
        inc=max(0,peeq-prev.get(eid,0));prev[eid]=peeq
        p=points[eid].update(inc,float(row['triaxiality']),float(row['characteristic_length']),fracture_curve,failure_displacement)
        result.append({'step':row['step'],'element_id':eid,'peeq':peeq,
                       'initiation_index':p.initiation_index,'damage':p.damage,'failed':int(p.failed)})
    with Path(output_csv).open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=['step','element_id','peeq','initiation_index','damage','failed'])
        w.writeheader();w.writerows(result)
    return result
