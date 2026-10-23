import argparse
import json
from pathlib import Path
from .geometry import core_polygon
from .mesh import mesh_wedges
from .solver import write_deck,run_ccx
from .advanced_solver import write_advanced_deck
from .paper_modes import write_paper_deck
from .plastic import from_reference

def main():
    p=argparse.ArgumentParser(description='Independent reentrant honeycomb FE pipeline')
    p.add_argument('--config',default='configs/dong2019_thin.json')
    p.add_argument('--plastic',required=True,help='Original Python setup containing PLASTIC_DATA; only numbers are read')
    p.add_argument('--output',default='results/reentrant.inp')
    p.add_argument('--run',action='store_true',help='Run installed CalculiX ccx')
    p.add_argument('--ccx',default='ccx')
    p.add_argument('--paper-modes',action='store_true',help='Dong2019 candidate: steel platens, friction contact, snapshots at 0.12/0.35/0.65')
    p.add_argument('--advanced',action='store_true',help='Experimental rigid bonded-grip model')
    p.add_argument('--self-contact',action='store_true',help='Experimental unvalidated self-contact, requires config acknowledgement')
    p.add_argument('--ductile-damage',action='store_true',help='Fail fast: coupled damage not yet implemented')
    a=p.parse_args()
    c=json.loads(Path(a.config).read_text())
    geom=core_polygon(c)
    nodes,elements,top,bottom=mesh_wedges(geom,c['depth'],c['mesh_size'],c['depth_layers'])
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    if a.paper_modes:
        deck=write_paper_deck(out,nodes,elements,top,bottom,c,from_reference(a.plastic),self_contact=a.self_contact,damage=a.ductile_damage)
    elif a.advanced or a.self_contact or a.ductile_damage:
        deck=write_advanced_deck(out,nodes,elements,top,bottom,c,from_reference(a.plastic),self_contact=a.self_contact,ductile_damage=a.ductile_damage)
    else:
        deck=write_deck(out,nodes,elements,top,bottom,c,from_reference(a.plastic))
    print('Geometry area %.4f mm2 | nodes %d | wedges %d | top %d | bottom %d'%(geom.area,len(nodes),len(elements),len(top),len(bottom)))
    print('Wrote',deck)
    if a.run:run_ccx(deck,a.ccx)
if __name__=='__main__':main()
