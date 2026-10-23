import argparse
import json
from pathlib import Path
from .geometry import core_polygon
from .mesh import mesh_wedges
from .solver import write_deck,run_ccx
from .plastic import from_reference

def main():
    p=argparse.ArgumentParser(description='Independent reentrant honeycomb FE pipeline')
    p.add_argument('--config',default='configs/dong2019_thin.json')
    p.add_argument('--plastic',required=True,help='Original Python setup containing PLASTIC_DATA; only numbers are read')
    p.add_argument('--output',default='results/reentrant.inp')
    p.add_argument('--run',action='store_true',help='Run installed CalculiX ccx')
    p.add_argument('--ccx',default='ccx')
    a=p.parse_args()
    c=json.loads(Path(a.config).read_text())
    geom=core_polygon(c)
    nodes,elements,top,bottom=mesh_wedges(geom,c['depth'],c['mesh_size'],c['depth_layers'])
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    deck=write_deck(out,nodes,elements,top,bottom,c,from_reference(a.plastic))
    print('Geometry area %.4f mm2 | nodes %d | wedges %d | top %d | bottom %d'%(geom.area,len(nodes),len(elements),len(top),len(bottom)))
    print('Wrote',deck)
    if a.run:run_ccx(deck,a.ccx)
if __name__=='__main__':main()
