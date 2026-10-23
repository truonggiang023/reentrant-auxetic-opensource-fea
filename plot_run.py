"""Plot committed static solver frames; never substitute an unconverged result."""
import argparse
from pathlib import Path
import numpy as np
from reentrant_fea.frd import committed_fields,read_core_deck


def main():
    p=argparse.ArgumentParser()
    p.add_argument('deck',type=Path)
    p.add_argument('--output',type=Path,default=Path('results/target_frames.png'))
    p.add_argument('--tolerance',type=float,default=.002)
    args=p.parse_args()
    xyz,elements,meta=read_core_deck(args.deck)
    nodes={i+1:p for i,p in enumerate(xyz)}
    elements=(elements+1).tolist()
    core=set(v for el in elements for v in el)
    end_strain=meta['end_nominal_strain']
    targets=[.12,.35,.65];best={x:(float('inf'),None,None) for x in targets}
    for time,name,data in committed_fields(args.deck.with_suffix('.frd')):
        if name!='DISP':continue
        strain=time/meta['period']*end_strain
        for target in targets:
            error=abs(strain-target)
            if error<best[target][0]:best[target]=error,strain,data
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    zmin=min(nodes[v][2] for v in core)
    faces=[el[:len(el)//2] for el in elements if all(abs(nodes[v][2]-zmin)<1e-7 for v in el[:len(el)//2])]
    fig,axs=plt.subplots(1,3,figsize=(15,5))
    for ax,target in zip(axs,targets):
        error,actual,data=best[target]
        if data is None or error>args.tolerance:
            ax.text(.5,.5,f'No converged frame near strain {target}',ha='center',va='center',transform=ax.transAxes)
            ax.set_axis_off();continue
        polygons=[np.array([nodes[v][:2]+np.array(data[v][:2]) for v in face]) for face in faces]
        ax.add_collection(PolyCollection(polygons,facecolors='#b9cbdc',edgecolors='#314b63',linewidths=.2))
        ax.autoscale();ax.set_aspect('equal');ax.set_xlabel('x (mm)');ax.set_ylabel('y (mm)')
        ax.set_title(f'Target {target}; solved {actual:.5f}')
    fig.suptitle('CalculiX solver frames; experimental mode agreement is not established')
    fig.tight_layout();args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=160)
    print(args.output)

if __name__=='__main__':main()
