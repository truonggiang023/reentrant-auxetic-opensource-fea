"""Generate/run Dong-2019-type CalculiX compression with sliding steel platens.

Run from repository root:
    python run_paper.py --run
    python run_paper.py --self-contact --run  # exploratory only
The solver binary is the repo's existing portable ccx213.exe by default.
"""
import argparse
import json
import os
import subprocess
from pathlib import Path
import numpy as np

from reentrant_fea.geometry import core_polygon
from reentrant_fea.mesh import mesh_wedges
from reentrant_fea.plastic import from_reference
from reentrant_fea.paper_contact_solver import write_paper_contact_deck


ROOT = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config",default="configs/dong2019_thin.json")
    ap.add_argument("--plastic",default="reference_setup_numeric_source.py")
    ap.add_argument("--output",default="results/dong2019_contact.inp")
    ap.add_argument("--ccx",default="tools/calculix/bin/ccx/ccx213.exe")
    ap.add_argument("--run",action="store_true")
    ap.add_argument("--threads",type=int,default=16)
    ap.add_argument("--self-contact",action="store_true",
                    help="EXPERIMENTAL: CCX213 single-surface self contact may not work")
    ap.add_argument("--imperfection-mm",type=float,default=0.0,
                    help="Explicit asymmetric X perturbation (mm) to explore mode selection; 0 by default")
    args=ap.parse_args()
    c=json.loads((ROOT/args.config).read_text(encoding="utf8"))
    poly=core_polygon(c)
    nodes,els,top,bottom=mesh_wedges(poly,c["depth"],c["mesh_size"],c["depth_layers"])
    imperfection=None
    if args.imperfection_mm:
        if args.imperfection_mm<0 or args.imperfection_mm>0.1*c["wall"]:
            raise ValueError("Limit proposed geometric imperfection to 0..10% of wall thickness")
        ymin,ymax=nodes[:,1].min(),nodes[:,1].max()
        yi=(nodes[:,1]-ymin)/(ymax-ymin)
        imperfection=np.zeros((len(nodes),2))
        imperfection[:,0]=args.imperfection_mm*np.sin(np.pi*yi)*np.cos(np.pi*yi)
    out=(ROOT/args.output).resolve()
    meta=write_paper_contact_deck(out,nodes,els,c,
                                  from_reference(ROOT/args.plastic),
                                  imperfections=imperfection,
                                  self_contact=args.self_contact)
    print(json.dumps(meta,indent=2))
    print("IMPORTANT: no ductile damage/deletion in this CCX213 deck.")
    print("IMPORTANT: compare snapshots at nominal strain 0.12, 0.35, 0.65.")
    print("IMPORTANT: compare plateau stress with Dong et al. thin-wall 0.28 MPa.")
    print("WARNING: paper used C3D8R; this mesh uses C3D6; mesh-equivalence is NOT established.")
    if args.run:
        exe=(ROOT/args.ccx).resolve()
        if not exe.is_file():
            raise FileNotFoundError("CalculiX executable missing: "+str(exe))
        env=os.environ.copy()
        for key in ("OMP_NUM_THREADS","CCX_NPROC_STIFFNESS",
                    "CCX_NPROC_RESULTS","CCX_NPROC_EQUATION_SOLVER"):
            env[key]=str(args.threads)
        subprocess.run([str(exe),"-i",out.stem],
                       cwd=out.parent,env=env,check=True)
        print("CalculiX exit code 0; inspect .sta/.dat for completed increments.")
        print("Results:",out.with_suffix(".frd"))


if __name__=="__main__":
    main()
