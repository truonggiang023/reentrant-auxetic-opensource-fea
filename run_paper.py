"""Generate/run Dong-2019-type CalculiX compression with sliding steel platens.

Run from repository root:
    python run_paper.py --run
    python run_paper.py --target-output-only --run  # exact morphology frames
Linux defaults to ccx on PATH; the tested Linux solver is CalculiX 2.20.
"""
import argparse
import json
import os
import subprocess
import shutil
from reentrant_fea.ccx_results import validate_run
from pathlib import Path
import numpy as np

from reentrant_fea.geometry import core_polygon
from reentrant_fea.mesh import mesh_wedges,mesh_bricks
from reentrant_fea.plastic import from_reference
from reentrant_fea.paper_contact_solver import write_paper_contact_deck
from reentrant_fea.provenance import record_generation


ROOT = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config",default="configs/dong2019_contact.json")
    ap.add_argument("--plastic",default="reference_setup_numeric_source.py")
    ap.add_argument("--output",default="results/dong2019_contact.inp")
    ap.add_argument("--ccx",default=None)
    ap.add_argument("--run",action="store_true")
    ap.add_argument("--target-output-only",action="store_true", help="Request exact strain frames .12/.35/.65; sparse history limits offline initiation diagnostics")
    ap.add_argument("--threads",type=int,default=2)
    ap.add_argument("--self-contact",action="store_true",
                    help="Disjoint wall patches with shared-node adjacency exclusions; potentially expensive")
    ap.add_argument("--imperfection-mm",type=float,default=0.0,
                    help="Explicit asymmetric X perturbation (mm) to explore mode selection; 0 by default")
    args=ap.parse_args()
    if args.threads < 1:
        ap.error('--threads must be positive')
    c=json.loads((ROOT/args.config).read_text(encoding="utf8"))
    c["target_output_only"] = args.target_output_only
    poly=core_polygon(c)
    mesher=mesh_bricks if c.get("element_type") in ("C3D8","C3D8R") else mesh_wedges
    kwargs={"strategy":c.get("brick_mesh","recombine")} if mesher==mesh_bricks else {}
    nodes,els,top,bottom=mesher(poly,c["depth"],c["mesh_size"],c["depth_layers"],**kwargs)
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
    generation=record_generation(
        out, ROOT/args.config, ROOT/args.plastic, config=c, metadata=meta,
        source_paths=[Path(__file__), ROOT/'src/reentrant_fea/geometry.py',
                      ROOT/'src/reentrant_fea/mesh.py',
                      ROOT/'src/reentrant_fea/paper_contact_solver.py',
                      ROOT/'src/reentrant_fea/contact.py', ROOT/'native/ductile.c'])
    # Overrides can change geometry/contact without changing the JSON file.
    generation['invocation'] = {'self_contact':args.self_contact,
                              'imperfection_mm':args.imperfection_mm,
                              'target_output_only':args.target_output_only,
                              'threads':args.threads}
    out.with_suffix('.generation.json').write_text(json.dumps(generation,indent=2)+'\n')
    print(json.dumps(meta,indent=2))
    print("Material:",meta["damage_model"])
    if c.get("ductile_damage"):
        print("WARNING: fracture-locus assumptions are not a recovered paper calibration; failed elements retain residual stiffness.")
    print("IMPORTANT: compare snapshots at nominal strain 0.12, 0.35, 0.65.")
    print("IMPORTANT: compare plateau stress with Dong et al. thin-wall 0.28 MPa.")
    print(f"Mesh element type: {meta['element_type']}; mesh convergence and paper equivalence are not established.")
    if args.run:
        choice = args.ccx or ("ccx-damage" if c.get("ductile_damage") else ("tools/calculix/bin/ccx/ccx213.exe" if os.name == "nt" else "ccx"))
        resolved = shutil.which(choice)
        exe=Path(resolved) if resolved else (ROOT/choice).resolve()
        if not exe.is_file():
            raise FileNotFoundError("CalculiX executable missing: "+str(exe))
        env=os.environ.copy()
        if c.get('ductile_damage'):
            env['REENTRANT_MATERIAL_DATA']=str(out.with_suffix('.material.dat'))
        for key in ("OMP_NUM_THREADS","CCX_NPROC_STIFFNESS",
                    "CCX_NPROC_RESULTS","CCX_NPROC_EQUATION_SOLVER"):
            env[key]=str(args.threads)
        # Quarantine prior outputs to prevent stale files from passing validation.
        import tempfile
        stale = Path(tempfile.mkdtemp(prefix=out.stem + '-previous-', dir=out.parent))
        for suffix in ('.sta', '.dat', '.frd', '.cvg', '.log'):
            prior = out.with_suffix(suffix)
            if prior.exists():
                prior.rename(stale/prior.name)
        log=out.with_suffix('.log')
        with log.open('w') as stream:
            subprocess.run([str(exe),"-i",out.stem], cwd=out.parent,env=env,
                           stdout=stream,stderr=subprocess.STDOUT,check=True)
        if c.get('ductile_damage') and 'REENTRANT_MATERIAL_ABI=3' not in log.read_text():
            raise RuntimeError('Solver did not load the compiled damage material; use the qualified ccx-damage build')
        increments=validate_run(out,log)
        print(f'CalculiX completed {len(increments)} converged increments to step time 1.')
        print("Results:",out.with_suffix(".frd"))


if __name__=="__main__":
    main()
