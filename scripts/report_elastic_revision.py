"""Reassess the no-AMS native control and compare its bound static reference."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.probe_courant_elastic_lattice import analyze


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',type=Path,required=True)
    p.add_argument('--static',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    native=analyze(a.native)
    static=json.loads(a.static.read_text())
    mesh=json.loads(Path(static['shared_mesh']).read_text())
    native_meta=json.loads(a.native.read_text())
    native_mesh=json.loads(Path(native_meta['diagnostic']['mesh_file']).read_text())
    if any(mesh[key]!=native_mesh[key] for key in ('nodes','elements','boundary','material')):
        raise ValueError('Static and native core or boundary conditions differ')
    case=next(c for c in static['cases'] if c['element_type']=='C3D8I')
    if not case.get('usable'):
        raise ValueError('Static C3D8I reference did not complete')
    for name,digest in case['output_sha256'].items():
        if sha(a.static.parent/name)!=digest:
            raise ValueError('Static output changed since execution: '+name)
    if sha(static['shared_mesh'])!=static['shared_mesh_sha256']:
        raise ValueError('Static mesh description changed')
    if native['ams']['enabled']:
        raise ValueError('This report requires the no-AMS control')
    old_path=ROOT/'validation/matched_elastic_native_audit.json'
    old=next(r for r in json.loads(old_path.read_text())['cases'] if r['solid_formulation']==24)
    k=native['work_equivalent_stiffness_N_per_mm'];ks=case['stiffness_N_per_mm'];ka=old['work_equivalent_stiffness_N_per_mm']
    report=dict(status='completed matched elastic controls; fracture calibration still blocked',
        calibrated=False,experimental_agreement=False,
        source_sha256={str(a.native):sha(a.native),str(a.static):sha(a.static),str(old_path):sha(old_path)},
        exact_core_mesh_and_boundary_match=True,
        native_noams=native,static_reference=case,
        archived_AMS_reference=dict(source=str(old_path),stiffness_N_per_mm=ka,
                                   scope='archived component diagnostic; not rerun by this report'),
        relative_native_noams_minus_static=k/ks-1,
        relative_native_noams_minus_archived_ams=k/ka-1,
        relative_final_internal_minus_external_work=(native['final_internal_energy_Nmm']/native['final_external_work_Nmm']-1),
        conclusion='Removing AMS does not resolve elastic stiffness excess on this retained mesh. Investigate native formulation/discretization before fracture-locus fitting. No full-collapse or paper agreement qualification follows.')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,4))
    bars=ax.bar(['Static C3D8I','Native HEPH\narchived AMS','Native HEPH\nno AMS'],[ks,ka,k],color=['#6b7280','#d97706','#2563eb'])
    ax.bar_label(bars,fmt='%.2f');ax.set_ylabel('Work-equivalent stiffness (N/mm)')
    ax.set_title('Matched elastic core diagnostic; no material fit');ax.set_ylim(0,1200)
    fig.tight_layout();fig.savefig(a.output.with_suffix('.png'),dpi=180);plt.close(fig)
    print(json.dumps({key:report[key] for key in ('relative_native_noams_minus_static','relative_native_noams_minus_archived_ams','relative_final_internal_minus_external_work')},indent=2))


if __name__=='__main__':
    main()
