"""Probe implicit availability without mistaking runtime failure for convergence.

The 8-brick elastic cube retains interior Y equilibrium DOFs. Histories and
animations are deliberately omitted: prior output requests crashed at cycle 0
before the useful native capability diagnostic. No lattice fit is performed.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import time

import numpy as np

from reentrant_fea.radioss import write_radioss_deck, _ints, _reals
from reentrant_fea.radioss_results import normal_termination, solver_errors


CASES = [(14, 1), (14, 2), (17, 1), (18, 1), (24, 1)]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def generate(directory, solid, solver):
    nodes = np.array([[x,y,z] for z in (0., .5, 1.) for y in (0., .5, 1.)
                      for x in (0., .5, 1.)])
    idx = lambda i,j,k: k*9+j*3+i
    elements = np.array([[idx(i,j,k),idx(i+1,j,k),idx(i+1,j+1,k),idx(i,j+1,k),
                          idx(i,j,k+1),idx(i+1,j,k+1),idx(i+1,j+1,k+1),idx(i,j+1,k+1)]
                         for k in range(2) for j in range(2) for i in range(2)])
    meta = write_radioss_deck(directory/'cube', nodes, elements,
                              dict(young=69000., poisson=.3, compression_strain=.001),
                              [(225.,0.),(225.,1.)],duration=.001,tensile_coupon=True)
    path = Path(meta['starter'])
    old = path.read_text().splitlines()
    begin, stop = old.index('/MAT/PLAS_TAB/1'), old.index('/NODE')
    old[begin:stop] = ['/IMPLICIT', '/MAT/ELAST/1', 'Elastic capability cube',
                       _reals(2.67e-9), _reals(69000., .3)]
    prop = old.index('/PROP/SOLID/1')
    old[prop+1] = f'Implicit availability probe solid={solid}'
    old[prop+2] = _ints(solid, 0)+old[prop+2][20:]
    lines, i = [], 0
    while i < len(old):
        if old[i].startswith('/TH/'):
            i += 1
            while i < len(old) and not old[i].startswith('/'):
                i += 1
        elif old[i] == '/BCS/1':
            lines += ['/BCS/1', 'Bottom Y only', f'{"010":>6} 000'+_ints(0,1)]
            i += 3
        elif old[i] == '/BCS/2':
            i += 3
        else:
            lines.append(old[i])
            i += 1
    if lines.pop() != '/END':
        raise ValueError('Unexpected writer layout')
    lines += ['/GRNOD/NODE/5', 'All XZ fixed']
    for start in range(1, len(nodes)+1, 10):
        lines.append(_ints(*range(start, min(start+10,len(nodes)+1))))
    lines += ['/BCS/5', 'All XZ', f'{"101":>6} 000'+_ints(0,5), '/END']
    path.write_text('\n'.join(lines)+'\n')
    engine = Path(meta['engine'])
    engine.write_text('\n'.join(['#RADIOSS ENGINE','/RUN/cube/1',_reals(.001),
        '/IMPL/LINEAR',f'/IMPL/SOLVER/{solver}',_ints(2,0,0)+_reals(1e-10),
        '/IMPL/DTINI',_reals(.001),'/IMPL/DT/STOP',_reals(.001,.001),
        '/IMPL/PRINT/LINEAR/1','/PRINT/-1'])+'\n')
    meta['sha256'] = {p.name:digest(p) for p in (path,engine)}
    meta.pop('hardening_extrapolation',None)
    meta['diagnostic'] = dict(kind='implicit runtime availability only', solid=solid,
                              solver=solver, all_xz_fixed=True, interior_y_dofs=9,
                              script_sha256=digest(__file__), material='ELAST')
    (directory/'cube.radioss.json').write_text(json.dumps(meta,indent=2)+'\n')
    return meta


def probe(directory, solid, solver, binaries, timeout):
    directory.mkdir()
    meta = generate(directory,solid,solver)
    metadata_hash = digest(directory/'cube.radioss.json')
    environment = os.environ.copy()
    environment['OMP_NUM_THREADS'] = '1'
    stages, started = [], time.monotonic()
    for binary, arguments in zip(binaries, [['-i','cube_0000.rad','-np','1'],['-i','cube_0001.rad']]):
        log = directory/(Path(binary).name+'.log')
        with log.open('w') as stream:
            try:
                result = subprocess.run([binary,*arguments],cwd=directory,env=environment,
                                        stdout=stream,stderr=subprocess.STDOUT,timeout=timeout)
                stages.append(dict(executable=str(binary),exit_code=result.returncode))
            except subprocess.TimeoutExpired:
                stages.append(dict(executable=str(binary),timeout=True))
                break
        if stages[-1]['exit_code']:
            break
    outputs = [directory/name for name in ('cube_0000.out','cube_0001.out',
               'starter_linux64_gf.log','engine_linux64_gf.log') if (directory/name).is_file()]
    engine_text = (directory/'cube_0001.out').read_text(errors='replace') if (directory/'cube_0001.out').exists() else ''
    combined = '\n'.join(p.read_text(errors='replace') for p in outputs)
    missing_mumps = 'Fatal error: MUMPS required' in combined
    normal = normal_termination(engine_text)
    return dict(solid=solid,solver=solver,stages=stages,elapsed_wall_seconds=time.monotonic()-started,
                source_inputs_sha256=meta['sha256'],metadata_sha256_before_execution=metadata_hash,
                metadata_unchanged=digest(directory/'cube.radioss.json')==metadata_hash,
                output_sha256={p.name:digest(p) for p in outputs},
                engine_normal_termination=normal,engine_error_lines=solver_errors(engine_text),
                explicit_mumps_required_failure=missing_mumps,
                usable_as_static_stiffness_evidence=False,
                reason=('Runtime explicitly reports MUMPS required before equilibrium solve' if missing_mumps else
                        'Availability or failure alone is not an independently checked numerical stiffness result'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--report',type=Path,required=True)
    p.add_argument('--timeout-seconds',type=float,default=30.)
    args = p.parse_args()
    if not 0 < args.timeout_seconds < float('inf'):
        p.error('Require finite positive timeout')
    if args.output.exists() and any(args.output.iterdir()):
        p.error('Use a fresh directory to preserve existing evidence')
    binaries = [shutil.which(n) for n in ('starter_linux64_gf','engine_linux64_gf')]
    if not all(binaries):
        p.error('Source the official OpenCourant env.sh first')
    args.output.mkdir(parents=True,exist_ok=True)
    rows = [probe(args.output/f'solid{s}_solver{v}',s,v,binaries,args.timeout_seconds) for s,v in CASES]
    report = dict(status='Native implicit capability probe; no static solution or calibration established',
        script='scripts/probe_courant_implicit.py',script_sha256=digest(__file__),
        runtime_executables_sha256={path:digest(path) for path in binaries},
        geometry=dict(cube_side_mm=1.,bricks=8,nodes=27,free_interior_y_dofs=9),
        material=dict(young_MPa=69000.,poisson=.3,density_Mg_per_mm3=2.67e-9),
        loading=dict(top_y_mm=.001,bottom_y_mm=0.,all_xz_fixed=True),
        analytical_check_if_a_future_run_succeeds=dict(constrained_y_modulus_MPa=69000.*.7/(1.3*.4),
                                                     expected_top_force_N=69000.*.7/(1.3*.4)*.001),
        cases=rows, all_cases_explicitly_require_mumps=all(r['explicit_mumps_required_failure'] for r in rows),
        implicit_static_qualified=False,
        qualification='This official GNU executable cannot supply the proposed implicit reference in these tested cases. A separate MUMPS-enabled build would need its own runtime and numerical qualification.')
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_cases_explicitly_require_mumps':report['all_cases_explicitly_require_mumps'],
                      'static_stiffness_qualified':False},indent=2))


if __name__=='__main__':
    main()
