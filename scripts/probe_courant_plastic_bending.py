"""Compare native solid formulations in an original plastic wall cantilever.

This is a component diagnostic with a constrained tip face, not a calibration
test or a rigorous plastic-limit benchmark. Source the qualified native env.sh
before running. Raw native evidence stays in the requested output directory.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np

from reentrant_fea.radioss import write_radioss_deck
from reentrant_fea.radioss_results import normal_termination, read_history_csv, solver_errors


E, NU, T, L, B = 69000., .3, .3, 5.66, .3
NY, DISPLACEMENT, YIELD = 16, 1., 225.
CASES = [('heph24_two_wall', 24, 2, .005), ('solid18_two_wall', 18, 2, .005),
         ('heph24_four_wall', 24, 4, .002), ('solid18_four_wall', 18, 4, .002)]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def generate(directory, solid, nx, duration):
    nodes = np.array([[x,y,z] for z in (0., B) for y in np.linspace(0,L,NY+1)
                      for x in np.linspace(0,T,nx+1)])
    index = lambda i,j,k: k*(nx+1)*(NY+1)+j*(nx+1)+i
    elements = np.array([[index(i,j,0),index(i+1,j,0),index(i+1,j+1,0),index(i,j+1,0),
                          index(i,j,1),index(i+1,j,1),index(i+1,j+1,1),index(i,j+1,1)]
                         for j in range(NY) for i in range(nx)])
    metadata = write_radioss_deck(directory/'wall', nodes, elements,
                                  {'young':E, 'poisson':NU, 'compression_strain':DISPLACEMENT/L},
                                  [(YIELD,0.),(YIELD,1.)], duration=duration, tensile_coupon=True)
    starter = Path(metadata['starter'])
    original = starter.read_text().splitlines()
    prop = original.index('/PROP/SOLID/1')
    # ISMSTR=0 lets each native formulation select its supported finite-geometry default.
    original[prop+2] = f'{solid:10d}{0:10d}'+original[prop+2][20:]
    original[prop+1] = f'Plastic bending diagnostic Isolid={solid}'
    lines, pos = [], 0
    while pos < len(original):
        if original[pos] == '/BCS/1':
            lines += ['/BCS/1','Bottom XY clamp', f'{"110":>6} 000'+f'{0:10d}{1:10d}']
            pos += 3
        elif original[pos] == '/BCS/2':
            pos += 3
        elif original[pos] == '/IMPDISP/1':
            lines += ['/IMPDISP/1','Tip face same prescribed X; Y free',
                      original[pos+2].replace('Y','X'),original[pos+3]]
            pos += 4
        else:
            lines.append(original[pos])
            pos += 1
    assert lines.pop() == '/END'
    lines += ['/GRNOD/NODE/5','All plane-strain Z nodes']
    for first in range(1,len(nodes)+1,10):
        lines.append(''.join(f'{n:10d}' for n in range(first,min(first+10,len(nodes)+1))))
    lines += ['/BCS/5','All Z fixed',f'{"001":>6} 000'+f'{0:10d}{5:10d}',
              '/TH/NODE/7','Bottom face nodal histories',
              ''.join(f'{v:<10s}' for v in ('DX','DY','DZ','REACX','REACY','REACZ'))]
    bottom = np.flatnonzero(nodes[:,1] == 0.)+1
    for n in bottom:
        lines.append(f'{n:10d}{0:10d}base-node-{n}')
    lines += ['/END']
    starter.write_text('\n'.join(lines)+'\n')
    # The writer's original hash preceded the diagnostic-only deck edits above.
    metadata['sha256'] = {path.name:digest(path) for path in (starter,Path(metadata['engine']))}
    metadata['diagnostic'] = {'formulation':solid,'elements_through_wall':nx,
                              'duration_seconds':duration,'bottom_face_nodes':bottom.tolist(),
                              'generated_deck_sha256':digest(starter)}
    (directory/'wall.radioss.json').write_text(json.dumps(metadata,indent=2)+'\n')
    return metadata, nodes


def execute(directory, executable, arguments, timeout):
    environment = os.environ.copy()
    environment['OMP_NUM_THREADS'] = '1'
    with (directory/(Path(executable).name+'.log')).open('w') as output:
        result = subprocess.run([executable,*arguments],cwd=directory,env=environment,
                                stdout=output,stderr=subprocess.STDOUT,timeout=timeout)
    if result.returncode:
        raise RuntimeError('Native executable failed: '+executable)


def analyze(directory, label, solid, nx, duration, nodes):
    starter = (directory/'wall_0000.out').read_text()
    engine = (directory/'wall_0001.out').read_text()
    history = read_history_csv(directory/'wallT01.csv')
    headers = list(history)
    data = np.column_stack(list(history.values()))
    errors = solver_errors(engine)
    valid = (bool(re.search(r'^\s*0\s+ERROR\(S\)\s*$',starter,re.M))
             and bool(re.search(r'^\s*0\s+WARNING\(S\)\s*$',starter,re.M))
             and not solver_errors(starter) and not errors
             and normal_termination(engine) and np.isfinite(data).all()
             and not re.search(r'\b(?:nan|infinity|inf)\b',engine,re.I)
             and 'WARNING ZERO OR NEGATIVE VOLUME' not in engine
             and data[-1,0] > .995*duration)
    if not valid:
        raise ValueError('Unqualified native result: '+label)
    topcols = [i for i,h in enumerate(headers) if h.startswith('Loading face nodal histories')]
    basecols = [i for i,h in enumerate(headers) if h.startswith('Bottom face nodal histories')]
    assert len(topcols) == len(basecols) == 6*2*(nx+1)
    dt = np.diff(data[:,0])
    avg = lambda a: .5*(a[:-1]+a[1:])
    topforce = np.diff(data[:,topcols[3::6]],axis=0)/dt[:,None]
    basefx = np.diff(data[:,basecols[3::6]],axis=0)/dt[:,None]
    basefy = np.diff(data[:,basecols[4::6]],axis=0)/dt[:,None]
    topfy = np.diff(data[:,topcols[4::6]],axis=0)/dt[:,None]
    topuy = avg(data[:,topcols[1::6]])
    displacement = avg(data[:,topcols[::6]]).mean(1)
    force = topforce.sum(1)
    # Reactions act at their individual CURRENT positions; the tip forces can form a couple.
    exact_moment = (topforce*(L+topuy)).sum(1)
    tip_couple = (topforce*(topuy-topuy.mean(1)[:,None])).sum(1)
    base_x = nodes[nodes[:,1] == 0.,0]
    base_moment = (basefy*(base_x[None,:]-T/2)).sum(1)
    momentum_rate = np.diff(data[:,headers.index('X-MOMENTUM')])/dt
    balance = force+basefx.sum(1)-momentum_rate
    equivalent_plane_strain_yield = 2/np.sqrt(3)*YIELD
    reference = equivalent_plane_strain_yield*B*T*T/4
    samples = []
    for target in (.05,.1,.2,.25,.35,.5,.75,.95):
        i = int(np.argmin(abs(displacement-target)))
        samples.append(dict(target_tip_displacement_mm=target,tip_displacement_mm=float(displacement[i]),
                            force_N=float(force[i]),exact_top_load_moment_about_base_centroid_Nmm=float(exact_moment[i]),
                            tip_constraint_couple_Nmm=float(tip_couple[i]),
                            base_reaction_moment_Nmm=float(base_moment[i]),
                            support_axial_resultant_N=float(basefy[i].sum()),
                            applied_tip_axial_resultant_N=float(topfy[i].sum()),
                            exact_top_moment_over_ideal_beam_Mp=float(exact_moment[i]/reference),
                            force_balance_residual_N=float(balance[i])))
    files = ['wall_0000.rad','wall_0001.rad','wall_0000.out','wall_0001.out','wallT01','wallT01.csv']
    banner = (directory/'starter_linux64_gf.log').read_text()
    commit = re.search(r'CommitID:\s*([0-9a-f]{40})',banner)
    return dict(case=label,Isolid=solid,elements_through_wall=nx,elements_along_wall=NY,
                element_count=nx*NY,depth_layers=1,duration_seconds=duration,
                starter_errors=0,starter_warnings=0,engine_error_lines=0,finite_history=True,
                normal_termination=True,history_samples=len(data),last_history_time_seconds=float(data[-1,0]),
                runtime_binary_commit=commit.group(1) if commit else None,
                final_internal_energy_Nmm=float(data[-1,1]),final_kinetic_energy_Nmm=float(data[-1,2]),
                final_external_work_Nmm=float(data[-1,9]),final_plastic_work_Nmm=float(data[-1,16]),
                maximum_absolute_force_balance_residual_N=float(abs(balance).max()),
                reconstructed_displacement_work_Nmm=float(np.sum(force*np.diff(data[:,topcols[::6]].mean(1)))),
                samples=samples,input_and_output_sha256={name:digest(directory/name) for name in files})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--timeout-seconds',type=float,default=240.)
    parser.add_argument('--analyze-existing',action='store_true',
                        help='Reanalyze retained inputs and outputs without claiming execution-time metadata binding')
    args = parser.parse_args()
    if not 0 < args.timeout_seconds < float('inf'):
        parser.error('Require finite positive timeout')
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()) and not args.analyze_existing:
        parser.error('Use an empty output directory to preserve prior evidence')
    binaries = [shutil.which(n) for n in ('starter_linux64_gf','engine_linux64_gf','th_to_csv_linux64_gf')]
    if not all(binaries) and not args.analyze_existing:
        parser.error('Source the qualified OpenCourant env.sh first')
    output.mkdir(parents=True,exist_ok=True)
    results = []
    for label,solid,nx,duration in CASES:
        directory = output/label
        if args.analyze_existing:
            nodes = np.array([[x,y,z] for z in (0.,B) for y in np.linspace(0,L,NY+1)
                              for x in np.linspace(0,T,nx+1)])
        else:
            directory.mkdir()
            metadata,nodes = generate(directory,solid,nx,duration)
            for binary,arguments in zip(binaries,[['-i','wall_0000.rad','-np','1'],['-i','wall_0001.rad'],['wallT01']]):
                execute(directory,binary,arguments,args.timeout_seconds)
        row = analyze(directory,label,solid,nx,duration,nodes)
        results.append(row)
        print(label,json.dumps(row['samples'][-1]),flush=True)
    report = dict(status='Plastic bending component diagnostic; no calibration, solver-defect or lattice-cause conclusion',
                  reproduction=dict(command='source /workspace/opencourant/env.sh && '+
                                    ' '.join([sys.executable,'scripts/probe_courant_plastic_bending.py','--output',str(output),'--report',str(args.report)]),
                                    retained_output_analysis_flag='--analyze-existing' if args.analyze_existing else None,
                                    script='scripts/probe_courant_plastic_bending.py',script_sha256=digest(__file__)),
                  evidence_binding='Hashes identify retained actual input/output files; this custom diagnostic does not establish execution-time metadata binding through the production runner.',
                  analysis_of_existing_outputs=args.analyze_existing,
                  geometry_mm=dict(thickness=T,length=L,depth=B),
                  material=dict(young_MPa=E,poisson=NU,density_Mg_per_mm3=2.67e-9,
                                cauchy_flow_stress_MPa=YIELD,hardening='perfectly plastic',damage=False),
                  boundary_conditions='Bottom XY fixed, all Z fixed; all tip-face nodes receive identical X displacement, Y free.',
                  loading=dict(tip_displacement_mm=DISPLACEMENT,amplitude='101-point piecewise-linear quintic smoothstep',
                               history_samples_requested=400,OMP_threads=1,mass_scaling=False,native_default_time_step_scale=.9),
                  force_definition='Native REACX is cumulative impulse: average force = Delta(REACX)/Delta(time).',
                  moment_definition='Sum over individual top nodes: Fx*(L+Uy); same-X tip restraint contributes a finite-rotation couple.',
                  beam_reference=dict(formula='(2/sqrt(3))*yield_stress*depth*thickness^2/4',
                                      ideal_plane_strain_fully_plastic_Mp_Nmm=2/np.sqrt(3)*YIELD*B*T*T/4),
                  cases=results,
                  limitations=['The ideal slender-beam plastic moment is a reference, not a rigorous bound for this finite-rotation constrained-tip solid model.',
                               'Identical tip X constraints suppress projected width change and create opposing nodal loads and a nonzero tip couple.',
                               'Four-through-wall cases also use a shorter loading duration; their comparison is not a pure mesh-convergence study.',
                               'Two and four wall elements do not establish mesh convergence or resolve length/depth, rate, integration and boundary sensitivities.',
                               'Formulation differences in this diagnostic do not establish the cause of the experimental lattice stress mismatch.',
                               'Zero reported hourglass energy for physical stabilization does not mean no stabilization contribution.'])
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(report,indent=2)+'\n')
    print('Wrote '+str(args.report),flush=True)


if __name__ == '__main__':
    main()
