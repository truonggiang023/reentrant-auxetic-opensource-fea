"""Verify TYPE25 contact-impulse signs against prescribed-brick momentum.

The two solids are deliberately kinematically prescribed and penetrate. This
is a force-bookkeeping diagnostic; its large energy error excludes mechanical,
quasi-static, material or experimental qualification. No lattice input changes.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np

from reentrant_fea.radioss import _ints, _reals, _group, _bcs, _contact
from reentrant_fea.radioss_results import normal_termination, solver_errors, read_history_csv, read_engine_out


DURATION = .001
CASES = [('plus30_face', 30., True, .35), ('minus30_face', -30., True, .35),
         ('plus30_all', 30., False, .35), ('plus30_face_no_friction', 30., True, 0.)]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def box(x0, x1, y0, y1, z0, z1):
    return np.array([[x0,y0,z0], [x1,y0,z0], [x1,y1,z0], [x0,y1,z0],
                     [x0,y0,z1], [x1,y0,z1], [x1,y1,z1], [x0,y1,z1]], float)


def generate(folder, theta, face_only, friction):
    xyz = np.vstack([box(-2,2,0,.2,-.25,.25), box(-.5,.5,.21,.41,-.25,.25)])
    angle = np.deg2rad(theta)
    rotation = np.array([[np.cos(angle),-np.sin(angle),0],
                         [np.sin(angle),np.cos(angle),0], [0,0,1]])
    xyz = xyz@rotation.T
    motion = rotation@np.array([.2,-.03,0.])
    lines = ['#RADIOSS STARTER','/BEGIN','probe',_ints(2022,0),
             f'{"Mg":>20}{"mm":>20}{"s":>20}',f'{"Mg":>20}{"mm":>20}{"s":>20}',
             '/TITLE','Original inclined sliding contact force audit',
             '/MAT/ELAST/1','Elastic prescribed solids',_reals(2.67e-9),_reals(69000.,.3),'/NODE']
    lines += [_ints(i)+_reals(*p) for i,p in enumerate(xyz,1)]
    lines += ['/PROP/SOLID/1','HEPH',_ints(24,0,None,0,None,0,0,0)+_reals(0),
              _reals(0,0,0),_reals(0)+_ints(1),'/PART/1','Fixed target',_ints(1,1),
              '/BRICK/1',_ints(1,*range(1,9)),'/PART/2','Prescribed mover',_ints(1,1),
              '/BRICK/2',_ints(2,*range(9,17))]
    _group(lines,1,'Fixed target nodes',range(1,9))
    _bcs(lines,1,'Target fixed',1,'111')
    _group(lines,2,'Moving brick nodes',range(9,17))
    _bcs(lines,2,'Mover Z fixed',2,'001')
    lines += ['/FUNCT/20','Smooth translation fraction']
    for f in np.linspace(0,1,101):
        lines += [_reals(DURATION*f,f**3*(10-15*f+6*f*f))]
    for i,axis,amplitude in [(1,'X',motion[0]),(2,'Y',motion[1])]:
        lines += [f'/IMPDISP/{i}','Moving brick prescribed translation',
                  _ints(20)+f'{axis:>10}'+_ints(0,0,2,0,0),_reals(1,amplitude,0,0)]
    if face_only:
        lines += ['/SURF/SEG/1','Target top face',_ints(1,4,8,7,3),
                  '/SURF/SEG/2','Mover bottom face',_ints(2,9,10,14,13)]
    else:
        lines += ['/SURF/PART/ALL/1','Target all faces',_ints(1),
                  '/SURF/PART/ALL/2','Mover all faces',_ints(2)]
    _contact(lines,1,2,1,friction)
    index = lines.index('/INTER/TYPE25/1')
    lines[index+3] = _ints(0,None)+_reals(0,0,.002,.002)
    lines += ['/TH/INTER/1','Interface','DEF',_ints(1),'/TH/NODE/2','Target reaction',
              ''.join(f'{v:<10}' for v in ('REACX','REACY','REACZ'))]
    lines += [_ints(n,0)+f'target-{n}' for n in range(1,9)]
    lines += ['/TH/NODE/3','Mover histories',
              ''.join(f'{v:<10}' for v in ('VX','VY','VZ','REACX','REACY','REACZ'))]
    lines += [_ints(n,0)+f'mover-{n}' for n in range(9,17)]
    lines += ['/TH/PART/4','Parts','DEF',_ints(1,2),'/END']
    starter,engine = folder/'probe_0000.rad',folder/'probe_0001.rad'
    starter.write_text('\n'.join(lines)+'\n')
    engine.write_text('\n'.join(['#RADIOSS ENGINE','/RUN/probe/1',_reals(DURATION),
                                '/PRINT/-10000','/TFILE/4',_reals(DURATION/500),
                                '/ANIM/DT',_reals(0,DURATION)])+'\n')
    record = dict(theta_deg=theta,face_only=face_only,friction=friction,
                  input_sha256={p.name:digest(p) for p in (starter,engine)},
                  generation_script_sha256=digest(__file__))
    (folder/'generation.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


def analyze(folder, theta, face_only, friction, csv_name='probeT01.csv'):
    folder = Path(folder)
    for name in ('probe_0000.out','probe_0001.out'):
        text = (folder/name).read_text()
        if not normal_termination(text) or solver_errors(text):
            raise ValueError('Native failure: '+str(folder/name))
    starter = (folder/'probe_0000.out').read_text()
    if not re.search(r'^\s*0\s+WARNING\(S\)\s*$',starter,re.M):
        raise ValueError('Inspect starter warnings: '+str(folder))
    h = read_history_csv(folder/csv_name)
    contact = np.column_stack([h[k] for k in h if k.startswith('Interface')])
    if contact.shape[1] != 6:
        raise ValueError('Require six verified native DEF interface components')
    target = np.column_stack([h[k] for k in h if k.startswith('Target reaction')]).reshape(-1,8,3).sum(1)
    moving = np.column_stack([h[k] for k in h if k.startswith('Mover histories')]).reshape(-1,8,6)
    mass = 2.67e-9*.1  # Prescribed mover volume: 1 x .2 x .5 mm3.
    velocity = moving[:,:,:3].mean(1)
    external_impulse = moving[:,:,3:].sum(1)
    secondary_contact = mass*velocity-external_impulse
    selected = h['time'] > .55*DURATION
    scale = np.sqrt(np.mean(secondary_contact[selected]**2))
    if scale <= 0 or h['time'][-1] < .995*DURATION:
        raise ValueError('Missing nontrivial completed force diagnostic')
    normal,tangent = contact[:,:3],contact[:,3:]
    alternatives = {}
    for sign in (1,-1):
        for tangential_sign in (1,-1):
            predicted = sign*(normal+tangential_sign*tangent)
            key = f'{sign}*(FN{tangential_sign:+d}FT)'
            alternatives[key] = float(np.sqrt(np.mean((predicted[selected]-secondary_contact[selected])**2))/scale)
    engine = read_engine_out(folder/'probe_0001.out')
    files = ['probe_0000.rad','probe_0001.rad','probe_0000.out','probe_0001.out','probeT01',csv_name]
    generation = folder/'generation.json'
    bound = generation.is_file()
    if bound:
        for name,expected in json.loads(generation.read_text())['input_sha256'].items():
            if digest(folder/name) != expected:
                raise ValueError('Input changed after generation')
    result = dict(theta_deg=theta,face_only=face_only,friction=friction,
                  tangent_displacement_mm=.2,sample_count=len(h['time']),
                  native_stages_normal=True,starter_warnings=0,
                  candidate_relative_RMS_to_secondary_contact_impulse=alternatives,
                  target_vs_secondary_contact_relative_RMS=float(np.sqrt(np.mean((target[selected]-secondary_contact[selected])**2))/scale),
                  final_normal_impulse_Ns=normal[-1].tolist(),
                  final_tangential_channel_impulse_Ns=tangent[-1].tolist(),
                  final_target_reaction_impulse_Ns=target[-1].tolist(),
                  final_secondary_contact_impulse_Ns=secondary_contact[-1].tolist(),
                  max_sampled_global_energy_error_percent=max(abs(r['energy_error_percent']) for r in engine['samples']),
                  generation_hash_recorded=bound,
                  analysis_script_sha256=digest(__file__),
                  retained_input_output_sha256={name:digest(folder/name) for name in files},
                  calibrated=False,quasi_static_validated=False,
                  qualification='Force bookkeeping only on prescribed penetrating solids; large global energy error excludes mechanical and experimental validation.')
    if alternatives['1*(FN-1FT)'] > 1e-4 or result['target_vs_secondary_contact_relative_RMS'] > 1e-4:
        raise ValueError('Independent impulse balances do not verify the proposed convention')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--analyze-existing',action='store_true')
    parser.add_argument('--csv-name',default='probeT01.csv')
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()) and not args.analyze_existing:
        parser.error('Use a fresh output directory to preserve previous evidence')
    output.mkdir(parents=True,exist_ok=True)
    bins = [shutil.which(n) for n in ('starter_linux64_gf','engine_linux64_gf','th_to_csv_linux64_gf')]
    if not args.analyze_existing and not all(bins):
        parser.error('Source the qualified native env.sh first')
    results = []
    for name,theta,face,mu in CASES:
        folder = output/name
        if not args.analyze_existing:
            folder.mkdir()
            generate(folder,theta,face,mu)
            env = os.environ.copy()
            env['OMP_NUM_THREADS'] = '1'
            commands = [['-i','probe_0000.rad','-np','1'],['-i','probe_0001.rad'],['probeT01']]
            for binary,options in zip(bins,commands):
                with (folder/(Path(binary).name+'.log')).open('w') as stream:
                    subprocess.run([binary,*options],cwd=folder,env=env,stdout=stream,
                                   stderr=subprocess.STDOUT,check=True,timeout=120)
        row = analyze(folder,theta,face,mu,args.csv_name)
        results.append(row)
        print(name,row['candidate_relative_RMS_to_secondary_contact_impulse'],flush=True)
    report = dict(status='TYPE25 native DEF impulse convention verified for the current writer/runtime; not a calibration',
                  definition='Contact impulse on SECONDARY = (FNX,FNY,FNZ) - (FTX,FTY,FTZ). Opposite resultant acts on MAIN.',
                  independent_balance='C_secondary = mass*velocity - sum(prescribed-mover REAC impulses) = sum(fixed-target REAC impulses)',
                  runtime_scope='Current TYPE25 writer settings, face contact Iedge=1000, official OpenCourant latest-20261011; no universal sign assertion for other interface algorithms',
                  analysis_script_sha256=digest(__file__),cases=results,
                  execution_of_new_cases=not args.analyze_existing,
                  limitations=['All solid translations are prescribed and solids penetrate; these are bookkeeping controls, not valid compression specimens.',
                               'Large global energy errors are retained; no mechanical, energy, quasi-static or experimental qualification follows.',
                               'A generation hash is not an execution-time metadata binding; retained historical cases without generation.json are explicitly marked.'])
    (output/'force_audit.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__ == '__main__':
    main()
