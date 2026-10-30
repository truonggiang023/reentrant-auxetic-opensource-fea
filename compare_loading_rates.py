"""Compare completed, bound native RP runs differing only in loading duration.

Uses observed relative RP displacement and signed REACY interval reactions.
A rate-independent material's dynamic loading sensitivity is not a calibrated
constitutive strain-rate law. No material inputs or curves are rescaled here.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from assess_courant import assess
from reentrant_fea.calibration import curve_objective
from reentrant_fea.paper_metrics import curve_metrics
from calibration_study import reference_curves


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def normalized_starter(meta):
    lines=Path(meta['starter']).read_text().splitlines()
    # Only model title/prefix and the time coordinate of FUNCT20 may differ.
    result=[];skip=False
    for line in lines:
        if line=='/FUNCT/20':skip=True;continue
        if skip and line.startswith('/'):skip=False
        if skip:continue
        result.append('<MODEL_PREFIX>' if line==meta['prefix'] else line)
    return result


def normalized_engine(meta):
    records=[];key=None;data=[]
    for line in Path(meta['engine']).read_text().splitlines()+['/<END>']:
        line=line.strip()
        if not line or line.startswith('#'):continue
        if line.startswith('/'):
            if key is not None:
                if key.startswith('/RUN/'):
                    if len(data)!=1 or abs(float(data[0])-meta['duration_seconds'])>meta['duration_seconds']*1e-10:
                        raise ValueError('Engine end time differs from declared duration')
                    records.append(('/RUN',[]))
                elif not key.startswith(('/ANIM','/PRINT','/TFILE','/RFILE')):
                    records.append((key,data))
            key,data=line,[]
        else:data.append(line)
    if sum(key=='/RUN' for key,_ in records)!=1:
        raise ValueError('Rate comparison requires one native run segment')
    return records


def compare(fast,slow,reference_dir):
    metas=[];assessments=[]
    for prefix in (fast,slow):
        mp=Path(str(prefix)+'.radioss.json');hp=Path(str(prefix)+'T01.csv')
        m=json.loads(mp.read_text());a=assess(mp,hp,force_source='rp')
        if not a['completed'] or not a['metadata_binding_verified']:
            raise ValueError('Both RP cases must have completed, bound observed motion')
        metas.append(m);assessments.append(a)
    fm,sm=metas
    if not fm['duration_seconds']<sm['duration_seconds']:
        raise ValueError('Slow duration must exceed fast duration')
    if normalized_starter(fm)!=normalized_starter(sm):
        raise ValueError('Rate comparison starter physics differ beyond amplitude timing')
    if normalized_engine(fm)!=normalized_engine(sm):
        raise ValueError('Rate comparison engine physics differ beyond duration/output')
    for field in ('height_mm','reference_area_mm2','target_nominal_strain','ams','platen_reference_points','fracture_locus'):
        if fm.get(field)!=sm.get(field):raise ValueError('Rate control differs in '+field)
    a,b=[np.asarray(m['displacement_fraction_samples']) for m in metas]
    if (a.shape!=b.shape or not np.allclose(a[:,1],b[:,1],rtol=0,atol=1e-12)
            or not np.allclose(a[:,0]/fm['duration_seconds'],b[:,0]/sm['duration_seconds'],rtol=0,atol=1e-12)):
        raise ValueError('Normalized displacement histories differ')
    curves=[np.asarray([[r['strain'],r['nominal_stress_MPa']] for r in x['reaction_history']]) for x in assessments]
    lo,hi=.005,.119
    for c in curves:
        if c[0,0]>lo or c[-1,0]<hi or np.any(np.diff(c[:,0])<=0):
            raise ValueError('Observed strain coverage must span .005-.119 monotonically')
    x=np.linspace(lo,hi,400);y=[np.interp(x,c[:,0],c[:,1]) for c in curves]
    references=reference_curves(Path(reference_dir)/'dong2019_fig12_thin.csv')
    report=dict(status='completed loading-duration comparison; material calibration not established',
                calibrated=False,constitutive_rate_law_fitted=False,
                same_starter_physics_verified=True,interval_strain=[lo,hi],
                relative_L2_fast_minus_slow=float(np.linalg.norm(y[0]-y[1])/np.linalg.norm(y[1])),
                mean_stress_fast_minus_slow_MPa=float(np.mean(y[0]-y[1])),cases=[])
    for prefix,m,a,c in zip((fast,slow),metas,assessments,curves):
        inside=c[(c[:,0]>lo)&(c[:,0]<hi),0]
        grid=np.r_[lo,inside,hi]
        bounded=np.column_stack((grid,np.interp(grid,c[:,0],c[:,1])))
        metrics=(curve_metrics(bounded) if np.all(bounded[:,1]>=0) else
                 dict(available=False,reason='Negative signed stress within fit window; no clipping'))
        mp=Path(str(prefix)+'.radioss.json');receipt=Path(str(prefix)+'.run.json')
        report['cases'].append(dict(prefix=str(prefix),duration_seconds=m['duration_seconds'],
            nominal_mean_strain_rate_per_s=m['target_nominal_strain']/m['duration_seconds'],
            peak_nominal_strain_rate_per_s=m.get('peak_imposed_speed_mm_per_s',
                (1. if m.get('loading_profile')=='linear' else 1.875)*m['height_mm']*m['target_nominal_strain']/m['duration_seconds'])/m['height_mm'],
            metadata_sha256=sha(mp),execution_receipt_sha256=sha(receipt),
            last_force_sample=a['reaction_history'][-1],
            global_energy_screen=a['global_energy_screen'],core_energy_screen=a['core_energy_screen'],
            metrics=metrics,metrics_scope='Observed .005-.119 fit window only; raw signed histories retained',
            negative_samples_outside_fit=int(np.count_nonzero((c[:,1]<0)&((c[:,0]<lo)|(c[:,0]>hi)))),
            paper_discrepancy=curve_objective(c,references,(.03,.119))))
    return report,curves


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fast',type=Path,required=True);p.add_argument('--slow',type=Path,required=True)
    p.add_argument('--reference-dir',type=Path,default=Path(__file__).resolve().parent/'validation/paper_reference');p.add_argument('--output',type=Path,required=True)
    p.add_argument('--wait-seconds',type=float,default=0,help='Wait for actual execution receipt; never treats elapsed time as completion')
    a=p.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
    status=a.output.with_suffix('.status.json');deadline=time.monotonic()+a.wait_seconds
    while not Path(str(a.slow)+'.run.json').exists():
        logfile=Path(str(a.slow)+'.log')
        if logfile.exists() and 'Traceback (most recent call last)' in logfile.read_text(errors='replace'):
            status.write_text(json.dumps({'state':'solver_failed','log':str(logfile),'calibrated':False},indent=2));raise RuntimeError('Slow run failed; inspect '+str(logfile))
        status.write_text(json.dumps({'state':'waiting_for_verified_solver_completion','slow_prefix':str(a.slow),'calibrated':False},indent=2))
        if time.monotonic()>=deadline:raise RuntimeError('Comparison pending actual solver completion')
        time.sleep(5)
    report,curves=compare(a.fast,a.slow,a.reference_dir)
    a.output.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,4.5))
    for c,case in zip(curves,report['cases']):ax.plot(c[:,0],c[:,1],label=f"T={case['duration_seconds']:g} s (native RP)")
    for i,c in enumerate(reference_curves(a.reference_dir/'dong2019_fig12_thin.csv'),1):ax.plot(*np.asarray(c).T,ls='--',label=f'Paper experiment {i}')
    ax.set_xlim(0,.12);ax.set_xlabel('Observed nominal compressive strain');ax.set_ylabel('Nominal stress (MPa)');ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(a.output.with_suffix('.png'),dpi=180);plt.close(fig)
    status.write_text(json.dumps({'state':'comparison_completed','report':str(a.output.with_suffix('.json')),'calibrated':False},indent=2))
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
