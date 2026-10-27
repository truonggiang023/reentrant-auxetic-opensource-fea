#!/usr/bin/env python3
"""Matched small-strain CCX elasticity on an exact retained native lattice.

No contact, platens, density, plasticity, damage or geometric nonlinearity is
included. This isolates the assembled elastic response of the same core mesh.
The exported mesh_and_bcs.json is the protocol for independent native controls.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
E, NU, DELTA = 69000., .3, .01


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def extract_mesh(path):
    nodes, elements = {}, []
    mode = ''
    for line in Path(path).read_text().splitlines():
        if line.startswith('/'):
            mode = line
            continue
        if mode == '/NODE':
            p = line.split()
            if p:
                nodes[int(p[0])] = list(map(float, p[1:4]))
        elif mode == '/BRICK/1':
            p = list(map(int, line.split()))
            if p:
                if len(p) != 9:
                    raise ValueError('Expected eight-node core bricks')
                elements.append(p)
    used = sorted({n for row in elements for n in row[1:]})
    if not elements or any(n not in nodes for n in used):
        raise ValueError('Missing core mesh in native starter')
    xyz = np.asarray([nodes[n] for n in used])
    bottom = [n for n in used if abs(nodes[n][1]-xyz[:, 1].min()) <= 1e-7]
    top = [n for n in used if abs(nodes[n][1]-xyz[:, 1].max()) <= 1e-7]
    center = .5*(xyz[:, 0].min()+xyz[:, 0].max())
    anchor = min(bottom, key=lambda n: (abs(nodes[n][0]-center), n))
    return dict(source_starter=str(Path(path).resolve()), source_sha256=sha(path),
                nodes=[[n, *nodes[n]] for n in used], elements=elements,
                boundary=dict(all_z_fixed=used, bottom_y_fixed=bottom, top_y_displaced=top,
                              one_bottom_x_pin=anchor, prescribed_top_y_mm=-DELTA),
                material=dict(young_MPa=E, poisson=NU),
                dimensions_mm=dict(width=float(np.ptp(xyz[:, 0])),
                                   height=float(np.ptp(xyz[:, 1])), depth=float(np.ptp(xyz[:, 2]))),
                protocol='Exact original core IDs/coordinates/connectivity; all Z fixed, bottom Y fixed, top Y=-.01mm; one bottom center-X pin; all remaining X/Y free. Linear elasticity, no NLGEOM/contact/plastic/damage/platens.',
                reference_area_mm2=1900.)


def nset(lines, name, values):
    lines.append('*NSET,NSET='+name)
    lines += [', '.join(map(str, values[i:i+12])) for i in range(0, len(values), 12)]


def make_deck(mesh, path, typ):
    lines = ['*HEADING', 'Matched linear-elastic lattice; exact retained native core mesh', '*NODE']
    lines += [f'{n}, {x:.15g}, {y:.15g}, {z:.15g}' for n, x, y, z in mesh['nodes']]
    lines.append(f'*ELEMENT,TYPE={typ},ELSET=CORE')
    lines += [', '.join(map(str, row)) for row in mesh['elements']]
    bcs = mesh['boundary']
    for name, ids in [('ALL', bcs['all_z_fixed']), ('BOTTOM', bcs['bottom_y_fixed']),
                      ('TOP', bcs['top_y_displaced']), ('XPIN', [bcs['one_bottom_x_pin']])]:
        nset(lines, name, ids)
    lines += ['*MATERIAL,NAME=CONFIRMED_ELASTIC', '*ELASTIC', f'{E:g}, {NU:g}',
              '*SOLID SECTION,ELSET=CORE,MATERIAL=CONFIRMED_ELASTIC', '*BOUNDARY',
              'ALL,3,3,0.', 'BOTTOM,2,2,0.', 'XPIN,1,1,0.', '*STEP', '*STATIC', '1., 1., 1.e-6, 1.',
              '*BOUNDARY', f'TOP,2,2,{-DELTA:.15g}',
              '*NODE PRINT,NSET=TOP,TOTALS=ONLY', 'RF',
              '*NODE PRINT,NSET=BOTTOM,TOTALS=ONLY', 'RF',
              '*EL PRINT,ELSET=CORE,TOTALS=ONLY', 'ELSE,EVOL',
              '*NODE FILE', 'U,RF', '*EL FILE', 'S', '*END STEP']
    Path(path).write_text('\n'.join(lines)+'\n')


def final_total(text, label):
    # Each requested total is followed by a blank line and its numeric row.
    pattern = re.compile(r'^\s*'+re.escape(label)+r'.*\n\s*\n\s*([^\n]+)', re.M)
    rows = pattern.findall(text)
    if not rows:
        raise ValueError('Missing actual DAT total: '+label)
    return [float(x.replace('D', 'E')) for x in rows[-1].split()]


def run_case(directory, mesh, typ, executable, threads, timeout):
    path = directory/(typ.lower()+'.inp')
    if path.exists():
        raise FileExistsError('Use a fresh output directory: '+str(path))
    make_deck(mesh, path, typ)
    log = path.with_suffix('.log')
    env = os.environ.copy()
    for name in ('OMP_NUM_THREADS', 'CCX_NPROC_STIFFNESS', 'CCX_NPROC_RESULTS', 'CCX_NPROC_EQUATION_SOLVER'):
        env[name] = str(threads)
    started = time.monotonic()
    with log.open('w') as stream:
        result = subprocess.run([executable, '-i', path.stem], cwd=directory, env=env,
                                stdout=stream, stderr=subprocess.STDOUT, timeout=timeout)
    record = dict(element_type=typ, exit_code=result.returncode,
                  elapsed_wall_seconds=time.monotonic()-started, input_sha256=sha(path))
    text = log.read_text(errors='replace')
    version = re.search(r'CalculiX Version ([^,\s]+)', text)
    record.update(solver_reported_version=version.group(1) if version else None,
                  solver_warning_count=text.upper().count('*WARNING'))
    dat = path.with_suffix('.dat')
    if result.returncode or '*ERROR' in text.upper() or 'Job finished' not in text or not dat.is_file():
        record.update(usable=False, reason='Solver failure or missing completion/output; inspect retained log')
        return record
    status = path.with_suffix('.sta').read_text().splitlines()
    accepted = [l.split() for l in status if len(l.split()) == 7 and l.split()[0].isdigit() and l.split()[2].isdigit()]
    if not accepted or abs(float(accepted[-1][5])-1.) > 1e-8:
        record.update(usable=False, reason='Did not reach full prescribed small load')
        return record
    data = dat.read_text()
    top = final_total(data, 'total force (fx,fy,fz) for set TOP')
    bottom = final_total(data, 'total force (fx,fy,fz) for set BOTTOM')
    energy = final_total(data, 'total internal energy for set CORE')[0]
    volume = final_total(data, 'total volume for set CORE')[0]
    stiffness = -top[1]/DELTA
    work = .5*(-top[1])*DELTA
    record.update(usable=True, top_total_reaction_N=top, bottom_total_reaction_N=bottom,
                  y_force_balance_N=top[1]+bottom[1], stiffness_N_per_mm=stiffness,
                  nominal_effective_modulus_MPa=stiffness*mesh['dimensions_mm']['height']/1900.,
                  reaction_work_Nmm=work, solver_internal_energy_Nmm=energy,
                  relative_reported_energy_minus_work=(energy-work)/work,
                  solver_total_volume_mm3=volume,
                  output_sha256={p.name: sha(p) for p in (dat, log, path.with_suffix('.sta'), path.with_suffix('.frd'))})
    return record


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, default=ROOT/'results/opencourant/thin_free_x_ams1e6_0000.rad')
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--report', type=Path)
    ap.add_argument('--types', nargs='+', choices=['C3D8R', 'C3D8', 'C3D8I'], default=['C3D8R', 'C3D8', 'C3D8I'])
    ap.add_argument('--ccx', default='/workspace/bin/ccx')
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--timeout', type=float, default=300.)
    ap.add_argument('--prepare-only', action='store_true')
    args = ap.parse_args()
    if args.threads < 1 or not np.isfinite(args.timeout) or args.timeout <= 0:
        raise ValueError('Require positive thread count and finite positive timeout')
    directory = args.output.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    mesh = extract_mesh(args.source)
    shared = directory/'mesh_and_bcs.json'
    if shared.exists():
        if json.loads(shared.read_text()) != mesh:
            raise ValueError('Existing shared mesh differs; preserve previous benchmark output')
    else:
        shared.write_text(json.dumps(mesh, separators=(',', ':'))+'\n')
    print(json.dumps(dict(shared_mesh=str(shared), shared_mesh_sha256=sha(shared),
                          nodes=len(mesh['nodes']), elements=len(mesh['elements']),
                          top_nodes=len(mesh['boundary']['top_y_displaced']),
                          bottom_nodes=len(mesh['boundary']['bottom_y_fixed']),
                          x_pin=mesh['boundary']['one_bottom_x_pin']), indent=2), flush=True)
    if args.prepare_only:
        return
    executable = shutil.which(args.ccx)
    if executable is None:
        raise FileNotFoundError(args.ccx)
    records = []
    for typ in args.types:
        record = run_case(directory, mesh, typ, executable, args.threads, args.timeout)
        records.append(record)
        print(json.dumps(record, indent=2), flush=True)
    report = dict(status='matched linear-elastic component diagnostic; no fracture calibration claim',
                  script='scripts/probe_matched_elastic_lattice.py', script_sha256=sha(__file__),
                  shared_mesh=str(shared), shared_mesh_sha256=sha(shared),
                  source_starter=mesh['source_starter'], source_starter_sha256=mesh['source_sha256'],
                  node_count=len(mesh['nodes']), element_count=len(mesh['elements']),
                  material=mesh['material'], dimensions_mm=mesh['dimensions_mm'],
                  boundary= {k: len(v) if isinstance(v,list) else v for k,v in mesh['boundary'].items()},
                  protocol=mesh['protocol'], geometric_nonlinearity=False,
                  ccx_executable=executable, ccx_wrapper_sha256=sha(executable),
                  threads=args.threads, cases=records,
                  energy_note='Reaction work is F*delta/2 for this linear load. Compare reported ELSE separately; do not assume reduced-integration stabilization energy is included.')
    target = args.report or directory/'ccx_elastic_report.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
