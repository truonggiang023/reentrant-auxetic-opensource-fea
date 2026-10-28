#!/usr/bin/env python3
"""Conforming in-plane refinement and static CCX elastic-response comparison.

Every original coordinate/ID and both extrusion planes are retained. Shared
edge midpoints and per-face centers split each coarse quad into four children;
no depth subdivision, contact, plasticity or fracture is introduced. One
refinement comparison is a sensitivity check, not a convergence proof.
"""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import shutil
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.probe_matched_elastic_lattice import run_case, sha


def geometry(mesh):
    nodes = {int(row[0]): np.asarray(row[1:], dtype=float) for row in mesh['nodes']}
    xyz = np.asarray([row[1:] for row in mesh['nodes']], dtype=float)
    cells = np.asarray([[nodes[n] for n in row[1:]] for row in mesh['elements']])
    if cells.ndim != 3 or cells.shape[1:] != (8, 3):
        raise ValueError('Expected eight-node extruded hexes')
    if not np.isfinite(cells).all():
        raise ValueError('Nonfinite mesh coordinates')
    front, back = cells[:, :4], cells[:, 4:]
    if (np.max(np.abs(front[:, :, :2]-back[:, :, :2])) > 1e-10
            or np.max(np.ptp(front[:, :, 2], axis=1)) > 1e-10
            or np.max(np.ptp(back[:, :, 2], axis=1)) > 1e-10
            or np.any(back[:, 0, 2] <= front[:, 0, 2])):
        raise ValueError('Expected matching front/back extruded quads in increasing Z order')
    edge = np.roll(front[:, :, :2], -1, axis=1)-front[:, :, :2]
    previous = -np.roll(edge, 1, axis=1)
    jac = edge[:, :, 0]*previous[:, :, 1]-edge[:, :, 1]*previous[:, :, 0]
    area = .5*np.sum(front[:, :, 0]*np.roll(front[:, :, 1], -1, axis=1)
                     - front[:, :, 1]*np.roll(front[:, :, 0], -1, axis=1), axis=1)
    if np.any(area <= 0) or np.any(jac <= 0):
        raise ValueError('Nonpositive quad area or corner Jacobian')
    depths = back[:, 0, 2]-front[:, 0, 2]
    planes = np.unique(xyz[:, 2])
    if len(planes) != 2:
        raise ValueError('This in-plane benchmark requires exactly two original Z planes')
    return dict(node_count=len(nodes), element_count=len(cells),
                original_z_planes_mm=planes.tolist(),
                summed_plan_area_mm2=float(area.sum()),
                summed_extruded_volume_mm3=float(np.sum(area*depths)),
                minimum_corner_jacobian_mm2=float(jac.min()),
                maximum_corner_jacobian_ratio=float((jac.max(axis=1)/jac.min(axis=1)).max()),
                all_corner_jacobians_positive=True,
                dimensions_mm=dict(width=float(np.ptp(xyz[:, 0])), height=float(np.ptp(xyz[:, 1])),
                                   depth=float(np.ptp(xyz[:, 2]))))


def refine_in_plane(mesh):
    """Split each extrusion quad into four, sharing mids across every old edge."""
    result = deepcopy(mesh)
    coordinates = {int(row[0]): np.asarray(row[1:], dtype=float) for row in mesh['nodes']}
    original = {n: xyz.copy() for n, xyz in coordinates.items()}
    next_node = max(coordinates)+1
    edge_nodes, face_nodes = {}, {}

    def new_node(xyz):
        nonlocal next_node
        ident = next_node
        coordinates[ident] = xyz
        next_node += 1
        return ident

    def split_face(face):
        a, b, c, d = map(int, face)
        mid = []
        for first, second in ((a, b), (b, c), (c, d), (d, a)):
            key = tuple(sorted((first, second)))
            if key not in edge_nodes:
                edge_nodes[key] = new_node(.5*(coordinates[first]+coordinates[second]))
            mid.append(edge_nodes[key])
        key = tuple(sorted(face))
        if key not in face_nodes:
            face_nodes[key] = new_node(np.mean([coordinates[n] for n in face], axis=0))
        ab, bc, cd, da = mid
        center = face_nodes[key]
        return [(a, ab, center, da), (ab, b, bc, center),
                (center, bc, c, cd), (da, center, cd, d)]

    children, parent_map = [], []
    for row in mesh['elements']:
        parent, *verts = row
        front, back = split_face(verts[:4]), split_face(verts[4:])
        for child_front, child_back in zip(front, back):
            ident = len(children)+1
            children.append([ident, *child_front, *child_back])
            parent_map.append([ident, parent])
    result['nodes'] = [[n, *coordinates[n].tolist()] for n in sorted(coordinates)]
    result['elements'] = children
    ymin = min(xyz[1] for xyz in coordinates.values())
    ymax = max(xyz[1] for xyz in coordinates.values())
    result['boundary'] = dict(
        all_z_fixed=sorted(coordinates),
        bottom_y_fixed=[n for n in sorted(coordinates) if abs(coordinates[n][1]-ymin) <= 1e-7],
        top_y_displaced=[n for n in sorted(coordinates) if abs(coordinates[n][1]-ymax) <= 1e-7],
        one_bottom_x_pin=mesh['boundary']['one_bottom_x_pin'],
        prescribed_top_y_mm=mesh['boundary']['prescribed_top_y_mm'])
    if not all(np.array_equal(coordinates[n], xyz) for n, xyz in original.items()):
        raise ValueError('Refinement changed original coordinates')
    before, after = geometry(mesh), geometry(result)
    for quantity in ('summed_plan_area_mm2', 'summed_extruded_volume_mm3'):
        if not np.isclose(before[quantity], after[quantity], rtol=1e-12, atol=1e-9):
            raise ValueError('Refinement changed mesh area/volume')
    if before['original_z_planes_mm'] != after['original_z_planes_mm']:
        raise ValueError('Refinement changed extrusion planes')
    result['refinement'] = dict(method='Shared edge midpoints and bilinear face centers; four children per old quad; no depth refinement',
                                parent_element_count=len(mesh['elements']), original_node_count=len(original),
                                original_nodes_preserved_exactly=True, parent_child_map=parent_map,
                                geometry_before=before, geometry_after=after)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh', type=Path,
                    default=Path('/workspace/fea-research/matched-elastic-lattice-final/mesh_and_bcs.json'))
    ap.add_argument('--coarse-report', type=Path, default=ROOT/'validation/matched_elastic_ccx_audit.json')
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--report', type=Path)
    ap.add_argument('--types', nargs='+', choices=['C3D8R', 'C3D8I', 'C3D8'], default=['C3D8R', 'C3D8I'])
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--timeout', type=float, default=300.)
    ap.add_argument('--ccx', default='/workspace/bin/ccx')
    ap.add_argument('--prepare-only', action='store_true')
    args = ap.parse_args()
    if args.threads < 1 or not np.isfinite(args.timeout) or args.timeout <= 0:
        raise ValueError('Require positive threads and finite positive timeout')
    mesh = json.loads(args.mesh.read_text())
    coarse = json.loads(args.coarse_report.read_text())
    if coarse['shared_mesh_sha256'] != sha(args.mesh):
        raise ValueError('Coarse report belongs to a different mesh')
    fine = refine_in_plane(mesh)
    directory = args.output.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    shared = directory/'mesh_and_bcs.json'
    if shared.exists():
        if json.loads(shared.read_text()) != fine:
            raise ValueError('Existing refined mesh differs; preserve previous output')
    else:
        shared.write_text(json.dumps(fine, separators=(',', ':'))+'\n')
    print(json.dumps(dict(shared_mesh=str(shared), shared_mesh_sha256=sha(shared),
                          geometry=geometry(fine), top_nodes=len(fine['boundary']['top_y_displaced']),
                          bottom_nodes=len(fine['boundary']['bottom_y_fixed'])), indent=2), flush=True)
    if args.prepare_only:
        return
    executable = shutil.which(args.ccx)
    if executable is None:
        raise FileNotFoundError(args.ccx)
    records = []
    for typ in args.types:
        record = run_case(directory, fine, typ, executable, args.threads, args.timeout)
        reference = next((c for c in coarse['cases'] if c['element_type'] == typ and c['usable']), None)
        if record['usable'] and reference:
            record['coarse_stiffness_N_per_mm'] = reference['stiffness_N_per_mm']
            record['relative_fine_minus_coarse_stiffness'] = record['stiffness_N_per_mm']/reference['stiffness_N_per_mm']-1
            record['coarse_reaction_work_Nmm'] = reference['reaction_work_Nmm']
        records.append(record)
        print(json.dumps(record, indent=2), flush=True)
    report = dict(status='single conforming refinement sensitivity check; not a mesh-convergence proof or fracture calibration',
                  script='scripts/probe_elastic_mesh_refinement.py', script_sha256=sha(__file__),
                  solver_driver_sha256=sha(ROOT/'scripts/probe_matched_elastic_lattice.py'),
                  coarse_mesh=str(args.mesh), coarse_mesh_sha256=sha(args.mesh),
                  coarse_report=str(args.coarse_report), coarse_report_sha256=sha(args.coarse_report),
                  refined_mesh=str(shared), refined_mesh_sha256=sha(shared),
                  refinement={k:v for k,v in fine['refinement'].items() if k != 'parent_child_map'},
                  node_counts={'top':len(fine['boundary']['top_y_displaced']), 'bottom':len(fine['boundary']['bottom_y_fixed'])},
                  unchanged_x_pin=fine['boundary']['one_bottom_x_pin'],
                  material=fine['material'], prescribed_top_y_mm=fine['boundary']['prescribed_top_y_mm'],
                  cases=records,
                  limitations=['Same original polygon and Z depth; finite-thickness discretization and element formulation sensitivity remain.',
                               'No paper curve was used to tune stiffness, E or the refinement.',
                               'A converged elastic reference requires checking subsequent refinements and integration/locking sensitivity; fracture/large-deformation/contact convergence is separate.'])
    target = args.report or directory/'refinement_report.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
