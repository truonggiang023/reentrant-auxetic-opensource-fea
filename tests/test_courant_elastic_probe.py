"""Regression for diagnostic history partition after shared writer extensions."""
import hashlib
import json
from pathlib import Path

from scripts.probe_courant_elastic_lattice import generate


def test_generated_probe_keeps_top_and_bottom_histories_separate(tmp_path):
    top, bottom = [3, 4, 7, 8], [1, 2, 5, 6]
    xyz = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],
           [0,0,1],[1,0,1],[1,1,1],[0,1,1]]
    shared = dict(nodes=[[i+1, *p] for i,p in enumerate(xyz)],
                  elements=[[1,1,2,3,4,5,6,7,8]],
                  boundary=dict(top_y_displaced=top, bottom_y_fixed=bottom,
                                prescribed_top_y_mm=-.001, all_z_fixed=list(range(1,9)),
                                one_bottom_x_pin=1),
                  material=dict(young_MPa=69000., poisson=.3),
                  reference_area_mm2=1.)
    mesh = tmp_path/'mesh.json'
    mesh.write_text(json.dumps(shared))
    meta = generate(mesh, tmp_path/'probe', duration=.004, ams_step=None, solid=24)
    starter = Path(meta['starter'])
    lines = starter.read_text().splitlines()
    def history_nodes(keyword):
        start = lines.index(keyword)+3
        end = next(i for i in range(start, len(lines)) if lines[i].startswith('/'))
        return [int(line[:10]) for line in lines[start:end]]
    assert history_nodes('/TH/NODE/3') == top
    assert history_nodes('/TH/NODE/4') == bottom
    assert not set(history_nodes('/TH/NODE/3')) & set(bottom)
    assert meta['loading_face_nodes'] == top
    assert meta['bottom_support_nodes'] == bottom
    assert meta['diagnostic']['bottom_nodes'] == bottom
    # The rewritten history is bound, not merely an untracked deck edit.
    assert meta['sha256'][starter.name] == hashlib.sha256(starter.read_bytes()).hexdigest()
