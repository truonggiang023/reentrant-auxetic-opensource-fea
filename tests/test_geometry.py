import json
from pathlib import Path
from reentrant_fea.geometry import core_polygon,parameters

def config():return json.loads(Path('configs/dong2019_thin.json').read_text())
def test_geometry_valid():
    p=core_polygon(config())
    assert p.is_valid and p.area>0
    assert p.bounds[3]>47.7

def test_morphology_shift():
    c=config(); c['lambda_m']=0
    assert parameters(c)[1]==0
