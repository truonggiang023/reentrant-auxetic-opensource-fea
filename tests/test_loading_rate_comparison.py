import json
from pathlib import Path
import pytest
import compare_loading_rates as mod


def pair(tmp_path,monkeypatch):
    prefixes=[];assessments={}
    for name,duration,scale in [('fast',.01,1.1),('slow',.04,1.)]:
        prefix=tmp_path/name;deck=tmp_path/(name+'_0000.rad')
        deck.write_text(f'/BEGIN\n{name}\n/FUNCT/20\nTime fraction\n0 0\n{duration} 1\n/NODE\n1 0 0 0\n/END\n')
        engine=tmp_path/(name+'_0001.rad')
        engine.write_text(f'/RUN/{name}/1\n{duration}\n/DT/AMS/2\n.67 1e-6\n')
        m=dict(prefix=name,starter=str(deck),engine=str(engine),duration_seconds=duration,height_mm=48.,reference_area_mm2=1900.,target_nominal_strain=.12,ams={'enabled':True},platen_reference_points=True,fracture_locus=None,displacement_fraction_samples=[[0.,0.],[duration,1.]])
        Path(str(prefix)+'.radioss.json').write_text(json.dumps(m));Path(str(prefix)+'.run.json').write_text('{}');prefixes.append(prefix)
        assessments[name]=dict(completed=True,metadata_binding_verified=True,global_energy_screen={},core_energy_screen={},reaction_history=[dict(strain=x,nominal_stress_MPa=y*scale) for x,y in [(0.,-.001),(.05,.3),(.12,.35)]])
    monkeypatch.setattr(mod,'assess',lambda mp,hp,**kw:assessments[Path(mp).name.split('.')[0]])
    return prefixes


def test_rate_comparison_keeps_signed_startup_and_matches_physics(tmp_path,monkeypatch):
    f,s=pair(tmp_path,monkeypatch)
    r,c=mod.compare(f,s,Path(__file__).resolve().parents[1]/'validation/paper_reference')
    assert r['relative_L2_fast_minus_slow']==pytest.approx(.1)
    assert r['same_starter_physics_verified']
    assert r['cases'][0]['negative_samples_outside_fit']==1
    assert c[0][0,1]<0
    assert not r['calibrated'] and not r['constitutive_rate_law_fitted']


def test_rate_comparison_rejects_changed_geometry(tmp_path,monkeypatch):
    f,s=pair(tmp_path,monkeypatch)
    p=tmp_path/'slow_0000.rad';p.write_text(p.read_text().replace('1 0 0 0','1 1 0 0'))
    with pytest.raises(ValueError,match='physics differ'):
        mod.compare(f,s,Path(__file__).resolve().parents[1]/'validation/paper_reference')


def test_rate_comparison_rejects_extra_dynamic_damping(tmp_path,monkeypatch):
    f,s=pair(tmp_path,monkeypatch)
    p=tmp_path/'slow_0001.rad';p.write_text(p.read_text()+'/DYREL\n1. 1.\n')
    with pytest.raises(ValueError,match='engine physics'):
        mod.compare(f,s,Path(__file__).resolve().parents[1]/'validation/paper_reference')
