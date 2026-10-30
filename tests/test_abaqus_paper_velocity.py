import ast
from pathlib import Path
from reentrant_fea.plastic import load_plastic


def test_paper_velocity_variant_preserves_material_and_requests_linear_physical_loading():
    root=Path(__file__).resolve().parents[1];p=root/'abaqus/run_paper_velocity012.py';tree=ast.parse(p.read_text())
    values={}
    for n in tree.body:
        if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name):
            try:values[n.targets[0].id]=ast.literal_eval(n.value)
            except (ValueError,TypeError):pass
    assert values['COMPRESSION_STRAIN']==.12
    assert values['PAPER_COMPRESSION_SPEED_MM_MIN']==2.
    assert values['USE_PAPER_PHYSICAL_SPEED'] is True
    assert values['ENABLE_MASS_SCALING'] is False
    assert values['SUBMIT_JOB'] is False
    assert load_plastic(p)==load_plastic(root/'reference_setup_numeric_source.py')
    amplitudes=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='TabularAmplitude']
    assert len(amplitudes)==1
    assert any(k.arg=='smooth' and ast.literal_eval(k.value)==0. for k in amplitudes[0].keywords)
    assert not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='SmoothStepAmplitude' for n in ast.walk(tree))
    assert "assembly.Set(name='CORE_ALL_U3', nodes=core_nodes)" in p.read_text()
