"""An otherwise plausible curve cannot be attached to a different solver run."""
import hashlib
import json

import pytest

from assess_courant import assess


def evidence(tmp_path):
    prefix = 'model'
    files = {name: tmp_path/(prefix+name) for name in
             ('_0000.rad', '_0001.rad', '_0000.out', '_0001.out', 'T01', 'T01.csv')}
    for suffix, path in files.items():
        path.write_text('input' if suffix.endswith('.rad') else 'output')
    files['_0001.out'].write_text('NORMAL TERMINATION\n')
    files['T01.csv'].write_text('time,INTERNAL ENERGY,KINETIC ENERGY\n0,0,0\n1,1,0\n')
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    inputs = {files[key].name: digest(files[key]) for key in ('_0000.rad', '_0001.rad')}
    metadata = tmp_path/'model.radioss.json'
    metadata.write_text(json.dumps(dict(prefix=prefix, backend='OpenCourant', sha256=inputs,
                                        duration_seconds=1., wall_mm=.3, time_history_interval_seconds=.01)))
    record = dict(input_sha256=inputs, output_sha256={files[key].name: digest(files[key])
                                                    for key in ('_0000.out', '_0001.out', 'T01', 'T01.csv')})
    return metadata, files, tmp_path/'model.run.json', record


def test_no_execution_binding_is_diagnostic_not_completed(tmp_path):
    metadata, files, manifest, record = evidence(tmp_path)
    result = assess(metadata, files['T01.csv'])
    assert result['engine_normal_termination']
    assert not result['run_binding_verified'] and not result['completed']
    manifest.write_text(json.dumps(record))
    result = assess(metadata, files['T01.csv'])
    assert result['run_binding_verified'] and result['completed']
    assert result['quasi_static_validated'] is False


def test_replacing_history_after_execution_is_rejected(tmp_path):
    metadata, files, manifest, record = evidence(tmp_path)
    manifest.write_text(json.dumps(record))
    # Same duration and plausible values still belong to different evidence.
    files['T01.csv'].write_text('time,INTERNAL ENERGY,KINETIC ENERGY\n0,0,0\n1,2,0\n')
    with pytest.raises(ValueError, match='changed since execution'):
        assess(metadata, files['T01.csv'])


@pytest.mark.parametrize('message', ['ABNORMAL TERMINATION', 'NOT NORMAL TERMINATION'])
def test_abnormal_termination_cannot_complete_even_with_bound_final_history(tmp_path, message):
    metadata, files, manifest, record = evidence(tmp_path)
    files['_0001.out'].write_text(message+'\n')
    record['output_sha256'][files['_0001.out'].name] = hashlib.sha256(files['_0001.out'].read_bytes()).hexdigest()
    manifest.write_text(json.dumps(record))
    result = assess(metadata, files['T01.csv'])
    assert result['run_binding_verified']
    assert not result['engine_normal_termination']
    assert not result['completed']


def test_error_line_cannot_be_overridden_by_normal_trailer(tmp_path):
    metadata, files, manifest, record = evidence(tmp_path)
    files['_0001.out'].write_text('**ERROR : TIME STEP LESS OR EQUAL DTMIN\nNORMAL TERMINATION\n')
    record['output_sha256'][files['_0001.out'].name] = hashlib.sha256(files['_0001.out'].read_bytes()).hexdigest()
    manifest.write_text(json.dumps(record))
    result = assess(metadata, files['T01.csv'])
    assert result['engine_normal_termination'] and result['errors']
    assert not result['completed']


def test_execution_metadata_cannot_be_relabelled(tmp_path):
    metadata, files, manifest, record = evidence(tmp_path)
    record['metadata_sha256'] = hashlib.sha256(metadata.read_bytes()).hexdigest()
    manifest.write_text(json.dumps(record))
    assert assess(metadata, files['T01.csv'])['metadata_binding_verified']
    data = json.loads(metadata.read_text())
    data['wall_mm'] = 1.
    metadata.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='metadata changed since execution'):
        assess(metadata, files['T01.csv'])


@pytest.mark.parametrize('top_displacement,completed', [(-.12,True),(-.06,False)])
def test_rp_completion_requires_observed_motion_and_bound_history(tmp_path,top_displacement,completed):
    import csv
    from reentrant_fea.native_rp_history import VARIABLES
    metadata, files, receipt, record=evidence(tmp_path)
    m=json.loads(metadata.read_text())
    m.update(platen_reference_points=True,loading_node_variables=VARIABLES,
             loading_face_nodes=[1],bottom_support_nodes=[2],height_mm=1.,
             reference_area_mm2=1900.,target_nominal_strain=.12)
    metadata.write_text(json.dumps(m))
    keys=['time']+[f'Loading face nodal histories {node} load-node-{node} var {j+6*(node-1)}'
                   for node in (1,2) for j in range(6)]
    values=[1.,0.,top_displacement,0.,0.,-10.,0.,0.,0.,0.,0.,10.,0.]
    with files['T01.csv'].open('w',newline='') as stream:
        w=csv.writer(stream);w.writerow(keys);w.writerow([0.]*13);w.writerow(values)
    record['metadata_sha256']=hashlib.sha256(metadata.read_bytes()).hexdigest()
    record['output_sha256'][files['T01.csv'].name]=hashlib.sha256(files['T01.csv'].read_bytes()).hexdigest()
    receipt.write_text(json.dumps(record))
    r=assess(metadata,files['T01.csv'],force_source='rp')
    assert r['completed'] is completed
    assert r['metadata_binding_verified']
    assert r['last_observed_nominal_strain']==pytest.approx(-top_displacement)
    assert r['reaction_history'][-1]['nominal_stress_MPa']==pytest.approx(10/1900.)
    with pytest.raises(ValueError,match='RP assessment or'):
        assess(metadata,files['T01.csv'],force_source='rp',force_column=keys[5])
