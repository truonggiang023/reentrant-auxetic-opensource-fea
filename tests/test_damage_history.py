from reentrant_fea.damage_history import evaluate_damage_csv

def test_damage_history(tmp_path):
    inp=tmp_path/'in.csv';out=tmp_path/'out.csv'
    inp.write_text('step,element_id,peeq,triaxiality,characteristic_length\n0,1,0,0.333,1\n1,1,0.047,0.333,1\n2,1,0.049,0.333,1\n')
    result=evaluate_damage_csv(inp,out,[(.333,.047)],.001)
    assert result[-1]['failed']==1
    assert out.exists()
