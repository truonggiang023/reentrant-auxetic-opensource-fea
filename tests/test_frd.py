from reentrant_fea.frd import committed_fields


def test_reject_best_unconverged_frame(tmp_path):
    path=tmp_path/'failed.frd'
    path.write_text('  100CL  101 1.00000E-01 1\n -4  DISP 4 1\n -1         1 1.00000E-01 0.00000E+00 0.00000E+00\n -3\n'
                    '  100CL  101 1.00100E-01 1\n -4  DISP 4 1\n -1         1 1.00000E+09 0.00000E+00 0.00000E+00\n -3\n')
    path.with_suffix('.sta').write_text('1 1 1 3 .1 .1 .1\n1 2 1U 50 .1 .1 .0001\n')
    rows=list(committed_fields(path))
    assert len(rows)==1
    assert rows[0][2][1][0]==.1


def test_projected_overlap_detects_missed_contact():
    import numpy as np
    from reentrant_fea.frd import projected_overlap
    n=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,0,1],[0,1,1],
                [2,0,0],[3,0,0],[2,1,0],[2,0,1],[3,0,1],[2,1,1]],float)
    e=np.array([[0,1,2,3,4,5],[6,7,8,9,10,11]])
    u={i+1:[0.,0.,0.] for i in range(12)}
    assert projected_overlap(n,e,u)['overlap_pairs']==0
    for i in range(7,13):u[i]=[-1.5,0.,0.]
    report=projected_overlap(n,e,u)
    assert report['overlap_pairs']==1
    assert abs(report['max_overlap_mm2']-.125)<1e-12
