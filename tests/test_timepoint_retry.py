"""Numerical regression for accepted output-point retry bookkeeping."""
import ctypes
from pathlib import Path
import subprocess
import numpy as np


def test_retried_completed_point_cannot_create_zero_increment(tmp_path):
    header=Path(__file__).resolve().parents[1]/'native'
    source=tmp_path/'driver.c'
    source.write_text('#define ITG int\n#include "timepoints.h"\n'
                      'int advance(int*n,int end,const double*p,double now,double period){'
                      'return reentrant_advance_completed_points(n,end,p,now,period);}\n')
    so=tmp_path/'driver.so'
    subprocess.run(['gcc','-shared','-fPIC','-I',str(header),str(source),'-lm','-o',str(so)],check=True)
    fun=ctypes.CDLL(str(so)).advance
    fun.argtypes=[ctypes.POINTER(ctypes.c_int),ctypes.c_int,
                  np.ctypeslib.ndpointer(dtype=np.float64),ctypes.c_double,ctypes.c_double]
    # The first scheduled interval failed, cleared itp and rewound inext;
    # a retry then accepted time .185. It must not request .185-.185 next.
    points=np.array([.184615384615,0.,.185,0.,.19,0.],dtype='d')
    nxt=ctypes.c_int(2)
    assert fun(ctypes.byref(nxt),3,points,.185,1.)==1
    assert nxt.value==3 and points[2*nxt.value-2]-.185>0
    # Future points must not be silently swallowed by a broad tolerance.
    nxt=ctypes.c_int(2)
    assert fun(ctypes.byref(nxt),3,points,.185-1e-9,1.)==0
    assert nxt.value==2
    nxt=ctypes.c_int(2)
    assert fun(ctypes.byref(nxt),3,points,.185-1e-17,1.)==1
