"""Fix CCX 2.20 retry bookkeeping that can schedule a zero time increment."""
import hashlib
from pathlib import Path
import sys

path=Path(sys.argv[1])/'checkconvergence.c'
expected='ce660d264f7d74e196a406ebde97aab8d14bb68bf5794512a4dc47405805de40'
old='\tif((*itpamp>0)&&(*idrct==0)){\n\t    if(*itp==1){'
new='''\tif((*itpamp>0)&&(*idrct==0)){
            /* REENTRANT: a retry may finish exactly at a point while itp
               was cleared and inext decremented during the cutback. */
            double tpnow=(namta[3**itpamp-1]<0)?*ttime+*time:*time;
            if(reentrant_advance_completed_points(inext,
                    namta[3**itpamp-2],amta,tpnow,*tper)>0) *itp=1;
\t    if(*itp==1){'''
include='#include "reentrant_timepoints.h"\n'
data=path.read_text();original=data.replace(include,'').replace(new,old)
if hashlib.sha256(original.encode()).hexdigest()!=expected:
    raise SystemExit('Unexpected checkconvergence.c checksum')
if original.count(old)!=1:raise SystemExit('Timepoint patch context is ambiguous')
patched=original.replace('#include "CalculiX.h"\n','#include "CalculiX.h"\n'+include).replace(old,new)
if data!=patched:path.write_text(patched)
