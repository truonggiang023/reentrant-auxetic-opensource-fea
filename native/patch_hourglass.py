"""Scoped CCX 2.20 patch for REENTRANT user-material hourglass control.

The stress-only interface permits no tangent update. Stock resultsmech then
uses elas(1) for hourglass forces, although it can be unset for a user material.
Use a deterministic positive elastic modulus degraded by SDV10 on BOTH force
and stiffness paths, exclusively for @REENTRANT materials. Other laws untouched.
"""
import hashlib
from pathlib import Path
import sys

FILES={
 'resultsmech.f':('94e7f6577af12a8fdb477c51e0fd087f472ad6067b03681f578b6fc266878250',
                  '              call hgforce(fn,elas,a,gs,vl,mi,konl)',
'''              if((amat(1:10).eq.'@REENTRANT').and.
     &           (nstate_.ge.10)) then
                reentrant_hg(1)=elconloc(1)*(1.d0-xstate(10,jj,i))
                call hgforce(fn,reentrant_hg,a,gs,vl,mi,konl)
              else
                call hgforce(fn,elas,a,gs,vl,mi,konl)
              endif'''),
 'e_c3d.f':('68d88d2ffbe73d875129f3c9a6c01530f0978a257d8e69603c927f65ccee2ab7',
              '          call hgstiffness(s,elas,a,gs)',
'''          if((amat(1:10).eq.'@REENTRANT').and.
     &       (nstate_.ge.10)) then
            reentrant_hg(1)=elcon(1,1,imat)*(1.d0-xstate(10,kk,nelem))
            call hgstiffness(s,reentrant_hg,a,gs)
          else
            call hgstiffness(s,elas,a,gs)
          endif''')}

for name,(expected,old,new) in FILES.items():
    path=Path(sys.argv[1])/name
    data=path.read_text()
    # Reconstruct and verify the authenticated original before writing or
    # accepting an already-applied patch. Unknown vendor changes are errors.
    original=data
    if 'real*8 reentrant_hg(1)' in data:
        previous=new.replace('elcon(1,1,imat)', 'elconloc(1)')
        applied=new if new in data else previous if previous in data else None
        if applied is None:raise SystemExit('Unrecognized hourglass patch: '+name)
        original=data.replace('      real*8 reentrant_hg(1)\n','').replace(applied,old)
    if hashlib.sha256(original.encode()).hexdigest()!=expected:
        raise SystemExit('Unexpected upstream source checksum: '+name)
    patched=original.replace('      implicit none\n','      implicit none\n      real*8 reentrant_hg(1)\n',1).replace(old,new)
    if patched!=data:path.write_text(patched)
