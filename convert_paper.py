"""Convert contact simulation FRD to VTU for ParaView (requires ccx2paraview)."""
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
frd=ROOT/"results"/"dong2019_contact.frd"
out=ROOT/"results"/"dong2019_contact_vtu"
if not frd.is_file():
    raise FileNotFoundError(frd)
out.mkdir(exist_ok=True,parents=True)
subprocess.run([sys.executable,"-m","ccx2paraview",str(frd),"vtu"],
               cwd=out,check=True)
files=sorted(out.glob("*.vtu"))
if not files:
    raise RuntimeError("Converter exited without VTU output")
print("Open in ParaView:",out)
print("Use Warp By Vector on displacement U and compare nominal strain snapshots")
print("Reference snapshots: 0.12, 0.35, 0.65; not guaranteed mode matches")
