
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
FRD = ROOT / "results" / "reentrant.frd"
OUT = ROOT / "results" / "vtk"

if not FRD.is_file():
    raise FileNotFoundError(f"FRD file not found: {FRD}")

OUT.mkdir(parents=True, exist_ok=True)

cmd = [
    sys.executable,
    "-m",
    "ccx2paraview",
    str(FRD),
    "vtu",
]

print("Converting:", FRD)

result = subprocess.run(
    cmd,
    cwd=OUT,
    capture_output=True,
    text=True,
)

print(result.stdout)

if result.returncode != 0:
    print(result.stderr)
    raise RuntimeError("FRD conversion failed")

files = sorted(OUT.glob("*.vtu"))

print(f"Conversion completed. Found {len(files)} VTU files.")

for file in files[:10]:
    print(file.name)
