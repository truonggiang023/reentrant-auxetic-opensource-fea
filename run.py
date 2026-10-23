
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent

CCX = (
    ROOT
    / "tools"
    / "calculix"
    / "bin"
    / "ccx"
    / "ccx213.exe"
)

RESULTS = ROOT / "results"
INPUT = RESULTS / "reentrant.inp"

if not CCX.is_file():
    raise FileNotFoundError(f"CalculiX not found: {CCX}")

if not INPUT.is_file():
    raise FileNotFoundError(f"FEA input not found: {INPUT}")

print(f"Using CalculiX: {CCX.name}")
print(f"Input: {INPUT}")

subprocess.run(
    [str(CCX), "-i", "reentrant"],
    cwd=RESULTS,
    check=True,
)

print("CalculiX process finished. Check .sta and .dat for solver status.")
