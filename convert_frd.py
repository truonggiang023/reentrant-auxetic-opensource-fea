"""Convert CalculiX FRD to VTU, with guarded symmetric-tensor eigensolver fallback.

The fallback only accepts finite, nearly symmetric 3x3 tensors. It does not
silently discard imaginary eigenvalue components from invalid results.
"""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
FRD = ROOT / "results" / "reentrant.frd"
OUT = ROOT / "results" / "vtk"

if not FRD.is_file():
    raise FileNotFoundError(f"FRD file not found: {FRD}")

OUT.mkdir(parents=True, exist_ok=True)

def convert(guarded=False):
    if not guarded:
        cmd = [sys.executable, "-m", "ccx2paraview", str(FRD), "vtu"]
    else:
        # ccx2paraview uses eigvals + sorted() for principal stresses.
        # For real symmetric stress tensors, eigvalsh is the stable algorithm.
        code = """
import sys
import numpy as np
import ccx2paraview.cli as cli
_original = np.linalg.eigvals
def checked_eigvals(a):
    a = np.asarray(a)
    if a.shape == (3, 3):
        if not np.isfinite(a).all():
            raise ValueError("Nonfinite FRD tensor: cannot calculate principal values")
        scale = max(1.0, float(np.max(np.abs(a))))
        asymmetry = float(np.max(np.abs(a - a.T)))
        if asymmetry > 1e-6 * scale:
            raise ValueError(
                f"FRD tensor is not symmetric (asymmetry={asymmetry:g}); "
                "inspect solver results, no eigenvalue coercion performed"
            )
        return np.linalg.eigvalsh(0.5 * (a + a.T))
    return _original(a)
np.linalg.eigvals = checked_eigvals
sys.argv = ["ccx2paraview", sys.argv[1], "vtu"]
cli.main()
"""
        cmd = [sys.executable, "-c", code, str(FRD)]
    return subprocess.run(cmd, cwd=OUT, capture_output=True, text=True)

print("Converting:", FRD)
result = convert()
if result.returncode and "not supported between instances of 'complex'" in result.stderr:
    print("Retrying with guarded symmetric-tensor principal-value calculation.")
    result = convert(guarded=True)

print(result.stdout)
if result.returncode:
    print(result.stderr, file=sys.stderr)
    raise RuntimeError(
        "FRD conversion failed. Inspect the solver .sta/.dat files and "
        "FRD tensor validity; do not treat a diverged run as a valid result."
    )

files = sorted(OUT.glob("*.vtu"))
if not files:
    raise RuntimeError("Conversion returned success but produced no VTU files")
print(f"Conversion completed. Found {len(files)} VTU files.")
for file in files[:10]:
    print(file.name)
