"""Convert CalculiX FRD to VTU; optionally skip derived principal tensors.

Usage:
    python convert_frd.py
    python convert_frd.py --frd results/dong2019_contact.frd --skip-principal

--skip-principal retains original FRD displacement, stress, and strain arrays.
It bypasses only ccx2paraview's derived principal-value calculation.
"""
import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
<<<<<<< HEAD
=======
FRD = ROOT / "results" / "dong2019_contact.frd"
OUT = ROOT / "results" / "vtk"
>>>>>>> 2f5d2f9 (Update CalculiX FEA)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--frd", type=Path, default=Path("results/reentrant.frd"))
    parser.add_argument("--out", type=Path, default=Path("results/vtk"))
    parser.add_argument("--skip-principal", action="store_true")
    args = parser.parse_args()
    frd = args.frd if args.frd.is_absolute() else ROOT / args.frd
    out = args.out if args.out.is_absolute() else ROOT / args.out
    if not frd.is_file():
        raise FileNotFoundError(f"FRD file not found: {frd}")
    out.mkdir(parents=True, exist_ok=True)

    # The converter calls calculate_principal for S/E even when only U is
    # needed. Skip derived principal arrays while retaining original fields.
    # Return original block rather than a fabricated stress tensor.
    script = """
import sys
import ccx2paraview.common as common
import ccx2paraview.cli as cli
common.FRD.calculate_principal = lambda self, block: block
sys.argv = ["ccx2paraview", sys.argv[1], "vtu"]
cli.main()
"""
    def invoke(skip):
        cmd = ([sys.executable, "-c", script, str(frd)] if skip
               else [sys.executable, "-m", "ccx2paraview", str(frd), "vtu"])
        return subprocess.run(cmd, cwd=out, capture_output=True, text=True)

    result = invoke(args.skip_principal)
    if result.returncode and not args.skip_principal and (
        "not supported between instances of 'complex'" in result.stderr
        or "FRD tensor is not symmetric" in result.stderr
    ):
        print("Principal-value conversion failed; retrying without derived principal arrays.")
        result = invoke(True)

    if result.stdout:
        print(result.stdout)
    if result.returncode:
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(
            "FRD conversion failed; see converter traceback and verify FRD/solver integrity."
        )
    files = sorted(out.glob("*.vtu"))
    if not files:
        raise RuntimeError("No VTU files were generated")
    print(f"Generated {len(files)} VTU files in {out}")
    print("Warning: VTU conversion success does not establish solver convergence.")

if __name__ == "__main__":
    main()
