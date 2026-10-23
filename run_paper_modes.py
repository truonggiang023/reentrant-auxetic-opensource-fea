"""Generate the paper-mode deck and solve it using the user's local CalculiX."""
from pathlib import Path
import argparse
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ccx", default=str(ROOT / "tools/calculix/bin/ccx/ccx213.exe"))
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--self-contact", action="store_true",
                   help="EXPERIMENTAL; same-surface self contact may not work in 2.13")
    p.add_argument("--generate-only", action="store_true")
    a = p.parse_args()
    cmd = [
        sys.executable, "-m", "reentrant_fea.cli",
        "--paper-modes",
        "--config", "configs/dong2019_paper_modes.json",
        "--plastic", "reference_setup_numeric_source.py",
        "--output", "results/reentrant.inp",
    ]
    if a.self_contact:
        cmd.append("--self-contact")
    subprocess.run(cmd, cwd=ROOT, check=True)
    if a.generate_only:
        return
    ccx = Path(a.ccx).resolve()
    if not ccx.is_file():
        raise FileNotFoundError(f"CalculiX executable is missing: {ccx}")
    env = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "CCX_NPROC_STIFFNESS",
                 "CCX_NPROC_EQUATION_SOLVER", "CCX_NPROC_RESULTS"):
        env[name] = str(a.threads)
    print(f"[CCX] {ccx} | requested threads={a.threads}", flush=True)
    subprocess.run([str(ccx), "-i", "reentrant"], cwd=ROOT / "results",
                   env=env, check=True)
    print("[END] Inspect results/reentrant.sta, .dat, .frd; successful exit is not validation.")


if __name__ == "__main__":
    main()
