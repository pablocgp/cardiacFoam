#!/usr/bin/env python3
"""Rollback audit over every result tree that feeds the paper.

Why the run summaries and not the nonlinearResiduals histories: the solver
writes "non-converged steps" and "rolled-back steps" into every
*_transient.dat itself, so the count does not depend on this script
reconstructing it from per-iteration traces. It is also two orders of magnitude
cheaper - the 2-D tree alone holds 151k residual files.

A run with a rolled-back step has a frozen field for that step, so its fitted
order means nothing. The audit exists to prove that no run behind a published
number had one.
"""
from __future__ import annotations
import sys
from pathlib import Path

TREES = {
    "2D nuevo (Distributed, np4 scotch)":
        "/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro"
        "/highOrderManufacturedFDAImplicitPETScDistributed/results_np4_threads1_scotch",
    "3D nuevo (barrido del paper)":
        "/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/mms3D_sweep_paper",
    "3D viejo (solver PETSc no distribuido)":
        "/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro"
        "/highOrderManufacturedFDAImplicitPETSc/results",
    "3D Tests3D (barridos de escalado)":
        "/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro"
        "/highOrderManufacturedFDAImplicitPETScDistributed_Tests3D",
}

def audit(root: Path):
    runs = steps = nonconv = rolled = 0
    sin_clave = []
    culpables = []
    for p in root.rglob("*_transient.dat"):
        txt = p.read_text(errors="replace")
        d = {}
        for line in txt.splitlines():
            k, sep, v = line.partition("=")
            if sep:
                d[k.strip()] = v.strip()
        if "non-converged steps" not in d:
            sin_clave.append(p)
            continue
        runs += 1
        steps += int(float(d.get("Number of steps", 0)))
        nc = int(float(d["non-converged steps"]))
        rb = int(float(d.get("rolled-back steps", 0)))
        nonconv += nc
        rolled += rb
        if nc or rb:
            culpables.append((p, nc, rb))
    return runs, steps, nonconv, rolled, sin_clave, culpables

def main() -> int:
    total_mal = 0
    print(f"{'arbol':42s} {'corridas':>9s} {'pasos':>10s} {'no conv':>8s} {'rollback':>9s} {'sin clave':>10s}")
    for nombre, ruta in TREES.items():
        root = Path(ruta)
        if not root.is_dir():
            print(f"{nombre:42s} {'(no existe)':>9s}")
            continue
        runs, steps, nc, rb, sin, culp = audit(root)
        total_mal += nc + rb
        print(f"{nombre:42s} {runs:9d} {steps:10d} {nc:8d} {rb:9d} {len(sin):10d}")
        for p, n, r in culp[:10]:
            print(f"    !! {p}  no_conv={n} rollback={r}")
    print()
    print("VEREDICTO:", "limpio, ningun paso no convergido ni revertido"
          if total_mal == 0 else f"**{total_mal} eventos**, ver arriba")
    return 0 if total_mal == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
