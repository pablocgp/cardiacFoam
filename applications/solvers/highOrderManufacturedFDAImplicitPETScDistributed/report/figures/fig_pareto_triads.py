#!/usr/bin/env python3
"""Accuracy-cost Pareto front over the (Vm, states, Iion) triad space.

Reads the harvested table pareto_2d_picard_dt1e-3.csv (one row per triad x N x
mesh x alpha, from results_np4_threads1_scotch) and draws a 2x2 grid of panels,
one per (mesh family, stabilisation alpha).

Form: emphasis. Every run is a faint gray dot; three families are highlighted
as their mesh-refinement curves; the Pareto front (non-dominated runs in cost
and error) is drawn in neutral ink, not a series color, because it is not a
series. Three highlighted families is the all-pairs cap for a scatter form.

Palette: categorical slots 1-3 of the reference palette. The validator could
not be run on this machine (no JavaScript runtime); palette.md documents these
three slots as passing all-pairs in light mode (worst CVD dE 9.2, normal-vision
dE 24.0).
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE / "pareto_2d_picard_dt1e-3.csv"

SURFACE = "#ffffff"          # print figure: plain white page
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
CONTEXT = "#c3c2b7"          # de-emphasised runs

# Fixed identity -> color. Never re-assigned by rank.
HIGHLIGHT = [
    ("p3/CCp2/p2", "#2a78d6", r"$(p3,\ \mathrm{CC}p2,\ p2)$  cheapest high-order"),
    ("p3/CCp3/p3", "#eb6834", r"$(p3,\ \mathrm{CC}p3,\ p3)$  diagonal"),
    ("NO/na/NO",   "#1baf7a", r"$(\mathrm{NO},\ -,\ \mathrm{NO})$  2nd-order FV"),
]

PANELS = [
    ("hexa", "0p0e00", r"hexahedral, $\alpha=0$"),
    ("hexa", "1p0em01", r"hexahedral, $\alpha=0.1$"),
    ("triangular_Unstr", "0p0e00", r"unstructured triangular, $\alpha=0$"),
    ("triangular_Unstr", "1p0em01", r"unstructured triangular, $\alpha=0.1$"),
]


def load():
    groups = defaultdict(list)
    with open(DATA) as fh:
        for r in csv.DictReader(fh):
            cost = float(r["timeLoop_s"]) if r["timeLoop_s"] else 0.0
            err = float(r["Vm_L2"])
            if cost <= 0.0 or err <= 0.0:
                continue
            groups[(r["mesh"], r["alpha"])].append(dict(
                triad=f'{r["vm"]}/{r["states"]}/{r["iion"]}',
                N=int(r["N"]), cost=cost, err=err))
    return groups


def pareto_front(points):
    pts = sorted(points, key=lambda p: (p["cost"], p["err"]))
    front, best = [], float("inf")
    for p in pts:
        if p["err"] < best:
            front.append(p)
            best = p["err"]
    return front


def style_axes(ax):
    ax.set_facecolor(SURFACE)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="major", color=GRID, linewidth=0.6, linestyle="-")
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(BASELINE)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelcolor=INK_SECONDARY, labelsize=8, length=3)


def main():
    groups = load()
    if not groups:
        print(f"no data in {DATA}")
        return 1

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Liberation Sans"],
        "mathtext.fontset": "dejavusans",
    })

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), sharex=True, sharey=True)
    fig.patch.set_facecolor(SURFACE)

    for ax, (mesh, alpha, title) in zip(axes.flat, PANELS):
        pts = groups[(mesh, alpha)]
        style_axes(ax)

        ax.scatter([p["cost"] for p in pts], [p["err"] for p in pts],
                   s=3, color=CONTEXT, alpha=0.35, linewidths=0, zorder=1,
                   rasterized=True)

        front = pareto_front(pts)
        ax.step([p["cost"] for p in front], [p["err"] for p in front],
                where="post", color=INK, linewidth=1.0, zorder=3)

        for triad, color, _ in HIGHLIGHT:
            fam = sorted((p for p in pts if p["triad"] == triad), key=lambda p: p["N"])
            if not fam:
                continue
            # Every 10th mesh level, so the refinement curve reads as a curve and
            # not as a band of overlapping markers.
            fam = fam[::10] + ([fam[-1]] if (len(fam) - 1) % 10 else [])
            ax.plot([p["cost"] for p in fam], [p["err"] for p in fam],
                    color=color, linewidth=1.5, solid_capstyle="round", zorder=4)
            ax.scatter([p["cost"] for p in fam], [p["err"] for p in fam],
                       s=22, color=color, edgecolors=SURFACE, linewidths=1.2, zorder=5)

        ax.set_title(title, fontsize=9, color=INK, loc="left", pad=4)

    for ax in axes[1, :]:
        ax.set_xlabel("time-loop wall time per run [s]", fontsize=8, color=INK_SECONDARY)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$V_m$ error, $L_2$ norm", fontsize=8, color=INK_SECONDARY)

    handles = [plt.Line2D([], [], color=c, linewidth=1.5, marker="o", markersize=5,
                          markeredgecolor=SURFACE, label=lab) for _, c, lab in HIGHLIGHT]
    handles.append(plt.Line2D([], [], color=INK, linewidth=1.0, label="Pareto front (all 52 triads)"))
    handles.append(plt.Line2D([], [], color=CONTEXT, marker="o", linestyle="none",
                              markersize=3, label="every run"))
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=8, labelcolor=INK_SECONDARY, bbox_to_anchor=(0.5, 0.0))

    fig.tight_layout(rect=(0, 0.11, 1, 1))

    for ext in ("pdf", "png"):
        out = HERE / f"fig_pareto_triads.{ext}"
        fig.savefig(out, dpi=200, facecolor=SURFACE)
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
