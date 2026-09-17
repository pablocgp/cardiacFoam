#!/usr/bin/env python3
"""Spatial-order figures for the paper, 2-D and 3-D, from the raw .dat files.

Two figures, each answering one claim:

  fig_order_vm      error of Vm against mesh size, for the four reconstruction
                    levels on the diagonal (NO, p1, p2, p3), in four panels:
                    2-D hexa, 2-D triangular, 3-D hexa, 3-D tetrahedral. This
                    is the odd/even statement and the loss of consistency of
                    the standard operator on simplices.

  fig_cc_vs_gp      error of the STATES at cell centres, cellCentredReconstruct
                    against gaussPointODE, same four panels. GP is pinned at
                    order 2 while CC follows the scheme; Vm converges
                    identically either way, which is drawn as context.

Both read the raw errors_vs_N_*.dat trees rather than any generated table, so
the numbers in the caption can be traced to a run directory.

Colour, and why it is not four hues: the four reconstruction levels are an
ORDERED scale (NO < p1 < p2 < p3), so they get one hue light-to-dark, which is
the rule for ordered categories. The CC/GP pair is genuinely categorical and
uses the two slots already validated for the Pareto figure. The palette
validator could not be run here - no JavaScript runtime on this machine - which
is a further reason not to invent a fourth categorical hue by eye.
"""

import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
RUN = Path("/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro")
TREE_2D = RUN / "highOrderManufacturedFDAImplicitPETScDistributed/results_np4_threads1_scotch/example0/2D_MMS_Picard"
TREE_3D = RUN / "mms3D_sweep_paper/barrido/results_np4_scotch/example0"

SURFACE = "#ffffff"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
GUIDE = "#c3c2b7"

# Ordered ramp, one hue, light -> dark, for the four reconstruction levels.
RAMP = ["#b9d4f2", "#7fb0e6", "#3f83d0", "#14416f"]
# Categorical pair, the two slots used by the Pareto figure.
BLUE, ORANGE = "#2a78d6", "#eb6834"

NIVELES = [
    ("NO", "hoVm_NO_states_na_hoIion_NO", r"$\mathrm{NO}$"),
    ("p1", "hoVm_p1_states_CCp1_hoIion_p1", r"$p1$"),
    ("p2", "hoVm_p2_states_CCp2_hoIion_p2", r"$p2$"),
    ("p3", "hoVm_p3_states_CCp3_hoIion_p3", r"$p3$"),
]


def leer(dat: Path, columna: str):
    """[(N, valor)] de un errors_vs_N_*.dat, ordenado por N."""
    lineas = dat.read_text().strip().split("\n")
    h = lineas[0].lstrip("# ").split()
    iN, iv = h.index("N"), h.index(columna)
    filas = []
    for linea in lineas[1:]:
        c = linea.split()
        filas.append((float(c[iN]), float(c[iv])))
    return sorted(filas)


_INDICE: dict[Path, list[Path]] = {}


def buscar(raiz: Path, malla: str, config: str):
    """El errors_vs_N de una configuración, o None si no existe.

    El árbol 2-D tiene ~151k archivos, así que el listado se hace UNA vez por
    raíz y se reutiliza; recorrerlo por configuración costaba 16 pasadas.
    """
    if raiz not in _INDICE:
        _INDICE[raiz] = list(raiz.rglob("errors_vs_N_*.dat"))
    for p in _INDICE[raiz]:
        s = str(p)
        if f"_mesh_{malla}_" in s and f"/{config}_" in s:
            return p
    return None


def pendiente(puntos):
    """Ajuste por mínimos cuadrados de log(error) contra log(h), h = 1/N."""
    import math
    xs = [math.log(1.0 / n) for n, _ in puntos]
    ys = [math.log(v) for _, v in puntos]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    return num / den if den else float("nan")


def submuestrear(puntos, maximo=12):
    if len(puntos) <= maximo:
        return puntos
    paso = len(puntos) // maximo
    fuera = puntos[::paso]
    if fuera[-1] != puntos[-1]:
        fuera.append(puntos[-1])
    return fuera


def tabla_pendientes(ax, entradas):
    """Las pendientes ajustadas, como tabla en el panel.

    No van pegadas al extremo de cada curva porque ahí se pisan: en 2-D hexa
    NO, p1 y p2 terminan a la misma altura, y en 3-D hexa dos series ajustan
    1.98. Una columna alineada se lee de una pasada y no colisiona nunca.
    """
    # Abajo a la derecha: con el error creciendo hacia h grande, ese cuadrante
    # queda libre en los cuatro paneles. El fondo es el de la superficie, sin
    # borde, para que el texto gane si alguna curva pasa por ahi.
    n = len(entradas)
    for i, (etiqueta, valor, color) in enumerate(entradas):
        ax.annotate(f"{etiqueta}  {valor:.2f}",
                    xy=(0.97, 0.03 + 0.075 * (n - 1 - i)),
                    xycoords="axes fraction",
                    fontsize=7.5, color=color, va="bottom", ha="right",
                    bbox=dict(facecolor=SURFACE, edgecolor="none",
                              boxstyle="square,pad=0.15", alpha=0.85))


def estilo(ax, titulo):
    ax.set_facecolor(SURFACE)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="major", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    for lado in ("left", "bottom"):
        ax.spines[lado].set_color(BASELINE)
        ax.spines[lado].set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelcolor=INK_SECONDARY, labelsize=8, length=3)
    ax.set_title(titulo, fontsize=9, color=INK, loc="left", pad=4)


PANELES = [
    (TREE_2D, "hexa", "2-D hexahedral"),
    (TREE_2D, "triangular_Unstr", "2-D unstructured triangular"),
    (TREE_3D, "hexa", "3-D hexahedral"),
    (TREE_3D, "triangular_Unstr", "3-D unstructured tetrahedral"),
]


def figura_orden_vm():
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2))
    fig.patch.set_facecolor(SURFACE)
    manejadores = {}

    for ax, (raiz, malla, titulo) in zip(axes.flat, PANELES):
        estilo(ax, titulo)
        entradas = []
        for (clave, config, etiqueta), color in zip(NIVELES, RAMP):
            dat = buscar(raiz, malla, config)
            if not dat:
                continue
            pts = leer(dat, "Vm_L2")
            p = pendiente(pts)
            m = submuestrear(pts)
            hs = [1.0 / n for n, _ in m]
            vs = [v for _, v in m]
            (ln,) = ax.plot(hs, vs, color=color, linewidth=1.6, marker="o",
                            markersize=4, markeredgecolor=SURFACE,
                            markeredgewidth=0.8)
            manejadores[etiqueta] = ln
            entradas.append((etiqueta, p, color))
        tabla_pendientes(ax, entradas)
        ax.set_xlabel(r"$h = 1/N$", fontsize=8, color=INK_SECONDARY)
        ax.set_ylabel(r"$V_m$ error, $L_2$ norm", fontsize=8, color=INK_SECONDARY)

    fig.legend(handles=list(manejadores.values()), labels=list(manejadores),
               loc="lower center", ncol=4, frameon=False, fontsize=8,
               labelcolor=INK_SECONDARY, bbox_to_anchor=(0.5, 0.0),
               title=r"$V_m$ reconstruction (fitted order shown per panel)",
               title_fontsize=8)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"fig_order_vm.{ext}", dpi=200, facecolor=SURFACE)
    print("escrito fig_order_vm.pdf/.png")


def figura_cc_vs_gp():
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2))
    fig.patch.set_facecolor(SURFACE)
    manejadores = {}

    for ax, (raiz, malla, titulo) in zip(axes.flat, PANELES):
        estilo(ax, titulo)
        entradas = []
        for config, color, etiqueta in (
            ("hoVm_p3_states_CCp3_hoIion_p3", BLUE, "states at cell centres (CC)"),
            ("hoVm_p3_states_GP_hoIion_p3", ORANGE, "states at Gauss points (GP)"),
        ):
            dat = buscar(raiz, malla, config)
            if not dat:
                continue
            pts = leer(dat, "u1_L2")
            p = pendiente(pts)
            m = submuestrear(pts)
            hs = [1.0 / n for n, _ in m]
            vs = [v for _, v in m]
            (ln,) = ax.plot(hs, vs, color=color, linewidth=1.6, marker="o",
                            markersize=4, markeredgecolor=SURFACE,
                            markeredgewidth=0.8)
            manejadores[etiqueta] = ln
            entradas.append(("CC" if color == BLUE else "GP", p, color))
        tabla_pendientes(ax, entradas)
        ax.set_xlabel(r"$h = 1/N$", fontsize=8, color=INK_SECONDARY)
        ax.set_ylabel(r"$u_1$ error at cell centres, $L_2$ norm", fontsize=8,
                      color=INK_SECONDARY)

    fig.legend(handles=list(manejadores.values()), labels=list(manejadores),
               loc="lower center", ncol=2, frameon=False, fontsize=8,
               labelcolor=INK_SECONDARY, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"fig_cc_vs_gp.{ext}", dpi=200, facecolor=SURFACE)
    print("escrito fig_cc_vs_gp.pdf/.png")


def main():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Liberation Sans"],
        "mathtext.fontset": "dejavusans",
    })
    faltan = [str(t) for t in (TREE_2D, TREE_3D) if not t.is_dir()]
    if faltan:
        print("faltan árboles de resultados:", ", ".join(faltan))
        return 1
    figura_orden_vm()
    figura_cc_vs_gp()
    return 0


if __name__ == "__main__":
    sys.exit(main())
