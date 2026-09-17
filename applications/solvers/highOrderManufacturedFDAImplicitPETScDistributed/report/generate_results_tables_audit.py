#!/usr/bin/env python3
"""Generate the audited LaTeX result tables (English and Spanish) for the
PETSc MMS report.

This script reuses the order/cost extraction of ``generate_results_tables.py``
and adds two nonlinear-robustness columns read from the per-run
``nonlinearResiduals_*.dat`` histories:

  RB      number of time steps that exhausted the nonlinear iteration budget,
          accumulated over every mesh size of the row.  For ``diagonalIion``
          and ``JFNK`` (``nonlinearAcceptUnconverged=false``) such a step is
          rolled back, i.e. Vm, Iion and the states are restored to their
          start-of-step values while the clock still advances.  For ``Picard``
          the same event means "last iterate accepted with a warning".

  n_nl    median/maximum number of nonlinear iterations per time step, again
          accumulated over every mesh size of the row.  This is the margin to
          the rollback condition: the budget is ``nonlinearIterations=2000``.

Files are written with an ``_en`` or ``_es`` suffix; the original
``generate_results_tables.py`` and the unsuffixed tables it produces are left
untouched.  Re-run this script whenever the result tree changes: it rescans
every ``nonlinearResiduals`` history, which takes about a minute.
"""

from __future__ import annotations

import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path

import generate_results_tables as G


REPORT_DIR = Path(__file__).resolve().parent
TABLE_DIR = REPORT_DIR / "tables"

MESH_RE = re.compile(r"^(?P<dim>[23]D)_mesh_(?P<mesh>.+)_alpha_(?P<alpha>.+)$")
CONFIG_RE = re.compile(
    r"^hoVm_(?P<vm>[^_]+)_states_(?P<states>[^_]+)_hoIion_(?P<iion>[^_]+)"
    r"_(?P<scheme>[^_]+)_(?P<mass>[^_]+)_(?P<method>[^_]+)_(?P<ode>[^_]+)$"
)

# Column index of the residual history files (header is a single '#' line):
#   0 time_s          1 step             2 nonlinearMethod  3 stateODESolver
#   4 iter            5 linearIterations 6 linearError      7 newtonResidual
#   8 Vm_relL2        9 u1_relL2        10 u2_relL2        11 u3_relL2
#  12 Iion_relL2     13 lineSearchIters 14 converged
COL_STEP = 1
COL_ITER = 4
COL_LS = 13
COL_CONVERGED = 14
N_COLS = 15


# Lines of the per-run summary that carry the AUTHORITATIVE step accounting.
# The solver computes them from its per-step record, which -- unlike the
# per-iteration residual trace -- covers every step of the run.
SUMMARY_STEPS = "Number of steps"
SUMMARY_NONCONV = "non-converged steps"
SUMMARY_ROLLED = "rolled-back steps"
SUMMARY_MAXIT = "max it used"


def read_run_summary(directory, dim: str, n_cells: str, dt_tag: str) -> dict:
    """The scalars the solver wrote for one run, or {} if there is no summary.

    The file sits next to the residual history it belongs to and is named
    ``{dim}_{N}_cells_dt_{dt}_transient.dat``.
    """
    path = directory / f"{dim}_{n_cells}_cells_dt_{dt_tag}_transient.dat"
    if not path.is_file():
        return {}

    wanted = (SUMMARY_STEPS, SUMMARY_NONCONV, SUMMARY_ROLLED, SUMMARY_MAXIT)
    out: dict[str, int] = {}
    for line in path.read_text().splitlines():
        if "=" not in line:
            continue
        label, _, value = line.partition("=")
        label = label.strip()
        if label in wanted:
            try:
                out[label] = int(float(value.strip()))
            except ValueError:
                pass
    return out


def expected_startup_steps(scheme_tag: str) -> int:
    """Steps a BDF run spends in its ESDIRK startup, from the directory name.

    The directory records the scheme as ``BDF3`` or, once a startup is in play,
    as ``BDF3-ESDIRK3``.  Only the ESDIRK startup takes steps: ``exact`` and
    ``constant`` fabricate the history before the time loop begins, so they
    contribute nothing here.  A BDF-k startup runs k-1 steps.
    """
    scheme, _, startup = scheme_tag.partition("-")
    if not startup.lower().startswith("esdirk"):
        return 0
    match = re.fullmatch(r"(?i)bdf(\d+)", scheme)
    return max(int(match.group(1)) - 1, 0) if match else 0


def scan_nonlinear_histories() -> dict[tuple[str, ...], dict]:
    """Walk every nonlinearResiduals file and aggregate per table row.

    The iteration DISTRIBUTION comes from the residual trace, which is the only
    place it exists.  The step COUNT and the rollback count do not: they are
    read from the per-run summary instead.

    That split is not cosmetic.  The ESDIRK startup writes no rows to the
    residual trace -- its stages run their own fixed-point iteration and never
    reach the residual writer -- so a BDF-k run is missing its first k-1 steps
    from that file.  Deriving RB by counting groups of rows would therefore
    report zero rollbacks for a startup step that actually exhausted its
    budget: a clean bill of health that was never checked.  The summary's
    "non-converged steps" is computed by the solver from its per-step record,
    which covers the startup, so it is the number that answers the question RB
    asks.  Runs predating the BDF work carry the same line, so nothing has to
    be regenerated.

    The missing rows are still counted, and a discrepancy larger than the k-1
    the startup explains is reported rather than absorbed.
    """
    stats: dict[tuple[str, ...], dict] = defaultdict(
        lambda: {
            "runs": 0,
            "steps": 0,
            "rollbacks": 0,
            "iters": [],
            "iter_max": 0,
            "linesearch": 0,
            "meshes": set(),
            # Runs whose summary could not be read: their RB is unknown rather
            # than zero, and saying so is the point.
            "no_summary": 0,
            # Steps present in the summary and absent from the trace, beyond
            # the k-1 the ESDIRK startup accounts for.
            "unexplained_gap": 0,
        }
    )

    for path in G.RESULT_ROOT.rglob("nonlinearResiduals_*.dat"):
        cfg_dir = path.parent.parent
        mesh_dir = cfg_dir.parent
        mesh_match = MESH_RE.match(mesh_dir.name)
        cfg_match = CONFIG_RE.match(cfg_dir.name)
        if not mesh_match or not cfg_match:
            continue

        dt_tag = path.parent.name.split("_dt_")[-1]
        n_match = re.search(r"_N(\d+)_", path.name)
        if not n_match:
            continue

        key = (
            mesh_match["dim"],
            cfg_match["method"],
            dt_tag,
            mesh_match["mesh"],
            mesh_match["alpha"],
            cfg_match["vm"],
            cfg_match["states"],
            cfg_match["iion"],
        )
        entry = stats[key]
        entry["runs"] += 1
        entry["meshes"].add(int(n_match.group(1)))

        # --- Iteration distribution, from the trace: the only source there is.
        step = None
        iters = 0
        trace_steps = 0
        with path.open() as handle:
            next(handle, None)
            for line in handle:
                parts = line.split()
                if len(parts) < N_COLS:
                    continue
                this_step = int(parts[COL_STEP])
                if this_step != step:
                    if step is not None:
                        trace_steps += 1
                        entry["iters"].append(iters)
                        entry["iter_max"] = max(entry["iter_max"], iters)
                    step = this_step
                    iters = 0
                iters = max(iters, int(parts[COL_ITER]))
                entry["linesearch"] += max(int(parts[COL_LS]), 0)

        if step is not None:
            trace_steps += 1
            entry["iters"].append(iters)
            entry["iter_max"] = max(entry["iter_max"], iters)

        # --- Step and rollback counts, from the summary: authoritative.
        summary = read_run_summary(
            path.parent, mesh_match["dim"], n_match.group(1), dt_tag
        )

        if not summary:
            # Without a summary the rollback count for this run is a lower
            # bound and not a measurement. Recorded as such instead of being
            # quietly treated as zero.
            entry["no_summary"] += 1
            entry["steps"] += trace_steps
            continue

        summary_steps = summary.get(SUMMARY_STEPS, trace_steps)
        entry["steps"] += summary_steps
        entry["rollbacks"] += summary.get(SUMMARY_NONCONV, 0)
        entry["iter_max"] = max(entry["iter_max"], summary.get(SUMMARY_MAXIT, 0))

        # Cross-check. The trace should be short by EXACTLY the startup steps.
        # Any other discrepancy means the two records disagree for a reason
        # nobody has accounted for, which is worth surfacing rather than
        # averaging away.
        gap = summary_steps - trace_steps
        entry["unexplained_gap"] += abs(
            gap - expected_startup_steps(cfg_match["scheme"])
        )

    return stats


def robustness_cells(entry: dict | None) -> tuple[str, str]:
    if not entry or not entry["iters"]:
        return "--", "--"

    median = statistics.median(entry["iters"])

    # RB is only a measurement where every run of the row had a summary to read
    # it from. A row missing one is marked rather than reported as a number:
    # "0 rollbacks" and "0 rollbacks observed among the runs we could check"
    # are different claims, and this table is cited as the first.
    rb = str(entry["rollbacks"])
    if entry.get("no_summary"):
        rb = f"{rb}$^{{\\dagger}}$"
    if entry.get("unexplained_gap"):
        rb = f"{rb}$^{{\\ddagger}}$"

    return rb, f"{median:.2f}/{entry['iter_max']}"


def dt_to_latex(tag: str) -> str:
    return G.dt_to_latex(tag)


def range_label(dim: str, key: str) -> str:
    n_min, n_max = G.N_RANGES[dim][key]
    return f"$N={n_min}$--${n_max}$"


# ----------------------------------------------------------------------- #
# Per-language strings.  Everything that differs between the English and the
# Spanish report lives here; the emission code below is shared.
# ----------------------------------------------------------------------- #

LANGS = {
    "en": {
        "mesh": lambda dim, mesh: (
            ("tetra unstr" if dim == "3D" else "triangular unstr")
            if mesh == "triangular_Unstr"
            else mesh.replace("_", " ")
        ),
        "fit": "Fit",
        "meshcol": "Mesh type",
        "caption": (
            "{dim} manufactured results for {method}, \\code{{crankNicolson}}, "
            "\\code{{consistent}} mass, \\code{{RKF45}} states, and "
            "$\\Delta t={dt}$. The first convergence block is fitted over "
            "{all_lbl}; the second over {fine_lbl}. The first column records "
            "the actual \\code{{Vm-states-Iion}} triplet used by the run. "
            "\\textbf{{RB}} is the number of time steps that exhausted the "
            "nonlinear iteration budget, accumulated over the {n_meshes} mesh "
            "sizes of the row ({n_steps} time steps per row in total). {note} "
            "\\textbf{{$n_{{\\mathrm{{nl}}}}$}} reports the median and the "
            "maximum number of nonlinear iterations per step over those same "
            "steps; the configured budget is "
            "\\code{{nonlinearIterations}}$=2000$."
        ),
        "note_picard": (
            "\\code{Picard} accepts its last iterate with a warning instead "
            "of rolling back, so for this method the column counts steps "
            "accepted without converging."
        ),
        "note_other": (
            "With \\code{nonlinearAcceptUnconverged=false} an unconverged "
            "step is reverted to its start-of-step values while the clock "
            "still advances."
        ),
        "empty": "No aggregated error file found for this case.",
        "subsection": lambda dim: (
            f"{'Two' if dim == '2D' else 'Three'}-dimensional manufactured problem"
        ),
        "intro": [
            "The tables in this subsection are generated directly from the",
            r"current \path{{{root}}} result tree. The method-specific roots are {roots}.",
            r"All entries use \code{{crankNicolson}}, the \code{{consistent}} mass matrix, and \code{{RKF45}} for the state ODEs.",
            r"The first column reports \code{{Vm-states-Iion}}: the $V_m$ reconstruction degree, the state treatment (\code{{CCp*}}, \code{{GP}}, or \code{{na}}), and the $\Iion$ quadrature/reconstruction setting.",
            "Fits use {all_lbl} for the broad mesh range and {fine_lbl} for the fine-mesh range.",
            r"The last two columns audit the nonlinear solver: \code{{RB}} counts rolled-back steps and $n_{{\mathrm{{nl}}}}$ summarises the nonlinear iterations per step (median/maximum).",
            "",
        ],
        "audit_caption": (
            r"Global nonlinear-convergence audit over the complete "
            r"\code{example0} sweep. \code{RB} counts the time steps that "
            r"exhausted the \code{nonlinearIterations}$=2000$ budget (and that "
            r"would therefore have been rolled back for "
            r"\code{diagonalIion}/\code{JFNK}, or accepted with a warning for "
            r"\code{Picard}). $\bar n_{\mathrm{nl}}$ and $n^{\max}_{\mathrm{nl}}$ "
            r"are the mean and maximum number of nonlinear iterations per step; "
            r"$\Sigma_{\mathrm{LS}}$ is the total number of Armijo line-search "
            r"backtracks."
        ),
        "audit_head": (
            r"Dim. & Method & $\Delta t$ & Runs & Steps & RB & "
            r"$\bar n_{\mathrm{nl}}$ & $n^{\max}_{\mathrm{nl}}$ & "
            r"$\Sigma_{\mathrm{LS}}$ \\"
        ),
        "audit_total": "Total",
    },
    "es": {
        "mesh": lambda dim, mesh: (
            ("tetra no estr." if dim == "3D" else "triang. no estr.")
            if mesh == "triangular_Unstr"
            else mesh.replace("_", " ")
        ),
        "fit": "Ajuste",
        "meshcol": "Tipo de malla",
        "caption": (
            "Resultados manufacturados {dim} para {method}, "
            "\\code{{crankNicolson}}, masa \\code{{consistent}}, estados "
            "\\code{{RKF45}} y $\\Delta t={dt}$. El primer bloque de "
            "convergencia se ajusta sobre {all_lbl}; el segundo sobre "
            "{fine_lbl}. La primera columna registra la terna "
            "\\code{{Vm-states-Iion}} efectivamente usada en la corrida. "
            "\\textbf{{RB}} es la cantidad de pasos de tiempo que agotaron el "
            "presupuesto de iteraciones no lineales, acumulada sobre los "
            "{n_meshes} tamaños de malla de la fila ({n_steps} pasos de tiempo "
            "en total por fila). {note} "
            "\\textbf{{$n_{{\\mathrm{{nl}}}}$}} reporta la mediana y el "
            "máximo de iteraciones no lineales por paso sobre esos mismos "
            "pasos; el presupuesto configurado es "
            "\\code{{nonlinearIterations}}$=2000$."
        ),
        "note_picard": (
            "\\code{Picard} acepta el último iterado con advertencia en vez "
            "de hacer rollback, de modo que para este método la columna cuenta "
            "pasos aceptados sin converger."
        ),
        "note_other": (
            "Con \\code{nonlinearAcceptUnconverged=false} un paso no "
            "convergido se revierte a los valores del inicio del paso mientras "
            "el reloj avanza igual."
        ),
        "empty": "No se encontró archivo de errores agregado para este caso.",
        "subsection": lambda dim: (
            f"Problema manufacturado {'bidimensional' if dim == '2D' else 'tridimensional'}"
        ),
        "intro": [
            "Las tablas de esta subsección se generan directamente desde el árbol de",
            r"resultados \path{{{root}}}. Las raíces por método son {roots}.",
            r"Todas las entradas usan \code{{crankNicolson}}, la matriz de masa",
            r"\code{{consistent}} y \code{{RKF45}} para las ODEs de estado.",
            r"La primera columna reporta \code{{Vm-states-Iion}}: el grado de reconstrucción",
            r"de $V_m$, el tratamiento de estados (\code{{CCp*}}, \code{{GP}} o \code{{na}}) y el",
            r"ajuste de cuadratura/reconstrucción de $\Iion$.",
            "Los ajustes usan {all_lbl} para el rango amplio de mallas y",
            "{fine_lbl} para el rango de mallas finas.",
            r"Las dos últimas columnas auditan el solver no lineal: \code{{RB}} cuenta los",
            r"pasos revertidos (rollback) y $n_{{\mathrm{{nl}}}}$ resume las iteraciones no",
            r"lineales por paso (mediana/máximo).",
            "",
        ],
        "audit_caption": (
            r"Auditoría global de convergencia no lineal sobre el barrido "
            r"completo de \code{example0}. \code{RB} cuenta los pasos de tiempo "
            r"que agotaron el presupuesto de \code{nonlinearIterations}$=2000$ "
            r"(y que por lo tanto habrían sido revertidos para "
            r"\code{diagonalIion}/\code{JFNK}, o aceptados con advertencia para "
            r"\code{Picard}). $\bar n_{\mathrm{nl}}$ y $n^{\max}_{\mathrm{nl}}$ "
            r"son la media y el máximo de iteraciones no lineales por paso; "
            r"$\Sigma_{\mathrm{LS}}$ es el total de retrocesos de la búsqueda de "
            r"línea de Armijo."
        ),
        "audit_head": (
            r"Dim. & Método & $\Delta t$ & Corridas & Pasos & RB & "
            r"$\bar n_{\mathrm{nl}}$ & $n^{\max}_{\mathrm{nl}}$ & "
            r"$\Sigma_{\mathrm{LS}}$ \\"
        ),
        "audit_total": "Total",
    },
}


def table_label(lang: str, dim: str, method: str, dt_tag: str) -> str:
    return f"tab:{lang}-results-{dim.lower()}-{method.lower()}-{dt_tag}"


def table_file_name(lang: str, dim: str, method: str, dt_tag: str) -> str:
    return f"results_{dim}_{method}_dt_{dt_tag}_{lang}.tex"


def write_table(
    lang: str,
    dim: str,
    method: str,
    dt_tag: str,
    rows: list,
    stats: dict[tuple[str, ...], dict],
) -> str:
    L = LANGS[lang]
    all_lbl = range_label(dim, "all")
    fine_lbl = range_label(dim, "fine")

    sample = next(
        (s for k, s in stats.items()
         if k[0] == dim and k[1] == method and k[2] == dt_tag),
        None,
    )
    n_meshes = len(sample["meshes"]) if sample else 0
    n_steps = sample["steps"] if sample else 0

    note = L["note_picard"] if method == "Picard" else L["note_other"]
    caption = L["caption"].format(
        dim=dim,
        method=G.method_label(method),
        dt=dt_to_latex(dt_tag),
        all_lbl=all_lbl,
        fine_lbl=fine_lbl,
        n_meshes=n_meshes,
        n_steps=n_steps,
        note=note,
    )

    header_group = (
        r"& & & \multicolumn{6}{c}{" + L["fit"] + " " + all_lbl + r"} & "
        r"\multicolumn{6}{c}{" + L["fit"] + " " + fine_lbl + r"} & & & & \\"
    )
    header_cols = (
        r"\code{Vm-states-Iion} & " + L["meshcol"] + r" & $\alpha$ & "
        r"$V_m$ & $V_m^G$ & $u_1$ & $u_1^G$ & $u_2$ & $u_2^G$ & "
        r"$V_m$ & $V_m^G$ & $u_1$ & $u_1^G$ & $u_2$ & $u_2^G$ & "
        r"$t_{\Sigma}$ [s] & RSS$_{\Sigma}$ [MB] & RB & $n_{\mathrm{nl}}$ \\"
    )

    lines = [
        r"\tiny",
        r"\setlength{\tabcolsep}{1.4pt}",
        r"\begin{longtable}{lllrrrrrrrrrrrrrrrr}",
        rf"\caption{{{caption}}}\label{{{table_label(lang, dim, method, dt_tag)}}}\\",
        r"\toprule",
        header_group,
        r"\cmidrule(lr){4-9}\cmidrule(lr){10-15}",
        header_cols,
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        header_group,
        r"\cmidrule(lr){4-9}\cmidrule(lr){10-15}",
        header_cols,
        r"\midrule",
        r"\endhead",
    ]

    last_triplet = None
    for row in rows:
        triplet = G.model_triplet(row)
        if last_triplet is not None and triplet != last_triplet:
            lines.append(r"\hline")
        last_triplet = triplet

        entry = stats.get(
            (dim, method, dt_tag, row.mesh, row.alpha_tag,
             row.vm, row.states, row.iion)
        )
        rollbacks, iter_cell = robustness_cells(entry)

        values = [
            G.tex_code(triplet),
            L["mesh"](dim, row.mesh),
            G.format_alpha(row.alpha_value, row.alpha_tag),
            *(G.format_order(row.order_all[m]) for m in G.METRICS),
            *(G.format_order(row.order_fine[m]) for m in G.METRICS),
            G.format_cost(row.total_time),
            G.format_cost(row.total_peak_rss),
            rollbacks,
            iter_cell,
        ]
        lines.append(" & ".join(values) + r" \\")

    if not rows:
        lines.append(r"\multicolumn{19}{c}{" + L["empty"] + r"} \\")

    lines.extend([r"\bottomrule", r"\end{longtable}", r"\normalsize"])
    file_name = table_file_name(lang, dim, method, dt_tag)
    (TABLE_DIR / file_name).write_text("\n".join(lines) + "\n")
    return file_name


def write_section_fragment(
    lang: str, dim: str, table_files: dict[tuple[str, str], str]
) -> None:
    L = LANGS[lang]
    roots = ", ".join(
        rf"\path{{{G.result_root(dim, m)}}}" for m in G.METHODS
    )
    lines = [
        rf"\subsection{{{L['subsection'](dim)}}}",
        rf"\label{{sec:{lang}-results-{dim.lower()}}}",
    ]
    for tpl in L["intro"]:
        lines.append(
            tpl.format(
                root=G.RESULT_ROOT,
                roots=roots,
                all_lbl=range_label(dim, "all"),
                fine_lbl=range_label(dim, "fine"),
            )
        )
    for method in G.METHODS:
        lines.append(rf"\subsubsection{{{G.method_label(method)}}}")
        for dt_tag in G.DT_TAGS:
            lines.append(rf"\paragraph{{$\Delta t={dt_to_latex(dt_tag)}$}}")
            lines.append(rf"\input{{tables/{table_files[(method, dt_tag)]}}}")
        lines.append("")
    (TABLE_DIR / f"results_{dim}_all_{lang}.tex").write_text(
        "\n".join(lines) + "\n"
    )


def write_audit_summary(lang: str, stats: dict[tuple[str, ...], dict]) -> None:
    """Emit the global audit table (one line per dim/method/dt)."""
    L = LANGS[lang]
    agg = defaultdict(
        lambda: {"runs": 0, "steps": 0, "rollbacks": 0,
                 "iters": [], "imax": 0, "ls": 0}
    )
    for key, entry in stats.items():
        a = agg[(key[0], key[1], key[2])]
        a["runs"] += entry["runs"]
        a["steps"] += entry["steps"]
        a["rollbacks"] += entry["rollbacks"]
        a["iters"].extend(entry["iters"])
        a["imax"] = max(a["imax"], entry["iter_max"])
        a["ls"] += entry["linesearch"]

    lines = [
        r"\begin{table}[h!]",
        r"\centering",
        r"\footnotesize",
        rf"\caption{{{L['audit_caption']}}}",
        rf"\label{{tab:{lang}-rollback-audit}}",
        r"\begin{tabular}{@{}llrrrrrrr@{}}",
        r"\toprule",
        L["audit_head"],
        r"\midrule",
    ]
    for dim in ("2D", "3D"):
        for method in G.METHODS:
            for dt_tag in G.DT_TAGS:
                a = agg.get((dim, method, dt_tag))
                if not a or not a["iters"]:
                    continue
                mean = sum(a["iters"]) / len(a["iters"])
                lines.append(
                    f"{dim} & {G.method_label(method)} & ${dt_to_latex(dt_tag)}$ & "
                    f"{a['runs']} & {a['steps']} & \\textbf{{{a['rollbacks']}}} & "
                    f"{mean:.3f} & {a['imax']} & {a['ls']} \\\\"
                )
        if dim == "2D":
            lines.append(r"\midrule")

    total_runs = sum(a["runs"] for a in agg.values())
    total_steps = sum(a["steps"] for a in agg.values())
    total_rb = sum(a["rollbacks"] for a in agg.values())
    lines.extend([
        r"\midrule",
        rf"\multicolumn{{3}}{{@{{}}l}}{{{L['audit_total']}}} & {total_runs} & "
        rf"{total_steps} & \textbf{{{total_rb}}} & & & \\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ])
    (TABLE_DIR / f"rollback_audit_{lang}.tex").write_text(
        "\n".join(lines) + "\n"
    )


def main() -> None:
    TABLE_DIR.mkdir(exist_ok=True)
    print("scanning nonlinearResiduals histories ...")
    stats = scan_nonlinear_histories()
    print(f"  {len(stats)} configurations")

    order_rows = {
        (dim, method, dt_tag): G.collect_rows(dim, method, dt_tag)
        for dim in ("2D", "3D")
        for method in G.METHODS
        for dt_tag in G.DT_TAGS
    }

    for lang in ("en", "es"):
        write_audit_summary(lang, stats)
        for dim in ("2D", "3D"):
            table_files = {}
            for method in G.METHODS:
                for dt_tag in G.DT_TAGS:
                    rows = order_rows[(dim, method, dt_tag)]
                    table_files[(method, dt_tag)] = write_table(
                        lang, dim, method, dt_tag, rows, stats
                    )
            write_section_fragment(lang, dim, table_files)
            print(f"  wrote {dim} tables ({lang})")

    summary = {
        "_".join(k): {
            "runs": v["runs"],
            "steps": v["steps"],
            "rollbacks": v["rollbacks"],
            "iter_max": v["iter_max"],
            "iter_mean": sum(v["iters"]) / len(v["iters"]) if v["iters"] else None,
            "linesearch": v["linesearch"],
        }
        for k, v in stats.items()
    }
    (TABLE_DIR / "rollback_audit.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
