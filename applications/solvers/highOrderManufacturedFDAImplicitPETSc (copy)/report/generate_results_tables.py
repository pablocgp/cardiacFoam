#!/usr/bin/env python3
"""Generate LaTeX result tables for the manufactured JFNK report (2D and 3D)."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path


REPORT_DIR = Path(__file__).resolve().parent
TABLE_DIR = REPORT_DIR / "tables"

TUTORIAL_ROOT = Path(
    "/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro"
)

# 2D: every nonlinear method lives under the same root.
RESULT_ROOT_2D = TUTORIAL_ROOT / "highOrderManufacturedFDAImplicitJfNK/results/example0/2D"

# 3D: each nonlinear method lives in its own tutorial root. The "_1" siblings
# only hold in-progress runs without an aggregated errors_vs_N file and are
# therefore harmless to include alongside their main counterparts.
RESULT_ROOTS_3D: dict[str, tuple[Path, ...]] = {
    "JFNK": (
        TUTORIAL_ROOT / "highOrderManufacturedFDAImplicitJfNK/results/example0/3D",
    ),
    "diagonalIion": (
        TUTORIAL_ROOT / "highOrderManufacturedFDAImplicitJfNK_diagonalIion/results/example0/3D",
        TUTORIAL_ROOT / "highOrderManufacturedFDAImplicitJfNK_diagonalIion_1/results/example0/3D",
    ),
    "Picard": (
        TUTORIAL_ROOT / "highOrderManufacturedFDAImplicitJfNK_Picard/results/example0/3D",
        TUTORIAL_ROOT / "highOrderManufacturedFDAImplicitJfNK_Picard_1/results/example0/3D",
    ),
}

METHODS = ("Picard", "diagonalIion", "JFNK")
DT_TAGS = ("1p0em02", "1p0em03")
SCHEMES = ("backwardEuler", "crankNicolson")
MASSES = ("consistent", "lumped")
LEVELS = ("NO", "p1", "p2", "p3")
METRICS = ("Vm_L2", "VmG_L2", "u1_L2", "u1G_L2", "u2_L2", "u2G_L2")

# Mesh ranges used for the two log-log fits, per dimension.
N_RANGES: dict[str, dict[str, tuple[int, int]]] = {
    "2D": {"all": (10, 100), "fine": (70, 100)},
    "3D": {"all": (10, 30), "fine": (20, 30)},
}


CONFIG_RE = re.compile(
    r"^hoVm_(?P<vm>.+)_hoStates_(?P<iion>.+)_scheme_(?P<scheme>.+)"
    r"_mass_(?P<mass>.+)_nonlin_(?P<method>.+)_ode_(?P<ode>.+)$"
)


@dataclass(frozen=True)
class ResultRow:
    vm: str
    iion: str
    scheme: str
    mass: str
    method: str
    dt_tag: str
    order_all: dict[str, float | None]
    order_fine: dict[str, float | None]
    total_time: float | None
    total_peak_rss: float | None


def dt_to_latex(tag: str) -> str:
    if tag == "1p0em02":
        return r"$10^{-2}$"
    if tag == "1p0em03":
        return r"$10^{-3}$"
    return r"\code{" + tag.replace("_", r"\_") + r"}"


def scheme_short(scheme: str) -> str:
    return {"backwardEuler": "bE", "crankNicolson": "cN"}.get(scheme, scheme)


def mass_short(mass: str) -> str:
    return {"consistent": "cons.", "lumped": "lump."}.get(mass, mass)


def method_label(method: str) -> str:
    if method == "diagonalIion":
        return r"\code{diagonalIion}"
    return r"\code{" + method + r"}"


def tex_code(value: str) -> str:
    return r"\code{" + value.replace("_", r"\_") + r"}"


def parse_error_file(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    if not lines:
        return []
    header = re.sub(r"^#\s*", "", lines[0].strip()).split()
    rows: list[dict[str, str]] = []
    for line in lines[1:]:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) < len(header):
            continue
        rows.append(dict(zip(header, parts)))
    return rows


def fit_order(
    rows: list[dict[str, str]],
    metric: str,
    n_min: int | None = None,
    n_max: int | None = None,
) -> float | None:
    xs: list[float] = []
    ys: list[float] = []
    for row in rows:
        try:
            n_val = int(float(row["N"]))
            err = float(row[metric])
        except (KeyError, ValueError):
            continue
        if n_min is not None and n_val < n_min:
            continue
        if n_max is not None and n_val > n_max:
            continue
        if not math.isfinite(err) or err <= 0:
            continue
        h = 1.0 / n_val
        xs.append(math.log(h))
        ys.append(math.log(err))
    if len(xs) < 2:
        return None
    xbar = sum(xs) / len(xs)
    ybar = sum(ys) / len(ys)
    den = sum((x - xbar) ** 2 for x in xs)
    if den <= 0:
        return None
    num = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys))
    return num / den


def read_scalar(pattern: str, text: str) -> float | None:
    match = re.search(pattern, text)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def accumulated_cost(error_file: Path, rows: list[dict[str, str]]) -> tuple[float | None, float | None]:
    total_time = 0.0
    total_rss = 0.0
    time_count = 0
    rss_count = 0
    for row in rows:
        dat_name = row.get("dat_file")
        if not dat_name:
            continue
        dat_path = error_file.parent / dat_name
        if not dat_path.exists():
            continue
        text = dat_path.read_text(errors="ignore")
        t_val = read_scalar(r"\btotal\s*=\s*([0-9.eE+-]+)", text)
        m_val = read_scalar(r"peakRSS_MB\s*=\s*([0-9.eE+-]+)", text)
        if t_val is not None:
            total_time += t_val
            time_count += 1
        if m_val is not None:
            total_rss += m_val
            rss_count += 1
    return (total_time if time_count else None, total_rss if rss_count else None)


def result_roots(dim: str, method: str) -> tuple[Path, ...]:
    if dim == "2D":
        return (RESULT_ROOT_2D,)
    return RESULT_ROOTS_3D.get(method, ())


def collect_rows(dim: str, method: str, dt_tag: str) -> list[ResultRow]:
    n_min_all, n_max_all = N_RANGES[dim]["all"]
    n_min_fine, n_max_fine = N_RANGES[dim]["fine"]

    rows: list[ResultRow] = []
    seen_configs: set[str] = set()
    for root in result_roots(dim, method):
        if not root.exists():
            continue
        for error_file in root.glob(f"*/endT_*_dt_{dt_tag}/errors_vs_N_dt_{dt_tag}.dat"):
            config_name = error_file.parent.parent.name
            match = CONFIG_RE.match(config_name)
            if not match or match.group("method") != method:
                continue
            if config_name in seen_configs:
                continue
            seen_configs.add(config_name)
            error_rows = parse_error_file(error_file)
            orders_all = {
                metric: fit_order(error_rows, metric, n_min=n_min_all, n_max=n_max_all)
                for metric in METRICS
            }
            orders_fine = {
                metric: fit_order(error_rows, metric, n_min=n_min_fine, n_max=n_max_fine)
                for metric in METRICS
            }
            total_time, total_rss = accumulated_cost(error_file, error_rows)
            rows.append(
                ResultRow(
                    vm=match.group("vm"),
                    iion=match.group("iion"),
                    scheme=match.group("scheme"),
                    mass=match.group("mass"),
                    method=method,
                    dt_tag=dt_tag,
                    order_all=orders_all,
                    order_fine=orders_fine,
                    total_time=total_time,
                    total_peak_rss=total_rss,
                )
            )
    return sorted(rows, key=sort_key)


def sort_key(row: ResultRow) -> tuple[int, int, int, int]:
    return (
        SCHEMES.index(row.scheme),
        MASSES.index(row.mass),
        LEVELS.index(row.vm),
        LEVELS.index(row.iion),
    )


def format_order(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "--"
    return f"{value:.3f}"


def format_cost(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "--"
    if abs(value) >= 1000.0:
        return f"{value:.1f}"
    if abs(value) >= 100.0:
        return f"{value:.2f}"
    return f"{value:.3f}"


def figure_label(dim: str, method: str, dt_tag: str, scheme: str, mass: str) -> str:
    return f"fig:results-{dim.lower()}-{method.lower()}-{dt_tag}-{scheme.lower()}-{mass}"


def table_label(dim: str, method: str, dt_tag: str) -> str:
    return f"tab:results-{dim.lower()}-{method.lower()}-{dt_tag}"


def table_file_name(dim: str, method: str, dt_tag: str) -> str:
    return f"results_{dim}_{method}_dt_{dt_tag}.tex"


def panel_file(dim: str, method: str, dt_tag: str, scheme: str, mass: str, field: str, metric: str) -> str:
    return (
        f"figures/implicit_{dim}_{method}_{scheme}_{mass}_{field}_{metric}"
        f"_dt_{dt_tag}.pdf"
    )


def range_label(dim: str, key: str) -> str:
    n_min, n_max = N_RANGES[dim][key]
    return f"$N={n_min}$--${n_max}$"


def write_table(dim: str, method: str, dt_tag: str, rows: list[ResultRow]) -> str:
    label = table_label(dim, method, dt_tag)
    fig_refs = [
        r"\ref{" + figure_label(dim, method, dt_tag, scheme, mass) + r"}"
        for scheme in SCHEMES
        for mass in MASSES
    ]
    all_lbl = range_label(dim, "all")
    fine_lbl = range_label(dim, "fine")
    caption = (
        f"{dim} manufactured results for {method_label(method)} with "
        f"$\\Delta t={dt_to_latex(dt_tag).strip('$')}$."
        f" The first six convergence columns are log--log fits over meshes "
        f"{all_lbl}; the next six columns are fits over {fine_lbl}."
        " The time column is the sum of the reported wall-clock total over all "
        "mesh runs in the row, and the memory column is the sum of the reported "
        "peak resident memory in MB over the same mesh runs."
        f" The corresponding convergence figures are Figs.~{', '.join(fig_refs)}."
    )
    lines = [
        r"\tiny",
        r"\setlength{\tabcolsep}{1.7pt}",
        r"\begin{longtable}{llllrrrrrrrrrrrrrr}",
        rf"\caption{{{caption}}}\label{{{label}}}\\",
        r"\toprule",
        r"& & & & \multicolumn{6}{c}{Fit " + all_lbl + r"} & "
        r"\multicolumn{6}{c}{Fit " + fine_lbl + r"} & & \\",
        r"\cmidrule(lr){5-10}\cmidrule(lr){11-16}",
        r"$V_m$ & $I_{\mathrm{ion}}$ & Scheme & Mass & "
        r"$V_m$ & $V_m^G$ & $u_1$ & $u_1^G$ & $u_2$ & $u_2^G$ & "
        r"$V_m$ & $V_m^G$ & $u_1$ & $u_1^G$ & $u_2$ & $u_2^G$ & "
        r"$t_{\Sigma}$ [s] & RSS$_{\Sigma}$ [MB] \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"& & & & \multicolumn{6}{c}{Fit " + all_lbl + r"} & "
        r"\multicolumn{6}{c}{Fit " + fine_lbl + r"} & & \\",
        r"\cmidrule(lr){5-10}\cmidrule(lr){11-16}",
        r"$V_m$ & $I_{\mathrm{ion}}$ & Scheme & Mass & "
        r"$V_m$ & $V_m^G$ & $u_1$ & $u_1^G$ & $u_2$ & $u_2^G$ & "
        r"$V_m$ & $V_m^G$ & $u_1$ & $u_1^G$ & $u_2$ & $u_2^G$ & "
        r"$t_{\Sigma}$ [s] & RSS$_{\Sigma}$ [MB] \\",
        r"\midrule",
        r"\endhead",
    ]
    last_group: tuple[str, str] | None = None
    for row in rows:
        group = (row.scheme, row.mass)
        if last_group is not None and group != last_group:
            lines.append(r"\hline")
        last_group = group
        values = [
            tex_code(row.vm),
            tex_code(row.iion),
            scheme_short(row.scheme),
            mass_short(row.mass),
            *(format_order(row.order_all[m]) for m in METRICS),
            *(format_order(row.order_fine[m]) for m in METRICS),
            format_cost(row.total_time),
            format_cost(row.total_peak_rss),
        ]
        lines.append(" & ".join(values) + r" \\")
    lines.extend([r"\bottomrule", r"\end{longtable}", r"\normalsize"])
    output = "\n".join(lines) + "\n"
    (TABLE_DIR / table_file_name(dim, method, dt_tag)).write_text(output)
    return table_file_name(dim, method, dt_tag)


def write_figure_block(dim: str, method: str, dt_tag: str, scheme: str, mass: str) -> str:
    fields = (("Vm", "V_m"), ("u1", "u_1"), ("u2", "u_2"))
    metrics = (("L2", "L_2"), ("G2", "G_2"))
    label = figure_label(dim, method, dt_tag, scheme, mass)
    tab_label = table_label(dim, method, dt_tag)
    lines = [
        r"\begin{figure}[p]",
        r"\centering",
    ]
    panel = 0
    for field, field_tex in fields:
        for metric, metric_tex in metrics:
            panel += 1
            suffix = "a" if panel == 1 else "b" if panel == 2 else "c" if panel == 3 else "d" if panel == 4 else "e" if panel == 5 else "f"
            lines.extend(
                [
                    r"\begin{subfigure}{0.485\linewidth}",
                    r"\centering",
                    rf"\maybeincludegraphics[width=\linewidth]{{{panel_file(dim, method, dt_tag, scheme, mass, field, metric)}}}",
                    rf"\caption{{$\|{field_tex}\|_{{{metric_tex}}}$}}",
                    rf"\label{{{label}-{suffix}}}",
                    r"\end{subfigure}",
                ]
            )
            if panel % 2 == 1:
                lines.append(r"\hfill")
            elif panel < 6:
                lines.append(r"\par\vspace{0.8em}")
    lines.extend(
        [
            rf"\caption{{{dim} convergence curves for {method_label(method)}, "
            rf"$\Delta t={dt_to_latex(dt_tag).strip('$')}$, {scheme_short(scheme)} "
            rf"with {mass_short(mass)} mass. The six panels correspond to the six "
            rf"order columns in Table~\ref{{{tab_label}}} for this scheme/mass "
            rf"block. The plot legends report the two fitted orders and the "
            rf"accumulated runtime/memory cost for each curve.}}",
            rf"\label{{{label}}}",
            r"\end{figure}",
            r"\FloatBarrier",
        ]
    )
    return "\n".join(lines) + "\n"


def section_intro(dim: str) -> list[str]:
    if dim == "2D":
        all_lbl = range_label(dim, "all")
        fine_lbl = range_label(dim, "fine")
        return [
            r"\subsection{Two-dimensional manufactured problem}",
            "The tables in this subsection are generated from",
            r"\path{/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/highOrderManufacturedFDAImplicitJfNK/results/example0/2D/}.",
            "For each nonlinear method and time step, one table reports the observed",
            "orders for all available high-order configurations, time schemes, and mass",
            "matrices. Each table is followed by the corresponding six-panel convergence",
            f"figures, one figure for each scheme/mass block. Fits use {all_lbl} for the",
            f"general order column and {fine_lbl} for the asymptotic column.",
            "",
        ]
    return [
        r"\subsection{Three-dimensional manufactured problem}",
        "The tables in this subsection collect the 3D results that are currently",
        "available. The Picard, diagonal-$\\Iion$, and JFNK sweeps live in separate",
        "tutorial roots:",
        r"\begin{itemize}",
        r"\item \path{highOrderManufacturedFDAImplicitJfNK/results/example0/3D/} (\code{JFNK}),",
        r"\item \path{highOrderManufacturedFDAImplicitJfNK_diagonalIion/results/example0/3D/} (and \path{_1}),",
        r"\item \path{highOrderManufacturedFDAImplicitJfNK_Picard/results/example0/3D/} (and \path{_1}).",
        r"\end{itemize}",
        f"Fits use {range_label('3D', 'all')} for the general order column and",
        f"{range_label('3D', 'fine')} for the asymptotic column, in line with the",
        "coarser 3D mesh sweep $N\\in[10,30]$. Missing rows correspond to",
        "configurations that have not yet been completed or whose error file is",
        "not yet aggregated.",
        "",
    ]


def write_section_fragment(dim: str, table_files: dict[tuple[str, str], str]) -> None:
    lines = list(section_intro(dim))
    for method in METHODS:
        lines.append(rf"\subsubsection{{{method_label(method)}}}")
        for dt_tag in DT_TAGS:
            lines.append(rf"\paragraph{{$\Delta t={dt_to_latex(dt_tag).strip('$')}$}}")
            lines.append(rf"\input{{tables/{table_files[(method, dt_tag)]}}}")
            for scheme in SCHEMES:
                for mass in MASSES:
                    lines.append(write_figure_block(dim, method, dt_tag, scheme, mass))
    (TABLE_DIR / f"results_{dim}_all.tex").write_text("\n".join(lines) + "\n")


def build_dim(dim: str) -> None:
    table_files: dict[tuple[str, str], str] = {}
    for method in METHODS:
        for dt_tag in DT_TAGS:
            rows = collect_rows(dim, method, dt_tag)
            table_files[(method, dt_tag)] = write_table(dim, method, dt_tag, rows)
    write_section_fragment(dim, table_files)


def main() -> None:
    TABLE_DIR.mkdir(exist_ok=True)
    build_dim("2D")
    build_dim("3D")


if __name__ == "__main__":
    main()
