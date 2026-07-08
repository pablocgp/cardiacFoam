#!/usr/bin/env python3
"""Generate LaTeX result tables for the PETSc manufactured-solution report."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path


REPORT_DIR = Path(__file__).resolve().parent
TABLE_DIR = REPORT_DIR / "tables"

RESULT_ROOT = Path(
    "/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro"
    "/highOrderManufacturedFDAImplicitPETSc/results/example0"
)

METHODS = ("Picard", "diagonalIion", "JFNK")
DT_TAGS = ("1p0em02", "1p0em03")
METRICS = ("Vm_L2", "VmG_L2", "u1_L2", "u1G_L2", "u2_L2", "u2G_L2")

N_RANGES: dict[str, dict[str, tuple[int, int]]] = {
    "2D": {"all": (10, 100), "fine": (70, 100)},
    "3D": {"all": (10, 30), "fine": (20, 30)},
}

LEVELS = ("NO", "p1", "p2", "p3")
STATE_LEVELS = ("na", "CCp1", "CCp2", "CCp3", "GP")
MESH_ORDER = ("hexa", "triangular_Unstr")

MESH_RE = re.compile(r"^(?P<dim>[23]D)_mesh_(?P<mesh>.+)_alpha_(?P<alpha>.+)$")
CONFIG_RE = re.compile(
    r"^hoVm_(?P<vm>[^_]+)_states_(?P<states>[^_]+)_hoIion_(?P<iion>[^_]+)"
    r"_(?P<scheme>[^_]+)_(?P<mass>[^_]+)_(?P<method>[^_]+)_(?P<ode>[^_]+)$"
)


@dataclass(frozen=True)
class ResultRow:
    vm: str
    states: str
    iion: str
    mesh: str
    alpha_tag: str
    alpha_value: float | None
    method: str
    dt_tag: str
    order_all: dict[str, float | None]
    order_fine: dict[str, float | None]
    total_time: float | None
    total_peak_rss: float | None


def dt_to_latex(tag: str) -> str:
    if tag == "1p0em02":
        return r"10^{-2}"
    if tag == "1p0em03":
        return r"10^{-3}"
    return r"\code{" + tag.replace("_", r"\_") + r"}"


def method_label(method: str) -> str:
    return r"\code{" + method + r"}"


def tex_code(value: str) -> str:
    return r"\code{" + value.replace("_", r"\_") + r"}"


def model_triplet(row: ResultRow) -> str:
    return f"{row.vm}--{row.states}--{row.iion}"


def mesh_label(dim: str, mesh: str) -> str:
    if mesh == "triangular_Unstr":
        return "tetra unstr" if dim == "3D" else "triangular unstr"
    return mesh.replace("_", " ")


def parse_alpha_tag(tag: str) -> float | None:
    try:
        return float(tag.replace("p", ".").replace("em", "e-").replace("e00", "e0"))
    except ValueError:
        return None


def format_alpha(value: float | None, tag: str) -> str:
    if value is None or not math.isfinite(value):
        return tex_code(tag)
    if abs(value) < 1.0:
        return f"{value:.1f}"
    return f"{value:g}"


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
    n_min: int,
    n_max: int,
) -> float | None:
    xs: list[float] = []
    ys: list[float] = []
    for row in rows:
        try:
            n_val = int(float(row["N"]))
            err = float(row[metric])
        except (KeyError, ValueError):
            continue
        if n_val < n_min or n_val > n_max:
            continue
        if not math.isfinite(err) or err <= 0:
            continue
        xs.append(math.log(1.0 / n_val))
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


def result_root(dim: str, method: str) -> Path:
    return RESULT_ROOT / f"{dim}_CN_cons_{method}_RKF45"


def collect_rows(dim: str, method: str, dt_tag: str) -> list[ResultRow]:
    n_min_all, n_max_all = N_RANGES[dim]["all"]
    n_min_fine, n_max_fine = N_RANGES[dim]["fine"]

    rows: list[ResultRow] = []
    root = result_root(dim, method)
    if not root.exists():
        return rows

    for error_file in root.glob(f"*/hoVm_*/endT_*_dt_{dt_tag}/errors_vs_N_dt_{dt_tag}.dat"):
        mesh_dir = error_file.parents[2].name
        config_name = error_file.parents[1].name

        mesh_match = MESH_RE.match(mesh_dir)
        config_match = CONFIG_RE.match(config_name)
        if not mesh_match or not config_match:
            continue
        if mesh_match.group("dim") != dim:
            continue
        if config_match.group("method") != method:
            continue

        error_rows = parse_error_file(error_file)
        if not error_rows:
            continue

        alpha_value = None
        try:
            alpha_value = float(error_rows[0]["stabilisation_alpha"])
        except (KeyError, ValueError):
            alpha_value = parse_alpha_tag(mesh_match.group("alpha"))

        total_time, total_peak_rss = accumulated_cost(error_file, error_rows)

        rows.append(
            ResultRow(
                vm=config_match.group("vm"),
                states=config_match.group("states"),
                iion=config_match.group("iion"),
                mesh=mesh_match.group("mesh"),
                alpha_tag=mesh_match.group("alpha"),
                alpha_value=alpha_value,
                method=method,
                dt_tag=dt_tag,
                order_all={
                    metric: fit_order(error_rows, metric, n_min_all, n_max_all)
                    for metric in METRICS
                },
                order_fine={
                    metric: fit_order(error_rows, metric, n_min_fine, n_max_fine)
                    for metric in METRICS
                },
                total_time=total_time,
                total_peak_rss=total_peak_rss,
            )
        )

    return sorted(rows, key=sort_key)


def index_or_end(values: tuple[str, ...], value: str) -> int:
    try:
        return values.index(value)
    except ValueError:
        return len(values)


def sort_key(row: ResultRow) -> tuple[int, int, int, int, int, float, str]:
    alpha = row.alpha_value if row.alpha_value is not None else math.inf
    return (
        index_or_end(LEVELS, row.vm),
        index_or_end(STATE_LEVELS, row.states),
        index_or_end(LEVELS, row.iion),
        index_or_end(MESH_ORDER, row.mesh),
        0 if row.alpha_tag == "0p0e00" else 1,
        alpha,
        row.alpha_tag,
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


def table_label(dim: str, method: str, dt_tag: str) -> str:
    return f"tab:results-{dim.lower()}-{method.lower()}-{dt_tag}"


def table_file_name(dim: str, method: str, dt_tag: str) -> str:
    return f"results_{dim}_{method}_dt_{dt_tag}.tex"


def range_label(dim: str, key: str) -> str:
    n_min, n_max = N_RANGES[dim][key]
    return f"$N={n_min}$--${n_max}$"


def write_table(dim: str, method: str, dt_tag: str, rows: list[ResultRow]) -> str:
    all_lbl = range_label(dim, "all")
    fine_lbl = range_label(dim, "fine")
    caption = (
        f"{dim} manufactured results for {method_label(method)}, "
        r"\code{crankNicolson}, \code{consistent} mass, \code{RKF45} states, "
        f"and $\\Delta t={dt_to_latex(dt_tag)}$. "
        f"The first convergence block is fitted over {all_lbl}; "
        f"the second over {fine_lbl}. "
        "The first column records the actual "
        r"\code{Vm-states-Iion} triplet used by the run."
    )

    lines = [
        r"\tiny",
        r"\setlength{\tabcolsep}{1.6pt}",
        r"\begin{longtable}{lllrrrrrrrrrrrrrr}",
        rf"\caption{{{caption}}}\label{{{table_label(dim, method, dt_tag)}}}\\",
        r"\toprule",
        r"& & & \multicolumn{6}{c}{Fit " + all_lbl + r"} & "
        r"\multicolumn{6}{c}{Fit " + fine_lbl + r"} & & \\",
        r"\cmidrule(lr){4-9}\cmidrule(lr){10-15}",
        r"\code{Vm-states-Iion} & Mesh type & $\alpha$ & "
        r"$V_m$ & $V_m^G$ & $u_1$ & $u_1^G$ & $u_2$ & $u_2^G$ & "
        r"$V_m$ & $V_m^G$ & $u_1$ & $u_1^G$ & $u_2$ & $u_2^G$ & "
        r"$t_{\Sigma}$ [s] & RSS$_{\Sigma}$ [MB] \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"& & & \multicolumn{6}{c}{Fit " + all_lbl + r"} & "
        r"\multicolumn{6}{c}{Fit " + fine_lbl + r"} & & \\",
        r"\cmidrule(lr){4-9}\cmidrule(lr){10-15}",
        r"\code{Vm-states-Iion} & Mesh type & $\alpha$ & "
        r"$V_m$ & $V_m^G$ & $u_1$ & $u_1^G$ & $u_2$ & $u_2^G$ & "
        r"$V_m$ & $V_m^G$ & $u_1$ & $u_1^G$ & $u_2$ & $u_2^G$ & "
        r"$t_{\Sigma}$ [s] & RSS$_{\Sigma}$ [MB] \\",
        r"\midrule",
        r"\endhead",
    ]

    last_triplet: str | None = None
    for row in rows:
        triplet = model_triplet(row)
        if last_triplet is not None and triplet != last_triplet:
            lines.append(r"\hline")
        last_triplet = triplet

        values = [
            tex_code(triplet),
            mesh_label(dim, row.mesh),
            format_alpha(row.alpha_value, row.alpha_tag),
            *(format_order(row.order_all[m]) for m in METRICS),
            *(format_order(row.order_fine[m]) for m in METRICS),
            format_cost(row.total_time),
            format_cost(row.total_peak_rss),
        ]
        lines.append(" & ".join(values) + r" \\")

    if not rows:
        lines.append(r"\multicolumn{17}{c}{No aggregated error file found for this case.} \\")

    lines.extend([r"\bottomrule", r"\end{longtable}", r"\normalsize"])
    output = "\n".join(lines) + "\n"
    file_name = table_file_name(dim, method, dt_tag)
    (TABLE_DIR / file_name).write_text(output)
    return file_name


def section_intro(dim: str) -> list[str]:
    roots = ", ".join(
        rf"\path{{{result_root(dim, method)}}}" for method in METHODS
    )
    return [
        rf"\subsection{{{'Two' if dim == '2D' else 'Three'}-dimensional manufactured problem}}",
        "The tables in this subsection are generated directly from the current",
        rf"\path{{{RESULT_ROOT}}} result tree. The method-specific roots are {roots}.",
        r"All entries use \code{crankNicolson}, the \code{consistent} mass matrix, and \code{RKF45} for the state ODEs.",
        r"The first column reports \code{Vm-states-Iion}: the $V_m$ reconstruction degree, the state treatment (\code{CCp*}, \code{GP}, or \code{na}), and the $\Iion$ quadrature/reconstruction setting.",
        f"Fits use {range_label(dim, 'all')} for the broad mesh range and {range_label(dim, 'fine')} for the fine-mesh range.",
        "",
    ]


def write_section_fragment(dim: str, table_files: dict[tuple[str, str], str]) -> None:
    lines = section_intro(dim)
    for method in METHODS:
        lines.append(rf"\subsubsection{{{method_label(method)}}}")
        for dt_tag in DT_TAGS:
            lines.append(rf"\paragraph{{$\Delta t={dt_to_latex(dt_tag)}$}}")
            lines.append(rf"\input{{tables/{table_files[(method, dt_tag)]}}}")
        lines.append("")
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
