# Niederer Benchmark Activation-Time Comparison

## Objective

This report documents the activation-time comparison workflow for the Niederer
slab benchmark in this repository. The goal is to compare activation-time
results across several simulation configurations and to track how the numerical
choices affect agreement with the reference benchmark.

The benchmark data and plotting assets used here come from the local tutorial
tree:

- `tutorials/NiedererEtAl2012`
- `tutorials/highOrderNiedererEtAl2012_p3`
- `tutorials/highOrderNiedererEtAl2012_p3/setupNiedererEtAl2012`

## Case Selection

The repository contains two closely related benchmark tutorials:

- `tutorials/NiedererEtAl2012`
- `tutorials/highOrderNiedererEtAl2012_p3`

For activation-time comparison, the high-order tutorial is the natural sweep
target because it already includes:

- an automated parameter sweep script
- ParaView-based line and point export
- plotting scripts that compare `activationTime` against reference data

The high-order setup code is intended to vary:

- `Δx = 0.5, 0.2, 0.1 mm`
- `Δt = 0.05, 0.01, 0.005 ms` in the configuration file, with local overrides
- `tissue = epicardialCells`
- `ionicModel = TNNP`
- `solver = explicit` in the current local configuration

The checked-in exported data currently covers only the `explicit` runs with
`Δt = 0.005 ms` for `Δx = 0.5 mm` and `Δx = 0.2 mm` under
`testSolvers/explicit_2nd_order/`.

The standard tutorial setup differs slightly: its config still enables both
`explicit` and `implicit`, but that path is not currently usable in this
workspace because the case does not contain a `0/` initial-condition directory.

## How To Run

The documented pipeline is:

1. Update the case dictionaries for mesh, time step, tissue model, ionic model,
   and solver choice.
2. Run the tutorial case via the helper shell script.
3. Touch `case.foam` for ParaView.
4. Export line and point activation data with ParaView.
5. Compare the exported CSVs against the reference benchmark plots.

The relevant entry points are:

- `tutorials/highOrderNiedererEtAl2012_p3/Allrun`
- `tutorials/highOrderNiedererEtAl2012_p3/setupNiedererEtAl2012/simulation/run_cases.sh`
- `tutorials/highOrderNiedererEtAl2012_p3/setupNiedererEtAl2012/main_setup_Niederer.py`
- `tutorials/highOrderNiedererEtAl2012_p3/setupNiedererEtAl2012/simulation/setup_multiple_simulations_Niederer.py`

In the current workspace, a fresh end-to-end rerun requires an OpenFOAM
environment and ParaView tooling that are not loaded in this shell by default.
The run agent verified that the required tools become available after sourcing
`/home/pablo/OpenFOAM-v2312/etc/bashrc`, and that `pvpython` resolves to
`/usr/bin/pvpython`.

## Post-Processing Pipeline

The post-processing flow is based on two ParaView exports and two plotting
scripts:

- `exportParaview_Niederer.py`
  - extracts activation data from the final timestep
  - writes line data through `PlotOverLine`
  - writes point data through `ResampleWithDataset`
  - produces CSV files named like:
    - `*_line_DT<dt>_DX<dx>.csv`
    - `*_points_DT<dt>_DX<dx>.csv`

- `line_postProcessing.py`
  - reads the line CSVs
  - removes nonessential columns
  - overlays cardiacFoam curves with the digitized reference data
  - supports a toggle between all runs and a single `ΔT` comparison view

- `points_postProcessing.py`
  - reads the point CSVs
  - builds 3D activation-time surfaces over the `(Δx, Δt)` sweep
  - identifies the earliest-activated point and omits it from the subplot grid

The point and line plots expect `activationTime` in the exported CSVs, which is
consistent with the ParaView export script in this tutorial.

The line post-processing script also writes two HTML artifacts in its working
directory:

- `cardiacFoam_allSimulations.html`
- `Niederer_vs_cardiacFoam.html`

The main comparison script for benchmark matching is `line_postProcessing.py`.
The point-surface script is useful for convergence visualization, but it does
not compare against the digitized Niederer reference data directly.

## Results

### Simulation Matrix

| Case | Solver | Ionic model | Tissue | `Δx` (mm) | `Δt` (ms) | Case folder | Run status |
| --- | --- | --- | --- | ---: | ---: | --- | --- |
| 1 | explicit | TNNP | epicardialCells | 0.5 | 0.005 | `testSolvers/explicit_2nd_order` | CSVs present in repo |
| 2 | explicit | TNNP | epicardialCells | 0.2 | 0.005 | `testSolvers/explicit_2nd_order` | CSVs present in repo |
| 3 | explicit | TNNP | epicardialCells | 0.1 | 0.005 | `testSolvers/explicit_2nd_order` | no checked-in CSV found |

### Activation-Time Summary

| Case | Reference source | Metric | Value | Units | Notes |
| --- | --- | --- | --- | --- | --- |
| 1 | exported point CSV | Earliest sampled activation time | `1.2495` | ms | from `explicit_TNNP_epicardialCells_points_DT005_DX05.csv` |
| 2 | exported point CSV | Earliest sampled activation time | `1.2484` | ms | from `explicit_TNNP_epicardialCells_points_DT005_DX02.csv` |
| 3 | exported point CSV | Earliest sampled activation time | `TODO` | ms | no `DX01` point CSV checked into the tutorial |

### Additional Observed Metrics

| Case | Metric | Value | Units | Source |
| --- | --- | --- | --- | --- |
| 1 | Latest sampled point activation time | `146.09` | ms | `explicit_TNNP_epicardialCells_points_DT005_DX05.csv` |
| 2 | Latest sampled point activation time | `57.723` | ms | `explicit_TNNP_epicardialCells_points_DT005_DX02.csv` |
| 1 | Diagonal line endpoint activation time | `146.09` | ms | `explicit_TNNP_epicardialCells_line_DT005_DX05.csv` |
| 2 | Diagonal line endpoint activation time | `57.723` | ms | `explicit_TNNP_epicardialCells_line_DT005_DX02.csv` |

### Output Artifacts

| Artifact | Expected location | Status |
| --- | --- | --- |
| Line CSV exports | `tutorials/highOrderNiedererEtAl2012_p3/testSolvers/explicit_2nd_order/` | present for `DX05` and `DX02` at `DT005` |
| Point CSV exports | `tutorials/highOrderNiedererEtAl2012_p3/testSolvers/explicit_2nd_order/` | present for `DX05` and `DX02` at `DT005` |
| Interactive comparison HTML | `tutorials/highOrderNiedererEtAl2012_p3/testSolvers/explicit_2nd_order/` | present as `127.0.0.1.html` and `3D_127.0.0.1.html` |
| Generated report | `tutorials/highOrderNiedererEtAl2012_p3/activation_time_comparison_report.md` | updated in this pass |

## Notes And Limitations

- The current local configuration contains hard-coded path assumptions for the
  setup scripts and should be checked before running on a new machine.
- `main_setup_Niederer.py` in the high-order tutorial refers to
  `setupNiedererEtAl2011`, while the actual folder in this repo is
  `setupNiedererEtAl2012`; that path mismatch should be fixed before using the
  CLI entry point for a full sweep.
- The high-order `line_postProcessing.py` currently uses `dt_target = 0.5`,
  which is inconsistent with both the active config (`0.05 ms`) and the
  checked-in exported CSVs (`0.005 ms`).
- `run_cases.sh` expects `$HOME/OpenFOAM-v2312/etc/bashrc`, and the current
  shell did not have `WM_PROJECT_DIR`, `blockMesh`, or the cardiacFoam solvers
  available by default when this report was assembled, even though those tools
  were confirmed to exist once the OpenFOAM bashrc was sourced.
- The standard tutorial is incomplete in this checkout because
  `tutorials/NiedererEtAl2012/0/` is missing, so its top-level `Allrun` cannot
  currently reproduce the benchmark from a clean start here.
- The checked-in CSV data is incomplete for the nominal three-resolution sweep:
  `DX=0.1 mm` outputs were not found under `testSolvers/explicit_2nd_order/`.
- If the sweep is expanded beyond the current `explicit` configuration, the
  comparison table should be updated to include the new solver family and
  parameter combinations.

## Next Fill-In Step

After the simulations are run, replace the remaining `TODO` entries above with:

- the actual case identifier
- the exported CSV filenames
- the activation-time metric values
- a short qualitative comparison against the Niederer reference curves
