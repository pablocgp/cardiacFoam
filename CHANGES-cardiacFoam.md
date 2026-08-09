# Port from LRE to solids4foam's movingLeastSquares

Log of every change made while moving the high-order cardiac solvers off the old
serial `LRE` library and onto `movingLeastSquares` from the solids4foam branch
`highOrder-nonLinGeomTotalLagTotalDispSolid-pablo`, so that they can run in MPI
parallel.

Why: `LRE` is serial by construction. `makeGlobalCellStencils` returns before the
(commented-out) MPI exchange and carries the note *"Implemented for serial run"*;
processor patches raise `NotImplemented`. `movingLeastSquares` stores global cell
IDs in its stencils and every evaluation entry point performs its own halo
exchange, so a caller needs no changes to be parallel-correct.

Convention used throughout: new code is marked `// Added for cardiacFoam: <why>`,
modified code `// Modified for cardiacFoam: <what changed and why>`.

---

## A. Changes inside solids4foam

This is the **complete** diff on the solids4foam tree. It is purely additive and
cannot change the behaviour of any existing caller.

### A.1 Three public accessors on `movingLeastSquares`

`src/solids4FoamModels/higherOrderHelpers/movingLeastSquares/movingLeastSquares.H`

Added `cellGradCoeffsData()`, `cellSecondGradCoeffsData()` and
`cellThirdGradCoeffsData()` next to the existing public forwarders
`stencilData()` and `faceGradCoeffsData()`. They expose the cell-centre
coefficient tables, which were private with no public route to them.

An external solver that assembles its own high-order operators (rather than going
through `hofvm`) needs these directly to build its matrix rows. The doc comment
also records the row layout, which was not written down anywhere: row `cellI`
follows `stencilData().cellsStencil()[cellI]` and carries one extra trailing
entry, at index `cellsStencil()[cellI].size()`, holding the cell's own
coefficient.

Verified non-breaking: `tutorials/solids/linearElasticity/cantilever2d
./Allrun highOrder` produces a solver log identical to the pre-change run except
for the wall-clock timestamp, the PID and the execution time. Same error norms
(`DDifference` 2.25762e-12 / 2.91684e-12 / 6.04668e-12), same 73 and 86 linear
iterations, same 2 nonlinear iterations.

---

## B. Bugs found in the solids4foam branch (not fixed here)

Both are independent of this port and are **not** needed for it. They are listed
so they can be raised separately.

1. `higherOrderHelpers/movingLeastSquares/movingLeastSquaresTemplates.C:166` —
   the `autoPtr`-returning `fGrad` overload calls `faceGrad(vf, ...)`, a member
   that does not exist anywhere in the repository. It is a compile error the
   moment that overload is instantiated; nothing instantiates it today. Callers
   should use the two-argument `fGrad` form, which is what this port does.

2. `higherOrderHelpers/movingLeastSquaresStencil/movingLeastSquaresStencil.C:1312`
   — `clear()` contains `remoteCentresMapPtr_();`, which *dereferences* the
   `autoPtr` instead of `.clear()`-ing it, leaving a stale remote-centres map
   after a topology change. Masked in the normal path because
   `movingLeastSquares::clear()` destroys the whole stencil.

---

## C. Changes inside cardiacFoam

### C.1 New: the compatibility adapter

`src/highOrderAdapter/highOrderInterp.H` (new, header-only)

Wraps a `movingLeastSquares` by composition and exposes the twelve `LRE` methods
the cardiac solvers actually used, so the ~7000-line solvers move across at one
place instead of at every call site. Deliberately temporary: once both cardiac
solvers are verified, the call sites can point straight at `movingLeastSquares`
and this header can be deleted.

Two design points worth knowing:

**The quadrature weight accessors are deliberately renamed.** `LRE` normalised
its weights (face weights summed to 1, cell weights to 1) so callers multiplied
by `|Sf|` or the cell volume themselves. `fvMeshQuadrature` returns PHYSICAL
weights: face weights sum to `|Sf|`, cell weights to the cell volume. Verified to
machine precision on both a 3-D tetrahedral mesh (error 0) and a 2-D mesh with
empty patches (3.6e-16 on faces, 4.5e-15 on cells). Keeping the old names would
have let the old arithmetic compile and silently rescale every flux by `|Sf|` —
an error that changes the answer without changing the convergence rate, so a
convergence study would not catch it. The adapter therefore exposes
`faceQuadWeightPhysical()` / `cellQuadWeightPhysical()`, turning each affected
site into a compile error.

**`mirrorPatches` is accepted and ignored, with a warning.** `LRE` could treat a
zeroGradient patch by mirroring the whole face stencil; `movingLeastSquares` has
the mirroring machinery (it uses it for symmetry patches) but no equivalent
opt-in. The cardiac solvers do not need one: they impose a homogeneous Neumann
condition by setting the boundary flux to zero, which is the correct
finite-volume treatment, and for the manufactured solution
`cos(pi x) cos(2 pi y) cos(3 pi z)` on the unit cube the normal derivative is
identically zero on all six faces, so it is exact.

Dictionary translation, all in one place: `N` -> `polynomialOrder`,
`Nn` -> `faceStencilExtraCells` and `cellStencilExtraCells` (same meaning in both
libraries: points added on top of the minimum for the polynomial order),
`weightFunction` -> `weightFunctionCoeffs/type`, `k` -> `weightFunctionCoeffs/k`.
Dropped because the new library has no equivalent and needs none:
`maxStencilSize` and `nLayers` (only sized a `DynamicList` in LRE's dead legacy
path), `useQRDecomposition` and `useGlobalStencils` (selected between LRE code
paths of which only one worked).

### C.2 `highOrderManufacturedFDAImplicitPETScParallel`

Mechanical renames: `LRE::symmTensor3Order` -> `highOrderInterp::` (16 sites),
`const LRE&` -> `const highOrderInterp&` (18), `LRE::cubicForm` -> `highOrderInterp::`
(5), object construction (3), `gradScalarFaceQuad` return type (1).

**The `|Sf|` factor removed at five sites**, per C.1: two in the stiffness matrix
(internal and boundary faces), one in the Dirichlet boundary vector, two in the
explicit flux path.

**`labelListList` -> `CompactListList<label>`** for the stencil accessors (5
declarations), following the new return type.

**`const labelList&` -> `const UList<label>` for the five `curStencil` bindings.**
This one is the important one, see D below.

`Make/files`: distinct binary name so the parallel port does not overwrite the
serial `highOrderManufacturedFDAImplicitPETSc`, which is kept as the convergence
reference.

`Make/options`:
- PETSc switched from the hardcoded system `petsc3.15` to `$(PETSC_DIR)/$(PETSC_ARCH)`,
  which is PETSc 3.23.2. `libsolids4FoamModels.so` links 3.23; having both in one
  process is an ABI hazard. Note the library there is `-lpetsc`, not `-lpetsc_real`.
  MPI is consistent: PETSc 3.23 and OpenFOAM's `sys-openmpi` both use the system
  `libmpi.so.40`.
- Eigen taken from the `solids4foam-ib-pc` tree actually being built against.
- **`sinclude $(SOLIDS4FOAM_DIR)/etc/wmake-options` added.** `$(VERSION_SPECIFIC_INC)`
  was already referenced but never defined, so it expanded to nothing. `LRE.H` did
  not care; `symmTensor3rdOrder.H` does, and without `-DOPENFOAM_COM` it takes the
  foam-extend branch and fails to compile. The same fix is needed in every
  cardiacFoam solver that includes the new headers.

**New permanent consistency check** after the stiffness assembly: every row of `K`
must sum to zero, because a constant field has zero gradient and the boundaries
carry no flux. Cheap (one pass over the non-zeros), globally reduced, and it
aborts if the matrix comes out empty. It is what caught the bug in D.

### C.3 `testHighOrderScalarGrad`

Ported to the adapter and given a permanent assertion that face quadrature
weights sum to `|Sf|` and cell weights to the cell volume, so a future change to
the weight convention upstream fails loudly here instead of silently rescaling a
flux somewhere else.

---

## D. The bug this port had to find, and how

Symptom: with the high-order path enabled the manufactured-solution error was
~100x too large **and flat** — 2.566e-2 at N=12 and 2.604e-2 at N=13, i.e. not
converging at all. A fixed error that does not reduce under refinement means the
operator is inconsistent, not merely inaccurate.

Cause:

```cpp
const labelList& curStencil = faceStencils[faceI];   // silently empty
```

`CompactListList<T>::operator[]` returns `const UList<T>` **by value** — a view,
not a `List`. `List<T>` derives from `UList<T>`, so binding the view to a
`const List<T>&` is a conversion the compiler accepts without any warning, and
the result is an **empty** list. `forAll(curStencil, cI)` therefore never
iterated, no triplets were generated, and the stiffness matrix was assembled with
**zero entries**: the solver was integrating the equation with no diffusion term
at all. `LRE` returned `labelListList`, i.e. `List<List<label>>`, so the binding
was natural and the problem did not exist before.

Five sites were affected: three in the mass and stabilisation assembly, two in
the stiffness assembly. Binding to `const UList<label>` fixes it and also avoids
a per-cell and per-face copy.

Worth recording because none of the ingredients was wrong. The bisection that
found it:

| step | result |
|---|---|
| low-order path only | exact, +0.0000% |
| high-order Iion/states only | exact, +0.0000% |
| high-order Vm only | fails |
| lumped instead of consistent mass | still fails -> stiffness, not mass |
| hexahedral instead of triangular | still fails -> not topology |
| 2-D quadrature weight assertion | 3.6e-16 -> not the `|Sf|` change |
| face gradients vs LRE baseline | +0.17% at p3 -> coefficients fine |
| **row sums and non-zeros of K** | **nnz = 0** |
| stencil size inside vs outside the loop | 0 vs 30 -> the binding |

Everything up to the second-to-last row confirmed correct ingredients. Only
measuring the assembled matrix showed the problem was in none of them. This is
the class of error worth adding a permanent check for, which is why C.2 keeps it.

**This bug class will recur when porting `highOrderElectroActivationFoamImplicitPETSc`:**
grep for `const labelList&` and `const List<` bound to anything returned by a
`CompactListList`.

---

## E. Verification so far

| what | result |
|---|---|
| solids4foam builds, `cantilever2d highOrder` | identical to pre-change |
| `testHighOrderScalarGrad` vs LRE, same mesh | p1 -3.6%, p2 +0.3%, p3 -1.7% |
| its convergence orders, two meshes | LRE 1.157/2.025/3.180 vs 1.135/2.026/3.172 |
| the same app, 1 rank vs 4 ranks | reduced norms identical to all printed digits |
| quadrature weights, 3-D and 2-D | 0 and 3.6e-16 |
| MMS `(NO,na,NO)` | +0.0000% |
| MMS `(p3,CCp3,p1)` triangular, alpha=0.1 | Vm_L2 -0.13%, order 2.129 vs 2.130 |

Serial reference for the MMS is captured in `/home/pablo/mms-baseline-serial/`
(208 configurations with fitted orders, plus the capture script), taken from the
archived results of the serial binary, which is still installed and still runs.

---

## F. Distributed linear algebra

Done in a separate solver, `highOrderManufacturedFDAImplicitPETScDistributed`, a
copy of the verified `...Parallel` taken once the port above was validated. The
two are kept side by side so the algebra change can be compared against a
known-good baseline rather than against memory.

The assembly loops were not touched. Instead the matrices are sized
`nLocalCells x nGlobalCells` - local rows, global columns - which is exactly the
MPIAIJ layout, so the existing Eigen assembly and the whole theta-scheme
arithmetic carry over unchanged. Two things made this cheap:

- **The columns were already global.** `movingLeastSquares` stencils hold global
  cell IDs, so `col = curStencil[cI]` needed nothing. Only the ~6 self/diagonal
  columns, which used the local index, needed `toGlobal()`.
- **Processor patches were already handled, by omission.** The boundary loop
  skips only `empty` and `zeroGradient`; anything else, including `processor`,
  falls through to the generic case, which contributes to the owner row using
  the face stencil that crosses the partition. That is the same pattern `hofvm`
  uses. No new branch was needed.

Changes: `MatCreateSeqAIJ`/`PETSC_COMM_SELF` -> `MatCreate` + `MatSetSizes` +
`MATAIJ` on `PETSC_COMM_WORLD` with `d_nnz`/`o_nnz` preallocation, factored into
`buildDistributedMat()`; `VecCreateSeq` -> `VecCreateMPI`; `KSPCreate` and
`MatCreateShell` on `COMM_WORLD`; five Eigen matrix-vector products replaced by
`MatMult` through a small `DistributedMatVec` helper; global reduction helpers
`gNorm`/`gDot`/`gSquaredNorm` applied to `relativeL2Norm`, the JFNK
finite-difference epsilon and the Armijo line search; `characteristicDx` moved
from `boundBox(mesh.points())` (local) to `mesh.bounds()` (reduced).

The row offset is a file-scope `gRowStart` set once in `main()` from
`globalIndex::localStart()`, rather than a parameter threaded through four
linear-solver signatures. It is taken from OpenFOAM rather than queried from
PETSc: the two distributions agree today, but relying on that would be a silent
trap if either changed.

Two things only surfaced when actually running in parallel:

1. **`ilu`/`lu` do not exist for MPIAIJ.** `KSPSetUp` fails with
   `PETSC_ERR_SUP` (92). They have to be wrapped in block-Jacobi, applying the
   requested factorisation to each rank's diagonal block. Handled by
   `petscParallelPcTypeName()`, which leaves serial behaviour untouched.

2. **`if (!Pstream::master()) return;` before writing a file deadlocks.**
   OpenFOAM's file handler can be collective, so the master enters a file
   operation the other ranks never reach. The symptom is nasty: the time loop
   completes, every nonlinear step converges, the timing block prints, and then
   the run simply stops - it looks like a solver problem, not I/O. The fix is
   that every rank constructs the stream and only the destination differs:
   master writes the canonical file to the global case directory, the others
   write a throwaway copy into their own `processorN/`. Causality was confirmed
   by reverting the guard and watching the run complete, not by inspection.

### Verification

Each step was arranged so that in serial it is a no-op (`toGlobal` is the
identity, `nGlobal == nLocal`, one rank), giving a check where any single
differing digit is a bug:

| after | serial result |
|---|---|
| global columns and `nLocal x nGlobal` sizing | 20/20 points identical bit for bit |
| MPIAIJ on `COMM_WORLD` | 20/20 identical |
| `MatMult` for all products | 20/20 identical |
| reductions, `characteristicDx`, output, block-Jacobi | 20/20 identical |

In parallel, on the full high-order configuration `(p3, CCp3, p3)`, 2-D
unstructured triangular, `N=20`, Picard, consistent mass, `dt=1e-3`, 50 steps:

| | 1 rank | 2 ranks | difference |
|---|---|---|---|
| Vm L1, alpha=0 | 6.18373376152e-04 | 6.18373356393e-04 | 3.2e-6 % |
| Vm L2, alpha=0 | 8.49348982299e-04 | 8.49348986105e-04 | 4.5e-7 % |
| Vm Linf, alpha=0 | 3.79586044420e-03 | 3.79586040315e-03 | 1.1e-6 % |
| Vm L1, alpha=0.1 | 5.66787570843e-05 | 5.68110947349e-05 | **0.23 %** |
| Vm L2, alpha=0.1 | 1.02623488473e-04 | 1.02724578089e-04 | 0.10 % |

The stiffness consistency check gives `max |row sum| = 1.39e-13` in both, i.e.
the operator annihilates constants **across the partition**. That was the one
genuine unknown of the plan - whether face stencils agree on both sides of a
processor face - and it is now answered: they do, or this would not be zero.

The alpha=0 agreement to seven or eight significant figures shows the
distributed algebra is exact to round-off. The remaining 0.23 % at alpha=0.1 is
therefore entirely the stabilisation, which is the one piece still missing (see
below). Running both alpha values is what isolates it; either alone would have
been ambiguous.

---

## G. Acceptance: convergence order under decomposition

The driver `run_convergence.py` gained `N_PROCS` and `DECOMPOSE_METHOD`: at 1 it
behaves exactly as before, above 1 it writes a `decomposeParDict`, runs
`decomposePar -force` and launches under `mpirun -np N ... -parallel`.
`clean_case()` now also removes `processor*`.

Sweep on `(p3, CCp3, p3)`, 2-D unstructured triangular, `alpha = 0.1`, Picard,
`dt = 1e-3`, N = 12, 16, 20, 24:

| N | 1 rank | 2 ranks | 4 ranks | 8 ranks |
|---|---|---|---|---|
| 12 | 1.299435e-04 | 1.299445e-04 | 1.300332e-04 | 1.302369e-04 |
| 16 | 4.924606e-05 | 4.925125e-05 | 4.928857e-05 | 4.931371e-05 |
| 20 | 2.238711e-05 | 2.239199e-05 | 2.239449e-05 | 2.242872e-05 |
| 24 | 1.066183e-05 | 1.066348e-05 | 1.066357e-05 | 1.067235e-05 |
| **fitted order** | **3.5855** | **3.5853** | **3.5864** | **3.5870** |

Maximum deviation in the fitted order: **0.0014**, against the 0.05 acceptance
tolerance. The raw errors agree to between 0.001 % and 0.09 %.

At 8 ranks these meshes leave roughly 100 cells per rank, which is where the
plan expected the p3 stencil search to run out of halo depth
(`haloDepthScale` defaults to `polynomialOrder*2.5`). It did not: no widening
was needed.

**JFNK also runs in parallel.** Same configuration with
`nonlinearMethod JFNK` and `jfnkPreconditioner diagonalIion`:

| N | 1 rank | 2 ranks | difference |
|---|---|---|---|
| 12 | 1.299434e-04 | 1.299445e-04 | 0.0009 % |
| 16 | 4.924602e-05 | 4.925122e-05 | 0.0106 % |

This is the path that depends on the reduced finite-difference epsilon and the
reduced Armijo line search; without those two the step size would have been
rank-dependent.

### Why the stabilisation gap turned out not to matter

> Superseded: the exchange sketched at the end of this subsection has now been
> implemented for the MMS solver. See section J.


`addCellTaylorExtrapolationCoeffs(row = own, cellI = nei)` reads the stencil and
the three coefficient tables of `nei`. On a processor-boundary face that cell
belongs to another rank, so the stabilisation term is dropped in the band along
the partition. On a single mesh this is visible: 0.23 % on the Vm L1 error at
`alpha = 0.1`, against 3e-6 % at `alpha = 0`, which is how the gap was isolated.

It does not change the convergence order, because the affected faces are a
surface within a volume: their share of the domain shrinks under refinement
faster than the discretisation error does. The planned fix - a one-off
`PstreamBuffers` exchange of cell centres, stencils and coefficient rows across
each `processorPolyPatch` - was therefore **not implemented**. Two supporting
measurements: in serial, switching the stabilisation off changes the fitted order
of the high-order configurations by at most 0.024 (4.323 -> 4.299 on hexa,
4.477 -> 4.460 on triangular), so the term is a small correction to begin with;
and the 1/2/4-rank agreement above is three orders of magnitude inside tolerance.

Worth writing down as a decision rather than an omission: if a future
configuration relies more heavily on the stabilisation - the standard low-order
operator on hexahedra is the obvious candidate, where alpha changes the order
from 2.000 to 2.294 - this needs revisiting, and the exchange sketched above is
the way to do it.

---

## I. The electro solver

`highOrderElectroActivationFoamImplicitPETScParallel` carries the port and
global matrix columns; `...Distributed` adds the distributed algebra and the
parallel-correctness fixes on top. The first is deliberately frozen once
verified, so every later stage has a serial-equivalent binary to compare with.

### I.1 What differed from the MMS

The LRE API surface is identical, so `highOrderInterp.H` needed no changes at
all. The mechanical port was 19 substitutions, each asserted on its expected
count before anything was written: 5 `curStencil` bindings to `const UList<label>`,
4 `|Sf|` factors dropped, 13 `const LRE&` signatures, and so on.

Two things the plan predicted did **not** hold, and both were worth the time
spent finding out:

**The Dirichlet gap does not exist.** No case in the benchmark uses a Dirichlet
condition on `Vm`; every tutorial and every branch of the driver's case writer
is `zeroGradient` (plus `empty` in 2-D), and `fixedVoltage` is not a class in
this repository. The scalar `patchFaceQuadValues` anticipated in section H is
not needed.

**The reference binary could not be rebuilt.** Rebuilding solids4foam-ib-pc
replaced `libsolids4FoamModels.so` in the shared platform directory, and the
serial electro binary needs 17 LRE symbols that library no longer has - it dies
at load with `undefined symbol: _ZTIN4Foam3LREE`. Relinking the old branch's
objects into a private directory got 95% of the way there and then hit a wall:
the LRE the binary was built against (`~/Desktop/LRE`, May) includes
`fixedVoltageFvPatchScalarField.H`, a class that exists nowhere on the system.

That turned out not to matter, because the archived study under
`tutorials_electro/highOrderNiedererEtAl2012ImplicitPETSc/results/` already IS
the pre-port reference, and a better one than a fresh capture: it stores every
input dictionary (`blockMeshDict.used`, `controlDict.used`, ...) beside a
`nonlinearResiduals.dat` of 29 columns at 12 significant digits, produced by the
reference binary itself.

### I.2 Two references, and why neither is bit-exact

**R0** is that archive. The `(NO,na,NO)` configuration - the only one that builds
no interpolator, so its result cannot depend on LRE vs MLS - reproduces it with
`step`, `iter`, `linearIterations` and `converged` **bit-identical** over 166
rows, `coupledResidual` agreeing to 1e-14 absolute on values of 4.7e-3, and no
state residual off by more than 1e-16. It is not byte-identical, for a reason
unrelated to the port: the reference binary linked PETSc 3.15 and this one links
3.23, so the Krylov solve takes a different rounding path.

**R1** (`electro-baseline-serial/capture_R1.py`, 12 configurations, reproducible
with zero differences on a second pass) is the post-port serial reference. It is
10 time steps, not 50: a deviation introduced by the distributed algebra shows up
in step 1, so what matters is the WIDTH of the fingerprint, not its length. The
whole sweep takes under two minutes, which matters because it is re-run after
every stage.

The distributed algebra cannot be bit-exact against R1 either, and it is worth
being explicit about why: replacing Eigen's sparse product with `MatMult` changes
the order in which each row is summed, and floating-point addition is not
associative. The deviations are 1e-12 to 1e-14 absolute - the last digit of the
output format. The criterion that does discriminate, and which all 12
configurations meet, is that the TRAJECTORY is identical: steps, Picard
iterations, linear iterations and the converged flag all match exactly.

### I.3 The serial bug fixed on the way in

`solveODE` addressed `algebraic_`, `rates_` and `stepMs_` by position in the
array it was handed, while the front-aware hybrid calls it twice per iteration -
once with `nF` compacted front Gauss points and once with `nS` compacted smooth
cells, both indexed from 0. Two different physical point sets sharing one bank.

The fix gives the scratch two disjoint banks (Gauss slots `[0, nPoints_)`, cell
slots `[nPoints_, nPoints_ + nCells)`) and an explicit slot map argument that
defaults to empty, so the other five call sites are untouched.

The plan predicted this would be byte-identical. **It is not**, and the reason is
more interesting than the prediction. Isolating it took three steps: a double
build showed the hybrid differing at ~4e-9 in the state residuals with
`coupledResidual` and `Vm_relL2` identical; a control build (scratch remapped,
`pointI` still passed as the logging ID) reproduced the fixed version exactly,
ruling out the logging argument; and instrumenting `stepMs_` showed what actually
happens. On the FIRST call every slot holds the constructor value 0.1 ms; from
the second on, every slot holds exactly `hMin = 1e-13`. So with the fix the value
is uniform. WITHOUT it, the front pass writes slots `0..nF-1` first, and the
smooth pass then reads a MIXTURE - 176 slots at 1e-13 and 442 still at 0.1 - so
the RKF45 sub-step sequence, and with it the rounding, depended on how many front
points there happened to be. Three orders of magnitude below the ODE's own
`relTol` of 1e-6, and the fixed behaviour is the defensible one.

Two earlier hypotheses of mine were wrong on the way here - FSAL carry of
`rates_` (`rkf45Trial` computes `k1` fresh and only writes `rates = k6`), and
"`stepMs` is uniform, therefore the aliasing swaps equal values" (true from the
second call onwards, false on the first). Both were settled by instrumenting
rather than by reading.

### I.4 Parallel correctness

Beyond the MMS's list:

**The orthogonal operator dropped every processor face.** Its boundary loop skips
`empty`/`zeroGradient` and then tests `fixedValue || fixedVoltage` with no
`else`, so a `processor` patch contributed nothing and each partition's block was
completely uncoupled. Confirmed by a differential test rather than by reading:
`(p3,GP,p3)`, which goes through `assembleHighOrderStiffnessMatrix` and handles
coupled patches in its generic case, ran on 2 ranks; `(NO,na,NO)` failed with
`KSP_DIVERGED_BREAKDOWN`.

The fix adds an `else if (patch.coupled())` branch using
`syncTools::swapBoundaryCellList` for the neighbouring conductivity and global
cell ID - not `conductivity.boundaryField()`, which on a processor patch holds
the interpolated FACE value rather than the neighbouring CELL value - and
`patch.delta()` for the cell-to-cell vector. Only the owner side is added; the
rank across the face adds its own.

The Rhie-Chow stabilisation jump is NOT added on processor faces: it needs the
neighbour's reconstruction stencil, which does not cross the partition. Inert at
the default `stabilisationAlpha = 0`; documented as an omission rather than
silently skipped.

**`nearestCellToPoint` searched locally**, so every rank believed it owned P8,
clipped `effectiveEndTime` with its own cell, and left the time loop at a
different step - a hang, preceded by corrupted physics since `effectiveEndTime`
also sets the last step's `dt`. Exactly one rank now returns a cell (min-reduce
on the distance, ties broken by lowest rank), and the activation time is read
with a max-reduction, which makes `stopPointTriggered` and `effectiveEndTime`
uniform by construction.

**`sampleIDW` took the k nearest LOCAL cells.** Each rank now builds its local
top-k as `(d2, globalCellID, value)`, the candidates are gathered, and every rank
merges the same set ordered by `(d2, globalCellID)`. That tie-break reproduces
serial exactly: the serial loop tests `d2 < bestD2` strictly, so among equal
distances it keeps the lowest index. The exact-hit early return moved AFTER the
gather - inside the search loop it would have had one rank return while the
others entered the collective.

**`dilateFrontCells` and the limiter stopped at the cut.** Both now exchange
through a `volScalarField` + `correctBoundaryConditions()`. Worth recording why
that works when it looks like it should not: constructing the field with
`zeroGradientFvPatchScalarField::typeName` does NOT put zeroGradient on the
processor patches, because `fvPatchField::New` overrides the requested type with
the patch's own constraint type - so they receive a real
`processorFvPatchScalarField` and the exchange happens. One exchange per dilation
layer.

**File output.** Six sites. Five used `runTime.path()`, which in parallel is
`<case>/processorN`, so results landed where the driver does not look; the sixth
was a bare relative filename that every rank opened and truncated. Every rank
constructs the stream and only the destination differs - `if (!Pstream::master())
return;` before opening deadlocks, as documented in section F.2.

### I.5 The bug the fingerprint could not see

`reportOperatorFingerprint` prints globally reduced `nnzGlobal`, `sum`, `sumAbs`
and `max|rowSum|` after each assembly. `nnzGlobal` matches the serial value
exactly at 1, 2, 4 and 8 ranks (3088 for the orthogonal operator, 23428 for p3,
22480 for the consistent mass), and `sumAbs` to 12 significant digits. That also
closes the halo-depth question: p3 at 8 ranks on 640 cells, ~80 cells per rank,
reproduces the serial pattern.

It is a strong check and it still missed a real bug. `assembleDiagonalMassMatrix`
wrote `triplets.emplace_back(cellI, cellI, coefficient)` - a LOCAL column. It
escaped the global-column sweep because it builds its triplet directly instead of
going through `addTripletIfNeeded`.

Putting an entry in the wrong COLUMN changes neither `nnz`, nor `sum`, nor
`sumAbs`, nor any row sum. The fingerprint of `M` read 640/640/640 in serial and
in parallel while being completely wrong. It only surfaced on combining `M` with
`L`, where the displaced diagonal no longer coincided with `L`'s: `AImplicit` had
3408 nonzeros at 2 ranks against 3088 in serial, exactly the per-rank row count
in excess.

The lesson is the same one as for the dropped processor faces: an operator
fingerprint is blind to anything that preserves the multiset of values. It has to
be paired with a solution comparison. Finding this took looking at a run that
appeared to hang, noticing it was actually advancing with Picard diverging, and
only then realising that A1 had been reported without A2 ever being run.

After the fix, the trajectory agrees:

| | coupledResidual | Vm_relL2 | Iion_relL2 |
|---|---|---|---|
| np=1 vs 2 | 1.29e-09 | 1.32e-11 | 1.76e-11 |
| np=1 vs 4 | 1.07e-09 | 1.32e-11 | 1.47e-11 |

and on the hybrid p3 configuration over 60 steps, 1.4e-09 at 2 ranks and 1.6e-08
at 4. The stop-point run gives an identical `activationStatus` at 1, 2 and 4
ranks: `OK 0.00175195581081 9 1`.

### I.6 A test that passed without testing anything

The front-cell probe is `returnReduce(sum(frontCell))` per step, an integer, so
it must match exactly. On the 10-step case it matched with AND without the
exchange - the front never reached a cut, so the test was vacuous. Extending to
60 steps produced the signal: without the exchange, 4 ranks gave
`75 76 77 77 81 81` against the serial `83 84 86 86 85 85`, the deficit being the
cells adjacent to the cut. With it, all three rank counts agree exactly.

A test that passes without exercising its condition proves nothing, and the only
way to notice is to break the fix on purpose and check the test fails.

### I.7 Results

**Speedup, fixed workload.** 640 cells, 200 steps, `endTime` fixed rather than
`stopAfterPointActivation` - a speedup measurement needs identical work at every
rank count, and stopping on activation could end different rank counts at
different steps. Sequential runs on an otherwise idle 14-physical-core machine.

| ranks | (NO,NO,NO) wall | speedup | (p3,CCp3,p3) wall | speedup |
|---|---|---|---|---|
| 1 | 120.5 s | 1.00x | 128.4 s | 1.00x |
| 2 | 71.3 s | 1.69x | 79.2 s | 1.62x |
| 4 | 39.9 s | 3.02x | 42.9 s | 2.99x |
| 8 | 42.2 s | 2.86x | 49.0 s | 2.62x |

Largest relative deviation of the residual trace against np=1, over
`coupledResidual`, `Vm_relL2`, `Iion_relL2` and `maxState_relL2`: 5.1e-07 /
2.8e-05 / 9.8e-07 for the first triad and 6.7e-08 / 5.0e-07 / 1.0e-05 for the
second. Both numbers belong in the same table: a speedup figure means nothing
unless the parallel run computes the same thing, and this solver already
produced one case where every operator fingerprint matched and the answer was
wrong (I.5).

It saturates at 8 ranks - 80 cells per rank, where the processor-face band stops
being negligible against the volume. The same shape the MMS solver showed, just
reached earlier because the mesh is smaller.

**A3, physical acceptance.** Full runs to P8 activation, hexa, dx = 0.5 mm,
dt = 1e-4, (p3, CCp3, p3). 8 ranks omitted deliberately: it had already
saturated, and 80 cells per rank would not be informative.

| ranks | status | cells activated | T_P8 [ms] | dev | CV [m/s] | dev | wall [s] | speedup |
|---|---|---|---|---|---|---|---|---|
| 1 | OK | 640 | 52.0795930112 | - | 0.38388 | - | 310.3 | 1.00x |
| 2 | OK | 640 | 52.0795930123 | 2.11e-11 | 0.38388 | 4.78e-11 | 177.0 | 1.75x |
| 4 | OK | 640 | 52.0795930133 | 4.03e-11 | 0.38388 | 7.74e-11 | 102.9 | 3.01x |

CV from a least-squares fit over 20 diagonal samples, R2 = 0.99649 in all three.

`OK` is the most important entry in each row: it is the binary conduction-block
detector, and it says the wavefront crosses every partition cut.

**The 0.5% tolerance proposed for T_P8 was wrong, and wrong in the safe
direction.** The argument for it was that `NONLINEAR_TOLERANCE = 1e-2` sets a
noise floor which a threshold-crossing time integrates over hundreds of steps.
That does not happen: the observed spread is 4e-11 relative, eight orders of
magnitude below the proposed bound. Picard converges to the same fixed point
deterministically, so a loose stopping residual does not introduce drift - it
stops the iteration at the same place on every rank, and the partition perturbs
only at rounding level.

This configuration is also the one that historically failed to ignite on hexa
meshes, which the earlier investigation traced to the state-reset bug. With that
fixed, (p3, CCp3, p3) ignites and propagates correctly.

**A4, decomposition invariance.** The same A3 run cut with `simple` (a regular
block split, n = (2 1 1) and (2 2 1)) instead of `scotch` (graph partitioning):

| run | status | T_P8 [ms] | dev vs serial | CV [m/s] |
|---|---|---|---|---|
| scotch np=2 | OK | 52.0795930123 | 2.11e-11 | 0.38388 |
| scotch np=4 | OK | 52.0795930133 | 4.03e-11 | 0.38388 |
| simple np=2 | OK | 52.0795930123 | 2.11e-11 | 0.38388 |
| simple np=4 | OK | 52.0795930064 | 9.22e-11 | 0.38388 |

Two genuinely different cut geometries agreeing is much stronger evidence than
several rank counts under one method.

**Unstructured triangular mesh** (the gmsh route, 1500 cells), which is the
coverage R1 could not provide:

| ranks | status | cells | T_P8 [ms] | dev | CV [m/s] | wall [s] | speedup |
|---|---|---|---|---|---|---|---|
| 1 | OK | 1500 | 51.9005456482 | - | 0.38259 | 669.0 | 1.00x |
| 2 | OK | 1500 | 51.9005456487 | 9.63e-12 | 0.38259 | 404.9 | 1.65x |
| 4 | OK | 1500 | 51.9005456516 | 6.55e-11 | 0.38259 | 235.9 | 2.84x |

CV differs from the hexa value by 0.34%, which is a discretisation difference
between the two meshes, not a parallel one - it is identical across rank counts
on each mesh.

Getting there required fixing the case, not the solver, and the fix corrected a
claim made earlier in this document. `decomposePar` - which, unlike the solver,
demands a patchField entry for every patch - failed with "inconsistent patch and
patchField types for patch type empty and patchField type fixedValue" on
`0/u1/boundaryField/zmin`. The gmsh mesh keeps the 3-D xmin..zmax patch names,
so those fields carried LITERAL entries declaring fixedValue on patches that are
`empty` in 2-D.

The blockMesh case had been fixed earlier by adding a `".*"` catch-all, and the
note written then said an empty patch has its type overridden by the patch
constraint. That is true only for a REGEX entry. A literal entry is not
overridden, and a mismatch is fatal. That asymmetry is precisely why the
blockMesh case worked (its walls/frontAndBack patches had no literal entry and
fell through to the regex) and the gmsh case did not. The obsolete literal
entries were removed.

**JFNK in parallel** (section H listed this as untested): np = 2 and 4 both run,
with the same number of nonlinear iterations as serial, agreement of 4e-11, and
`lineSearchIters` identical exactly - which is what validates the global
reductions added inside the Armijo line search.

## J. Stabilisation across processor faces (MMS solver)

Section G closed the stabilisation gap as a documented omission. It is now
implemented in `highOrderManufacturedFDAImplicitPETScDistributed`. The electro
solver is deliberately left as it was.

### What was missing

`assembleHighOrderStiffnessMatrix`'s boundary loop applied its face-jump
stabilisation only under `bcType == fixedValue | fixedVoltage`. A processor
patch passed the empty/zeroGradient filter, contributed its flux term, and then
fell through with no stabilisation at all - so at `alpha > 0` the operator was
missing a term in the band along every partition cut.

`assembleStandardOrthogonalStiffnessMatrix` (the `useHighOrder_Vm = false` path)
was worse: it had no coupled branch of any kind, and wrote LOCAL cell indices as
matrix columns for its orthogonal terms while `addCellGradientDotCoeffs` wrote
GLOBAL ones for its stabilisation terms. Both are invisible in serial, where
`gRowStart == 0`. That path could not have run correctly at more than one rank.

### How the far cell's row crosses the partition

The jump needs the neighbouring cell's reconstruction row - up to ~60
(global column, coefficient) pairs for a p3 stencil in 3-D - which no field swap
can carry. The rank that owns the cell builds the row itself and ships it
(`exchangeCoupledFaceRows`, one `PstreamBuffers` round per assembly):

- high-order path: `buildCellTaylorExtrapolationRow`, evaluated at the face
  centre, which both ranks have;
- orthogonal path: `buildCellGradientDotRow`, evaluated at MINUS the sending
  rank's `patch.delta()`, so the row arrives already expressed in the receiving
  rank's orientation.

Consequently no cell centres are exchanged and `mesh.C()`'s coupled values - a
`slicedVolVectorField`, valid only after an evaluation nothing here performs -
are never touched. `patch.delta()` supplies `dPN`; the far cell's conductivity
comes from `syncTools::swapBoundaryCellList`, since a processor patch field
holds the interpolated FACE value, not the neighbouring CELL value.

Only the local row is written on each side. Writing the jump from the point of
view of the local cell,

    a/V[c] * ( T_far(xf) - T_c(xf) )

removes the owner/neighbour convention entirely: the `+` the internal-face loop
puts in the owner row and the `-` it puts in the neighbour row are the same
expression seen from either cell, and `a` is orientation-independent (`|dPN . n|`
and `n . D . n`), so both ranks compute it identically without agreeing on
anything.

The Taylor coefficient of a single stencil entry was factored out
(`cellTaylorEntryCoeff`) so the shipped row and the locally assembled one cannot
drift apart.

### Acceptance

`np=1` vs `np=2` deviation of `Vm_L2`, alpha = 0.1, diagonalIion, against the
alpha = 0 control:

| dim | mesh | before | after | alpha=0 |
|---|---|---|---|---|
| 2D | hexa | 3.34e-05 | 4.55e-06 | 2.40e-09 |
| 2D | triangular | 1.99e-04 | 1.30e-07 | 2.04e-07 |
| 3D | hexa | 6.80e-03 | 7.38e-10 | 1.81e-10 |
| 3D | triangular | 7.40e-04 | 9.55e-04 | 9.70e-04 |

Serial results are unchanged to every digit. 3-D triangular collapsed onto its
alpha = 0 control, which is the acceptance criterion; that control is itself
wrong, see below.

2-D hexa did not reach its alpha = 0 control, and that turned out not to be a
discretisation issue: with `petscLinearPcType jacobi` - partition-independent,
unlike ilu, which PETSc wraps in block-Jacobi in parallel - np=1 and np=2 agree
to **all 13 printed digits**, deviation exactly 0. The 4.55e-06 is preconditioner
path dependence in a solve whose error norm is 2.7e-05 to begin with.

### The check that made this decidable

Setting `CF_DUMP_K` in the environment writes K as (global row, global column,
value), one file per rank. Mapping the parallel indices back through
`processorN/constant/polyMesh/cellProcAddressing` gives an entry-by-entry
comparison against the serial assembly. For 2-D hexa at alpha = 0.1 the fixed
high-order operator matches serial **exactly** (47012 nonzeros, identical
sparsity, no entry differing by more than 1e-12 relative), which is what
separated the residual 4.55e-06 above from anything the assembly does.

A row-sum or nnz fingerprint cannot do this: an omitted face contributes nothing
to either of its rows, so row sums still vanish, and a wrong COLUMN changes
neither nnz nor any sum.

### Two operator discrepancies this exposed, both in the MLS stencils

Both are independent of the stabilisation and are NOT fixed here.

Neither affects a convergence result. Comparing the `Vm` fields cell by cell in
3-D tet at alpha = 0.1, the np=1 vs np=2 difference is 2.2 % of the error
against the analytical solution (RMS 1.51e-05 against 6.84e-04) and 7.5e-06 of
the solution amplitude; over a refinement sweep N = 8, 10, 12, 14 the fitted
order is **2.90722821 at one rank against 2.90722635 at two**, a difference of
1.9e-06 against the 0.05 tolerance. The deviation does not grow under
refinement, which is what two consistent schemes of the same formal order
predict.

**3-D unstructured, any alpha.** At alpha = 0 the operator already differs
between np=1 and np=2: 552 of 4971 rows, up to 4.17 against a maximum |K| of
256, with essentially every entry of an affected row different - the signature of
a different stencil, not of a different coefficient. 164 of those rows own a face
on the cut; the rest sit within a stencil radius of it. Meshes are byte-identical
and `petscLinearPcType jacobi` changes nothing (9.548e-04 either way). This is
the 9.70e-04 that section G's successor left unexplained; it is now located, in
the flux term's face stencils, not in the solver.

**Low order (`useHighOrder_Vm = false`, p1: N=1, Nn=13, maxStencilSize 33), 2-D
hexa.** At alpha = 0 the operator matches serial exactly. At alpha = 0.1, 12
rows differ - all where the partition cut meets the domain boundary - and one
serial row has 26 columns against the parallel 25. Again a stencil-membership
difference, this time in the cell gradient rows, and again on the far side of the
cut. Worth keeping because it is a 900-cell reproducer of the same phenomenon as
the 3-D tetrahedral case, and because the p3 stencil on the same mesh reproduces
serial exactly - so whatever selects the stencil is partition-dependent only when
the requested stencil is small relative to the candidate set.

## H. Not done yet

**The JFNK path has not been exercised in parallel.** Done for the electro
solver (section I.7). Not yet done for the MMS solver.

**`run_convergence.py` has no parallel support.** Done for the electro solver in
`tutorials_electro/electroDistValidation/` (section I); the MMS driver already
had it.

**The electro solver has not been ported.** Done - see section I. Note that the
Dirichlet piece anticipated here turned out not to be needed: no case in the
benchmark uses a Dirichlet condition on `Vm`.

**The RKF45 adaptive-step carry is still discarded per call.** Noted while
fixing the hybrid scratch (I.3): every call restarts from ~1e-13 ms and climbs
by a factor 4 per sub-step, wasting roughly 20 sub-steps per point per step.
Making the carry real is a large results change and deserves its own decision.

**The triangular_Unstr meshes are not in R1.** They archive a gmsh `.msh` plus a
boundary-type rewrite performed by the tutorial driver, so reproducing them
exactly means going through the driver; they are covered by the acceptance sweep,
which does exactly that.

## K. The same stabilisation fix, ported to the electro solver

Section J did this for the MMS solver and deliberately left the electro one
alone. This is that follow-up. Same mechanism, ported verbatim; the differences
are all in how the target solver is written.

### What was missing

`assembleStandardOrthogonalStiffnessMatrix` added the orthogonal flux on a
processor face and then `continue`d past the Rhie-Chow jump, with a comment
saying why. `assembleHighOrderStiffnessMatrix` was worse: no coupled branch at
all, the stabilisation gated on `fixedValue`/`fixedVoltage`, so a processor
patch contributed flux and nothing else - and a `cyclic` patch would have been
treated as if it were Dirichlet, silently.

### What differed from the MMS port

- **Column convention.** The MMS threads `const globalIndex& gc =
  LREInterp.globalCells()` into every helper; this solver has file-scope
  `gCol()` / `gNCols()`. Every ported line was rewritten accordingly.
- **`LREInterp` is a nullable pointer** in the orthogonal assembly, not a
  reference. The row-building lambda captures it, so it is constructed only
  under `stabilisationAlpha > SMALL` - at alpha = 0 a null pointer is legal and
  must not be dereferenced.
- **Chunked assembly.** The high-order path assembles in blocks of 50000 faces.
  `exchangeCoupledFaceRows` is COLLECTIVE, so it sits before the boundary loop,
  outside any chunk loop: inside one it would deadlock as soon as two ranks had
  a different number of chunks. NOTE: at 640 cells there is only one chunk, so
  the runs below do NOT exercise this. It needs a case with >50000 faces.

### A regression the verification caught

Extracting `cellTaylorEntryCoeff` (so the row written locally and the row
shipped over the wire share one piece of arithmetic) changed the SERIAL result.
The serial check failed in exactly the informative way: bit-identical at
alpha = 0 in all six configurations, different at alpha = 0.1 in all six - which
points straight at the only code that runs only when alpha > 0.

Cause: the association order of the self coefficient. The original computed
`((1 + grad.d) + 0.5*quad) + cubic/6`; folding the `1` in at the end instead
gives `1 + ((grad.d + 0.5*quad) + cubic/6)`. Floating-point addition is not
associative, and over 50 steps of a stiff system that grew to 4.1e-07 relative.

Fixed by passing the zeroth-order term as a `base` PARAMETER rather than adding
it afterwards, so the arithmetic is bit-identical to what the solver produced
before the refactor and the shipped row still shares it.

Worth recording that the fix "worked" by every other measure while this
regression was present. Only a check that was purely defensive - serial must not
move at all - found it.

### Acceptance

Electro 2-D, dx = 0.5 mm, 50 steps, (p3,CCp3,p3), `pcType = jacobi`, np=1 vs
np=2 deviation of the residual trace. Two DISTINCT binaries were built
(`electroPreStab` = pre-port, `...Distributed` = ported) because a first attempt
had rebuilds landing while the baseline sweep was still running, which would
have made a null result indistinguishable from a working fix.

| mesh | method | before | after | alpha=0 control | gain |
|---|---|---|---|---|---|
| hexa | Picard | 3.23e-09 | 9.06e-12 | 8.53e-12 | 356x |
| hexa | JFNK | 2.58e-05 | 2.59e-05 | 2.48e-05 | 1x |
| hexa | diagonalIion | 1.27e-08 | 4.51e-10 | 2.28e-10 | 28x |
| triangular | Picard | 2.45e-07 | 4.68e-11 | 9.55e-12 | 5240x |
| triangular | JFNK | 4.44e-06 | 4.43e-06 | 4.64e-06 | 1x |
| triangular | diagonalIion | 1.59e-07 | 6.17e-09 | 5.59e-10 | 26x |

Serial bit-identical (SHA-256) in all 12 configurations, alpha = 0 and 0.1.

The assembled paths collapse essentially onto the alpha = 0 control, which is
the achievable floor. **JFNK does not move, and that is expected**: its
deviation is dominated by its own preconditioner - the Eigen ILUT of the
diagonal block, which is partition-dependent by construction - not by the
assembled stabilisation. The same 2.48e-05 appears with `ilu` and with `jacobi`,
which is what identified it.

### Not done

- The collective-inside-chunks hazard is reasoned about but NOT exercised.
- `CF_DUMP_K` / `CF_DUMP_STENCILS` were not ported, so there is no
  entry-by-entry comparison of K against serial for this solver.
- Cyclic patches now `FatalError` in both assemblies rather than being silently
  mistreated. No benchmark case uses one.

