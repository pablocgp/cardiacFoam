# Handoff: stabilisation across processor faces

## Done (this session)

Rhie-Chow / face-jump stabilisation on coupled faces is implemented in
`applications/solvers/highOrderManufacturedFDAImplicitPETScDistributed/`, for
BOTH assemblies:

- `assembleHighOrderStiffnessMatrix` - the term that was missing;
- `assembleStandardOrthogonalStiffnessMatrix` - which had no coupled branch at
  all, not even the orthogonal flux, and mixed local and global matrix columns.
  That path could not have run at more than one rank.

Full write-up in `CHANGES-cardiacFoam.md`, section J. The mechanism in one line:
the rank that owns a cell builds its reconstruction row itself
(`exchangeCoupledFaceRows` + `buildCellTaylorExtrapolationRow` /
`buildCellGradientDotRow`) and ships (global column, coefficient) pairs, so
nothing but the shared face geometry has to agree across the cut.

Acceptance, `np=1` vs `np=2` deviation of `Vm_L2` at alpha = 0.1:

| dim | mesh | before | after | alpha=0 control |
|---|---|---|---|---|
| 2D | hexa | 3.34e-05 | 4.55e-06 | 2.40e-09 |
| 2D | triangular | 1.99e-04 | 1.30e-07 | 2.04e-07 |
| 3D | hexa | 6.80e-03 | 7.38e-10 | 1.81e-10 |
| 3D | triangular | 7.40e-04 | 9.55e-04 | 9.70e-04 |

Serial unchanged to every digit. The 2-D hexa residual is not discretisation:
with `PETSC_LINEAR_PC_TYPE = "jacobi"` np=1 and np=2 agree to all printed
digits, and the assembled K matches serial entry by entry.

## Done: the electro solver too

Ported in a later session. Full write-up in `CHANGES-cardiacFoam.md`, section K.
Serial bit-identical in all 12 configurations; np=1 vs np=2 deviation at
alpha = 0.1 improved 356x / 28x (2-D hexa, Picard / diagonalIion) and
5240x / 26x (2-D triangular), collapsing onto the alpha = 0 control. JFNK is
unchanged by design - its deviation is its own preconditioner, not the
assembled stabilisation.

Two things that bit, worth knowing before the next port of this kind:

- Extracting `cellTaylorEntryCoeff` changed the SERIAL answer by 4.1e-07,
  because folding the zeroth-order `1` in at the end re-associates the sum.
  Pass it as a `base` parameter instead.
- Build the pre-change binary under its OWN name before measuring a baseline.
  A rebuild landing mid-sweep makes a null result look identical to a working
  fix.

Still open for the electro solver: the collective-inside-chunked-assembly
hazard is reasoned about but not exercised (needs >50000 faces), and
`CF_DUMP_K` was not ported, so there is no entry-by-entry K comparison.

## Open problem: the MLS stencil is not a function of the geometry

Diagnosed. It is NOT the stabilisation, and the work is in solids4foam's
`movingLeastSquaresStencil::buildFacesStencil`, not in the solver.

Two dumps do the measuring, both env-gated in the MMS solver: `CF_DUMP_K`
(global row / global column / value) and `CF_DUMP_STENCILS` (cell stencils,
face stencils, cell centres). Map the parallel indices back through
`processorN/constant/polyMesh/{cellProcAddressing,faceProcAddressing}`.

### What was measured

3-D tet, 4971 cells, alpha = 0, np=1 vs np=2:

- 552 of 4971 rows of K differ, up to 4.17 against max |K| = 256.
- **536 of 10676 face stencils differ**, all at exactly the same size (80): the
  members are SUBSTITUTED, not truncated - 929 swaps in all.
- 98 % of the swapped members sit in the last 10 % of the distance-ordered
  stencil, i.e. right at the N-th cut. Median relative position 0.988.
- Parallel GAINS 443 members that live on the other rank and LOSES 691 that are
  local.
- **9 cut faces get two different stencils from their two ranks** - so the two
  rows of one face are built from different flux stencils, and the scheme is not
  conservative across the cut. This one is a defect on its own terms, with no
  serial comparison needed.
- Cell centres differ between the serial and the decomposed mesh by up to
  3.3e-16 on a scale of 1 (157 of 4971 cells): the centroid sum runs over a
  renumbered face list. Harmless at `relTol = 1e-8`, but it rules the
  centroid out as a stable sort key.

### How much it actually matters

Measured, not inferred. 3-D tet, alpha = 0.1, `jacobi`, comparing the `Vm`
FIELDS cell by cell (np2 reconstructed) rather than the deviation of a norm:

- `||Vm_np1 - Vm_np2||` (RMS) = 1.51e-05, max 2.49e-04
- `||Vm - exact||` (RMS) = 6.84e-04, max 8.20e-03
- so the partition-to-partition difference is **2.2 % of the discretisation
  error**, and 7.5e-06 of the solution amplitude.

Refinement sweep N = 8, 10, 12, 14, np=1 against np=2:

| N | Vm_L2 np1 | Vm_L2 np2 | rel. dev |
|---|---|---|---|
| 8 | 1.187049e-03 | 1.186528e-03 | 4.39e-04 |
| 10 | 6.429373e-04 | 6.423235e-04 | 9.55e-04 |
| 12 | 3.975660e-04 | 3.975437e-04 | 5.60e-05 |
| 14 | 2.275209e-04 | 2.273524e-04 | 7.41e-04 |

**Fitted order: 2.90722821 (np1) against 2.90722635 (np2), a difference of
1.9e-06** against the 0.05 acceptance tolerance. Every other metric agrees to
1e-5 or better; `VmG_L2` gives 3.65710763 against 3.65717835.

The deviation does not grow under refinement - it stays at the 1e-4..1e-3 level,
i.e. the perturbation shrinks at the same rate as the error. That is what the
two stencils being consistent schemes of the same formal order predicts, and it
is why the fitted order is untouched.

So: real, worth fixing, and it does NOT invalidate any convergence result.

### The cause

`buildFacesStencil` builds its candidate set by growing complete layers of
`mesh_.cellCells()` until it holds `ceil(1.5*N)` cells, then takes the nearest N
by distance. That candidate set is TOPOLOGICAL. On a tet mesh a topological
neighbourhood is not a geometric ball, so with only 0.5*N cells of slack the
nearest-N of the candidates is not the nearest-N of the mesh.

In serial that is all that happens: `procCells()` is empty, so Phase 2 always
takes the early return and every stencil is the nearest N of a layer set.

In parallel the layer growth stops at the cut, and the far side is supplied by
`remoteCandidates()`, which returns every cell of the neighbour rank whose
centre lies in this rank's face bounding box grown by the halo depth - a
GEOMETRIC set, and a generous one. Near a cut the selection therefore mixes a
topological candidate set on one side with a geometric one on the other. That is
why parallel gains remote members and loses local ones, and why
`haloDepthScale` 7.5 -> 16 changed nothing: the box was never the binding
constraint.

The 2-D hexa p1 case adds a second, independent defect: the parallel branch
filters candidates by `d2 <= R2` BEFORE the tie expansion (`cut*(1+relTol)`)
runs, so a tie group straddling R2 is clipped - serial keeps 17 members, parallel
16. The serial early return never clips because it expands over the full local
list.

### How to fix it

1. **Make the candidate set geometric on both sides.** Query cell centres by
   radius (`indexedOctree`) instead of growing `cellCells` layers, with the
   radius taken from the layer estimate as an upper bound, then nearest N plus
   the tie group. The stencil becomes a function of (face centre, cell centres,
   N) alone - invariant under partitioning, and no longer hostage to tet
   connectivity.
2. **Stop the tie clipping**: inflate R2 by `(1 + relTol)` before filtering
   remote candidates, and route the serial early return through the same
   selection code so the two branches cannot differ.
3. **Assert completeness** instead of failing silently: after selecting, check
   that the N-th distance is strictly inside the queried radius. That check
   would have turned this into an immediate error rather than a 1e-3 drift.
4. **One face, one stencil.** Have the lower-ranked side of a cut face own the
   stencil and send it, or make the selection deterministic enough that both
   sides agree by construction. Without this the flux leaving a cell is not the
   flux entering its neighbour.

`relTol` (1e-8), `overSampleFactor` (1.5) and `searchExpansionFactor` (1.3) are
hard-coded default arguments, not dictionary entries, so none of this can be
explored from a case - it has to be a library patch.

**Note before touching it**: serial is the truncated one, so the fix changes
serial results too, and with them every convergence table in the study. Worth
deciding deliberately rather than as a side effect.

## Watch out - five bugs of one family found this way

A local index used where PETSc/Eigen expects a global one is INVISIBLE in serial
(`gRowStart == 0`, `nLocal == nGlobal`). In the electro solver: diagonal mass
column, `addDiagonalToMatrix` column, `updateValues` size check, `updateValues`
row offset. In the MMS solver: the whole orthogonal assembly (this session).

The operator fingerprint (nnz / sum / sumAbs / row sums) CANNOT see a wrong
column, and cannot see an omitted face either - it must be paired with a
solution comparison, or with the `CF_DUMP_K` entry-by-entry diff.

Also: to compare np=1 against np>1 iteration-by-iteration, set
`PETSC_LINEAR_PC_TYPE = "jacobi"`. ilu/lu are wrapped in block-Jacobi in
parallel and are therefore partition-dependent by construction - that alone
accounts for a 4.5e-06 deviation in the 2-D hexa MMS error norm.
