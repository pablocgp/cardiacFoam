/*---------------------------------------------------------------------------*\
    cardiacFoam - shared code of the monodomain solvers
    (highOrderElectroActivationFoamImplicitPETScDistributed and
    highOrderManufacturedFDAImplicitPETScDistributed)

    Diffusion-operator assembly: face coefficients, the
    non-orthogonal correction, the processor-face row exchange.

    INCLUDED INSIDE THE SOLVER'S ANONYMOUS NAMESPACE, after its system,
    OpenFOAM, Eigen and PETSc includes, through monodomainCommon.H. It
    therefore includes nothing itself and has internal linkage in each
    executable. Code here is used verbatim by both solvers: a change made here
    changes both, which is the point.
\*---------------------------------------------------------------------------*/

    // Added for cardiacFoam: exchange one variable-length matrix row per
    // coupled face.
    //
    // The rank that OWNS a cell builds that cell's reconstruction row itself
    // and ships (global column, coefficient) pairs. Nothing but the shared face
    // geometry has to agree across the cut - in particular no stencil, no cell
    // centre and no owner/neighbour convention.
    //
    // COLLECTIVE. Call it exactly once per assembly, outside any chunk loop:
    // inside one it would deadlock as soon as two ranks had a different number
    // of chunks. The decision to call it must also be uniform across ranks.
    //
    // Results are indexed by boundary-face position (faceI - nInternalFaces);
    // non-processor boundary faces are left empty.
    void exchangeCoupledFaceRows
    (
        const fvMesh& mesh,
        const std::function
        <
            void(const label patchI, const label faceI, labelList&, scalarField&)
        >& buildRow,
        List<labelList>& nbrCols,
        List<scalarField>& nbrCoeffs
    )
    {
        nbrCols.setSize(mesh.nBoundaryFaces());
        nbrCoeffs.setSize(mesh.nBoundaryFaces());

        if (!Pstream::parRun())
        {
            return;
        }

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            const processorFvPatch* ppPtr = isA<processorFvPatch>(patch);

            if (!ppPtr)
            {
                continue;
            }

            List<labelList> sendCols(patch.size());
            List<scalarField> sendCoeffs(patch.size());

            forAll(patch, faceI)
            {
                buildRow(patchI, faceI, sendCols[faceI], sendCoeffs[faceI]);
            }

            UOPstream toNbr(ppPtr->neighbProcNo(), pBufs);
            toNbr << sendCols << sendCoeffs;
        }

        pBufs.finishedSends();

        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            const processorFvPatch* ppPtr = isA<processorFvPatch>(patch);

            if (!ppPtr)
            {
                continue;
            }

            List<labelList> recvCols;
            List<scalarField> recvCoeffs;

            UIPstream fromNbr(ppPtr->neighbProcNo(), pBufs);
            fromNbr >> recvCols >> recvCoeffs;

            if (recvCols.size() != patch.size() || recvCoeffs.size() != patch.size())
            {
                FatalErrorInFunction
                    << "Received " << recvCols.size() << " rows for patch "
                    << patch.name() << " of size " << patch.size()
                    << exit(FatalError);
            }

            const label bStart = patch.start() - mesh.nInternalFaces();

            forAll(patch, faceI)
            {
                nbrCols[bStart + faceI].transfer(recvCols[faceI]);
                nbrCoeffs[bStart + faceI].transfer(recvCoeffs[faceI]);
            }
        }
    }

    // Two-point orthogonal diffusion coefficient of a face,
    // |Sf| (n . D . e)/|d|, with e the unit owner-to-neighbour direction. The
    // standard low-order finite-volume approximation of the diffusive flux; it
    // neglects the non-orthogonal part of d relative to Sf.
    scalar orthogonalDiffusionCoeff
    (
        const vector& Sf,
        const vector& d,
        const tensor& D
    )
    {
        const scalar area = mag(Sf) + VSMALL;
        const vector n = Sf/area;
        const scalar dMag = mag(d) + VSMALL;
        const vector e = d/dMag;

        return area*(n & (D & e))/dMag;
    }

    // Diffusivity resolved along a face normal, n . D . n, floored at SMALL so
    // that a degenerate tensor cannot produce a zero stabilisation scale.
    scalar normalDiffusivity(const tensor& D, const vector& n)
    {
        return max(mag(n & (D & n)), SMALL);
    }

    // Scale of the Rhie-Chow-like face-jump stabilisation:
    // alpha |Sf| (n.D.n)/|d.n|. Zero when alpha is zero, which switches the
    // whole stabilisation off. Note it depends on the face orientation only
    // through |d.n|, so the two sides of a face - including the two ranks of a
    // processor face - compute the same value.
    scalar stabilisationFaceCoeff
    (
        const scalar alpha,
        const scalar area,
        const tensor& D,
        const vector& n,
        const vector& d
    )
    {
        if (alpha <= SMALL)
        {
            return 0.0;
        }

        const scalar dn = max(mag(d & n), VSMALL);
        return alpha*area*normalDiffusivity(D, n)/dn;
    }

    // Added for cardiacFoam, ported from the MMS solver: the part of the
    // diffusive flux that a two-point flux cannot see.
    //
    // The exact flux through a face is |Sf| n . D . grad(V)|_f. Splitting the
    // gradient along the owner-to-neighbour direction e = d/|d| and the rest,
    //
    //     |Sf| n.D.grad V = |Sf| (n.D.e)(grad V . e)  +  |Sf| (grad V . w)
    //
    // and approximating (grad V . e) by (V_nei - V_own)/|d| turns the FIRST term
    // into exactly orthogonalDiffusionCoeff() times the two-point difference.
    // This function returns the vector w of the SECOND term, which the two-point
    // flux drops:
    //
    //     g = D . n,    w = g - (g . e) e
    //
    // i.e. the part of D.n orthogonal to the line joining the cell centres.
    //
    // Two cases make it non-zero, and both are real here:
    //   - a NON-ORTHOGONAL mesh, where e is not parallel to n. This is what
    //     makes the uncorrected two-point flux inconsistent on triangles and
    //     tetrahedra: measured on the MMS, order 0.81 instead of 2.
    //   - an ANISOTROPIC D, where D.n is not parallel to n even on an orthogonal
    //     mesh. With fibre-oriented conductivity that is every face whose
    //     normal is not aligned with the fibres, i.e. most of a heart mesh.
    //
    // With D aligned with an orthogonal mesh (the Niederer slab: diagonal D,
    // axis-aligned hexahedra), g = lambda n and e = n, so w = 0 EXACTLY. That
    // is why the term carries no dictionary key: where it would cost something
    // it is also required for consistency, and where it is not required it
    // costs nothing.
    //
    // Same decomposition as OpenFOAM's gaussLaplacianScheme (nonOrthDeltaCoeffs
    // plus nonOrthCorrectionVectors), except that here the correction is
    // assembled IMPLICITLY, through the reconstruction rows, instead of being
    // deferred to the right-hand side.
    vector nonOrthFluxVector(const vector& Sf, const vector& d, const tensor& D)
    {
        const vector n = Sf/(mag(Sf) + VSMALL);
        const vector e = d/(mag(d) + VSMALL);
        const vector g = D & n;
        const vector w = g - (g & e)*e;

        // w is the difference of two nearly equal vectors when e is nearly
        // parallel to D.n, so where it should vanish it comes out at ROUNDOFF
        // rather than at zero - face normals and cell centres are themselves
        // computed by summation. That residue is not harmless: multiplied by
        // |Sf|/V and the reconstruction coefficients it lands above the ABSOLUTE
        // SMALL that addTripletIfNeeded filters on and survives into K. Measured
        // on the hexahedral MMS at N=20: 3043 extra entries in K, all at 8.0e-16
        // of the largest, moving Vm_L2 by 2.8e-07 relative.
        //
        // So the test is RELATIVE to |g|, the scale w is a remainder of, with a
        // threshold inside a measured gap: on the unstructured MMS mesh the
        // term's largest entry is 2.3e-02 of the largest entry of K, fourteen
        // orders above the hexahedral residue. 1e-12 sits ~300x above the noise
        // and ~1e11 below the signal.
        const scalar wRelativeFloor = 1.0e-12;

        if (mag(w) <= wRelativeFloor*max(mag(g), VSMALL))
        {
            return vector::zero;
        }

        return w;
    }

    // EXPLICIT high-order diffusion: evaluates the face fluxes
    // n . (D grad Vm) by LRE quadrature over each face and takes their
    // divergence, writing both the face flux field and the cell Laplacian.
    //
    // This is the diagnostic/post-processing counterpart of
    // assembleHighOrderStiffnessMatrix, which builds the same operator as a
    // matrix. Boundary faces with zeroGradient or empty carry no flux;
    // gradScalarFaceQuad performs its own halo exchange, so this is
    // parallel-correct as written.
    void computeHighOrderLaplacian
    (
        const volScalarField& Vm,
        const volTensorField& conductivity,
        const highOrderInterp& LREInterp_Vm,
        surfaceScalarField& fluxVm_HO,
        volScalarField& lapVm
    )
    {
        const fvMesh& mesh = Vm.mesh();

        autoPtr<CompactListList<vector>> gradVmQuadPtr =
            LREInterp_Vm.gradScalarFaceQuad(Vm);
        // Modified for cardiacFoam: CompactListList, not List<List>. It is
        // indexed identically but stores contiguously.
        CompactListList<vector>& gradVmQuad = gradVmQuadPtr.ref();
        const CompactListList<scalar>& faceQuadW = LREInterp_Vm.faceQuadWeightPhysical();

        const surfaceVectorField nHat(mesh.Sf()/mesh.magSf());

        for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
        {
            const label owner = mesh.owner()[faceI];
            const vector& faceNormal = nHat[faceI];

            fluxVm_HO[faceI] = 0.0;

            forAll(gradVmQuad[faceI], pI)
            {
                const vector Dg = conductivity[owner] & gradVmQuad[faceI][pI];
                // Modified for cardiacFoam: no |Sf| factor - the
                // physical quadrature weights already sum to it.
                fluxVm_HO[faceI] +=
                    (faceNormal & Dg)*faceQuadW[faceI][pI];
            }
        }

        forAll(fluxVm_HO.boundaryField(), patchI)
        {
            scalarField& patchFlux = fluxVm_HO.boundaryFieldRef()[patchI];

            if (patchFlux.size() == 0)
            {
                continue;
            }

            const word bcType = Vm.boundaryField()[patchI].type();
            if
            (
                bcType == zeroGradientFvPatchScalarField::typeName
             || bcType == "empty"
            )
            {
                patchFlux = 0.0;
                continue;
            }

            const label start = mesh.boundaryMesh()[patchI].start();
            const vectorField& pNormals = nHat.boundaryField()[patchI];

            forAll(patchFlux, faceI)
            {
                const label globalFaceI = start + faceI;
                const label owner = mesh.owner()[globalFaceI];
                const vector& faceNormal = pNormals[faceI];

                patchFlux[faceI] = 0.0;

                forAll(gradVmQuad[globalFaceI], pI)
                {
                    const vector Dg = conductivity[owner] & gradVmQuad[globalFaceI][pI];
                    // Modified for cardiacFoam: no |Sf| factor.
                    patchFlux[faceI] +=
                        (faceNormal & Dg)*faceQuadW[globalFaceI][pI];
                }
            }
        }

        lapVm = fvc::div(fluxVm_HO);
    }

    // Add coefficient*diagonal[cell] to the diagonal of an assembled matrix,
    // in place. Used by the diagonalIion linearisation, which folds a fresh
    // -theta*dIion/dVm into A on every nonlinear iteration.
    void addDiagonalToMatrix
    (
        const scalarField& diagonal,
        const scalar coefficient,
        SpMat& A
    )
    {
        // Modified for cardiacFoam: the column must be GLOBAL, like every
        // other column in these matrices. Same trap as the diagonal mass
        // matrix: a wrong column changes neither nnz nor any row sum, so the
        // operator fingerprint cannot see it.
        forAll(diagonal, cellI)
        {
            A.coeffRef(cellI, gCol(cellI)) += coefficient*diagonal[cellI];
        }

        A.makeCompressed();
    }

    // Lumped mass matrix: coefficient (chi*Cm, or 1) on the diagonal. Rows are
    // local, the column is global.
    void assembleDiagonalMassMatrix
    (
        const fvMesh& mesh,
        const scalar coefficient,
        SpMat& M
    )
    {
        std::vector<Triplet> triplets;
        triplets.reserve(mesh.nCells());

        forAll(mesh.C(), cellI)
        {
            // Modified for cardiacFoam: the column must be GLOBAL. This site
            // was missed by the first sweep because it builds its triplet
            // directly instead of going through addTripletIfNeeded.
            //
            // It is worth recording why the operator fingerprint did not catch
            // it: putting an entry in the wrong COLUMN changes neither nnz, nor
            // sum, nor sumAbs, nor any row sum. Only combining M with L, where
            // the misplaced diagonal no longer coincided with L's, exposed it -
            // as 320 extra nonzeros in AImplicit at two ranks, exactly the
            // per-rank row count.
            triplets.emplace_back(cellI, gCol(cellI), coefficient);
        }

        M.resize(mesh.nCells(), gNCols());
        M.setFromTriplets(triplets.begin(), triplets.end());
        M.makeCompressed();

        reportOperatorFingerprint("M (lumped)", M);
    }

    // Expand grad Vm|_cellI . d into matrix coefficients: one entry per cell of
    // that cell's LRE stencil plus one for the cell itself, all with GLOBAL
    // column indices. The building block of the stabilisation jump.
    void addCellGradientDotCoeffs
    (
        std::vector<Triplet>& triplets,
        const label row,
        const scalar scale,
        const label cellI,
        const vector& d,
        const highOrderInterp& LREInterp
    )
    {
        if (mag(scale) <= SMALL)
        {
            return;
        }

        const CompactListList<label>& stencils =
            LREInterp.globalCellStencils();
        const CompactListList<vector>& gradCoeffs = LREInterp.QRGradCoeffs();
        // Modified for cardiacFoam: must be UList<label> BY VALUE, not
                // a labelList reference - see the note on
                // CompactListList in the header comment.
                const UList<label> curStencil = stencils[cellI];
        const label selfCoeffI = curStencil.size();

            // Added for cardiacFoam: the coefficient rows carry one
            // trailing entry for the cell itself, so this index must
            // be the last one. If the stencil silently came back
            // empty, selfCoeffI is 0, which is still in range and
            // reads a neighbour's coefficient instead.
            if (selfCoeffI != gradCoeffs[cellI].size() - 1)
            {
                FatalErrorInFunction
                    << "Stencil/coefficient mismatch at cell " << cellI
                    << ": stencil size " << selfCoeffI
                    << ", coefficient row " << gradCoeffs[cellI].size()
                    << " (expected row = size + 1)."
                    << abort(FatalError);
            }

        forAll(curStencil, cI)
        {
            addTripletIfNeeded
            (
                triplets,
                row,
                curStencil[cI],
                scale*(gradCoeffs[cellI][cI] & d)
            );
        }

        addTripletIfNeeded
        (
            triplets,
            row,
            gCol(cellI),
            scale*(gradCoeffs[cellI][selfCoeffI] & d)
        );
    }

    // Added for cardiacFoam: grad|_cellI . d as (global column, coefficient)
    // pairs. The wire-format twin of addCellGradientDotCoeffs; no 1.0 term.
    void buildCellGradientDotRow
    (
        labelList& cols,
        scalarField& coeffs,
        const label cellI,
        const vector& d,
        const highOrderInterp& LREInterp
    )
    {
        const CompactListList<vector>& gradCoeffs = LREInterp.QRGradCoeffs();

        const UList<label> curStencil = LREInterp.globalCellStencils()[cellI];
        const label selfCoeffI = curStencil.size();

        cols.setSize(selfCoeffI + 1);
        coeffs.setSize(selfCoeffI + 1);

        forAll(curStencil, cI)
        {
            cols[cI] = curStencil[cI];
            coeffs[cI] = gradCoeffs[cellI][cI] & d;
        }

        cols[selfCoeffI] = gCol(cellI);
        coeffs[selfCoeffI] = gradCoeffs[cellI][selfCoeffI] & d;
    }

    // Added for cardiacFoam: the Taylor coefficient of ONE entry of a cell's
    // reconstruction row, to the polynomial order of the interpolator
    // (gradient, plus Hessian at p2+, plus third derivatives at p3).
    //
    // Factored out so the row written straight into the matrix
    // (addCellTaylorExtrapolationCoeffs) and the row SHIPPED to another rank
    // (buildCellTaylorExtrapolationRow) cannot drift apart. Two copies of this
    // arithmetic would fail as an inconsistency between ranks, which is one of
    // the more expensive things to diagnose.
    // "base" is the zeroth-order term: 0 for a stencil entry, 1 for the self
    // entry. It is a PARAMETER rather than something the caller adds
    // afterwards because floating-point addition is not associative: the
    // original code computed ((1 + grad.d) + 0.5*quad) + cubic/6, and folding
    // the 1 in at the end instead changed the serial result by 4e-07 relative
    // over 50 steps. Threading it through keeps the arithmetic bit-identical
    // to what this solver produced before the refactor.
    scalar cellTaylorEntryCoeff
    (
        const label cellI,
        const label entryI,
        const vector& d,
        const highOrderInterp& LREInterp,
        const bool twoD,
        const scalar base = 0.0
    )
    {
        scalar coeff = base + (LREInterp.QRGradCoeffs()[cellI][entryI] & d);

        if (LREInterp.order() >= 2)
        {
            coeff +=
                0.5*quadraticForm(LREInterp.cellHessianCoeffs()[cellI][entryI], d);
        }

        if (LREInterp.order() >= 3)
        {
            coeff +=
                (1.0/6.0)
               *highOrderInterp::cubicForm
                (
                    LREInterp.cellThirdDerivCoeffs()[cellI][entryI],
                    d,
                    twoD
                );
        }

        return coeff;
    }

    // Added for cardiacFoam: the same Taylor row as
    // addCellTaylorExtrapolationCoeffs, but as (global column, coefficient)
    // pairs instead of triplets, so it can be sent to another rank.
    void buildCellTaylorExtrapolationRow
    (
        labelList& cols,
        scalarField& coeffs,
        const label cellI,
        const point& evalPoint,
        const highOrderInterp& LREInterp,
        const vectorField& C,
        const bool twoD
    )
    {
        const vector d = evalPoint - C[cellI];

        const UList<label> curStencil = LREInterp.globalCellStencils()[cellI];
        const label selfCoeffI = curStencil.size();

        cols.setSize(selfCoeffI + 1);
        coeffs.setSize(selfCoeffI + 1);

        forAll(curStencil, cI)
        {
            cols[cI] = curStencil[cI];
            coeffs[cI] = cellTaylorEntryCoeff(cellI, cI, d, LREInterp, twoD);
        }

        cols[selfCoeffI] = gCol(cellI);
        coeffs[selfCoeffI] =
            cellTaylorEntryCoeff(cellI, selfCoeffI, d, LREInterp, twoD, 1.0);
    }

    // Expand the Taylor reconstruction of cell cellI evaluated at evalPoint
    // into matrix coefficients, to the polynomial order of the interpolator.
    // Same row layout as addCellGradientDotCoeffs; the self entry additionally
    // carries the 1.0 of the zeroth-order term.
    //
    // Modified for cardiacFoam: the per-entry arithmetic now lives in
    // cellTaylorEntryCoeff, shared with the wire-format row builder.
    void addCellTaylorExtrapolationCoeffs
    (
        std::vector<Triplet>& triplets,
        const label row,
        const scalar scale,
        const label cellI,
        const point& evalPoint,
        const highOrderInterp& LREInterp,
        const vectorField& C,
        const bool twoD
    )
    {
        if (mag(scale) <= SMALL)
        {
            return;
        }

        const vector d = evalPoint - C[cellI];

        // Modified for cardiacFoam: must be UList<label> BY VALUE, not a
        // labelList reference - CompactListList::operator[] returns a view.
        const UList<label> curStencil = LREInterp.globalCellStencils()[cellI];
        const label selfCoeffI = curStencil.size();

        // Added for cardiacFoam: the coefficient rows carry one trailing entry
        // for the cell itself, so this index must be the last one. If the
        // stencil silently came back empty, selfCoeffI is 0, which is still in
        // range and reads a neighbour's coefficient instead.
        if (selfCoeffI != LREInterp.QRGradCoeffs()[cellI].size() - 1)
        {
            FatalErrorInFunction
                << "Stencil/coefficient mismatch at cell " << cellI
                << ": stencil size " << selfCoeffI
                << ", coefficient row " << LREInterp.QRGradCoeffs()[cellI].size()
                << " (expected row = size + 1)."
                << abort(FatalError);
        }

        forAll(curStencil, cI)
        {
            addTripletIfNeeded
            (
                triplets,
                row,
                curStencil[cI],
                scale*cellTaylorEntryCoeff(cellI, cI, d, LREInterp, twoD)
            );
        }

        addTripletIfNeeded
        (
            triplets,
            row,
            gCol(cellI),
            scale*cellTaylorEntryCoeff(cellI, selfCoeffI, d, LREInterp, twoD, 1.0)
        );
    }

    // Stabilisation jump on an internal face, added to one row:
    //
    //     scale*(Vm[nei] - Vm[own])
    //       - 0.5*scale*(grad Vm|_own . dPN)
    //       - 0.5*scale*(grad Vm|_nei . dPN)
    //
    // i.e. the difference between the two cell values and the value the two
    // reconstructions predict, which vanishes for a field the reconstruction
    // represents exactly and therefore does not alter the formal order. The
    // caller passes +a/V[own] for the owner row and -a/V[nei] for the
    // neighbour's.
    void addRhieChowInternalJumpCoeffs
    (
        std::vector<Triplet>& triplets,
        const label row,
        const scalar scale,
        const label own,
        const label nei,
        const vector& dPN,
        const highOrderInterp& LREInterp
    )
    {
        addTripletIfNeeded(triplets, row, gCol(nei), scale);
        addTripletIfNeeded(triplets, row, gCol(own), -scale);

        addCellGradientDotCoeffs
        (
            triplets,
            row,
            -0.5*scale,
            own,
            dPN,
            LREInterp
        );

        addCellGradientDotCoeffs
        (
            triplets,
            row,
            -0.5*scale,
            nei,
            dPN,
            LREInterp
        );
    }

    // Added for cardiacFoam: the same jump across a PROCESSOR face.
    //
    // Only the local row is written. The far rank writes its own, exactly as
    // the internal-face loop writes both rows at once when both cells are
    // local - so no owner/neighbour convention has to be agreed across the
    // partition.
    void addRhieChowCoupledJumpCoeffs
    (
        std::vector<Triplet>& triplets,
        const label row,
        const scalar scale,
        const label own,
        const label nbrGlobalCell,
        const vector& dPN,
        const labelList& nbrGradCols,
        const scalarField& nbrGradCoeffs,
        const highOrderInterp& LREInterp
    )
    {
        addTripletIfNeeded(triplets, row, nbrGlobalCell, scale);
        addTripletIfNeeded(triplets, row, gCol(own), -scale);

        addCellGradientDotCoeffs
        (
            triplets,
            row,
            -0.5*scale,
            own,
            dPN,
            LREInterp
        );

        // The remote row was built as grad|_nbr . dPN with dPN already
        // expressed in THIS rank's orientation (the sender used minus its own
        // delta), so it carries the same -0.5 factor as the local one.
        addRemoteRowCoeffs
        (
            triplets,
            row,
            -0.5*scale,
            nbrGradCols,
            nbrGradCoeffs
        );
    }

    // High-order consistent mass matrix: M_ij is the integral over cell i of
    // the LRE reconstruction basis of cell j, evaluated by cell quadrature.
    // Rows are local, columns global.
    //
    // Unlike the lumped alternative it is not diagonal - each row spans the
    // cell's whole stencil - which is what lets the temporal accuracy match the
    // spatial one, at the cost of a real linear solve.
    void assembleConsistentMassMatrixHO
    (
        const fvMesh& mesh,
        const highOrderInterp& LREInterp,
        const scalar coefficient,
        const bool compactRows,
        const label cellChunk,
        SpMat& M
    )
    {
        // Chunked assembly mirroring assembleHighOrderStiffnessMatrix.
        // For ~290K tet cells with p3 + Nn=60, the monolithic triplet vector
        // would peak at ~270 MB plus reallocation overhead. Chunking caps
        // the per-chunk triplet buffer at ~50 MB.
        const bool twoD = mesh.nGeometricD() == 2;
        const vectorField& C = mesh.C();

        const CompactListList<label>& stencils = LREInterp.globalCellStencils();
        const CompactListList<point>& cellQP = LREInterp.cellQuadPoints();
        const CompactListList<scalar>& cellQW = LREInterp.cellQuadWeightPhysical();
        const CompactListList<vector>& gradCoeffs = LREInterp.QRGradCoeffs();
        const CompactListList<symmTensor>& hessCoeffs =
            LREInterp.cellHessianCoeffs();
        const CompactListList<highOrderInterp::symmTensor3Order>& thirdCoeffs =
            LREInterp.cellThirdDerivCoeffs();

        const label nCells = mesh.nCells();
        // Local rows, global columns (MPIAIJ layout).
        const label nGlobalCells = gNCols();

        std::vector<Triplet> triplets;
        if (compactRows)
        {
            label compactReserve = 0;
            forAll(stencils, cellI)
            {
                compactReserve += stencils[cellI].size() + 1;
            }
            triplets.reserve(compactReserve);
        }
        else
        {
            triplets.reserve(cellChunk * 80);
        }

        SpMat Mlocal(nCells, nGlobalCells);
        bool mInitialised = false;

        auto flushChunk =
        [&]()
        {
            if (triplets.empty())
            {
                return;
            }
            Mlocal.setZero();
            Mlocal.setFromTriplets(triplets.begin(), triplets.end());
            Mlocal.makeCompressed();

            if (!mInitialised)
            {
                M.swap(Mlocal);
                Mlocal.resize(nCells, nGlobalCells);
                mInitialised = true;
            }
            else
            {
                M += Mlocal;
            }
            triplets.clear();
        };

        const label nStencils = stencils.size();

        for
        (
            label cellStart = 0;
            cellStart < nStencils;
            cellStart += cellChunk
        )
        {
            const label cellEnd = min(cellStart + cellChunk, nStencils);

            std::vector<scalar> rowValues;

            for (label cellI = cellStart; cellI < cellEnd; ++cellI)
            {
                // CompactListList::operator[] returns a
        // UList VIEW by value, not a List. Binding it to a const labelList&
        // compiles but yields an EMPTY list, silently dropping the whole
        // stencil. Bind to the returned view type instead.
        const UList<label> curStencil = stencils[cellI];
                const label selfCoeffI = curStencil.size();

                // The coefficient row carries one entry per stencil member
                // plus a trailing one for the cell itself; anything else is a
                // library/adapter mismatch, reported here rather than read out
                // of bounds.
                if (selfCoeffI != gradCoeffs[cellI].size() - 1)
                {
                    FatalErrorInFunction
                        << "Stencil/coefficient mismatch at cell " << cellI
                        << ": stencil size " << selfCoeffI
                        << ", coefficient row " << gradCoeffs[cellI].size()
                        << " (expected row = size + 1)."
                        << abort(FatalError);
                }
                if (compactRows)
                {
                    rowValues.assign(selfCoeffI + 1, 0.0);
                }

                scalar wSum = 0.0;
                forAll(cellQW[cellI], qpI)
                {
                    wSum += cellQW[cellI][qpI];
                }
                wSum = max(wSum, SMALL);

                forAll(cellQP[cellI], qpI)
                {
                    const scalar w = cellQW[cellI][qpI]/wSum;
                    const vector d = cellQP[cellI][qpI] - C[cellI];

                    forAll(curStencil, cI)
                    {
                        scalar coeff = (gradCoeffs[cellI][cI] & d);

                        if (LREInterp.order() >= 2)
                        {
                            coeff +=
                                0.5*quadraticForm(hessCoeffs[cellI][cI], d);
                        }

                        if (LREInterp.order() >= 3)
                        {
                            coeff +=
                                (1.0/6.0)
                               *highOrderInterp::cubicForm
                                (
                                    thirdCoeffs[cellI][cI],
                                    d,
                                    twoD
                                );
                        }

                        if (compactRows)
                        {
                            rowValues[cI] += coefficient*w*coeff;
                        }
                        else
                        {
                            addTripletIfNeeded
                            (
                                triplets,
                                cellI,
                                curStencil[cI],
                                coefficient*w*coeff
                            );
                        }
                    }

                    scalar selfCoeff =
                        1.0 + (gradCoeffs[cellI][selfCoeffI] & d);

                    if (LREInterp.order() >= 2)
                    {
                        selfCoeff +=
                            0.5
                           *quadraticForm(hessCoeffs[cellI][selfCoeffI], d);
                    }

                    if (LREInterp.order() >= 3)
                    {
                        selfCoeff +=
                            (1.0/6.0)
                           *highOrderInterp::cubicForm
                            (
                                thirdCoeffs[cellI][selfCoeffI],
                                d,
                                twoD
                            );
                    }

                    if (compactRows)
                    {
                        rowValues[selfCoeffI] += coefficient*w*selfCoeff;
                    }
                    else
                    {
                        addTripletIfNeeded
                        (
                            triplets,
                            cellI,
                            gCol(cellI),
                            coefficient*w*selfCoeff
                        );
                    }
                }

                if (compactRows)
                {
                    forAll(curStencil, cI)
                    {
                        addTripletIfNeeded
                        (
                            triplets,
                            cellI,
                            curStencil[cI],
                            rowValues[cI]
                        );
                    }

                    addTripletIfNeeded
                    (
                        triplets,
                        cellI,
                        gCol(cellI),
                        rowValues[selfCoeffI]
                    );
                }
            }

            flushChunk();
        }

        if (!mInitialised)
        {
            M.resize(nCells, nGlobalCells);
            M.setZero();
        }

        M.makeCompressed();
        std::vector<Triplet>().swap(triplets);

        reportOperatorFingerprint("M (consistent)", M);
    }

    // Low-order diffusion operator (useHighOrder_Vm = false): two-point
    // orthogonal fluxes over internal faces, Dirichlet boundary faces, and
    // coupled (processor) faces, divided by the cell volume so that K
    // represents div(D grad .) rather than its integral.
    //
    // The coupled branch is the one the serial ancestor did not have: without
    // it every partition's block of the matrix is uncoupled from its
    // neighbours', which does not produce a wrong answer but a KSP breakdown -
    // and the row-sum check cannot see it, because an omitted face contributes
    // nothing to either of its rows. The far cell's data comes from
    // syncTools::swapBoundaryCellList, never from a patch field.
    //
    // Modified for cardiacFoam, ported from the MMS solver: the flux carries
    // BOTH parts of the decomposition in nonOrthFluxVector() - the two-point
    // term along the line joining the cell centres, and the tangential
    // remainder, expanded implicitly through the reconstruction rows. The
    // second one is not optional. A bare two-point flux is INCONSISTENT on
    // triangles and tetrahedra (order 0.81 instead of 2 on the MMS, with an
    // error 230x larger) and with conductivity not aligned with the faces, so
    // it is not selectable. Where D is aligned with an orthogonal mesh the
    // tangential term is identically zero and costs nothing.
    //
    // The reconstruction is taken by reference: createFields.H now always
    // builds it, because this function always needs it.
    void assembleStandardOrthogonalStiffnessMatrix
    (
        const fvMesh& mesh,
        const volTensorField& conductivity,
        const highOrderInterp& LREInterp,
        const scalar stabilisationAlpha,
        SpMat& K
    )
    {
        const vectorField& C = mesh.C();
        const scalarField& V = mesh.V();
        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();

        std::vector<Triplet> triplets;
        triplets.reserve(4*mesh.nFaces());

        forAll(neighbour, faceI)
        {
            const label own = owner[faceI];
            const label nei = neighbour[faceI];

            const tensor Df = 0.5*(conductivity[own] + conductivity[nei]);
            const vector dPN = C[nei] - C[own];
            const scalar a = orthogonalDiffusionCoeff(mesh.Sf()[faceI], dPN, Df);

            addTripletIfNeeded(triplets, own, gCol(own), -a/max(V[own], SMALL));
            addTripletIfNeeded(triplets, own, gCol(nei),  a/max(V[own], SMALL));
            addTripletIfNeeded(triplets, nei, gCol(own),  a/max(V[nei], SMALL));
            addTripletIfNeeded(triplets, nei, gCol(nei), -a/max(V[nei], SMALL));

            // The tangential flux a two-point form would drop. grad V|_f is
            // taken as the average of the two cell gradients, which is what
            // OpenFOAM's correction does; here each of them is expanded into
            // the matrix through its reconstruction row, so the term is
            // IMPLICIT. Skipped where w is zero, which on an aligned orthogonal
            // mesh is every face.
            const vector w = nonOrthFluxVector(mesh.Sf()[faceI], dPN, Df);

            if (w != vector::zero)
            {
                const scalar areaHalf = 0.5*(mag(mesh.Sf()[faceI]) + VSMALL);

                for (const label cellI : {own, nei})
                {
                    addCellGradientDotCoeffs
                    (
                        triplets, own, areaHalf/max(V[own], SMALL), cellI, w,
                        LREInterp
                    );
                    addCellGradientDotCoeffs
                    (
                        triplets, nei, -areaHalf/max(V[nei], SMALL), cellI, w,
                        LREInterp
                    );
                }
            }
        }

        if (stabilisationAlpha > SMALL)
        {
            forAll(neighbour, faceI)
            {
                const label own = owner[faceI];
                const label nei = neighbour[faceI];

                const vector Sf = mesh.Sf()[faceI];
                const scalar area = mag(Sf) + VSMALL;
                const vector n = Sf/area;
                const vector dPN = C[nei] - C[own];
                const tensor Df = 0.5*(conductivity[own] + conductivity[nei]);
                const scalar a =
                    stabilisationFaceCoeff
                    (
                        stabilisationAlpha,
                        area,
                        Df,
                        n,
                        dPN
                    );

                addRhieChowInternalJumpCoeffs
                (
                    triplets,
                    own,
                    a/max(V[own], SMALL),
                    own,
                    nei,
                    dPN,
                    LREInterp
                );

                addRhieChowInternalJumpCoeffs
                (
                    triplets,
                    nei,
                    -a/max(V[nei], SMALL),
                    own,
                    nei,
                    dPN,
                    LREInterp
                );
            }
        }

        // Added for cardiacFoam: data from the cell on the FAR side of every
        // coupled (processor) face. Both are indexed by boundary-face position,
        // i.e. faceI - mesh.nInternalFaces().
        //
        // syncTools::swapBoundaryCellList is used rather than
        // conductivity.boundaryField() or mesh.C().boundaryField(): the latter
        // is a slicedVolVectorField whose coupled values are only correct after
        // an evaluation nothing here performs, and a processor patch field
        // holds the interpolated FACE value, not the neighbouring CELL value.
        List<tensor> nbrConductivity;
        List<label> nbrGlobalCell;
        List<labelList> nbrGradCols;
        List<scalarField> nbrGradCoeffs;
        // Added for cardiacFoam, ported from the MMS solver: far-side
        // reconstruction rows already dotted with the tangential flux vector w
        // of the non-orthogonal correction. A SECOND exchange is needed because
        // w is not dPN: the stabilisation jump ships grad|_far . dPN and this
        // ships grad|_far . w.
        List<labelList> nbrWCols;
        List<scalarField> nbrWCoeffs;

        if (Pstream::parRun())
        {
            syncTools::swapBoundaryCellList
            (
                mesh, conductivity.primitiveField(), nbrConductivity
            );

            labelList myGlobalCell(mesh.nCells());
            forAll(myGlobalCell, cellI)
            {
                myGlobalCell[cellI] = gCol(cellI);
            }
            syncTools::swapBoundaryCellList(mesh, myGlobalCell, nbrGlobalCell);

            // w flips sign across a processor face: the far rank sees -Sf and
            // -d, so its g and e both flip and w = g - (g.e)e flips with g.
            // Building the row with MINUS the local w therefore hands the
            // receiver exactly grad|_there . w_here, usable with the same sign
            // as its own gradient term - the same trick the stabilisation
            // exchange uses with -delta.
            //
            // Unconditional, unlike the stabilisation exchange below: the
            // correction does not depend on alpha. COLLECTIVE, and every rank
            // takes it. A face whose w is zero ships an empty row.
            {
                label lastWPatch = -1;
                vectorField wDelta;

                exchangeCoupledFaceRows
                (
                    mesh,
                    [&]
                    (
                        const label patchI,
                        const label faceI,
                        labelList& cols,
                        scalarField& coeffs
                    )
                    {
                        const fvPatch& patch = mesh.boundary()[patchI];
                        const label own = owner[patch.start() + faceI];

                        // Cached per patch: delta() rebuilds the whole field.
                        if (patchI != lastWPatch)
                        {
                            wDelta = patch.delta();
                            lastWPatch = patchI;
                        }

                        const label bFaceI =
                            patch.start() - mesh.nInternalFaces() + faceI;
                        const tensor Df =
                            0.5*(conductivity[own] + nbrConductivity[bFaceI]);
                        const vector w =
                            nonOrthFluxVector
                            (
                                mesh.Sf().boundaryField()[patchI][faceI],
                                wDelta[faceI],
                                Df
                            );

                        if (w == vector::zero)
                        {
                            cols.clear();
                            coeffs.clear();
                            return;
                        }

                        buildCellGradientDotRow
                        (
                            cols, coeffs, own, -w, LREInterp
                        );
                    },
                    nbrWCols,
                    nbrWCoeffs
                );
            }

            // Added for cardiacFoam: the far cell's gradient row, for the
            // stabilisation jump. Only under alpha > 0.
            //
            // The call is COLLECTIVE and the condition is uniform across ranks:
            // stabilisationAlpha comes from a dictionary.
            if (stabilisationAlpha > SMALL)
            {
                label lastDeltaPatch = -1;
                vectorField sendDelta;

                exchangeCoupledFaceRows
                (
                    mesh,
                    [&]
                    (
                        const label patchI,
                        const label faceI,
                        labelList& cols,
                        scalarField& coeffs
                    )
                    {
                        const fvPatch& patch = mesh.boundary()[patchI];
                        const label own = owner[patch.start() + faceI];

                        // Cached per patch: delta() rebuilds the whole field.
                        if (patchI != lastDeltaPatch)
                        {
                            sendDelta = patch.delta();
                            lastDeltaPatch = patchI;
                        }

                        // Built with MINUS this rank's delta, i.e. with the
                        // RECEIVER's owner-to-neighbour vector, so the row
                        // arrives as grad|_here . dPN_there and can be used
                        // with the same sign as the local gradient term.
                        buildCellGradientDotRow
                        (
                            cols,
                            coeffs,
                            own,
                            -sendDelta[faceI],
                            LREInterp
                        );
                    },
                    nbrGradCols,
                    nbrGradCoeffs
                );
            }
        }

        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            const word bcType =
                patch.lookupPatchField<volScalarField, scalar>("Vm").type();

            if
            (
                bcType == "empty"
             || bcType == zeroGradientFvPatchScalarField::typeName
            )
            {
                continue;
            }

            // Added for cardiacFoam: coupling across a processor boundary.
            //
            // This branch did not exist. The test below is
            // "fixedValue or fixedVoltage" with no else, so a processor patch
            // passed the empty/zeroGradient filter, failed that one and
            // contributed NOTHING - leaving every partition's block of the
            // matrix uncoupled from its neighbours. The symptom is not a wrong
            // answer but a KSP breakdown, and the usual sanity check does not
            // see it: omitting a face adds nothing, so row sums still vanish.
            //
            // Only the owner side is added; the rank on the other side of the
            // face adds its own, exactly as the internal-face loop adds both
            // sides when both cells are local.
            if (patch.coupled())
            {
                // Added for cardiacFoam: exchangeCoupledFaceRows can reach the
                // far rank of a processorFvPatch and nothing else. A cyclic is
                // coupled too, and taking this branch with no row to add would
                // be worse than the omission it replaces.
                if (!isA<processorFvPatch>(patch))
                {
                    FatalErrorInFunction
                        << "Coupled patch " << patch.name() << " of type "
                        << patch.type() << " is not a processor patch; the"
                        << " coupling across it is not implemented."
                        << exit(FatalError);
                }

                const vectorField patchDelta(patch.delta());
                const label bStart = patch.start() - mesh.nInternalFaces();

                forAll(patch, faceI)
                {
                    const label gf = patch.start() + faceI;
                    const label own = owner[gf];
                    const label bFaceI = bStart + faceI;

                    const tensor Df =
                        0.5*(conductivity[own] + nbrConductivity[bFaceI]);

                    const scalar a =
                        orthogonalDiffusionCoeff
                        (
                            mesh.Sf().boundaryField()[patchI][faceI],
                            patchDelta[faceI],
                            Df
                        );

                    addTripletIfNeeded
                    (
                        triplets, own, gCol(own), -a/max(V[own], SMALL)
                    );
                    addTripletIfNeeded
                    (
                        triplets, own, nbrGlobalCell[bFaceI],
                        a/max(V[own], SMALL)
                    );

                    // Added for cardiacFoam, ported from the MMS solver: the
                    // tangential flux. Own half from the local reconstruction
                    // row, far half from the row that arrived already dotted
                    // with w. The far rank computed the same |w| and shipped an
                    // empty row exactly when this w is zero.
                    const vector w =
                        nonOrthFluxVector
                        (
                            mesh.Sf().boundaryField()[patchI][faceI],
                            patchDelta[faceI],
                            Df
                        );

                    if (w != vector::zero)
                    {
                        const scalar areaHalf =
                            0.5
                           *(mag(mesh.Sf().boundaryField()[patchI][faceI])
                           + VSMALL);
                        const scalar scale = areaHalf/max(V[own], SMALL);

                        addCellGradientDotCoeffs
                        (
                            triplets, own, scale, own, w, LREInterp
                        );
                        addRemoteRowCoeffs
                        (
                            triplets, own, scale,
                            nbrWCols[bFaceI], nbrWCoeffs[bFaceI]
                        );
                    }

                    // Modified for cardiacFoam: the stabilisation jump used to
                    // be skipped here, because it needs the neighbouring cell's
                    // reconstruction stencil and that does not cross the
                    // partition. It now does: the far rank builds its own
                    // gradient row and ships (global column, coefficient)
                    // pairs, so nothing but the shared face geometry has to
                    // agree. Inert at the default stabilisationAlpha = 0.
                    if (stabilisationAlpha > SMALL)
                    {
                        const vector Sf =
                            mesh.Sf().boundaryField()[patchI][faceI];
                        const scalar area = mag(Sf) + VSMALL;
                        const vector n = Sf/area;

                        // Orientation-independent (|d.n| and n.D.n), so both
                        // ranks compute the same value with no agreement.
                        const scalar aStab =
                            stabilisationFaceCoeff
                            (
                                stabilisationAlpha,
                                area,
                                Df,
                                n,
                                patchDelta[faceI]
                            );

                        addRhieChowCoupledJumpCoeffs
                        (
                            triplets,
                            own,
                            aStab/max(V[own], SMALL),
                            own,
                            nbrGlobalCell[bFaceI],
                            patchDelta[faceI],
                            nbrGradCols[bFaceI],
                            nbrGradCoeffs[bFaceI],
                            LREInterp
                        );
                    }
                }

                continue;
            }

        }

        K.resize(mesh.nCells(), gNCols());
        K.setFromTriplets(triplets.begin(), triplets.end());
        K.makeCompressed();

        reportOperatorFingerprint("K (orthogonal)", K);
        std::vector<Triplet>().swap(triplets);
    }

    // High-order stiffness matrix K representing the discrete diffusion
    // operator  div(D * grad Vm)  via face-quadrature LRE reconstruction.
    //
    // For each internal face the flux through the face is integrated using
    // Gauss quadrature on the face:
    //
    //     Phi_f = sum_{qp on f}  area * w_qp * ( n . D . grad Vm |_qp )
    //
    // with grad Vm at the quadrature point expressed as a linear combination
    // of the LRE face-stencil values:
    //
    //     grad Vm |_qp = sum_j  gCoeff[face, qp, j] * Vm_j
    //
    // The contribution +Phi_f is added to the owner row (outward flux) and
    // -Phi_f to the neighbour row (inward flux from neighbour's viewpoint).
    // Because gCoeff for the owner cell itself gives a negative contribution
    // in the outward-normal direction, the resulting K is negative semi-
    // definite (negative diagonal, positive off-diagonals), as expected for
    // a Laplacian discretisation.
    void assembleHighOrderStiffnessMatrix
    (
        const fvMesh& mesh,
        const volTensorField& conductivity,
        const highOrderInterp& LREInterp,
        const scalar stabilisationAlpha,
        const label tripletsPerFaceReserve,
        const label faceChunk,
        SpMat& K
    )
    {
        const bool twoD = mesh.nGeometricD() == 2;
        const vectorField& C = mesh.C();
        const surfaceVectorField& Cf = mesh.Cf();
        const CompactListList<point>& faceQP = LREInterp.faceQuadPoints();
        const CompactListList<scalar>& faceQW = LREInterp.faceQuadWeightPhysical();
        const CompactListList<label>& faceStencils =
            LREInterp.globalFaceStencils();
        const List<CompactListList<vector>>& faceGradCoeffs =
            LREInterp.QRGradFaceGPCoeffs();

        const scalarField& V = mesh.V();
        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();

        const label nCells = mesh.nCells();
        const label nInternalFaces = neighbour.size();

        std::vector<Triplet> triplets;
        const label reserveFaces = min(max(mesh.nFaces(), label(1)), faceChunk);
        triplets.reserve(reserveFaces*tripletsPerFaceReserve);

        SpMat Klocal(nCells, gNCols());
        bool kInitialised = false;

        auto flushChunk =
        [&]()
        {
            if (triplets.empty())
            {
                return;
            }

            Klocal.setZero();
            Klocal.setFromTriplets(triplets.begin(), triplets.end());
            Klocal.makeCompressed();

            if (!kInitialised)
            {
                K.swap(Klocal);
                Klocal.resize(nCells, gNCols());
                kInitialised = true;
            }
            else
            {
                K += Klocal;
            }

            triplets.clear();
        };

        // -------- Internal faces -------------------------------------------
        for
        (
            label faceStart = 0;
            faceStart < nInternalFaces;
            faceStart += faceChunk
        )
        {
            const label faceEnd = min(faceStart + faceChunk, nInternalFaces);

            for (label faceI = faceStart; faceI < faceEnd; ++faceI)
            {
                const label own = owner[faceI];
                const label nei = neighbour[faceI];
                const vector Sf = mesh.Sf()[faceI];
                const scalar area = mag(Sf) + VSMALL;
                const vector n = Sf/area;            // unit normal owner -> nei
                // Modified for cardiacFoam: must be UList<label> BY VALUE, not
                // a labelList reference - see the note on
                // CompactListList in the header comment.
                const UList<label> curStencil = faceStencils[faceI];

                forAll(faceQP[faceI], qpI)
                {
                    const scalar w = faceQW[faceI][qpI];

                    forAll(curStencil, cI)
                    {
                        const label col = curStencil[cI];
                        const vector gCoeff = faceGradCoeffs[faceI][qpI][cI];
                        // Modified for cardiacFoam: the FACE tensor is the
                        // average of the two cells, not the owner's. The same
                        // fluxCoeff goes into both rows below, so using
                        // conductivity[own] made the face's diffusivity depend
                        // on which cell the mesh happens to call the owner.
                        // Invisible with a uniform D, where the two tensors are
                        // equal - and with a non-uniform D it made K depend on
                        // the PARTITION, because a face that is internal in
                        // serial becomes a coupled face whose owner is
                        // rank-local. Measured before the fix on an alternating
                        // D: max|dK|/max|K| = 2.1e-01 at np=2 and np=4, against
                        // 4.8e-15 with a uniform D (test M15).
                        //
                        // With a uniform D this is bit-identical to the old
                        // expression: 0.5*(x + x) is exact in floating point.
                        const tensor Df =
                            0.5*(conductivity[own] + conductivity[nei]);

                        const scalar fluxCoeff =
                            // Modified for cardiacFoam: no area
                            // factor, w is a physical weight.
                            w*(n & (Df & gCoeff));

                        addTripletIfNeeded
                        (
                            triplets,
                            own,
                            col,
                            fluxCoeff/max(V[own], SMALL)
                        );

                        addTripletIfNeeded
                        (
                            triplets,
                            nei,
                            col,
                            -fluxCoeff/max(V[nei], SMALL)
                        );
                    }
                }
            }

            if (stabilisationAlpha > SMALL)
            {
                for (label faceI = faceStart; faceI < faceEnd; ++faceI)
                {
                    const label own = owner[faceI];
                    const label nei = neighbour[faceI];

                    const vector Sf = mesh.Sf()[faceI];
                    const scalar area = mag(Sf) + VSMALL;
                    const vector n = Sf/area;
                    const vector dPN = C[nei] - C[own];
                    const tensor Df = 0.5*(conductivity[own] + conductivity[nei]);
                    const scalar a =
                        stabilisationFaceCoeff
                        (
                            stabilisationAlpha,
                            area,
                            Df,
                            n,
                            dPN
                        );
                    const point xf = Cf[faceI];

                    addCellTaylorExtrapolationCoeffs
                    (
                        triplets,
                        own,
                        a/max(V[own], SMALL),
                        nei,
                        xf,
                        LREInterp,
                        C,
                        twoD
                    );
                    addCellTaylorExtrapolationCoeffs
                    (
                        triplets,
                        own,
                        -a/max(V[own], SMALL),
                        own,
                        xf,
                        LREInterp,
                        C,
                        twoD
                    );
                    addCellTaylorExtrapolationCoeffs
                    (
                        triplets,
                        nei,
                        -a/max(V[nei], SMALL),
                        nei,
                        xf,
                        LREInterp,
                        C,
                        twoD
                    );
                    addCellTaylorExtrapolationCoeffs
                    (
                        triplets,
                        nei,
                        a/max(V[nei], SMALL),
                        own,
                        xf,
                        LREInterp,
                        C,
                        twoD
                    );
                }
            }

            flushChunk();
        }

        // Added for cardiacFoam: data needed by the face-jump stabilisation on
        // coupled (processor) faces, which the boundary loop below did not
        // handle at all - the stabilisation was gated on fixedValue/
        // fixedVoltage, so a processor patch contributed the flux term and
        // nothing else. Invisible in serial, and invisible in the operator
        // fingerprint too: an omitted term contributes nothing to either of its
        // rows, so the row sums still vanish.
        //
        // Two pieces cross the boundary:
        //   - the far cell's conductivity, so both ranks form the same face
        //     coefficient (syncTools gives the CELL value; a processor patch
        //     field would give the interpolated FACE value);
        //   - the far cell's Taylor reconstruction row, evaluated at the shared
        //     face centre by the rank that owns the cell.
        //
        // COLLECTIVE, and deliberately here: the internal-face loop above
        // assembles in chunks, and this call inside a chunk loop would deadlock
        // as soon as two ranks had a different number of chunks.
        List<tensor> nbrConductivity;
        List<labelList> nbrTaylorCols;
        List<scalarField> nbrTaylorCoeffs;

        // Modified for cardiacFoam: the conductivity swap is NOT gated on the
        // stabilisation any more. The coupled-face flux below needs the far
        // cell's tensor to form the face average, whatever alpha is - see the
        // note at the internal-face flux. The low-order assembly swaps it on
        // Pstream::parRun() alone for the same reason.
        //
        // Still here rather than inside the chunk loop, and still collective:
        // every rank takes this branch or none does.
        if (Pstream::parRun())
        {
            syncTools::swapBoundaryCellList
            (
                mesh, conductivity.primitiveField(), nbrConductivity
            );
        }

        if (stabilisationAlpha > SMALL && Pstream::parRun())
        {
            exchangeCoupledFaceRows
            (
                mesh,
                [&]
                (
                    const label patchI,
                    const label faceI,
                    labelList& cols,
                    scalarField& coeffs
                )
                {
                    const fvPatch& patch = mesh.boundary()[patchI];
                    const label own = owner[patch.start() + faceI];

                    buildCellTaylorExtrapolationRow
                    (
                        cols,
                        coeffs,
                        own,
                        Cf.boundaryField()[patchI][faceI],
                        LREInterp,
                        C,
                        twoD
                    );
                },
                nbrTaylorCols,
                nbrTaylorCoeffs
            );
        }

        // -------- Boundary faces -------------------------------------------
        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            const word bcType =
                patch.lookupPatchField<volScalarField, scalar>("Vm").type();

            // Insulating (Neumann) and empty patches contribute zero flux.
            // Modified for cardiacFoam: Dirichlet patches were removed from
            // both solvers (refused at start-up), so what reaches the loop
            // below is processor patches only.
            if
            (
                bcType == "empty"
             || bcType == zeroGradientFvPatchScalarField::typeName
            )
            {
                continue;
            }

            // Added for cardiacFoam: processor patches only. A cyclic is
            // coupled too and would otherwise fall through into the flux loop
            // below and be treated as if it were Dirichlet - silently wrong.
            if (patch.coupled() && !isA<processorFvPatch>(patch))
            {
                FatalErrorInFunction
                    << "Coupled patch " << patch.name() << " of type "
                    << patch.type() << " is not a processor patch; the"
                    << " stabilisation jump across it is not implemented."
                    << exit(FatalError);
            }

            const bool stabiliseCoupled =
                stabilisationAlpha > SMALL
             && Pstream::parRun()
             && patch.coupled();

            // Owner-to-neighbour vector across a coupled face. patch.delta() is
            // used rather than C[nei] - C[own] because mesh.C() is a
            // slicedVolVectorField whose coupled values are only valid after an
            // evaluation nothing here performs.
            vectorField patchDelta;
            if (stabiliseCoupled)
            {
                patchDelta = patch.delta();
            }

            const label bStart = patch.start() - mesh.nInternalFaces();

            forAll(patch, faceI)
            {
                const label globalFaceI = patch.start() + faceI;
                const label own = owner[globalFaceI];
                const vector Sf = mesh.Sf().boundaryField()[patchI][faceI];
                const scalar area = mag(Sf) + VSMALL;
                const vector n = Sf/area;
                // Modified for cardiacFoam: must be UList<label> BY VALUE, not
                // a labelList reference - see the note on
                // CompactListList in the header comment.
                // Modified for cardiacFoam: the face tensor, matching the
                // internal-face loop. This loop is reached by BOTH processor
                // patches and Dirichlet patches, and the two are not the same
                // face:
                //
                //   - across a processor cut there IS a cell on the far side, so
                //     the face tensor is the average of the two. This is the half
                //     of the fix that closes the partition dependence; changing
                //     only the internal-face loop fixes the SERIAL operator and
                //     leaves parallel building the same physical face with two
                //     different tensors, one per rank.
                //   - on a physical boundary there is no neighbour cell, so the
                //     owner's tensor IS the face tensor. Indexing nbrConductivity
                //     here would read an entry swapBoundaryCellList never filled.
                //
                // patch.coupled() implies processorFvPatch (anything else
                // FatalErrors above), which implies parRun, so nbrConductivity is
                // populated whenever this branch is taken.
                const tensor Dface =
                    patch.coupled()
                  ? 0.5*(conductivity[own] + nbrConductivity[bStart + faceI])
                  : conductivity[own];

                const UList<label> curStencil = faceStencils[globalFaceI];

                forAll(faceQP[globalFaceI], qpI)
                {
                    const scalar w = faceQW[globalFaceI][qpI];

                    forAll(curStencil, cI)
                    {
                        const label col = curStencil[cI];
                        const vector gCoeff =
                            faceGradCoeffs[globalFaceI][qpI][cI];
                        const scalar fluxCoeff =
                            // Modified for cardiacFoam: no area
                            // factor, w is a physical weight.
                            w*(n & (Dface & gCoeff));

                        addTripletIfNeeded(triplets, own, col, fluxCoeff/max(V[own], SMALL));
                    }
                }

                // Added for cardiacFoam: the same jump as the internal-face
                // loop, written from the LOCAL cell's point of view:
                //
                //     a/V[own] * ( T_far(xf) - T_own(xf) )
                //
                // where T_c is the Taylor reconstruction of cell c. The
                // internal loop writes this with a + sign into the owner row
                // and a - sign into the neighbour row; those two signs are the
                // same expression seen from either cell, so no owner/neighbour
                // convention has to be agreed across the partition. The far
                // rank adds its own copy.
                if (stabiliseCoupled)
                {
                    const label bFaceI = bStart + faceI;
                    const point xf = Cf.boundaryField()[patchI][faceI];
                    const tensor Df =
                        0.5*(conductivity[own] + nbrConductivity[bFaceI]);

                    // Orientation-independent (|d.n|, n.D.n): both ranks get
                    // the same value with no agreement needed.
                    const scalar a =
                        stabilisationFaceCoeff
                        (
                            stabilisationAlpha,
                            area,
                            Df,
                            n,
                            patchDelta[faceI]
                        );

                    addRemoteRowCoeffs
                    (
                        triplets,
                        own,
                        a/max(V[own], SMALL),
                        nbrTaylorCols[bFaceI],
                        nbrTaylorCoeffs[bFaceI]
                    );

                    addCellTaylorExtrapolationCoeffs
                    (
                        triplets,
                        own,
                        -a/max(V[own], SMALL),
                        own,
                        xf,
                        LREInterp,
                        C,
                        twoD
                    );
                }
            }
        }

        flushChunk();

        if (!kInitialised)
        {
            K.resize(nCells, gNCols());
            K.setZero();
        }

        K.makeCompressed();

        reportOperatorFingerprint("K (high order)", K);
        std::vector<Triplet>().swap(triplets);
    }
