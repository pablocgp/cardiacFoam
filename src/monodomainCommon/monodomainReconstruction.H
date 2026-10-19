/*---------------------------------------------------------------------------*\
    cardiacFoam - shared code of the monodomain solvers
    (highOrderElectroActivationFoamImplicitPETScDistributed and
    highOrderManufacturedFDAImplicitPETScDistributed)

    Taylor reconstruction from the MLS derivatives.

    INCLUDED INSIDE THE SOLVER'S ANONYMOUS NAMESPACE, after its system,
    OpenFOAM, Eigen and PETSc includes, through monodomainCommon.H. It
    therefore includes nothing itself and has internal linkage in each
    executable. Code here is used verbatim by both solvers: a change made here
    changes both, which is the point.
\*---------------------------------------------------------------------------*/

    // d . H . d for a symmetric tensor, written out to exploit the symmetry.
    // The second-order term of a Taylor reconstruction.
    scalar quadraticForm(const symmTensor& H, const vector& d)
    {
        return
            H.xx()*d.x()*d.x()
          + 2.0*H.xy()*d.x()*d.y()
          + 2.0*H.xz()*d.x()*d.z()
          + H.yy()*d.y()*d.y()
          + 2.0*H.yz()*d.y()*d.z()
          + H.zz()*d.z()*d.z();
    }

    // Evaluate a Taylor reconstruction at offset d from a cell centre:
    //
    //     c + grad.d + (1/2) d.H.d + (1/6) T3(d,d,d)
    //
    // The Hessian and third-derivative pointers are null below the corresponding
    // polynomial order, which is how the same routine serves p1, p2 and p3. In
    // 2-D the cubic form drops its out-of-plane components (twoD).
    scalar reconstructFromTaylor
    (
        const scalar c,
        const vector& grad,
        const symmTensor* H,
        const highOrderInterp::symmTensor3Order* T3,
        const vector& d,
        const bool twoD
    )
    {
        scalar val = c + (grad & d);

        if (H)
        {
            val += 0.5*quadraticForm(*H, d);
        }

        if (T3)
        {
            val += (1.0/6.0)*highOrderInterp::cubicForm(*T3, d, twoD);
        }

        return val;
    }

    // Reconstruct K cell-centred fields at the ionic quadrature points, at the
    // order of the interpolator - the shared core of every quadrature-point
    // reconstruction in both solvers (Vm, the three MMS states, the TNNP
    // states).
    //
    // Fused and threaded (ported from the MMS solver). The library's
    // grad()/hessian()/thirdDeriv() are three passes over the mesh, each with
    // its own halo exchange, field allocation and serial stencil walk, plus a
    // fourth loop for the Taylor expansion; here the halo is exchanged once
    // per field and one threaded pass gathers the stencil, forms every
    // derivative and expands. Measured on the MMS at p3 on 300x300, the
    // unfused state reconstruction was 70% of the ionic work and capped OpenMP
    // at 1.20x. The arithmetic is identical to movingLeastSquares::grad/
    // secondGrad/thirdGrad, including the trailing coefficient that carries
    // the cell's own contribution at index stencil.size().
    //
    // All K fields go through ONE walk of the stencil: each entry is located
    // (local index, or the processor/slot of a remote cell, a hash lookup - as
    // movingLeastSquares::fieldValue does, which is private to the library)
    // once and its coefficients are loaded once, for every field. The first
    // shared version took one field per call, which repeated the walk K times
    // and made the MMS state reconstruction 30% slower than the solver's own
    // three-field loop it replaced (3-D M18: 3.54 s -> 4.61 s). Per field the
    // accumulation order is unchanged, so the results are bit-identical to
    // reconstructing the fields one at a time.
    //
    // write(ip, k, value) stores each result, so each solver keeps its own
    // data layout without a copy; ip is the flat quadrature-point index and k
    // the field. Only the cells with writeCellMask true are written (all when
    // it is null).
    //
    // limitToNeighbourRange clamps each value to the [min, max] of the cell
    // and its face neighbours, including those across a processor face (a
    // Barth-Jespersen-style limiter), field by field. The physiological solver
    // uses it for the stiff ionic states, whose p3 reconstruction across an
    // under-resolved front overshoots and can drive the ionic algebra to
    // overflow; the MMS solver never does, since on a smooth field it only
    // clips legitimate extrema.
    template<class Writer>
    void reconstructCellFieldsAtIntegrationPoints
    (
        const fvMesh& mesh,
        const UList<const scalarField*>& fields,
        const highOrderInterp& LREInterp,
        const CompactListList<point>& cellQuadPoints,
        const Writer& write,
        const boolList* writeCellMask,
        const bool limitToNeighbourRange,
        const bool useOpenMP,
        const label openMPThreshold,
        scalar* haloWallTime = nullptr,
        scalar* loopWallTime = nullptr
    )
    {
        const vectorField& C = mesh.C();
        const bool twoD = mesh.nGeometricD() == 2;
        const label nCells = mesh.nCells();
        const label order = LREInterp.order();
        const label nFields = fields.size();

        // Where each cell's quadrature points start in the flat array, taken
        // from the container so that no running counter ties the iterations.
        const labelList& ipOffsets = cellQuadPoints.offsets();

        const auto tHaloStart = std::chrono::steady_clock::now();

        // Halo exchange: MPI, once per field, outside the threaded region.
        const auto& stencilData = LREInterp.mls().stencilData();
        const globalIndex& globalCells = stencilData.globalCells();
        const Map<FixedList<label, 2>>& remoteLoc =
            stencilData.remoteCellLocation();
        const CompactListList<label>& stencils = stencilData.cellsStencil();
        List<List<Field<scalar>>> remote(nFields);
        for (label k = 0; k < nFields; ++k)
        {
            remote[k] = stencilData.remoteFieldPerProc(*fields[k]);
        }

        // The limiter's admissible range must include neighbours on the far
        // side of a processor face: cellCells() stops at the cut, and a cell
        // next to one would otherwise be clamped harder than in serial.
        List<scalarField> coupledLo;
        List<scalarField> coupledHi;
        if (limitToNeighbourRange)
        {
            coupledLo.setSize(nFields, scalarField(nCells, GREAT));
            coupledHi.setSize(nFields, scalarField(nCells, -GREAT));

            if (Pstream::parRun())
            {
                for (label k = 0; k < nFields; ++k)
                {
                    List<scalar> nbrValues;
                    syncTools::swapBoundaryCellList(mesh, *fields[k], nbrValues);

                    forAll(mesh.boundary(), patchI)
                    {
                        const fvPatch& p = mesh.boundary()[patchI];
                        if (!p.coupled())
                        {
                            continue;
                        }

                        const labelUList& faceCells = p.faceCells();
                        const label bStart = p.start() - mesh.nInternalFaces();

                        forAll(faceCells, i)
                        {
                            const label cellI = faceCells[i];
                            const scalar v = nbrValues[bStart + i];
                            coupledLo[k][cellI] = min(coupledLo[k][cellI], v);
                            coupledHi[k][cellI] = max(coupledHi[k][cellI], v);
                        }
                    }
                }
            }
        }
        const labelListList& cellCells =
            limitToNeighbourRange ? mesh.cellCells() : labelListList::null();

        const CompactListList<vector>& gradCoeffs = LREInterp.QRGradCoeffs();
        const CompactListList<symmTensor>& hessCoeffs =
            LREInterp.cellHessianCoeffs();
        const CompactListList<highOrderInterp::symmTensor3Order>& thirdCoeffs =
            LREInterp.cellThirdDerivCoeffs();

        const bool haveHess = (order >= 2);
        const bool haveThird = (order >= 3);

        const auto tLoopStart = std::chrono::steady_clock::now();

        #pragma omp parallel if(useOpenMP && nCells >= openMPThreshold)
        {
            // Per-thread accumulators, allocated once per call, not per cell.
            List<vector> g(nFields);
            List<symmTensor> H(nFields);
            List<highOrderInterp::symmTensor3Order> T(nFields);
            scalarList v(nFields);
            scalarList vc(nFields);
            scalarList lo(nFields);
            scalarList hi(nFields);

            // Each iteration writes only its own cell's quadrature points.
            #pragma omp for schedule(static)
            for (label cellI = 0; cellI < nCells; ++cellI)
            {
                if (writeCellMask && !(*writeCellMask)[cellI])
                {
                    continue;
                }

                const UList<label>& stencil = stencils[cellI];
                const UList<vector>& gc = gradCoeffs[cellI];

                for (label k = 0; k < nFields; ++k)
                {
                    g[k] = vector::zero;
                    H[k] = symmTensor::zero;
                    T[k] = highOrderInterp::symmTensor3Order::zero;
                }

                forAll(stencil, cI)
                {
                    // Locate the entry once, then read it in every field.
                    const label gid = stencil[cI];
                    if (globalCells.isLocal(gid))
                    {
                        const label l = globalCells.toLocal(gid);
                        for (label k = 0; k < nFields; ++k)
                        {
                            v[k] = (*fields[k])[l];
                        }
                    }
                    else
                    {
                        const FixedList<label, 2>& loc = remoteLoc[gid];
                        for (label k = 0; k < nFields; ++k)
                        {
                            v[k] = remote[k][loc[0]][loc[1]];
                        }
                    }

                    const vector& gcI = gc[cI];
                    for (label k = 0; k < nFields; ++k)
                    {
                        g[k] += gcI*v[k];
                    }

                    if (haveHess)
                    {
                        const symmTensor& hcI = hessCoeffs[cellI][cI];
                        for (label k = 0; k < nFields; ++k)
                        {
                            H[k] += hcI*v[k];
                        }
                    }

                    if (haveThird)
                    {
                        const highOrderInterp::symmTensor3Order& tcI =
                            thirdCoeffs[cellI][cI];
                        for (label k = 0; k < nFields; ++k)
                        {
                            T[k] += tcI*v[k];
                        }
                    }
                }

                // Trailing entry: the cell's own contribution.
                const label last = stencil.size();
                for (label k = 0; k < nFields; ++k)
                {
                    vc[k] = (*fields[k])[cellI];
                    g[k] += gc[last]*vc[k];

                    if (haveHess)
                    {
                        H[k] += hessCoeffs[cellI][last]*vc[k];
                    }

                    if (haveThird)
                    {
                        T[k] += thirdCoeffs[cellI][last]*vc[k];
                    }
                }

                if (limitToNeighbourRange)
                {
                    const labelList& nbr = cellCells[cellI];
                    for (label k = 0; k < nFields; ++k)
                    {
                        const scalarField& f = *fields[k];
                        lo[k] = vc[k];
                        hi[k] = vc[k];
                        forAll(nbr, j)
                        {
                            lo[k] = min(lo[k], f[nbr[j]]);
                            hi[k] = max(hi[k], f[nbr[j]]);
                        }
                        lo[k] = min(lo[k], coupledLo[k][cellI]);
                        hi[k] = max(hi[k], coupledHi[k][cellI]);
                    }
                }

                const vector& xc = C[cellI];
                const label base = ipOffsets[cellI];
                const UList<point>& qp = cellQuadPoints[cellI];

                forAll(qp, qI)
                {
                    const vector d = qp[qI] - xc;

                    for (label k = 0; k < nFields; ++k)
                    {
                        scalar value =
                            reconstructFromTaylor
                            (
                                vc[k],
                                g[k],
                                haveHess ? &H[k] : nullptr,
                                haveThird ? &T[k] : nullptr,
                                d,
                                twoD
                            );

                        if (limitToNeighbourRange)
                        {
                            value = min(max(value, lo[k]), hi[k]);
                        }

                        write(base + qI, k, value);
                    }
                }
            }
        }

        const auto tEnd = std::chrono::steady_clock::now();
        if (haloWallTime)
        {
            *haloWallTime +=
                std::chrono::duration<scalar>(tLoopStart - tHaloStart).count();
        }
        if (loopWallTime)
        {
            *loopWallTime +=
                std::chrono::duration<scalar>(tEnd - tLoopStart).count();
        }
    }

    // Vm at every ionic quadrature point. With a high-order Vm, the shared
    // reconstruction core above at the interpolator's order; otherwise a
    // linear extrapolation from the finite-volume cell gradient. This is what
    // makes the source term high-order: a nonlinear function evaluated at the
    // cell centre and multiplied by the volume is first-order however good the
    // diffusion operator is.
    void reconstructVmAtIionIntegrationPoints
    (
        const volScalarField& Vm,
        const Switch useHighOrderVm,
        const highOrderInterp& LREInterp_Vm,
        const highOrderInterp& LREInterp_Iion,
        scalarField& VmIntegrationPoints,
        const Switch useOpenMP,
        const label openMPThreshold,
        const Switch profileTimings,
        FineGrainedTimings& fineTimings
    )
    {
        const fvMesh& mesh = Vm.mesh();
        const CompactListList<point>& cellIionQuadP =
            LREInterp_Iion.cellQuadPoints();

        if (useHighOrderVm)
        {
            auto write = [&](const label ip, const label, const scalar v)
            {
                VmIntegrationPoints[ip] = v;
            };

            const scalarField* field = &Vm.primitiveField();

            reconstructCellFieldsAtIntegrationPoints
            (
                mesh,
                UList<const scalarField*>(&field, 1),
                LREInterp_Vm,
                cellIionQuadP,
                write,
                nullptr,
                false,
                useOpenMP,
                openMPThreshold,
                profileTimings ? &fineTimings.vmReconDerivWallTime : nullptr,
                profileTimings ? &fineTimings.vmReconExpandWallTime : nullptr
            );

            if (profileTimings)
            {
                ++fineTimings.vmReconCalls;
            }
            return;
        }

        // Low order: fvc::grad is a library call; the expansion is threaded.
        const vectorField& C = mesh.C();
        const label nCells = mesh.nCells();
        const labelList& ipOffsets = cellIionQuadP.offsets();
        const scalarField& VmI = Vm.primitiveField();

        tmp<volVectorField> tGradVm = fvc::grad(Vm);
        const vectorField& gradVm = tGradVm->internalField();

        #pragma omp parallel for schedule(static) \
                if(useOpenMP && nCells >= openMPThreshold)
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            const scalar Vc = VmI[cellI];
            const vector& gradVc = gradVm[cellI];
            const vector& xc = C[cellI];
            const label base = ipOffsets[cellI];
            const UList<point>& qp = cellIionQuadP[cellI];

            forAll(qp, qI)
            {
                VmIntegrationPoints[base + qI] = Vc + (gradVc & (qp[qI] - xc));
            }
        }
    }

    // A stored cell level of Vm - a BDF history level or an ESDIRK stage
    // value - at the ionic quadrature points, by the same reconstruction that
    // gives V^n there, so each level is exactly what its own step had at those
    // points. Lets the Gauss-point state ODE paths drive the ODEs with the
    // history polynomial without keeping any level at the quadrature points
    // between steps (stateODEGaussPointHistory). scratch is a field with Vm's
    // mesh and patch types; it is overwritten.
    void reconstructVmLevelAtIionIntegrationPoints
    (
        const scalarField& cellValues,
        volScalarField& scratch,
        const Switch useHighOrderVm,
        const highOrderInterp& LREInterp_Vm,
        const highOrderInterp& LREInterp_Iion,
        scalarField& VmLevelIntegrationPoints,
        const Switch useOpenMP,
        const label openMPThreshold,
        const Switch profileTimings,
        FineGrainedTimings& fineTimings
    )
    {
        scratch.primitiveFieldRef() = cellValues;
        scratch.correctBoundaryConditions();
        reconstructVmAtIionIntegrationPoints
        (
            scratch,
            useHighOrderVm,
            LREInterp_Vm,
            LREInterp_Iion,
            VmLevelIntegrationPoints,
            useOpenMP,
            openMPThreshold,
            profileTimings,
            fineTimings
        );
    }
