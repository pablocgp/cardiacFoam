/*---------------------------------------------------------------------------*\
    cardiacFoam - shared code of the monodomain solvers
    (highOrderElectroActivationFoamImplicitPETScDistributed and
    highOrderManufacturedFDAImplicitPETScDistributed)

    Types, global cell addressing, option validation, global
    norms, field<->vector conversions and triplet helpers.

    INCLUDED INSIDE THE SOLVER'S ANONYMOUS NAMESPACE, after its system,
    OpenFOAM, Eigen and PETSc includes, through monodomainCommon.H. It
    therefore includes nothing itself and has internal linkage in each
    executable. Code here is used verbatim by both solvers: a change made here
    changes both, which is the point.
\*---------------------------------------------------------------------------*/

    using SpMat = Eigen::SparseMatrix<scalar, Eigen::RowMajor>;
    using Triplet = Eigen::Triplet<scalar>;
    using EigVec = Eigen::Matrix<scalar, Eigen::Dynamic, 1>;

    // Global row offset of this rank's block of cells, i.e.
    // globalIndex::localStart(). Constant for the whole run (one mesh, one
    // decomposition), so it is set once, by initialiseGlobalCells(), rather
    // than threaded through every linear-solver signature. Zero in serial.
    label gRowStart = 0;
    bool gRowStartSet = false;

    // Read through an accessor that refuses to serve the initial 0. Every
    // PETSc row insertion reads it; a code path that touched a matrix before
    // initialiseGlobalCells() would offset by 0 on EVERY rank, so ranks > 0
    // would write their rows on top of rank 0's: a corrupt operator and a
    // plausible wrong answer, never an error. This makes it impossible rather
    // than merely unlikely.
    label globalRowStart()
    {
        if (!gRowStartSet)
        {
            FatalErrorInFunction
                << "The global row offset was read before it was set."
                << " It is assigned once in main by initialiseGlobalCells();"
                << " reaching a PETSc matrix before that offsets every rank's"
                << " rows by 0 and corrupts the operator."
                << exit(FatalError);
        }

        return gRowStart;
    }

    // Global cell addressing used for matrix COLUMNS. The assembly routines are
    // free functions, so a file-scope handle avoids threading an extra argument
    // through all of them. It is the same numbering the reconstructions use
    // (a globalIndex over the local cell counts, rank by rank).
    autoPtr<globalIndex> gCellsPtr;

    // Local cell index -> global column index.
    inline label gCol(const label localCellI)
    {
        return gCellsPtr->toGlobal(localCellI);
    }

    // Total number of columns (global cell count).
    inline label gNCols()
    {
        return gCellsPtr->totalSize();
    }

    // Build the global cell addressing and publish this rank's row offset.
    // Called once in main, after the mesh is read and before anything
    // assembles a matrix. In serial toGlobal() is the identity, totalSize() is
    // mesh.nCells() and the offset is zero, so the serial path is unchanged.
    void initialiseGlobalCells(const fvMesh& mesh)
    {
        gCellsPtr.reset(new globalIndex(mesh.nCells()));
        gRowStart = gCellsPtr->localStart();
        gRowStartSet = true;
    }

    // Modified for cardiacFoam, ported from the MMS solver: every dictionary
    // option has exactly ONE accepted spelling, compared case-sensitively, and
    // anything else stops the run at start-up with the list of valid values.
    //
    // This replaces lowerWord, which matched user input case-insensitively and
    // had accumulated aliases on top of it: SparseLU/sparselu/lu, Picard/picard,
    // diagonalIion/diagonal/localDiagonal, RKF45/RKT45/rkf45, off/false/nopc,
    // ilut, jakobi, amg, boomeramg and more. Two spellings for one thing is how
    // the SparseLU bug hid in the MMS: one test compared case-sensitively and
    // another did not, so the same word took a different path in each place.
    void requireOneOf
    (
        const char* key,
        const word& value,
        std::initializer_list<const char*> valid
    )
    {
        for (const char* v : valid)
        {
            if (value == v)
            {
                return;
            }
        }

        std::string options;
        for (const char* v : valid)
        {
            if (!options.empty())
            {
                options += ", ";
            }
            options += v;
        }

        FatalErrorInFunction
            << "Unknown " << key << ": " << value << nl
            << "Valid options are " << options.c_str()
            << " (exact spelling, case-sensitive)."
            << exit(FatalError);
    }

    // Added for cardiacFoam: globally reduced vector norms and inner product.
    //
    // An EigVec holds only this rank's block of cells, so Eigen's norm() and
    // dot() are local quantities. Used unreduced in a convergence test or in a
    // finite-difference step size they make the result depend on how the mesh
    // was decomposed, which looks like a discretisation problem rather than a
    // parallel bug. In serial these reduce to the plain Eigen calls.
    scalar gSquaredNorm(const EigVec& v)
    {
        scalar s = v.squaredNorm();
        reduce(s, sumOp<scalar>());
        return s;
    }

    // Global Euclidean norm of a distributed vector (each rank holds its own
    // block of rows).
    scalar gNorm(const EigVec& v)
    {
        return std::sqrt(gSquaredNorm(v));
    }

    // Global inner product of two distributed vectors. Needed wherever a
    // convergence test or an orthogonalisation would otherwise see only this
    // rank's block.
    scalar gDot(const EigVec& a, const EigVec& b)
    {
        scalar s = a.dot(b);
        reduce(s, sumOp<scalar>());
        return s;
    }

    // Copy the internal field of an OpenFOAM volScalarField into an Eigen
    // vector. In parallel this is this rank's block only, which is exactly the
    // row range its matrices own.
    EigVec fieldToEigVec(const volScalarField& fld)
    {
        const scalarField& f = fld.primitiveField();
        EigVec v(f.size());

        forAll(f, cellI)
        {
            v[cellI] = f[cellI];
        }

        return v;
    }

    // Inverse of fieldToEigVec: write a solution vector back into a field.
    // Boundary conditions are not evaluated here; the caller does it when
    // needed.
    void eigVecToField(const EigVec& v, volScalarField& fld)
    {
        scalarField& f = fld.primitiveFieldRef();

        forAll(f, cellI)
        {
            f[cellI] = v[cellI];
        }
    }

    // ||r|| / ||reference||, both global. The nonlinear stopping test for Vm
    // and the linear one for the Krylov solvers.
    scalar relativeL2Norm(const EigVec& r, const EigVec& reference)
    {
        // Modified for cardiacFoam: both norms are over LOCAL blocks only.
        return gNorm(r)/(gNorm(reference) + SMALL);
    }

    // Append a matrix entry, dropping values below SMALL. The filter is on the
    // INDIVIDUAL contribution, not on the accumulated one, so it is independent
    // of the order in which contributions arrive and of any chunking of the
    // assembly.
    void addTripletIfNeeded
    (
        std::vector<Triplet>& triplets,
        const label row,
        const label col,
        const scalar value
    )
    {
        if (mag(value) > SMALL)
        {
            triplets.emplace_back(row, col, value);
        }
    }

    // Added for cardiacFoam: pour a row received from another rank into the
    // triplets. The columns are already global, so nothing is translated here.
    void addRemoteRowCoeffs
    (
        std::vector<Triplet>& triplets,
        const label row,
        const scalar scale,
        const labelList& cols,
        const scalarField& coeffs
    )
    {
        if (mag(scale) <= SMALL)
        {
            return;
        }

        forAll(cols, i)
        {
            addTripletIfNeeded(triplets, row, cols[i], scale*coeffs[i]);
        }
    }

    // Project the per-integration-point ionic currents back to a cell-
    // centred field via the quadrature-weighted average:
    //
    //     Iion_cell = sum_q ( w_q * Iion_q )  /  sum_q w_q
    //
    // In-place clamp of a scalar field to a closed interval [vMin, vMax].
    //
    // Used as a protective barrier on the Vm vector passed into the ionic
    // ODE solver: when JFNK probes a candidate iterate far outside the
    // physiological range (~-100 mV to +60 mV), TNNP's gating dynamics
    // become extremely stiff and RKF45 stalls. Clamping the input keeps
    // the integrator in a regime where its time-step controller behaves.
    inline void clampPhysical
    (
        scalarField& v,
        const scalar vMin,
        const scalar vMax
    )
    {
        forAll(v, i)
        {
            if (v[i] < vMin) v[i] = vMin;
            else if (v[i] > vMax) v[i] = vMax;
        }
    }

    // Relative L2 change between two iterates, ||current - previous|| / ||current||,
    // reduced over all ranks. The two solvers used to normalise differently
    // (the MMS by the PREVIOUS iterate); the current one is the agreed
    // reference. Two fields of different sizes cannot be compared: that is a
    // programming error, reported rather than read out of bounds.
    // Modified for cardiacFoam: the denominator can be held at or above
    // denominatorFloor, a global L2 norm (so the result stays uniform across
    // ranks). A field that crosses zero - the net ionic current does, on the
    // plateau of the action potential - drives ||current|| to zero while the
    // iteration noise in it stays put, and the plain relative increment then
    // meets no tolerance at all. Found in a 0-D TNNP run with tolerances of
    // 1e-8: every scheme with dt <= 3.125e-6 stopped converging at the same
    // instant, t = 0.21157 s. With the default floor of 0 the result is
    // unchanged bit for bit.
    scalar relativeL2Difference
    (
        const scalarField& current,
        const scalarField& previous,
        const scalar denominatorFloor = 0.0
    )
    {
        if (current.size() != previous.size())
        {
            FatalErrorInFunction
                << "Cannot compare fields with different sizes: "
                << current.size() << " and " << previous.size()
                << abort(FatalError);
        }

        scalar num = 0.0;
        scalar den = 0.0;

        forAll(current, cellI)
        {
            const scalar diff = current[cellI] - previous[cellI];
            num += diff*diff;
            den += current[cellI]*current[cellI];
        }

        reduce(num, sumOp<scalar>());
        reduce(den, sumOp<scalar>());

        return std::sqrt(num)/(max(std::sqrt(den), denominatorFloor) + SMALL);
    }

    // Added for cardiacFoam: globally reduced fingerprint of an assembled
    // operator. This is the sharpest parallel check available, and it costs one
    // reduction on a run that need not advance in time at all.
    //
    // nnzGlobal must match the serial value EXACTLY: movingLeastSquares builds
    // its stencils geometrically, so they are invariant to how the mesh was
    // partitioned. A dropped processor face, an insufficient halo for a high
    // polynomial order, or a stencil that disagrees across a cut all change it.
    //
    // Note what maxAbsRowSum does NOT catch: omitting a face contributes
    // nothing to either of its rows, so the rows still sum to zero. That is why
    // the missing processor-face branch in the orthogonal assembly showed up as
    // a KSP breakdown rather than as a failed sanity check.
    void reportOperatorFingerprint(const word& name, const SpMat& A)
    {
        label nnz = A.nonZeros();
        scalar sumAll = 0.0;
        scalar sumAbs = 0.0;
        scalar maxAbsRowSum = 0.0;

        for (label row = 0; row < label(A.rows()); ++row)
        {
            scalar rowSum = 0.0;
            for (SpMat::InnerIterator it(A, row); it; ++it)
            {
                sumAll += it.value();
                sumAbs += mag(it.value());
                rowSum += it.value();
            }
            maxAbsRowSum = max(maxAbsRowSum, mag(rowSum));
        }

        reduce(nnz, sumOp<label>());
        reduce(sumAll, sumOp<scalar>());
        reduce(sumAbs, sumOp<scalar>());
        reduce(maxAbsRowSum, maxOp<scalar>());

        Info<< "    [operator " << name << "] nnzGlobal = " << nnz
            << "  sum = " << sumAll
            << "  sumAbs = " << sumAbs
            << "  maxAbsRowSum = " << maxAbsRowSum << endl;

        if (nnz == 0)
        {
            FatalErrorInFunction
                << "Operator " << name << " is empty." << abort(FatalError);
        }
    }

    // The stopping test of every nonlinear iteration - Picard, diagonalIion,
    // JFNK and the ESDIRK stages - in both solvers. The coupled residual of
    // the step equation, the Vm and Iion increments, and the states when
    // nonlinearRequireStatesConvergence asks for them. The MMS solver used to
    // stop Picard on the Vm and state increments alone; one test for every
    // method means a method comparison compares methods, not stopping rules.
    //
    // How a nonlinear iteration ends. The coupled residual can stop at a floor
    // above nonlinearTolerance - the true residual a preconditioned Krylov
    // solve leaves, or plain round-off on a large stiff system - and then
    // every further iteration repeats the same solve up to the iteration cap
    // (found on the MMS: Picard, ilu, dt = 2e-3, 2000 identical iterations per
    // step at a residual of 3.66e-12 against a tolerance of 1e-12). When the
    // increments already pass and the residual did not fall since the last
    // iteration, the iterate is as converged as this operator allows. The
    // increment tests are never loosened, so the worst case is the old
    // increments-only criterion of the MMS Picard; the caller counts these
    // exits and reports them.
    enum class NonlinearExit { Continue, Converged, ResidualFloor };

    NonlinearExit nonlinearIterationExit
    (
        const label corr,
        const label minIterations,
        const scalar coupledResidual,
        const scalar coupledResidualPrevIter,
        const scalar coupledTolerance,
        const scalar VmIncrement,
        const scalar VmTolerance,
        const scalar IionIncrement,
        const scalar IionTolerance,
        const bool requireStates,
        const scalar maxStateIncrement,
        const scalar statesTolerance
    )
    {
        const bool incrementsPass =
            corr + 1 >= minIterations
         && VmIncrement <= VmTolerance
         && IionIncrement <= IionTolerance
         && (!requireStates || maxStateIncrement <= statesTolerance);

        if (!incrementsPass)
        {
            return NonlinearExit::Continue;
        }
        if (coupledResidual <= coupledTolerance)
        {
            return NonlinearExit::Converged;
        }
        if (coupledResidual >= coupledResidualPrevIter)
        {
            return NonlinearExit::ResidualFloor;
        }
        return NonlinearExit::Continue;
    }

    // The caller's side of nonlinearIterationExit: true when the iteration
    // stops; a floor exit is counted and printed.
    bool nonlinearIterationStops
    (
        const NonlinearExit how,
        const scalar coupledResidual,
        const scalar coupledTolerance,
        label& nResidualFloorExits
    )
    {
        if (how == NonlinearExit::ResidualFloor)
        {
            ++nResidualFloorExits;
            Info<< "       Coupled residual at its floor (" << coupledResidual
                << " > nonlinearTolerance " << coupledTolerance
                << ", no longer falling) with the increments converged: "
                << "stopping" << nl;
        }
        return how != NonlinearExit::Continue;
    }

    // Wall time and call counts of the phases INSIDE the time loop, kept as one
    // struct rather than ten loose scalars so that writeSummary takes one more
    // argument instead of ten more positional ones.
    //
    // The pair that matters is nonlinearEval (the ionic ODEs and the source,
    // which is what the OpenMP threads split) against sparseLinearSolve (PETSc,
    // which is what the MPI ranks split): their ratio is the parallelisable
    // fraction of the loop MEASURED, rather than inferred backwards from an
    // observed speedup. The JFNK entries stay zero on a Picard run.
    //
    // Only accumulated when profileTimings is on. They are still written out
    // when it is off - as zeros, which a reader can tell apart from a key that
    // is not there at all.
    struct FineGrainedTimings
    {
        scalar nonlinearEvalWallTime = 0.0;
        label nonlinearEvalCalls = 0;
        scalar sparseLinearSolveWallTime = 0.0;
        label sparseLinearSolveCalls = 0;
        scalar gmresWallTime = 0.0;
        label gmresCalls = 0;
        scalar preconditionerSetupWallTime = 0.0;
        label preconditionerSetups = 0;
        scalar preconditionerApplyWallTime = 0.0;
        label preconditionerApplications = 0;

        // Las dos mitades de reconstructStatesAtIionIntegrationPoints, que es
        // lo unico caro del camino cellCentredReconstruct que NO tiene hilos.
        //
        // Separarlas decide donde vale la pena trabajar y donde no:
        //   stateReconDeriv  - grad/hessian/thirdDeriv, dentro de
        //                      movingLeastSquares (solids4foam). Hilarlo
        //                      significa tocar biblioteca compartida.
        //   stateReconExpand - la expansion de Taylor en los puntos de Gauss,
        //                      en este archivo. Se hila quitando el contador
        //                      corrido, que es un cambio chico y local.
        // Si la primera domina, el cambio chico no sirve de nada.
        scalar stateReconDerivWallTime = 0.0;
        scalar stateReconExpandWallTime = 0.0;
        label stateReconCalls = 0;

        // Lo mismo para reconstructVmAtIionIntegrationPoints, que recibio la
        // misma fusion. Se miden por separado y no se suman a las de arriba
        // porque son dos caminos independientes: el de Vm corre SIEMPRE que
        // useHighOrder_Iion este activo, mientras el de estados corre solo en
        // cellCentredReconstruct. Sumarlos esconderia cual de los dos paga.
        scalar vmReconDerivWallTime = 0.0;
        scalar vmReconExpandWallTime = 0.0;
        label vmReconCalls = 0;
    };

    // Peak resident set size of this process in MB, from getrusage. Linux
    // reports ru_maxrss in KiB and macOS in bytes, hence the platform split.
    scalar currentPeakRSSMB()
    {
        struct rusage usage;
        if (getrusage(RUSAGE_SELF, &usage) != 0)
        {
            return -1.0;
        }

#if defined(__APPLE__)
        return scalar(usage.ru_maxrss)/(1024.0*1024.0);
#else
        return scalar(usage.ru_maxrss)/1024.0;
#endif
    }

    // Print the peak RSS at a named point of the run. The high-order setup
    // (stencils, per-quadrature-point QR factorisations) is the memory peak of
    // these cases, so the checkpoints bracket it.
    void logMemoryCheckpoint(const word& stage)
    {
        Info<< "Memory checkpoint [" << stage << "]: peakRSS_MB = "
            << currentPeakRSSMB() << nl;
    }

    // Modified for cardiacFoam: one spelling per mode, lumped and consistent.
    // The aliases diagonal and consistentHO are gone (see requireOneOf).
    word normalizedMassMatrixType(const word& massMatrix)
    {
        requireOneOf("massMatrix", massMatrix, {"lumped", "consistent"});
        return massMatrix;
    }

    // Whether the memory-optimised assembly is in force. "auto" turns it on for
    // a 3-D high-order Vm above memoryOptimizationCellThreshold cells - counted
    // GLOBALLY. Both solvers used to compare the LOCAL count, so the mode
    // flipped with the rank count and an RSS or wall-time comparison across
    // ranks did not compare the same configuration.
    bool memoryOptimizationIsEffective
    (
        const word& mode,
        const bool useHighOrderVm,
        const label dim,
        const fvMesh& mesh,
        const label cellThreshold
    )
    {
        requireOneOf("memoryOptimization", mode, {"auto", "on", "off"});

        const label nCellsGlobal = returnReduce(mesh.nCells(), sumOp<label>());

        return
            mode == "on"
         || (
                mode == "auto"
             && useHighOrderVm
             && dim == 3
             && nCellsGlobal >= cellThreshold
            );
    }
