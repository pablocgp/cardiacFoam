/*---------------------------------------------------------------------------*\
    cardiacFoam - shared code of the monodomain solvers
    (highOrderElectroActivationFoamImplicitPETScDistributed and
    highOrderManufacturedFDAImplicitPETScDistributed)

    Krylov solvers: the in-house restarted GMRES (globally
    reduced, right-preconditioned), the JFNK preconditioner and the JFNK
    linear solve.

    INCLUDED INSIDE THE SOLVER'S ANONYMOUS NAMESPACE, after its system,
    OpenFOAM, Eigen and PETSc includes, through monodomainCommon.H. It
    therefore includes nothing itself and has internal linkage in each
    executable. Code here is used verbatim by both solvers: a change made here
    changes both, which is the point.
\*---------------------------------------------------------------------------*/

    // Added for cardiacFoam: the square diagonal block of a local-rows /
    // global-columns matrix, i.e. the columns this rank owns, renumbered to
    // start at zero.
    //
    // Eigen's factorisations need a SQUARE matrix. In parallel AImplicit is
    // nLocal x nGlobal - 320 x 640 on two ranks - and handing that to
    // IncompleteLUT aborts inside Eigen on a dimension assertion. Factorising
    // the diagonal block instead is precisely block-Jacobi ILU, the same
    // substitution petscParallelPcTypeName() makes for the assembled path.
    //
    // In serial gRowStart is 0 and the block is the whole matrix, so the
    // factorisation is bit-identical to what it was.
    SpMat diagonalBlock(const SpMat& A)
    {
        const label nLocal = A.rows();
        std::vector<Triplet> triplets;
        triplets.reserve(A.nonZeros());

        for (label row = 0; row < nLocal; ++row)
        {
            for (SpMat::InnerIterator it(A, row); it; ++it)
            {
                const label col = it.col() - globalRowStart();
                if (col >= 0 && col < nLocal)
                {
                    triplets.emplace_back(row, col, it.value());
                }
            }
        }

        SpMat block(nLocal, nLocal);
        block.setFromTriplets(triplets.begin(), triplets.end());
        block.makeCompressed();
        return block;
    }

    // Restarted GMRES(m) with Modified Gram-Schmidt orthogonalisation.
    //
    //   Solves  A * x = b  iteratively, where A is provided as a matrix-vector
    //   product functor (matVec). The Krylov subspace is rebuilt up to
    //   maxRestarts times to avoid storing arbitrarily many basis vectors when
    //   convergence is slow.
    //
    //   - krylovDim   : maximum dimension m of each inner subspace (GMRES(m))
    //   - maxRestarts : maximum number of outer restarts (0 = no restart)
    //   - tolerance   : stop when ||r_k||/||b|| < tolerance
    //   - iterations  : total number of inner Arnoldi steps performed (out)
    //   - estimatedError : final relative residual ||r_k||/||b|| (out)
    template<class MatrixVectorProduct>
    EigVec solveGMRES
    (
        const MatrixVectorProduct& matVec,
        const EigVec& b,
        const label krylovDim,
        const label maxRestarts,
        const scalar tolerance,
        label& iterations,
        scalar& estimatedError,
        const std::function<EigVec(const EigVec&)>* applyPreconditioner = nullptr
    )
    {
        // Right preconditioning (Phase A3): solve  A M^{-1} u = b,  x = M^{-1} u.
        // The Arnoldi operator becomes A M^{-1} and the correction V*y is mapped
        // back through M^{-1}; the least-squares residual still equals the true
        // residual ||b - A x||, and the restart residual uses the raw A. When no
        // preconditioner is supplied, applyPrec is the identity (unchanged GMRES).
        auto applyPrec = [&](const EigVec& v) -> EigVec
        {
            return applyPreconditioner ? (*applyPreconditioner)(v) : v;
        };
        const label n = b.size();

        // Modified for cardiacFoam: b.size() is this rank's block, so a Krylov
        // dimension capped by it would differ per rank - and the restart/break
        // logic below runs inside a COLLECTIVE solve, so ranks would leave the
        // loop at different iterations and deadlock. Cap by the global count.
        const label nGlobal = returnReduce(n, sumOp<label>());
        const label m = min(max(krylovDim, label(1)), label(nGlobal));
        const label maxOuter = max(maxRestarts, label(0));

        // Initial iterate x_0 = 0  =>  initial residual r_0 = b - A*0 = b
        EigVec x = EigVec::Zero(n);
        EigVec r = b;
        const scalar bNorm = gNorm(b);
        scalar beta = bNorm;
        iterations = 0;
        estimatedError = 1.0;

        if (bNorm <= SMALL)
        {
            // Right-hand side is essentially zero; trivial solution
            estimatedError = 0.0;
            return x;
        }

        // Storage for Arnoldi basis V and upper-Hessenberg matrix H.
        // We allocate once and reuse across restarts to avoid reallocation.
        std::vector<EigVec> V(m + 1, EigVec::Zero(n));
        Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> H =
            Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero
            (m + 1, m);

        for (label restart = 0; restart <= maxOuter; ++restart)
        {
            // Seed the Krylov subspace with the (normalised) current residual.
            V[0] = r/beta;
            H.setZero();
            bool toleranceMet = false;
            EigVec bestX = x;

            for (label j = 0; j < m; ++j)
            {
                ++iterations;

                // Arnoldi step: form w = A * M^{-1} * V[j]
                EigVec w = matVec(applyPrec(V[j]));

                // Modified Gram-Schmidt orthogonalisation.
                //
                //   Each H(i,j) is computed against the *currently updated* w,
                //   so loss-of-orthogonality is bounded by O(eps * kappa(A_m))
                //   instead of O(eps * kappa(A_m)^2) as in classical GS.
                //   This is important for stiff ionic systems (TNNP) where
                //   the Jacobian can be poorly conditioned.
                for (label i = 0; i <= j; ++i)
                {
                    H(i, j) = gDot(V[i], w);
                    w -= H(i, j)*V[i];
                }

                // Sub-diagonal entry of the Hessenberg matrix.
                // Happy breakdown (H(j+1,j) ~= 0) means the Krylov space is
                // already A-invariant and the current iterate is exact.
                H(j + 1, j) = gNorm(w);
                if (H(j + 1, j) > SMALL && j + 1 < m + 1)
                {
                    V[j + 1] = w/H(j + 1, j);
                }

                // Least-squares solve for the projected problem:
                //     min || beta*e_1 - H_j * y ||
                // We use Eigen's column-pivoted Householder QR. For small
                // Hessenberg blocks (j+2 x j+1) this is fast and robust.
                Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> Hj =
                    H.block(0, 0, j + 2, j + 1);
                EigVec g = EigVec::Zero(j + 2);
                g[0] = beta;
                const EigVec y = Hj.colPivHouseholderQr().solve(g);

                // Reconstruct candidate iterate: x_j = x_restart + M^{-1}(V_j * y)
                // The correction lives in the preconditioned space and is mapped
                // back with M^{-1} (identity when unpreconditioned).
                EigVec Vy = EigVec::Zero(n);
                for (label i = 0; i <= j; ++i)
                {
                    Vy += y[i]*V[i];
                }
                EigVec xj = x + applyPrec(Vy);

                // In GMRES theory the residual norm of the least-squares
                // problem equals the true residual norm:
                //     ||g - H_j*y|| = ||b - A*x_j||
                // Using it avoids an extra matVec call per inner iteration.
                const scalar relResidual =
                    (g - Hj*y).norm()/max(bNorm, SMALL);
                bestX = xj;
                estimatedError = relResidual;

                if (relResidual <= tolerance || H(j + 1, j) <= SMALL)
                {
                    toleranceMet = true;
                    break;
                }
            }

            // Commit the best iterate found in this inner cycle.
            x = bestX;
            if (toleranceMet) break;
            if (restart == maxOuter) break;

            // Restart: recompute the *true* residual r = b - A*x and re-seed
            // the next Krylov cycle. This costs one extra matrix-vector
            // product but recovers from subspace exhaustion.
            r = b - matVec(x);
            beta = gNorm(r);
            estimatedError = beta/max(bNorm, SMALL);
            if (beta <= SMALL*bNorm) break;
        }

        return x;
    }

    // The JFNK preconditioner P ~ J, shared by both solvers.
    //
    //   none          no preconditioning
    //   AImplicit     P = AImplicit, the linear part of the Jacobian
    //   diagonalIion  P = AImplicit - theta*diag(d source/d Vm), which adds the
    //                 local reaction Jacobian - the same linearisation as the
    //                 diagonalIion nonlinear method. The default: on the
    //                 upstroke it is the better approximation of J.
    //
    // Rebuilt at the first nonlinear iteration of every time step and then
    // every updateFrequency iterations. The sparsity pattern is that of
    // AImplicit and never changes, so after the first build only the values
    // are copied and only the numerical factorisation is redone.
    //
    // PETSc backend: factorised by its own KSP (preonly + ilu by default;
    // block-Jacobi in parallel, see petscParallelPcTypeName). Eigen backend:
    // ILUT of the LOCAL DIAGONAL BLOCK of P (see diagonalBlock), i.e.
    // block-Jacobi in parallel and the whole matrix in serial. The MMS solver
    // used to factorise the whole local matrix, which in parallel is
    // rectangular and made its Eigen JFNK path serial only; the electro
    // solver used AImplicit alone, refreshed only when dt changed. Both now do
    // this.
    class JfnkPreconditioner
    {
        const word type_;
        const bool petsc_;
        const word kspType_;
        const word pcType_;
        const scalar tolerance_;
        const label maxIterations_;
        const label restart_;
        const word optionsPrefix_;
        const bool useOptions_;
        const label fillFactor_;
        const scalar dropTolerance_;
        const label updateFrequency_;

        SpMat P_;
        bool pInitialised_ = false;
        PetscKspMatrixSolver petscSolver_;
        Eigen::IncompleteLUT<scalar> ilu_;

        bool readyThisStep_ = false;
        label setupCorr_ = -1;

        std::function<EigVec(const EigVec&)> applyFn_;

        // Timings, for the solvers' own summaries.
        scalar setupWallTime_ = 0.0;
        label setups_ = 0;
        scalar applyWallTime_ = 0.0;
        label applications_ = 0;

    public:

        JfnkPreconditioner
        (
            const word& type,
            const bool petsc,
            const word& kspType,
            const word& pcType,
            const scalar tolerance,
            const label maxIterations,
            const label restart,
            const word& optionsPrefix,
            const bool useOptions,
            const label fillFactor,
            const scalar dropTolerance,
            const label updateFrequency
        )
        :
            type_(type),
            petsc_(petsc),
            kspType_(kspType),
            pcType_(pcType),
            tolerance_(tolerance),
            maxIterations_(maxIterations),
            restart_(restart),
            optionsPrefix_(optionsPrefix),
            useOptions_(useOptions),
            fillFactor_(fillFactor),
            dropTolerance_(dropTolerance),
            updateFrequency_(max(updateFrequency, label(1)))
        {
            requireOneOf
            (
                "jfnkPreconditioner", type_, {"none", "AImplicit", "diagonalIion"}
            );
            applyFn_ = [this](const EigVec& r) -> EigVec { return apply(r); };
        }

        JfnkPreconditioner(const JfnkPreconditioner&) = delete;
        void operator=(const JfnkPreconditioner&) = delete;

        bool active() const
        {
            return type_ != "none";
        }

        //- True when update() needs the source derivative (diagonalIion)
        bool usesSourceDerivative() const
        {
            return type_ == "diagonalIion";
        }

        //- Call once at the start of every time step's nonlinear loop
        void beginStep()
        {
            readyThisStep_ = false;
        }

        bool needsUpdate(const label corr) const
        {
            return
                active()
             && (!readyThisStep_ || corr - setupCorr_ >= updateFrequency_);
        }

        //- Rebuild P from AImplicit and, for diagonalIion, the source
        //  derivative at the current iterate (null otherwise)
        void update
        (
            const SpMat& AImplicit,
            const scalarField* sourceDerivative,
            const scalar theta,
            const label corr
        )
        {
            const auto t0 = std::chrono::steady_clock::now();

            // A changed pattern cannot happen today (M and K are constant and
            // dt only scales them), but if it ever did, copying values into
            // the old pattern would be silently wrong: rebuild instead.
            const bool samePattern =
                pInitialised_
             && P_.nonZeros() == AImplicit.nonZeros()
             && P_.rows() == AImplicit.rows()
             && P_.cols() == AImplicit.cols();

            if (!samePattern)
            {
                P_ = AImplicit;
                pInitialised_ = true;
            }
            else
            {
                std::copy
                (
                    AImplicit.valuePtr(),
                    AImplicit.valuePtr() + AImplicit.nonZeros(),
                    P_.valuePtr()
                );
            }

            if (usesSourceDerivative())
            {
                if (!sourceDerivative)
                {
                    FatalErrorInFunction
                        << "jfnkPreconditioner diagonalIion needs the source"
                        << " derivative." << abort(FatalError);
                }
                addDiagonalToMatrix(*sourceDerivative, -theta, P_);
            }

            if (petsc_)
            {
                if (!petscSolver_.isInitialised() || !samePattern)
                {
                    petscSolver_.reset
                    (
                        P_,
                        kspType_,
                        pcType_,
                        tolerance_,
                        maxIterations_,
                        restart_,
                        optionsPrefix_,
                        useOptions_,
                        scalar(fillFactor_),
                        dropTolerance_
                    );
                }
                else
                {
                    petscSolver_.updateValues(P_);
                }
            }
            else
            {
                ilu_.setDroptol(dropTolerance_);
                ilu_.setFillfactor(fillFactor_);
                ilu_.compute(diagonalBlock(P_));

                if (ilu_.info() != Eigen::Success)
                {
                    FatalErrorInFunction
                        << "JFNK ILUT preconditioner setup failed"
                        << exit(FatalError);
                }
            }

            readyThisStep_ = true;
            setupCorr_ = corr;

            setupWallTime_ +=
                std::chrono::duration<scalar>
                (
                    std::chrono::steady_clock::now() - t0
                ).count();
            ++setups_;
        }

        //- z = P^{-1} r
        EigVec apply(const EigVec& r)
        {
            const auto t0 = std::chrono::steady_clock::now();

            EigVec z(r.size());
            if (petsc_)
            {
                label iterations = 0;
                scalar error = GREAT;
                z = petscSolver_.solve(r, iterations, error);
            }
            else
            {
                z = ilu_.solve(r);
            }

            applyWallTime_ +=
                std::chrono::duration<scalar>
                (
                    std::chrono::steady_clock::now() - t0
                ).count();
            ++applications_;

            return z;
        }

        //- The application as a function, or null when inactive
        const std::function<EigVec(const EigVec&)>* applyPtr() const
        {
            return active() ? &applyFn_ : nullptr;
        }

        scalar setupWallTime() const { return setupWallTime_; }
        label setups() const { return setups_; }
        scalar applyWallTime() const { return applyWallTime_; }
        label applications() const { return applications_; }
    };


    // Controls of the inner Krylov solve of JFNK, J delta = -R.
    struct JfnkKrylovControls
    {
        bool petsc = true;
        word kspType;
        word pcType;
        label restart = 30;
        label maxKrylovIterations = 30;
        label maxRestarts = 0;
        scalar tolerance = 1.0e-6;
        word optionsPrefix;
        bool useOptions = true;
    };


    // Solve J delta = -R for the Newton step, shared by both solvers.
    //
    // PETSc backend: KSP over a matrix-free MatShell, with the preconditioner
    // as a PCShell applied on the RIGHT (see PetscShellKspSolver). Eigen
    // backend: the in-house restarted GMRES, also right-preconditioned and
    // with globally reduced inner products, so it runs at any rank count.
    EigVec solveJfnkLinearSystem
    (
        const std::function<EigVec(const EigVec&)>& matVec,
        const EigVec& R,
        JfnkPreconditioner& pc,
        PetscShellKspSolver& shell,
        const JfnkKrylovControls& c,
        label& iterations,
        scalar& error
    )
    {
        if (c.petsc)
        {
            if (!shell.isInitialised())
            {
                shell.initialise
                (
                    static_cast<label>(R.size()),
                    c.kspType,
                    c.pcType,
                    c.restart,
                    max(c.maxKrylovIterations*(c.maxRestarts + 1), label(1)),
                    c.tolerance,
                    c.optionsPrefix,
                    c.useOptions,
                    pc.active()
                );
            }

            return shell.solve(matVec, pc.applyPtr(), -R, iterations, error);
        }

        return solveGMRES
        (
            matVec,
            -R,
            c.maxKrylovIterations,
            c.maxRestarts,
            c.tolerance,
            iterations,
            error,
            pc.applyPtr()
        );
    }

    // How a persistent linear solver should treat the matrix on a given solve
    // (Phase B: reuse factorisations instead of rebuilding them every call).
    //   Rebuild    : first use, or the sparsity pattern / dt changed -> full
    //                allocation + symbolic factorisation + numeric setup.
    //   ValuesOnly : same pattern, values changed (e.g. diagonalIion adds a
    //                fresh diagonal each nonlinear iteration) -> refresh values
    //                and re-factorise, but reuse the allocation and pattern.
    //   Reuse      : matrix identical to the previous solve (Picard with a
    //                constant operator) -> no matrix touch at all, just solve.
    enum class MatrixUpdate
    {
        Rebuild,
        ValuesOnly,
        Reuse
    };

    // Persistent Eigen factorisations reused across time steps / nonlinear
    // iterations (mirror of the persistent PetscKspMatrixSolver for the Eigen
    // backend). The symbolic analysis (analyzePattern) is done once; only the
    // numeric factorisation is repeated when the values change.
    struct PersistentEigenSolvers
    {
        Eigen::SparseLU<SpMat> lu;
        bool luAnalyzed = false;
        Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<scalar>> bicgstab;

        // Incomplete-LU preconditioner for the matrix-free JFNK Jacobian
        // (Phase A3). An exact LU of the wide-stencil high-order operator has
        // impractical fill-in; ILUT gives a cheap, bounded-fill approximation
        // of M^{-1} that is very effective here since AImplicit ~ M/dt is
        // strongly diagonally dominant for small dt.
    };

    // Solve A x = b with the Eigen backend: SparseLU (direct) or BiCGSTAB with
    // an ILUT preconditioner. Serial by construction - Eigen has no notion of a
    // distributed matrix - so this path is available only at one rank.
    //
    // updatePolicy is what makes repeated solves cheap: the symbolic analysis
    // depends only on the sparsity pattern and is done once, the numeric
    // factorisation is redone when the values change, and a wholly unchanged
    // matrix skips both.
    EigVec solveSparseSystemEigen
    (
        const SpMat& A,
        const EigVec& b,
        const word& linearSolver,
        const scalar tol,
        const label maxIter,
        label& linearIterations,
        scalar& linearError,
        PersistentEigenSolvers& persistent,
        const MatrixUpdate updatePolicy
    )
    {
        if (linearSolver == "SparseLU")
        {
            Eigen::SparseLU<SpMat>& solver = persistent.lu;

            // Symbolic analysis depends only on the sparsity pattern, which is
            // constant for the whole run: do it once. The numeric factorisation
            // is redone whenever the values change (Rebuild / ValuesOnly) and
            // skipped entirely when the matrix is unchanged (Reuse).
            if (updatePolicy != MatrixUpdate::Reuse)
            {
                if (!persistent.luAnalyzed)
                {
                    solver.analyzePattern(A);
                    persistent.luAnalyzed = true;
                }

                solver.factorize(A);

                if (solver.info() != Eigen::Success)
                {
                    FatalErrorInFunction
                        << "SparseLU factorization failed"
                        << exit(FatalError);
                }
            }

            EigVec x = solver.solve(b);
            linearIterations = 1;
            linearError = 0.0;

            if (solver.info() != Eigen::Success)
            {
                FatalErrorInFunction
                    << "SparseLU solve failed"
                    << exit(FatalError);
            }

            return x;
        }
        else if (linearSolver == "BiCGSTAB")
        {
            Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<scalar>>& solver =
                persistent.bicgstab;
            solver.setTolerance(tol);
            solver.setMaxIterations(maxIter);

            // compute() builds the incomplete-LU preconditioner; reuse it when
            // the matrix is unchanged.
            if (updatePolicy != MatrixUpdate::Reuse)
            {
                solver.compute(A);

                if (solver.info() != Eigen::Success)
                {
                    FatalErrorInFunction
                        << "BiCGSTAB setup failed"
                        << exit(FatalError);
                }
            }

            EigVec x = solver.solve(b);
            linearIterations = solver.iterations();
            linearError = solver.error();

            if (solver.info() != Eigen::Success)
            {
                FatalErrorInFunction
                    << "BiCGSTAB solve failed"
                    << exit(FatalError);
            }

            return x;
        }

        FatalErrorInFunction
            << "Unknown implicit linear solver: " << linearSolver << nl
            << "Valid options are SparseLU or BiCGSTAB"
            << exit(FatalError);

        return EigVec();
    }

    // Front door for the assembled linear solves (Picard and diagonalIion):
    // dispatches to PETSc or to Eigen according to linearSolverBackend, and
    // returns the iteration count and final relative residual through its
    // out-parameters.
    //
    // The two backends take the same convergence controls so that a case file
    // produces an algorithmically equivalent run either way; they differ only
    // in how the system is solved. In parallel only the PETSc path is usable.
    EigVec solveSparseSystem
    (
        const SpMat& A,
        const EigVec& b,
        const word& linearSolverBackend,
        const word& linearSolver,
        const word& petscKspType,
        const word& petscPcType,
        const label petscRestart,
        const word& petscOptionsPrefix,
        const bool petscUseOptions,
        const scalar tol,
        const label maxIter,
        label& linearIterations,
        scalar& linearError,
        PetscKspMatrixSolver& petscSolver,
        PersistentEigenSolvers& eigenSolvers,
        const MatrixUpdate updatePolicy,
        const bool reusePreconditioner
    )
    {
        if (usesPetscBackend(linearSolverBackend))
        {
            // Reuse the persistent KSP/factorisation across time steps and
            // nonlinear iterations (Phase B). The matrix is only reassembled
            // into PETSc when the pattern changes (Rebuild) or its values
            // change (ValuesOnly); a constant operator (Reuse) skips both.
            switch (updatePolicy)
            {
                case MatrixUpdate::Rebuild:
                {
                    word kspType;
                    word pcType;
                    resolvePetscKspPc
                    (
                        linearSolver, petscKspType, petscPcType, kspType, pcType
                    );

                    petscSolver.reset
                    (
                        A,
                        kspType,
                        pcType,
                        tol,
                        maxIter,
                        petscRestart,
                        petscOptionsPrefix,
                        petscUseOptions
                    );
                    break;
                }

                case MatrixUpdate::ValuesOnly:
                {
                    petscSolver.updateValues(A, reusePreconditioner);
                    break;
                }

                case MatrixUpdate::Reuse:
                {
                    // Nothing to do: the cached operator is still valid.
                    break;
                }
            }

            return petscSolver.solve(b, linearIterations, linearError);
        }

        // Added for cardiacFoam, ported from the MMS solver: the Eigen backend
        // is serial only, and the comments above said so without enforcing it.
        // A carries LOCAL rows and GLOBAL columns (the MPIAIJ layout), so in
        // parallel it is rectangular by design and Eigen's SparseLU/BiCGSTAB
        // die on an eigen_assert naming a file inside the library - not this
        // solver, not the key that selected them, and not the word "parallel".
        if (Pstream::parRun())
        {
            FatalErrorInFunction
                << "linearSolverBackend = " << linearSolverBackend
                << " selects the Eigen backend, which is serial only." << nl
                << "The assembled matrix has " << A.rows()
                << " local rows and " << A.cols()
                << " global columns, so Eigen would be handed a rectangular"
                << " system." << nl
                << "Set 'linearSolverBackend PETSc;' to run on "
                << Pstream::nProcs() << " ranks."
                << exit(FatalError);
        }

        return solveSparseSystemEigen
        (
            A,
            b,
            linearSolver,
            tol,
            maxIter,
            linearIterations,
            linearError,
            eigenSolvers,
            updatePolicy
        );
    }
