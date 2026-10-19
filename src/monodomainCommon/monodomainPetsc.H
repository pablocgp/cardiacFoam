/*---------------------------------------------------------------------------*\
    cardiacFoam - shared code of the monodomain solvers
    (highOrderElectroActivationFoamImplicitPETScDistributed and
    highOrderManufacturedFDAImplicitPETScDistributed)

    PETSc session, distributed matrices and products, and the
    resolution of KSP/PC names.

    INCLUDED INSIDE THE SOLVER'S ANONYMOUS NAMESPACE, after its system,
    OpenFOAM, Eigen and PETSc includes, through monodomainCommon.H. It
    therefore includes nothing itself and has internal linkage in each
    executable. Code here is used verbatim by both solvers: a change made here
    changes both, which is the point.
\*---------------------------------------------------------------------------*/

    // Abort with the calling context if a PETSc call returned an error code.
    // Every PETSc entry point below goes through this: PETSc reports failures by
    // return value, so an unchecked call fails silently and the symptom appears
    // much later, in the solution rather than at the call site.
    void checkPetscError(const PetscErrorCode ierr, const char* context)
    {
        if (ierr)
        {
            FatalErrorInFunction
                << "PETSc call failed in " << context
                << " with error code " << ierr
                << exit(FatalError);
        }
    }

    // Copy a PETSc Vec back into this rank's block of an Eigen vector.
    void copyPetscVecToEigVec(Vec src, EigVec& dst)
    {
        const PetscScalar* values = nullptr;
        checkPetscError(VecGetArrayRead(src, &values), "VecGetArrayRead");

        for (PetscInt i = 0; i < dst.size(); ++i)
        {
            dst[i] = static_cast<scalar>(values[i]);
        }

        checkPetscError
        (
            VecRestoreArrayRead(src, &values),
            "VecRestoreArrayRead"
        );
    }

    // Added for cardiacFoam: build a distributed PETSc matrix from an Eigen
    // matrix that holds LOCAL rows and GLOBAL columns.
    //
    // That layout is exactly MPIAIJ's, so the only translation needed is to
    // offset the row index by this rank's global row start. Preallocation is
    // split into the diagonal block (columns owned by this rank) and the
    // off-diagonal one; getting that split right is what stops assembly from
    // degenerating into repeated reallocation.
    //
    // In serial gRowStart is 0 and the global column count equals the local row
    // count, so this produces exactly the sequential matrix it replaces.
    void buildDistributedMat(const SpMat& A, Mat& out)
    {
        const PetscInt nLocalRows = static_cast<PetscInt>(A.rows());
        const PetscInt nGlobalCols = static_cast<PetscInt>(A.cols());
        const PetscInt rStart = static_cast<PetscInt>(globalRowStart());
        const PetscInt rEnd = rStart + nLocalRows;

        std::vector<PetscInt> dNnz(nLocalRows, 0);
        std::vector<PetscInt> oNnz(nLocalRows, 0);
        for (PetscInt row = 0; row < nLocalRows; ++row)
        {
            for (SpMat::InnerIterator it(A, row); it; ++it)
            {
                const PetscInt c = static_cast<PetscInt>(it.col());
                if (c >= rStart && c < rEnd)
                {
                    ++dNnz[row];
                }
                else
                {
                    ++oNnz[row];
                }
            }
        }

        checkPetscError(MatCreate(PETSC_COMM_WORLD, &out), "MatCreate");
        checkPetscError
        (
            MatSetSizes
            (
                out,
                nLocalRows,
                nLocalRows,     // local size of the vector this matrix multiplies
                PETSC_DETERMINE,
                nGlobalCols
            ),
            "MatSetSizes"
        );
        checkPetscError(MatSetType(out, MATAIJ), "MatSetType");
        checkPetscError
        (
            MatSeqAIJSetPreallocation
            (
                out, 0, dNnz.empty() ? nullptr : dNnz.data()
            ),
            "MatSeqAIJSetPreallocation"
        );
        checkPetscError
        (
            MatMPIAIJSetPreallocation
            (
                out,
                0, dNnz.empty() ? nullptr : dNnz.data(),
                0, oNnz.empty() ? nullptr : oNnz.data()
            ),
            "MatMPIAIJSetPreallocation"
        );
        checkPetscError
        (
            MatSetOption(out, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE),
            "MatSetOption"
        );

        for (PetscInt row = 0; row < nLocalRows; ++row)
        {
            for (SpMat::InnerIterator it(A, row); it; ++it)
            {
                checkPetscError
                (
                    MatSetValue
                    (
                        out,
                        rStart + static_cast<PetscInt>(it.row()),
                        static_cast<PetscInt>(it.col()),
                        static_cast<PetscScalar>(it.value()),
                        INSERT_VALUES
                    ),
                    "MatSetValue(buildDistributedMat)"
                );
            }
        }

        checkPetscError(MatAssemblyBegin(out, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin");
        checkPetscError(MatAssemblyEnd(out, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd");
    }

    // Added for cardiacFoam: distributed matrix-vector product.
    //
    // With local rows and global columns, Eigen can no longer perform A*x: it
    // would need the whole global vector, and each rank only holds its own
    // block. PETSc's MatMult does the halo exchange internally, so the product
    // becomes correct in parallel and unchanged in serial. The vectors are
    // wrapped around the caller's arrays, so there is no copy.
    class DistributedMatVec
    {
        Mat A_ = nullptr;
        mutable Vec x_ = nullptr;
        mutable Vec y_ = nullptr;
        label nLocal_ = 0;

    public:

        DistributedMatVec() = default;

        ~DistributedMatVec()
        {
            clear();
        }

        DistributedMatVec(const DistributedMatVec&) = delete;
        void operator=(const DistributedMatVec&) = delete;

        bool isInitialised() const
        {
            return A_ != nullptr;
        }

        void clear()
        {
            PetscBool finalized = PETSC_FALSE;
            PetscFinalized(&finalized);
            if (!finalized)
            {
                if (y_) VecDestroy(&y_);
                if (x_) VecDestroy(&x_);
                if (A_) MatDestroy(&A_);
            }
            y_ = nullptr;
            x_ = nullptr;
            A_ = nullptr;
            nLocal_ = 0;
        }

        void reset(const SpMat& A)
        {
            clear();
            nLocal_ = A.rows();
            buildDistributedMat(A, A_);

            checkPetscError
            (
                VecCreateMPI
                (
                    PETSC_COMM_WORLD,
                    static_cast<PetscInt>(nLocal_),
                    PETSC_DETERMINE,
                    &x_
                ),
                "VecCreateMPI(DistributedMatVec x)"
            );
            checkPetscError(VecDuplicate(x_, &y_), "VecDuplicate(DistributedMatVec y)");
        }

        EigVec operator()(const EigVec& x) const
        {
            if (!isInitialised())
            {
                FatalErrorInFunction
                    << "DistributedMatVec used before reset()"
                    << exit(FatalError);
            }
            if (x.size() != nLocal_)
            {
                FatalErrorInFunction
                    << "DistributedMatVec: vector has " << x.size()
                    << " entries, expected " << nLocal_
                    << exit(FatalError);
            }

            PetscScalar* xa = nullptr;
            checkPetscError(VecGetArray(x_, &xa), "VecGetArray(x)");
            for (label i = 0; i < nLocal_; ++i)
            {
                xa[i] = static_cast<PetscScalar>(x[i]);
            }
            checkPetscError(VecRestoreArray(x_, &xa), "VecRestoreArray(x)");

            checkPetscError(MatMult(A_, x_, y_), "MatMult");

            EigVec result(nLocal_);
            const PetscScalar* ya = nullptr;
            checkPetscError(VecGetArrayRead(y_, &ya), "VecGetArrayRead(y)");
            for (label i = 0; i < nLocal_; ++i)
            {
                result[i] = static_cast<scalar>(ya[i]);
            }
            checkPetscError(VecRestoreArrayRead(y_, &ya), "VecRestoreArrayRead(y)");

            return result;
        }
    };

    // RAII guard around PetscInitialize/PetscFinalize.
    //
    // Only finalises if it was this object that initialised PETSc: OpenFOAM's
    // own PETSc-based components, or a library linked into the same process,
    // may have initialised it first, and finalising someone else's session
    // would break them.
    class PetscSession
    {
        bool ownsSession_;

    public:

        PetscSession(int& argc, char**& argv)
        :
            ownsSession_(false)
        {
            PetscBool initialized = PETSC_FALSE;
            checkPetscError(PetscInitialized(&initialized), "PetscInitialized");

            if (!initialized)
            {
                checkPetscError
                (
                    PetscInitialize(&argc, &argv, nullptr, nullptr),
                    "PetscInitialize"
                );
                ownsSession_ = true;
            }
        }

        ~PetscSession()
        {
            if (ownsSession_)
            {
                PetscBool finalized = PETSC_FALSE;
                PetscBool initialized = PETSC_FALSE;

                PetscFinalized(&finalized);
                PetscInitialized(&initialized);

                if (initialized && !finalized)
                {
                    PetscFinalize();
                }
            }
        }

        PetscSession(const PetscSession&) = delete;
        void operator=(const PetscSession&) = delete;
    };

    // Give a KSP its own options prefix, so that several solvers living in the
    // same run (the PDE solve, the JFNK inner Krylov, its preconditioner) can
    // be tuned independently from the command line or from PETSC_OPTIONS.
    void setPetscOptionsPrefix(KSP ksp, const word& prefix)
    {
        if (prefix.size())
        {
            checkPetscError
            (
                KSPSetOptionsPrefix(ksp, prefix.c_str()),
                "KSPSetOptionsPrefix"
            );
        }
    }

    // True when a linearSolverBackend / jfnkLinearSolverBackend key selects
    // PETSc rather than the Eigen path.
    //
    // Modified for cardiacFoam, ported from the MMS solver: the name is
    // VALIDATED. It used to be tested against a two-entry whitelist (petsc,
    // ksp, any case) whose else-branch was Eigen, so a typo moved the whole run
    // onto the serial backend without a word in the log - the same failure mode
    // as "nonlinearMethod Picrad", which collapsed a method matrix onto one
    // method. One spelling each, PETSc and Eigen.
    bool usesPetscBackend(const word& backend)
    {
        requireOneOf("linear solver backend", backend, {"PETSc", "Eigen"});
        return backend == "PETSc";
    }

    // Added for cardiacFoam, ported from the MMS solver: the (KSP, PC) pair a
    // linearSolver request really means for PETSc, resolved in ONE place.
    //
    // A direct solve is preonly + lu whatever petscLinearPcType says, so that
    // is what SparseLU resolves to. On the MMS side the PC used to come from its
    // default (ilu) on two of the three PETSc paths: ONE application of ILU
    // standing in for the inverse of A, and preonly never checks a tolerance, so
    // Picard converged - to the fixed point of the wrong operator. BiCGSTAB
    // needs no case of its own: the default of petscLinearKspType is already
    // its translation (petscKspTypeForLinearSolver), and a petscLinearKspType
    // written explicitly is an explicit choice that should win.
    void resolvePetscKspPc
    (
        const word& linearSolver,
        const word& petscKspType,
        const word& petscPcType,
        word& kspType,
        word& pcType
    )
    {
        kspType = petscKspType;
        pcType = petscPcType;

        if (linearSolver == "SparseLU")
        {
            kspType = "preonly";
            pcType = "lu";
        }
    }

    // Added for cardiacFoam, ported from the MMS solver: the PETSc KSP a
    // backend-neutral linearSolver name stands for, used as the default of
    // petscLinearKspType. SparseLU becomes preonly because in PETSc a direct
    // solve is PCLU applied once (resolvePetscKspPc sets that PC).
    //
    // This replaces petscKspTypeName and petscPcTypeName, which accepted any
    // capitalisation plus a list of aliases and passed anything else through
    // to PETSc unchecked. The KSP and PC keys are now validated at start-up
    // against closed lists of PETSc's own names, so nothing is left to
    // translate.
    word petscKspTypeForLinearSolver(const word& linearSolver)
    {
        if (linearSolver == "GMRES")
        {
            return "gmres";
        }
        if (linearSolver == "BiCGSTAB")
        {
            return "bcgs";
        }
        return "preonly";
    }

    // Added for cardiacFoam: incomplete and complete LU are sequential-only
    // preconditioners; PETSc has no MPIAIJ implementation of either and
    // KSPSetUp fails with PETSC_ERR_SUP. In parallel they have to be wrapped in
    // a block-Jacobi preconditioner, which applies the requested factorisation
    // independently on each rank's diagonal block. That is the standard
    // substitution and is what "ilu in parallel" means in practice.
    //
    // Returns the outer PC name; subPc receives the factorisation to configure
    // on each block, or stays empty when no wrapping is needed.
    std::string petscParallelPcTypeName(const word& pcType, std::string& subPc)
    {
        const std::string s(pcType);
        subPc.clear();

        if (Pstream::parRun() && (s == "ilu" || s == "lu"))
        {
            subPc = s;
            return "bjacobi";
        }

        return s;
    }

    // True when the requested KSP is restarted GMRES, which is the only one
    // that needs a restart length to be configured.
    bool isGmresType(const word& kspType)
    {
        return kspType == "gmres";
    }

    // Payload for a PCShell: the preconditioner application and the local row
    // count. Lets an arbitrary C++ operator (here an ILU of the diffusion
    // matrix) act as a PETSc preconditioner.
    struct PetscShellPCContext
    {
        std::function<EigVec(const EigVec&)> apply;
        label n;
    };

    // Copy this rank's block of an Eigen vector into a PETSc Vec. Assembly and
    // the nonlinear layer work in Eigen; only the linear solve is PETSc's, so
    // each solve crosses this boundary twice.
    void copyEigVecToPetscVec(const EigVec& src, Vec dst)
    {
        PetscScalar* values = nullptr;
        checkPetscError(VecGetArray(dst, &values), "VecGetArray");

        for (PetscInt i = 0; i < src.size(); ++i)
        {
            values[i] = static_cast<PetscScalar>(src[i]);
        }

        checkPetscError(VecRestoreArray(dst, &values), "VecRestoreArray");
    }

    // Persistent PETSc KSP for an ASSEMBLED matrix, kept alive across time
    // steps and nonlinear iterations.
    //
    // The point is the MatrixUpdate policy: rebuilding the matrix and its
    // preconditioner on every solve dominates the cost when, as in Picard with
    // a constant operator, the matrix does not change at all. Holding the Mat,
    // the KSP and the work vectors lets a solve reduce to Reuse (solve only) or
    // ValuesOnly (refresh values and re-factorise, keep the pattern).
    class PetscKspMatrixSolver
    {
        Mat A_;
        KSP ksp_;
        Vec b_;
        Vec x_;
        label n_;

        // This rank's global row offset and the global
        // column count, so updateValues() rebuilds the same distributed matrix
        // that reset() created.
        PetscInt rowStart_ = 0;
        PetscInt nGlobalCols_ = 0;

    public:

        PetscKspMatrixSolver()
        :
            A_(nullptr),
            ksp_(nullptr),
            b_(nullptr),
            x_(nullptr),
            n_(0)
        {}

        ~PetscKspMatrixSolver()
        {
            clear();
        }

        PetscKspMatrixSolver(const PetscKspMatrixSolver&) = delete;
        void operator=(const PetscKspMatrixSolver&) = delete;

        void clear()
        {
            if (ksp_) checkPetscError(KSPDestroy(&ksp_), "KSPDestroy");
            if (A_) checkPetscError(MatDestroy(&A_), "MatDestroy");
            if (b_) checkPetscError(VecDestroy(&b_), "VecDestroy");
            if (x_) checkPetscError(VecDestroy(&x_), "VecDestroy");

            ksp_ = nullptr;
            A_ = nullptr;
            b_ = nullptr;
            x_ = nullptr;
            n_ = 0;
        }

        // Build a distributed MPIAIJ matrix on
        // PETSC_COMM_WORLD instead of a sequential one on PETSC_COMM_SELF.
        //
        // The Eigen matrix handed in has LOCAL rows and GLOBAL columns, which is
        // exactly the MPIAIJ layout, so the conversion is a straight copy once
        // the row indices are offset by this rank's global row start. That
        // offset comes from the file-scope gRowStart rather than from PETSc:
        // OpenFOAM's
        // globalIndex and PETSc's default row distribution happen to agree
        // (contiguous, in rank order), but relying on that coincidence would be
        // a silent trap if either ever changed.
        //
        // In serial rowStart is 0 and A.cols() == A.rows(), so this reduces to
        // the previous sequential behaviour.
        void reset
        (
            const SpMat& A,
            const word& kspType,
            const word& pcType,
            const scalar tolerance,
            const label maxIterations,
            const label restart,
            const word& optionsPrefix,
            const bool useOptions,
            const scalar factorFill = 1.0,
            const scalar dropTolerance = 0.0
        )
        {
            clear();

            n_ = A.rows();
            rowStart_ = static_cast<PetscInt>(globalRowStart());
            nGlobalCols_ = static_cast<PetscInt>(A.cols());
            const PetscInt nLocalRows = static_cast<PetscInt>(A.rows());

            buildDistributedMat(A, A_);

            checkPetscError
            (
                VecCreateMPI(PETSC_COMM_WORLD, nLocalRows, PETSC_DETERMINE, &b_),
                "VecCreateMPI(b)"
            );
            checkPetscError(VecDuplicate(b_, &x_), "VecDuplicate(x)");

            checkPetscError(KSPCreate(PETSC_COMM_WORLD, &ksp_), "KSPCreate");
            checkPetscError(KSPSetOperators(ksp_, A_, A_), "KSPSetOperators");

            const std::string kspName = std::string(kspType);
            // In parallel an ilu/lu request becomes
            // block-Jacobi with that factorisation on each block; see
            // petscParallelPcTypeName. subPcName is empty otherwise, and the
            // options below still key off the factorisation the user asked for.
            std::string subPcName;
            const std::string outerPcName =
                petscParallelPcTypeName(pcType, subPcName);
            const std::string pcName = std::string(pcType);
            checkPetscError(KSPSetType(ksp_, kspName.c_str()), "KSPSetType");

            PC pc = nullptr;
            checkPetscError(KSPGetPC(ksp_, &pc), "KSPGetPC");
            checkPetscError(PCSetType(pc, outerPcName.c_str()), "PCSetType");

            if (!subPcName.empty())
            {
                // The sub-KSPs only exist once the outer PC has been set up.
                checkPetscError(PCSetUp(pc), "PCSetUp(bjacobi)");

                PetscInt nLocalBlocks = 0;
                KSP* subKsp = nullptr;
                checkPetscError
                (
                    PCBJacobiGetSubKSP(pc, &nLocalBlocks, nullptr, &subKsp),
                    "PCBJacobiGetSubKSP"
                );
                for (PetscInt b = 0; b < nLocalBlocks; ++b)
                {
                    PC subPc = nullptr;
                    checkPetscError
                    (
                        KSPSetType(subKsp[b], "preonly"),
                        "KSPSetType(bjacobi sub)"
                    );
                    checkPetscError(KSPGetPC(subKsp[b], &subPc), "KSPGetPC(sub)");
                    checkPetscError
                    (
                        PCSetType(subPc, subPcName.c_str()),
                        "PCSetType(bjacobi sub)"
                    );
                    checkPetscError
                    (
                        PCFactorSetFill(subPc, max(factorFill, scalar(1.0))),
                        "PCFactorSetFill(bjacobi sub)"
                    );

                    // The drop tolerance applies to each block's ILU as it
                    // does to a serial one; it used to be skipped here, so
                    // the same dictionary gave a different factorisation in
                    // parallel.
                    if (subPcName == "ilu" && dropTolerance > 0.0)
                    {
                        checkPetscError
                        (
                            PCFactorSetDropTolerance
                            (
                                subPc,
                                dropTolerance,
                                dropTolerance,
                                1000
                            ),
                            "PCFactorSetDropTolerance(bjacobi sub)"
                        );
                    }
                }
            }
            if (subPcName.empty() && (pcName == "ilu" || pcName == "lu"))
            {
                checkPetscError
                (
                    PCFactorSetFill(pc, max(factorFill, scalar(1.0))),
                    "PCFactorSetFill"
                );
            }
            if (subPcName.empty() && pcName == "ilu" && dropTolerance > 0.0)
            {
                checkPetscError
                (
                    PCFactorSetDropTolerance
                    (
                        pc,
                        dropTolerance,
                        dropTolerance,
                        1000
                    ),
                    "PCFactorSetDropTolerance"
                );
            }

            checkPetscError
            (
                KSPSetTolerances
                (
                    ksp_,
                    tolerance,
                    PETSC_DEFAULT,
                    PETSC_DEFAULT,
                    max(maxIterations, label(1))
                ),
                "KSPSetTolerances"
            );

            if (isGmresType(kspType))
            {
                checkPetscError
                (
                    KSPGMRESSetRestart(ksp_, max(restart, label(1))),
                    "KSPGMRESSetRestart"
                );
            }

            // Stop on the TRUE residual, the one the nonlinear stopping test
            // measures. GMRES and BiCGStab default to left preconditioning,
            // whose tolerance applies to the preconditioned residual: with ilu
            // at a large dt the two differed by ~4x, the linear solve
            // "converged" to 1e-12 and left a true residual of 3.66e-12, and
            // Picard repeated that solve up to the iteration cap. Right
            // preconditioning, as the JFNK shell solver already uses; cg has no
            // right side, so it monitors the unpreconditioned norm instead.
            // Before KSPSetFromOptions, so PETSc options can still override.
            if (isGmresType(kspType) || kspName == "bcgs")
            {
                checkPetscError(KSPSetPCSide(ksp_, PC_RIGHT), "KSPSetPCSide");
            }
            else if (kspName == "cg")
            {
                checkPetscError
                (
                    KSPSetNormType(ksp_, KSP_NORM_UNPRECONDITIONED),
                    "KSPSetNormType"
                );
            }

            setPetscOptionsPrefix(ksp_, optionsPrefix);

            if (useOptions)
            {
                checkPetscError(KSPSetFromOptions(ksp_), "KSPSetFromOptions");
            }

            checkPetscError(KSPSetUp(ksp_), "KSPSetUp");
        }

        bool isInitialised() const
        {
            return ksp_ != nullptr && A_ != nullptr;
        }

        label size() const
        {
            return n_;
        }

        // Refresh the numerical values of the cached PETSc matrix without
        // destroying the KSP/PC/factorisation. Assumes the sparsity pattern
        // of A matches the pattern used in the previous reset(). The flag
        // MAT_NEW_NONZERO_ALLOCATION_ERR was set to PETSC_FALSE in reset(),
        // so any (unexpected) new entries will be tolerated, but the fast
        // path requires the pattern to be invariant.
        // Rows are local and columns global, so the
        // squareness check no longer applies and the row index must carry the
        // same global offset reset() used.
        void updateValues(const SpMat& A, const bool reusePreconditioner = false)
        {
            if (!isInitialised())
            {
                FatalErrorInFunction
                    << "PetscKspMatrixSolver::updateValues requires the "
                    << "solver to have been initialised via reset() first."
                    << exit(FatalError);
            }

            if (A.rows() != n_ || A.cols() != nGlobalCols_)
            {
                FatalErrorInFunction
                    << "PetscKspMatrixSolver::updateValues called with "
                    << "matrix of size " << A.rows() << "x" << A.cols()
                    << " but the cached PETSc Mat is " << n_
                    << "x" << nGlobalCols_
                    << exit(FatalError);
            }

            checkPetscError
            (
                MatZeroEntries(A_),
                "MatZeroEntries(updateValues)"
            );

            const PetscInt nRows = static_cast<PetscInt>(A.rows());
            for (PetscInt row = 0; row < nRows; ++row)
            {
                for (SpMat::InnerIterator it(A, row); it; ++it)
                {
                    const PetscInt r =
                        rowStart_ + static_cast<PetscInt>(it.row());
                    const PetscInt c = static_cast<PetscInt>(it.col());
                    const PetscScalar v =
                        static_cast<PetscScalar>(it.value());

                    checkPetscError
                    (
                        MatSetValue(A_, r, c, v, INSERT_VALUES),
                        "MatSetValue(updateValues)"
                    );
                }
            }

            checkPetscError
            (
                MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY),
                "MatAssemblyBegin(updateValues)"
            );
            checkPetscError
            (
                MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY),
                "MatAssemblyEnd(updateValues)"
            );

            // Force the PC to recompute with the new numerical values while
            // keeping the existing symbolic factorisation cached.
            checkPetscError
            (
                KSPSetReusePreconditioner
                (
                    ksp_,
                    reusePreconditioner ? PETSC_TRUE : PETSC_FALSE
                ),
                "KSPSetReusePreconditioner(updateValues)"
            );
            checkPetscError(KSPSetUp(ksp_), "KSPSetUp(updateValues)");
        }

        EigVec solve
        (
            const EigVec& rhs,
            label& iterations,
            scalar& estimatedError
        )
        {
            if (!ksp_)
            {
                FatalErrorInFunction
                    << "Attempted to use an unset PETSc KSP solver"
                    << exit(FatalError);
            }

            copyEigVecToPetscVec(rhs, b_);
            checkPetscError(VecSet(x_, 0.0), "VecSet(x)");
            checkPetscError(KSPSolve(ksp_, b_, x_), "KSPSolve");

            PetscInt petscIterations = 0;
            PetscReal residualNorm = 0.0;
            KSPConvergedReason reason;

            checkPetscError
            (
                KSPGetIterationNumber(ksp_, &petscIterations),
                "KSPGetIterationNumber"
            );
            checkPetscError
            (
                KSPGetResidualNorm(ksp_, &residualNorm),
                "KSPGetResidualNorm"
            );
            checkPetscError
            (
                KSPGetConvergedReason(ksp_, &reason),
                "KSPGetConvergedReason"
            );

            if (reason < 0)
            {
                FatalErrorInFunction
                    << "PETSc KSP diverged with reason " << reason
                    << exit(FatalError);
            }

            EigVec result(rhs.size());
            copyPetscVecToEigVec(x_, result);

            iterations = static_cast<label>(petscIterations);
            // Relative to the GLOBAL rhs norm. rhs.norm() is this rank's
            // block only, so the reported error used to depend on the rank
            // count (the convergence test itself is PETSc's and was right).
            estimatedError = static_cast<scalar>(residualNorm)
                /max(gNorm(rhs), scalar(SMALL));

            return result;
        }
    };

    // Payload handed to a PETSc MatShell: the matrix-free product to apply and
    // the local row count. Used by JFNK, whose Jacobian is never assembled -
    // it exists only as a directional finite difference of the residual.
    struct PetscShellMatVecContext
    {
        std::function<EigVec(const EigVec&)> matVec;
        label n;
    };

    // MatShell callback: PETSc asks for y = A*x, this unpacks the context and
    // delegates to the C++ functor, converting Vec <-> EigVec on the way.
    PetscErrorCode petscShellMatMult(Mat A, Vec x, Vec y)
    {
        void* rawContext = nullptr;
        PetscErrorCode ierr = MatShellGetContext(A, &rawContext);
        if (ierr) return ierr;

        PetscShellMatVecContext* context =
            static_cast<PetscShellMatVecContext*>(rawContext);

        EigVec xEigen(context->n);
        copyPetscVecToEigVec(x, xEigen);

        const EigVec yEigen = context->matVec(xEigen);
        copyEigVecToPetscVec(yEigen, y);

        return 0;
    }

    // PCShell callback: PETSc asks for y = M^{-1} x, this delegates to the
    // C++ functor.
    PetscErrorCode petscShellPCApply(PC pc, Vec x, Vec y)
    {
        void* rawContext = nullptr;
        PetscErrorCode ierr = PCShellGetContext(pc, &rawContext);
        if (ierr) return ierr;

        PetscShellPCContext* context =
            static_cast<PetscShellPCContext*>(rawContext);

        EigVec xEigen(context->n);
        copyPetscVecToEigVec(x, xEigen);

        const EigVec yEigen = context->apply(xEigen);
        copyEigVecToPetscVec(yEigen, y);

        return 0;
    }

    // Persistent shell-Krylov solver: keeps Mat (shell), Vec b/x
    // and KSP alive across solves. The mat-vec / apply-PC lambdas captured
    // by the JFNK loop change every Newton iteration, but the callback
    // structs live as members of this class so the Mat/KSP keep a stable
    // pointer to them; on each solve we just refresh the std::function
    // payload. Eliminates the per-call alloc/destroy of the Krylov subspace
    // (~25-40 MB) and reduces heap fragmentation across long runs.
    class PetscShellKspSolver
    {
        Mat A_;
        KSP ksp_;
        Vec b_;
        Vec x_;
        PetscShellMatVecContext matContext_;
        PetscShellPCContext pcContext_;
        label n_;

        // This rank's global row offset and the global
        // column count, so updateValues() rebuilds the same distributed matrix
        // that reset() created.
        PetscInt rowStart_ = 0;
        PetscInt nGlobalCols_ = 0;
        bool initialised_;
        bool hasShellPC_;

    public:

        PetscShellKspSolver()
        :
            A_(nullptr),
            ksp_(nullptr),
            b_(nullptr),
            x_(nullptr),
            matContext_{},
            pcContext_{},
            n_(0),
            initialised_(false),
            hasShellPC_(false)
        {
            matContext_.n = 0;
            pcContext_.n = 0;
        }

        ~PetscShellKspSolver()
        {
            clear();
        }

        PetscShellKspSolver(const PetscShellKspSolver&) = delete;
        void operator=(const PetscShellKspSolver&) = delete;

        void clear()
        {
            if (ksp_) checkPetscError(KSPDestroy(&ksp_), "KSPDestroy(cached shell)");
            if (A_) checkPetscError(MatDestroy(&A_), "MatDestroy(cached shell)");
            if (b_) checkPetscError(VecDestroy(&b_), "VecDestroy(cached shell b)");
            if (x_) checkPetscError(VecDestroy(&x_), "VecDestroy(cached shell x)");
            ksp_ = nullptr;
            A_ = nullptr;
            b_ = nullptr;
            x_ = nullptr;
            n_ = 0;
            initialised_ = false;
            hasShellPC_ = false;
            matContext_.matVec = nullptr;
            pcContext_.apply = nullptr;
        }

        bool isInitialised() const { return initialised_; }

        void initialise
        (
            const label n,
            const word& kspType,
            const word& pcType,
            const label restart,
            const label maxIterations,
            const scalar tolerance,
            const word& optionsPrefix,
            const bool useOptions,
            const bool withShellPC
        )
        {
            clear();

            n_ = n;
            const PetscInt nP = static_cast<PetscInt>(n);

            matContext_.n = n;
            pcContext_.n = n;

            // The matrix-free Jacobian becomes a
            // distributed shell operator. Local sizes are this rank's cell
            // count; the global sizes are left to PETSc, which sums them. The
            // shell's MatMult callback is elementwise over the local block, so
            // it needs no halo of its own: the coupling arrives through the
            // residual evaluation, which does its own exchange.
            checkPetscError
            (
                MatCreateShell
                (
                    PETSC_COMM_WORLD,
                    nP,
                    nP,
                    PETSC_DETERMINE,
                    PETSC_DETERMINE,
                    &matContext_,
                    &A_
                ),
                "MatCreateShell(cached)"
            );
            checkPetscError
            (
                MatShellSetOperation
                (
                    A_,
                    MATOP_MULT,
                    reinterpret_cast<void(*)(void)>(petscShellMatMult)
                ),
                "MatShellSetOperation(cached MATOP_MULT)"
            );

            checkPetscError
            (
                VecCreateMPI(PETSC_COMM_WORLD, nP, PETSC_DETERMINE, &b_),
                "VecCreateMPI(cached shell b)"
            );
            checkPetscError
            (
                VecDuplicate(b_, &x_),
                "VecDuplicate(cached shell x)"
            );

            checkPetscError
            (
                KSPCreate(PETSC_COMM_WORLD, &ksp_),
                "KSPCreate(cached shell)"
            );
            checkPetscError
            (
                KSPSetOperators(ksp_, A_, A_),
                "KSPSetOperators(cached shell)"
            );

            const std::string kspName = std::string(kspType);
            checkPetscError
            (
                KSPSetType(ksp_, kspName.c_str()),
                "KSPSetType(cached shell)"
            );
            checkPetscError
            (
                KSPSetTolerances
                (
                    ksp_,
                    tolerance,
                    PETSC_DEFAULT,
                    PETSC_DEFAULT,
                    max(maxIterations, label(1))
                ),
                "KSPSetTolerances(cached shell)"
            );

            if (isGmresType(kspType))
            {
                checkPetscError
                (
                    KSPGMRESSetRestart(ksp_, max(restart, label(1))),
                    "KSPGMRESSetRestart(cached shell)"
                );
            }

            PC pc = nullptr;
            checkPetscError(KSPGetPC(ksp_, &pc), "KSPGetPC(cached shell)");

            setPetscOptionsPrefix(ksp_, optionsPrefix);

            if (useOptions)
            {
                checkPetscError
                (
                    KSPSetFromOptions(ksp_),
                    "KSPSetFromOptions(cached shell)"
                );
                checkPetscError
                (
                    KSPGetPC(ksp_, &pc),
                    "KSPGetPC(cached shell, after options)"
                );
            }

            if (withShellPC)
            {
                checkPetscError
                (
                    PCSetType(pc, PCSHELL),
                    "PCSetType(cached shell PCSHELL)"
                );
                checkPetscError
                (
                    PCShellSetContext(pc, &pcContext_),
                    "PCShellSetContext(cached)"
                );
                checkPetscError
                (
                    PCShellSetApply(pc, petscShellPCApply),
                    "PCShellSetApply(cached)"
                );

                // RIGHT preconditioning: GMRES then minimises the TRUE
                // residual ||F + J delta||, which is what inexact-Newton
                // theory bounds and what the forcing term jfnkLinearTolerance
                // is meant to control, and the finite-difference matvec is
                // applied to the actual Jacobian action. Left preconditioning
                // (PETSc's default for GMRES) would test the preconditioned
                // residual instead, whose size depends on how good the
                // preconditioner happens to be. Standard for JFNK (Knoll &
                // Keyes 2004); both solvers use it.
                checkPetscError
                (
                    KSPSetPCSide(ksp_, PC_RIGHT),
                    "KSPSetPCSide(cached shell)"
                );
                hasShellPC_ = true;
            }
            else
            {
                // ilu/lu have no MPIAIJ implementation; in parallel they are
                // wrapped in block-Jacobi, as for the assembled solves.
                std::string subPcNameShell;
                const std::string pcName =
                    petscParallelPcTypeName(pcType, subPcNameShell);
                checkPetscError
                (
                    PCSetType(pc, pcName.c_str()),
                    "PCSetType(cached shell scalar PC)"
                );
                hasShellPC_ = false;
            }

            initialised_ = true;
        }

        EigVec solve
        (
            const std::function<EigVec(const EigVec&)>& matVec,
            const std::function<EigVec(const EigVec&)>* applyPreconditioner,
            const EigVec& rhs,
            label& iterations,
            scalar& estimatedError
        )
        {
            if (!initialised_)
            {
                FatalErrorInFunction
                    << "PetscShellKspSolver::solve called before initialise()"
                    << exit(FatalError);
            }

            if (rhs.size() != n_)
            {
                FatalErrorInFunction
                    << "PetscShellKspSolver::solve rhs size " << rhs.size()
                    << " differs from cached dimension " << n_
                    << exit(FatalError);
            }

            if (hasShellPC_ != static_cast<bool>(applyPreconditioner))
            {
                FatalErrorInFunction
                    << "PetscShellKspSolver::solve: shell-PC configuration "
                    << "changed between initialise() and solve()."
                    << exit(FatalError);
            }

            // Refresh callback payloads in place; the persistent Mat / KSP
            // hold stable pointers to matContext_ / pcContext_.
            matContext_.matVec = matVec;
            if (applyPreconditioner)
            {
                pcContext_.apply = *applyPreconditioner;
            }

            copyEigVecToPetscVec(rhs, b_);
            checkPetscError(VecSet(x_, 0.0), "VecSet(cached shell x)");
            checkPetscError(KSPSolve(ksp_, b_, x_), "KSPSolve(cached shell)");

            PetscInt petscIterations = 0;
            PetscReal residualNorm = 0.0;
            KSPConvergedReason reason;

            checkPetscError
            (
                KSPGetIterationNumber(ksp_, &petscIterations),
                "KSPGetIterationNumber(cached shell)"
            );
            checkPetscError
            (
                KSPGetResidualNorm(ksp_, &residualNorm),
                "KSPGetResidualNorm(cached shell)"
            );
            checkPetscError
            (
                KSPGetConvergedReason(ksp_, &reason),
                "KSPGetConvergedReason(cached shell)"
            );

            if (reason < 0)
            {
                FatalErrorInFunction
                    << "PETSc cached shell KSP diverged with reason "
                    << reason
                    << exit(FatalError);
            }

            EigVec result(rhs.size());
            copyPetscVecToEigVec(x_, result);

            iterations = static_cast<label>(petscIterations);
            // Relative to the GLOBAL rhs norm (see PetscKspMatrixSolver).
            estimatedError =
                static_cast<scalar>(residualNorm)
               /max(gNorm(rhs), scalar(SMALL));

            return result;
        }
    };
