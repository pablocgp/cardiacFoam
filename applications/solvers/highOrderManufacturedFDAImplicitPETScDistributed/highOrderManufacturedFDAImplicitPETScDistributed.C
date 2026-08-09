/*---------------------------------------------------------------------------*\
Added for cardiacFoam: distributed-linear-algebra variant.

    Copy of highOrderManufacturedFDAImplicitPETScParallel, taken once the
    LRE -> movingLeastSquares port had been verified against the serial
    reference. That sibling is the known-good baseline: it uses the parallel
    high-order library, and its assembly already produces matrices with local
    rows and GLOBAL columns (the MPIAIJ layout), but it still builds a single
    Eigen sparse matrix and hands the finished CSR to PETSc purely as a linear
    solver, on PETSC_COMM_SELF. That is serial by construction.

    Here the linear algebra becomes genuinely distributed: the matrices are
    created as MPIAIJ on PETSC_COMM_WORLD and the matrix-vector products go
    through MatMult instead of Eigen. Kept as a separate solver so the two
    paths can be run against each other rather than compared from memory.
\*---------------------------------------------------------------------------*/

#include <cmath>
#include <petscksp.h>
#include "fvCFD.H"
// Modified for cardiacFoam: LRE (serial) replaced by highOrderInterp, the
// adapter over solids4foam's parallel movingLeastSquares.
#include "highOrderInterp.H"
// Added for cardiacFoam: needed by the coupled-patch (processor face)
// stabilisation exchange in the stiffness assembly.
#include "globalIndex.H"
#include "syncTools.H"
#include "PstreamBuffers.H"
#include "processorFvPatch.H"
#include <chrono>
#include <algorithm>
#include <cctype>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>
#include <functional>
#include <sys/resource.h>
#include <string>
#include <vector>
#ifdef __GLIBC__
#include <malloc.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
    using SpMat = Eigen::SparseMatrix<scalar, Eigen::RowMajor>;
    // Added for cardiacFoam: global row offset of this rank's block of cells,
    // i.e. globalIndex::localStart(). Constant for the whole run (one mesh, one
    // decomposition), so it is set once in main() rather than threaded through
    // every linear-solver signature. Zero in serial.
    label gRowStart = 0;

    using Triplet = Eigen::Triplet<scalar>;
    using EigVec = Eigen::Matrix<scalar, Eigen::Dynamic, 1>;

    // Abort with the calling context if a PETSc call returned an error code.
    // PETSc reports failures by return value, so an unchecked call fails
    // silently and the symptom appears later, in the solution rather than at
    // the call site.
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
        const PetscInt rStart = static_cast<PetscInt>(gRowStart);
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


    // Lower-case copy of a dictionary word, so user input is matched
    // case-insensitively throughout.
    std::string lowerWord(const word& value)
    {
        std::string result(value.c_str());
        std::transform
        (
            result.begin(),
            result.end(),
            result.begin(),
            [](unsigned char c){ return std::tolower(c); }
        );
        return result;
    }

    // True when a linearSolverBackend / jfnkLinearSolverBackend key selects
    // PETSc rather than the Eigen path.
    bool usesPetscBackend(const word& backend)
    {
        const std::string b = lowerWord(backend);
        return b == "petsc" || b == "ksp";
    }

    // Map the backend-neutral solver names used by the dictionaries onto PETSc
    // KSP type names, so one case file drives either backend. Names PETSc
    // already knows pass through. The direct solvers become "preonly" - it is
    // PCLU that does the factorisation.
    std::string petscKspTypeName(const word& kspType)
    {
        const std::string s = lowerWord(kspType);

        if
        (
            s == "petsc"
         || s == "gmres"
         || s == "kspgmres"
         || s == "jfnk"
        )
        {
            return "gmres";
        }
        if (s == "bicgstab" || s == "bcgs")
        {
            return "bcgs";
        }
        if (s == "sparselu" || s == "lu")
        {
            return "preonly";
        }

        return s;
    }

    // Map preconditioner names onto PETSc PC type names, absorbing the
    // spellings the dictionaries have accumulated and the AMG aliases. Anything
    // else passes through, so any PC PETSc supports can be named directly.
    std::string petscPcTypeName(const word& pcType)
    {
        const std::string s = lowerWord(pcType);

        if
        (
            s == "off"
         || s == "false"
         || s == "none"
         || s == "nopc"
        )
        {
            return "none";
        }
        if (s == "ilut")
        {
            return "ilu";
        }
        if (s == "jakobi")
        {
            return "jacobi";
        }
        if (s == "lu" || s == "sparselu")
        {
            return "lu";
        }

        return s;
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
        const std::string s = petscPcTypeName(pcType);
        subPc.clear();

        if (Pstream::parRun() && (s == "ilu" || s == "lu"))
        {
            subPc = s;
            return "bjacobi";
        }

        return s;
    }

    // True when the requested KSP is restarted GMRES, the only one that needs
    // a restart length configured.
    bool isGmresType(const word& kspType)
    {
        return petscKspTypeName(kspType) == "gmres";
    }

    // RAII guard around PetscInitialize/PetscFinalize. Finalises only if it was
    // this object that initialised PETSc: another component may have
    // initialised it first, and finalising someone else's session breaks them.
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

    // Copy this rank's block of an Eigen vector into a PETSc Vec. Assembly and
    // the nonlinear layer work in Eigen and only the linear solve is PETSc's,
    // so every solve crosses this boundary twice.
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

    // Give a KSP its own options prefix, so the several solvers in one run
    // (PDE solve, JFNK inner Krylov, its preconditioner) can be tuned
    // independently from the command line or PETSC_OPTIONS.
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

        // Added for cardiacFoam: this rank's global row offset and the global
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

        // Modified for cardiacFoam: build a distributed MPIAIJ matrix on
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
            rowStart_ = static_cast<PetscInt>(gRowStart);
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

            const std::string kspName = petscKspTypeName(kspType);
            // Modified for cardiacFoam: in parallel an ilu/lu request becomes
            // block-Jacobi with that factorisation on each block; see
            // petscParallelPcTypeName. subPcName is empty otherwise, and the
            // options below still key off the factorisation the user asked for.
            std::string subPcName;
            const std::string outerPcName =
                petscParallelPcTypeName(pcType, subPcName);
            const std::string pcName = petscPcTypeName(pcType);
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
        // Modified for cardiacFoam: rows are local and columns global, so the
        // squareness check no longer applies and the row index must carry the
        // same global offset reset() used.
        void updateValues(const SpMat& A)
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
                KSPSetReusePreconditioner(ksp_, PETSC_FALSE),
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
            estimatedError = static_cast<scalar>(residualNorm)
                /max(rhs.norm(), scalar(SMALL));

            return result;
        }
    };

    // Payload for a PETSc MatShell: the matrix-free product to apply and the
    // local row count. Used by JFNK, whose Jacobian is never assembled - it
    // exists only as a directional finite difference of the residual.
    struct PetscShellMatVecContext
    {
        std::function<EigVec(const EigVec&)> matVec;
        label n;
    };

    // MatShell callback: PETSc asks for y = A*x; unpack the context and
    // delegate to the C++ functor, converting Vec <-> EigVec on the way.
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

    // Payload for a PCShell: the preconditioner application and the local row
    // count, which lets an arbitrary C++ operator act as a PETSc PC.
    struct PetscShellPCContext
    {
        std::function<EigVec(const EigVec&)> apply;
        label n;
    };

    // PCShell callback: PETSc asks for y = M^{-1} x; delegate to the functor.
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

        // Added for cardiacFoam: this rank's global row offset and the global
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

            // Modified for cardiacFoam: the matrix-free Jacobian becomes a
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

            const std::string kspName = petscKspTypeName(kspType);
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
                hasShellPC_ = true;
            }
            else
            {
                const std::string pcName = petscPcTypeName(pcType);
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
            estimatedError =
                static_cast<scalar>(residualNorm)
               /max(rhs.norm(), scalar(SMALL));

            return result;
        }
    };

    // Which manufactured field a generic routine is acting on. u3 is absent
    // because its exact solution is identically zero.
    enum ManufacturedFieldID
    {
        mfVm = 0,
        mfU1 = 1,
        mfU2 = 2
    };

    // Error norms of one field against its exact solution, at the CELL CENTRES:
    // volume-weighted L1 and L2, Linf, and the unweighted cell-error norm that
    // the convergence driver fits its rates against.
    struct FieldErrorSummary
    {
        scalar L1;
        scalar L2;
        scalar Linf;
        scalar errorCell;
    };

    // Error norms of one field at the QUADRATURE POINTS, i.e. of the high-order
    // reconstruction rather than of the cell average. normL2 is the norm of the
    // exact field itself and relL2 the resulting percentage, which is what
    // makes the p-refinement comparable across fields of different magnitude.
    struct ReconstructionErrorSummary
    {
        scalar L1;
        scalar L2;
        scalar Linf;
        scalar normL2;
        scalar relL2;
    };

    // One row of the per-time-step nonlinear history: iteration counts, the
    // coupled and per-field residuals, whether the step converged or was rolled
    // back, and the linear-solver and line-search work it took. Written to
    // nonlinearResiduals.dat when writeNonlinearResiduals is on, and the
    // primary evidence when a run has to be compared against a serial one
    // iteration by iteration.
    struct NonlinearConvergenceRecord
    {
        scalar time;
        label step;
        label iterations;
        label maxIterations;
        bool converged;
        bool rolledBack;
        scalar coupledResidual;
        scalar VmResidual;
        scalar u1Residual;
        scalar u2Residual;
        scalar u3Residual;
        scalar maxStateResidual;
        scalar IionResidual;
        label linearIterations;
        scalar linearError;
        label lineSearchIterations;
    };

    // Spatial factor of the manufactured Vm: a product of sines whose
    // wavenumbers increase with the coordinate index, so that no direction is
    // degenerate and the operator is exercised anisotropically. Chosen smooth
    // and non-polynomial, so a p3 reconstruction cannot be exact and the
    // measured order is the scheme's, not the test function's.
    scalar computeF(const point& p, const label dim)
    {
        const scalar pi = constant::mathematical::pi;

        if (dim == 1)
        {
            return std::cos(pi*p.x());
        }
        else if (dim == 2)
        {
            return std::cos(pi*p.x())*std::cos(2.0*pi*p.y());
        }

        return
            std::cos(pi*p.x())
           *std::cos(2.0*pi*p.y())
           *std::cos(3.0*pi*p.z());
    }

    // Spatial factor of the manufactured states, polynomial and bounded away
    // from zero (so the square root in u2 is safe).
    scalar computeG(const point& p, const label dim)
    {
        if (dim == 1)
        {
            return 1.0 + p.x();
        }
        else if (dim == 2)
        {
            return 1.0 + p.x()*sqr(p.y());
        }

        return 1.0 + p.x()*sqr(p.y())*pow3(p.z());
    }

    // Manufactured solution for Vm: sqrt(1+t) * F(x). The square root gives a
    // time dependence that no theta scheme integrates exactly, so the temporal
    // order is measurable too.
    scalar exactVm(const point& p, const scalar t, const label dim)
    {
        return std::sqrt(1.0 + t)*computeF(p, dim);
    }

    // Manufactured u1, deliberately coupled to Vm so the state ODEs cannot be
    // satisfied independently of the PDE.
    scalar exactU1(const point& p, const scalar t, const label dim)
    {
        return (1.0 + t)*computeG(p, dim) + exactVm(p, t, dim);
    }

    // Manufactured u2. Decays like 1/(1+t) and stays positive.
    scalar exactU2(const point& p, const scalar t, const label dim)
    {
        return 1.0/((1.0 + t)*std::sqrt(computeG(p, dim)));
    }

    // Manufactured u3 is identically zero: the third state is carried so that
    // the code path matches the physiological solver's, not to be verified.
    scalar exactU3(const point&, const scalar, const label)
    {
        return 0.0;
    }

    // Exact value of the field named by a ManufacturedFieldID, so the error and
    // boundary routines can be written once for all fields.
    scalar exactFieldValue
    (
        const ManufacturedFieldID fieldID,
        const point& p,
        const scalar t,
        const label dim
    )
    {
        switch (fieldID)
        {
            case mfVm:
                return exactVm(p, t, dim);
            case mfU1:
                return exactU1(p, t, dim);
            case mfU2:
                return exactU2(p, t, dim);
            default:
                return 0.0;
        }
    }

    // Right-hand side of the manufactured state ODEs, du/dt = R(V, u). These
    // are the rates the SOLVER integrates; the manufactured source that makes
    // the exact solution satisfy them is added separately.
    void reactionRates
    (
        const scalar V,
        const scalar u1,
        const scalar u2,
        const scalar u3,
        scalar& du1dt,
        scalar& du2dt,
        scalar& du3dt
    )
    {
        const scalar a = u1 + u3 - V;
        const scalar b = V - u3;

        du1dt = sqr(a)*sqr(u2) + 0.5*a*sqr(u2)*b;
        du2dt = -a*pow(u2, 3);
        du3dt = 0.0;
    }

    // Manufactured ionic current: a nonlinear surrogate of the physiological
    // Iion, cubic in the states and linear in (V - u3). It is nonlinear in Vm,
    // which is what makes the Picard / diagonalIion / JFNK distinction
    // meaningful on this problem.
    scalar ionicCurrentPDE
    (
        const scalar V,
        const scalar u1,
        const scalar u2,
        const scalar u3,
        const scalar beta,
        const scalar chiVal,
        const scalar CmVal
    )
    {
        return
            -0.5*(u1 + u3 - V)*sqr(u2)*(V - u3)
          + (beta/(chiVal*CmVal))*(V - u3);
    }

    // Source term of the Vm equation, i.e. -Iion. Kept as its own function
    // because the sign convention between "ionic current" and "source" is the
    // easiest thing to get wrong when comparing against the physiological
    // solver.
    scalar vmSourcePDE
    (
        const scalar V,
        const scalar u1,
        const scalar u2,
        const scalar u3,
        const scalar beta,
        const scalar chiVal,
        const scalar CmVal
    )
    {
        return -ionicCurrentPDE(V, u1, u2, u3, beta, chiVal, CmVal);
    }

    // Coefficient that makes the manufactured Vm an exact solution of the
    // diffusion operator: for the sine product of computeF, the Laplacian
    // contributes -pi^2 (D_xx + 4 D_yy + 9 D_zz), the wavenumbers squared. Read
    // from cell 0 because the manufactured setup uses a uniform conductivity.
    scalar computeBeta(const volTensorField& conductivity, const label dim)
    {
        const tensor& D = conductivity[0];
        const scalar pi2 = sqr(constant::mathematical::pi);

        if (dim == 1)
        {
            return -pi2*D.xx();
        }
        else if (dim == 2)
        {
            return -pi2*(D.xx() + 4.0*D.yy());
        }

        return -pi2*(D.xx() + 4.0*D.yy() + 9.0*D.zz());
    }

    // Mesh resolution N recovered from the global cell count, only for
    // labelling output. Reduced, so it does not depend on the decomposition.
    label estimatedN(const fvMesh& mesh)
    {
        const label nCellsGlobal =
            returnReduce(mesh.nCells(), sumOp<label>());

        const scalar dim = max(scalar(mesh.nGeometricD()), 1.0);
        return label(std::round(std::pow(scalar(nCellsGlobal), 1.0/dim)));
    }

    // Representative cell size h = (active domain measure / global cell
    // count)^(1/dim). This is the h the convergence driver fits its orders
    // against, so it must be a property of the mesh and not of the
    // decomposition - hence the reduced bounding box and cell count below.
    scalar characteristicDx(const fvMesh& mesh)
    {
        // Modified for cardiacFoam: mesh.bounds() is the globally reduced
        // bounding box; boundBox(mesh.points()) is this rank's own points only.
        // The cell count below was already reduced, so mixing the two gave a
        // domain measure that shrank with the rank count and an h that did not
        // match the mesh. run_convergence.py fits the convergence order against
        // this h, so the error would have surfaced as a wrong order rather than
        // as an obviously wrong number.
        const boundBox& bb = mesh.bounds();
        const vector span = bb.max() - bb.min();

        scalar activeMeasure = 1.0;
        if (mesh.nGeometricD() >= 1)
        {
            activeMeasure *= max(span.x(), SMALL);
        }
        if (mesh.nGeometricD() >= 2)
        {
            activeMeasure *= max(span.y(), SMALL);
        }
        if (mesh.nGeometricD() >= 3)
        {
            activeMeasure *= max(span.z(), SMALL);
        }

        const scalar nCellsGlobal =
            returnReduce(mesh.nCells(), sumOp<label>());

        return std::pow
        (
            activeMeasure/max(nCellsGlobal, scalar(1.0)),
            1.0/max(scalar(mesh.nGeometricD()), scalar(1.0))
        );
    }

    // Volume-weighted L1 norm, sum(|f| V)/sum(V), globally reduced. Weighting
    // by volume is what makes the norms comparable between meshes of different
    // cell size, and reduction is what makes them independent of the
    // decomposition.
    scalar volumeWeightedL1(const volScalarField& fld)
    {
        const scalarField& f = fld.primitiveField();
        const scalarField& V = fld.mesh().V();

        return gSum(mag(f)*V)/gSum(V);
    }

    // Volume-weighted L2 norm, sqrt(sum(f^2 V)/sum(V)), globally reduced. This
    // is the Vm_L2 the convergence tables fit their orders against.
    scalar volumeWeightedL2(const volScalarField& fld)
    {
        const scalarField& f = fld.primitiveField();
        const scalarField& V = fld.mesh().V();

        return std::sqrt(gSum(sqr(f)*V)/gSum(V));
    }

    // Global maximum of |f|.
    scalar linfNorm(const volScalarField& fld)
    {
        const scalarField magFld(mag(fld.primitiveField()));
        return gMax(magFld);
    }

    // Unweighted L2 over cells, sqrt(sum(f^2)), globally reduced. Reported
    // alongside the weighted norms because it is the one that does not hide a
    // localised error behind a small cell volume.
    scalar cellErrorNorm(const volScalarField& fld)
    {
        return std::sqrt(gSum(sqr(fld.primitiveField())));
    }

    // theta of the time scheme: 1 for backward Euler, 1/2 for Crank-Nicolson.
    // Anything else is a fatal dictionary error rather than a silent default.
    scalar thetaFromScheme(const word& scheme)
    {
        if (scheme == "backwardEuler")
        {
            return 1.0;
        }
        else if (scheme == "crankNicolson")
        {
            return 0.5;
        }

        FatalErrorInFunction
            << "Unknown implicitScheme: " << scheme << nl
            << "Valid options are backwardEuler or crankNicolson"
            << exit(FatalError);

        return 1.0;
    }

    // Collapse the accepted spellings of the massMatrix key onto the two modes
    // the code branches on, lumped and consistent.
    word normalizedMassMatrixType(const word& massMatrix)
    {
        if (massMatrix == "lumped" || massMatrix == "diagonal")
        {
            return "lumped";
        }
        else if (massMatrix == "consistent" || massMatrix == "consistentHO")
        {
            return "consistent";
        }

        FatalErrorInFunction
            << "Unknown massMatrix: " << massMatrix << nl
            << "Valid options are lumped (or diagonal) and consistent (or consistentHO)"
            << exit(FatalError);

        return "lumped";
    }

    FieldErrorSummary computeFieldErrorSummary(const volScalarField& fld)
    {
        FieldErrorSummary s;
        s.L1 = volumeWeightedL1(fld);
        s.L2 = volumeWeightedL2(fld);
        s.Linf = linfNorm(fld);
        s.errorCell = cellErrorNorm(fld);
        return s;
    }

    // d . H . d for a symmetric tensor, written out to exploit the symmetry:
    // the second-order term of a Taylor reconstruction.
    scalar quadraticForm(const symmTensor& H, const vector& d)
    {
        const scalar dx = d.x();
        const scalar dy = d.y();
        const scalar dz = d.z();

        return
            H.xx()*dx*dx
          + 2.0*H.xy()*dx*dy
          + 2.0*H.xz()*dx*dz
          + H.yy()*dy*dy
          + 2.0*H.yz()*dy*dz
          + H.zz()*dz*dz;
    }

    // Evaluate a Taylor reconstruction at offset d from a cell centre:
    //
    //     c + grad.d + (1/2) d.H.d + (1/6) T3(d,d,d)
    //
    // The Hessian and third-derivative pointers are null below the
    // corresponding polynomial order, which is how one routine serves p1, p2
    // and p3. In 2-D the cubic form drops its out-of-plane components.
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

    // Copy the internal field of a volScalarField into an Eigen vector. In
    // parallel this is this rank's block only, which is exactly the row range
    // its matrices own.
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

    // Same conversion for a source field, kept separate so the two uses read
    // differently at the call site.
    EigVec sourceToEigVec(const volScalarField& fld)
    {
        return fieldToEigVec(fld);
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

    // ||current - previous|| / ||current||, globally reduced: the nonlinear
    // stopping test.
    scalar relativeL2Difference
    (
        const scalarField& current,
        const scalarField& previous
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
            den += previous[cellI]*previous[cellI];
        }

        reduce(num, sumOp<scalar>());
        reduce(den, sumOp<scalar>());

        return std::sqrt(num)/(std::sqrt(den) + SMALL);
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

    // Global inner product of two distributed vectors, so a convergence test or
    // an orthogonalisation does not see only this rank's block.
    scalar gDot(const EigVec& a, const EigVec& b)
    {
        scalar s = a.dot(b);
        reduce(s, sumOp<scalar>());
        return s;
    }

    // Modified for cardiacFoam: globally reduced. Its neighbour
    // relativeL2Difference() already reduced; this one did not.
    scalar relativeL2Norm(const EigVec& r, const EigVec& reference)
    {
        return gNorm(r)/(gNorm(reference) + SMALL);
    }

    // Peak resident set size of this process in KB, from getrusage.
    scalar peakResidentSetSizeKB()
    {
        struct rusage usage;
        if (getrusage(RUSAGE_SELF, &usage) == 0)
        {
            return scalar(usage.ru_maxrss);
        }

        return -1.0;
    }

    // Print the peak RSS at a named point of the run. The high-order setup
    // (stencils, per-quadrature-point QR factorisations) is the memory peak of
    // these cases, so the checkpoints bracket it.
    void logMemoryCheckpoint(const char* stage)
    {
        const scalar peakMemoryKB = peakResidentSetSizeKB();
        Info<< "Memory checkpoint [" << stage << "]: peakRSS_MB = "
            << peakMemoryKB/1024.0 << endl;
    }

    // Modified for cardiacFoam: column globalised, see assembleDiagonalMassMatrix.
    void addDiagonalToMatrix
    (
        const scalarField& diagonal,
        const globalIndex& gc,
        const scalar coefficient,
        SpMat& A
    )
    {
        forAll(diagonal, cellI)
        {
            A.coeffRef(cellI, gc.toGlobal(cellI)) += coefficient*diagonal[cellI];
        }

        A.makeCompressed();
    }

    template<class MatrixVectorProduct>
    EigVec solveGMRES
    (
        const MatrixVectorProduct& matVec,
        const EigVec& b,
        const label krylovDim,
        const label maxRestarts,
        const scalar tolerance,
        label& iterations,
        scalar& estimatedError
    )
    {
        const label n = b.size();
        const label m = min(max(krylovDim, label(1)), label(n));
        const label maxOuter = max(maxRestarts, label(0));

        EigVec x = EigVec::Zero(n);
        EigVec r = b;
        const scalar bNorm = b.norm();
        scalar beta = bNorm;
        iterations = 0;
        estimatedError = 1.0;

        if (bNorm <= SMALL)
        {
            estimatedError = 0.0;
            return x;
        }

        std::vector<EigVec> V(m + 1, EigVec::Zero(n));
        Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> H =
            Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero
            (m + 1, m);

        for (label restart = 0; restart <= maxOuter; ++restart)
        {
            V[0] = r/beta;
            H.setZero();
            bool toleranceMet = false;
            EigVec bestX = x;

            for (label j = 0; j < m; ++j)
            {
                ++iterations;

                EigVec w = matVec(V[j]);

                for (label i = 0; i <= j; ++i)
                {
                    H(i, j) = V[i].dot(w);
                    w -= H(i, j)*V[i];
                }

                H(j + 1, j) = w.norm();
                if (H(j + 1, j) > SMALL && j + 1 < m + 1)
                {
                    V[j + 1] = w/H(j + 1, j);
                }

                Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> Hj =
                    H.block(0, 0, j + 2, j + 1);
                EigVec g = EigVec::Zero(j + 2);
                g[0] = beta;

                const EigVec y = Hj.colPivHouseholderQr().solve(g);

                EigVec xj = x;
                for (label i = 0; i <= j; ++i)
                {
                    xj += y[i]*V[i];
                }

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

            x = bestX;
            if (toleranceMet) break;
            if (restart == maxOuter) break;

            r = b - matVec(x);
            beta = r.norm();
            estimatedError = beta/max(bNorm, SMALL);
            if (beta <= SMALL*bNorm) break;
        }

        return x;
    }

    template<class MatrixVectorProduct, class Preconditioner>
    EigVec solveLeftPreconditionedGMRES
    (
        const MatrixVectorProduct& matVec,
        const Preconditioner& applyPreconditioner,
        const EigVec& b,
        const label krylovDim,
        const label maxRestarts,
        const scalar tolerance,
        label& iterations,
        scalar& estimatedError
    )
    {
        const label n = b.size();
        const label m = min(max(krylovDim, label(1)), label(n));
        const label maxOuter = max(maxRestarts, label(0));

        EigVec x = EigVec::Zero(n);
        EigVec bPrec = applyPreconditioner(b);
        EigVec r = bPrec;
        const scalar bNorm = bPrec.norm();
        scalar beta = bNorm;
        iterations = 0;
        estimatedError = 1.0;

        if (bNorm <= SMALL)
        {
            estimatedError = 0.0;
            return x;
        }

        std::vector<EigVec> V(m + 1, EigVec::Zero(n));
        Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> H =
            Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero
            (m + 1, m);

        for (label restart = 0; restart <= maxOuter; ++restart)
        {
            V[0] = r/beta;
            H.setZero();
            bool toleranceMet = false;
            EigVec bestX = x;

            for (label j = 0; j < m; ++j)
            {
                ++iterations;

                EigVec w = applyPreconditioner(matVec(V[j]));

                for (label i = 0; i <= j; ++i)
                {
                    H(i, j) = V[i].dot(w);
                    w -= H(i, j)*V[i];
                }

                H(j + 1, j) = w.norm();
                if (H(j + 1, j) > SMALL && j + 1 < m + 1)
                {
                    V[j + 1] = w/H(j + 1, j);
                }

                Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> Hj =
                    H.block(0, 0, j + 2, j + 1);
                EigVec g = EigVec::Zero(j + 2);
                g[0] = beta;

                const EigVec y = Hj.colPivHouseholderQr().solve(g);

                EigVec xj = x;
                for (label i = 0; i <= j; ++i)
                {
                    xj += y[i]*V[i];
                }

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

            x = bestX;
            if (toleranceMet) break;
            if (restart == maxOuter) break;

            r = applyPreconditioner(b - matVec(x));
            beta = r.norm();
            estimatedError = beta/max(bNorm, SMALL);
            if (beta <= SMALL*bNorm) break;
        }

        return x;
    }

    // Append a matrix entry, dropping values below SMALL. The filter is on the
    // INDIVIDUAL contribution, not the accumulated one, so it is independent of
    // the order contributions arrive in and of any chunking of the assembly.
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

    // Modified for cardiacFoam: matrices are now nLocalCells x nGlobalCells
    // (local rows, global columns), which is the MPIAIJ layout PETSc expects.
    // In serial nGlobal == nLocal and toGlobal() is the identity, so nothing
    // changes there.
    void assembleDiagonalMassMatrix
    (
        const fvMesh& mesh,
        const globalIndex& gc,
        const scalar coefficient,
        SpMat& M
    )
    {
        std::vector<Triplet> triplets;
        triplets.reserve(mesh.nCells());

        forAll(mesh.C(), cellI)
        {
            triplets.emplace_back(cellI, gc.toGlobal(cellI), coefficient);
        }

        M.resize(mesh.nCells(), gc.totalSize());
        M.setFromTriplets(triplets.begin(), triplets.end());
        M.makeCompressed();
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
        // Modified for cardiacFoam: local rows, global columns (MPIAIJ layout).
        const globalIndex& gc = LREInterp.globalCells();
        const label nGlobalCells = gc.totalSize();
        const label cellChunk = 50000;

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
                // Modified for cardiacFoam: CompactListList::operator[] returns a
        // UList VIEW by value, not a List. Binding it to a const labelList&
        // compiles but yields an EMPTY list, silently dropping the whole
        // stencil. Bind to the returned view type instead.
        const UList<label> curStencil = stencils[cellI];
                const label selfCoeffI = curStencil.size();
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
                            gc.toGlobal(cellI),
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
                        gc.toGlobal(cellI),
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
    }

    // Two-point orthogonal diffusion coefficient of a face,
    // |Sf| (n . D . e)/|d| with e the unit owner-to-neighbour direction: the
    // standard low-order finite-volume flux, which neglects the non-orthogonal
    // part of d relative to Sf.
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
    // a degenerate tensor cannot give a zero stabilisation scale.
    scalar normalDiffusivity(const tensor& D, const vector& n)
    {
        return max(mag(n & (D & n)), SMALL);
    }

    // Scale of the face-jump stabilisation, alpha |Sf| (n.D.n)/|d.n|. Zero when
    // alpha is zero, which switches the stabilisation off entirely. It depends
    // on the face orientation only through |d.n|, so the two sides of a face -
    // including the two ranks of a processor face - compute the same value.
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

        // Modified for cardiacFoam: matrix columns are global cell indices.
        const globalIndex& gc = LREInterp.globalCells();

        const CompactListList<label>& stencils = LREInterp.globalCellStencils();
        const CompactListList<vector>& gradCoeffs = LREInterp.QRGradCoeffs();
        // Modified for cardiacFoam: CompactListList::operator[] returns a
        // UList VIEW by value, not a List. Binding it to a const labelList&
        // compiles but yields an EMPTY list, silently dropping the whole
        // stencil. Bind to the returned view type instead.
        const UList<label> curStencil = stencils[cellI];
        const label selfCoeffI = curStencil.size();

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
            gc.toGlobal(cellI),
            scale*(gradCoeffs[cellI][selfCoeffI] & d)
        );
    }

    // Added for cardiacFoam: Taylor extrapolation coefficient of ONE entry of a
    // cell's reconstruction row, i.e. the weight with which the unknown in
    // column entryI enters the value reconstructed at cellI + d.
    //
    // Factored out so that the triplet writer below and the row builder used to
    // ship a cell's row across a processor boundary cannot drift apart: the
    // remote contribution has to be the SAME polynomial the local one uses.
    // The self entry (index == stencil size) additionally carries the 1.0 of
    // the zeroth-order term; that is added by the callers, not here.
    scalar cellTaylorEntryCoeff
    (
        const label cellI,
        const label entryI,
        const vector& d,
        const highOrderInterp& LREInterp,
        const bool twoD
    )
    {
        scalar coeff = LREInterp.QRGradCoeffs()[cellI][entryI] & d;

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

    // Added for cardiacFoam: the same row as addCellTaylorExtrapolationCoeffs
    // writes, materialised as (global column, coefficient) pairs instead of
    // being pushed straight into the triplet list.
    //
    // This is what crosses a processor boundary. The row is built by the rank
    // that OWNS the cell - the only rank that has its stencil and its
    // reconstruction coefficients - and evaluated at the shared face centre, so
    // the receiving rank needs neither the remote cell centre nor its stencil.
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
        const globalIndex& gc = LREInterp.globalCells();
        const vector d = evalPoint - C[cellI];

        // Modified for cardiacFoam: CompactListList::operator[] returns a
        // UList VIEW by value, not a List (see addCellGradientDotCoeffs).
        const UList<label> curStencil = LREInterp.globalCellStencils()[cellI];
        const label selfCoeffI = curStencil.size();

        cols.setSize(selfCoeffI + 1);
        coeffs.setSize(selfCoeffI + 1);

        forAll(curStencil, cI)
        {
            cols[cI] = curStencil[cI];
            coeffs[cI] = cellTaylorEntryCoeff(cellI, cI, d, LREInterp, twoD);
        }

        cols[selfCoeffI] = gc.toGlobal(cellI);
        coeffs[selfCoeffI] =
            1.0 + cellTaylorEntryCoeff(cellI, selfCoeffI, d, LREInterp, twoD);
    }

    // Added for cardiacFoam: add a pre-built remote row (global columns) to one
    // matrix row. Used only for rows that arrived from another rank.
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

    // Added for cardiacFoam: (global column, coefficient) form of the row that
    // addCellGradientDotCoeffs writes, i.e. of grad Vm|_cellI . d.
    //
    // Same purpose as buildCellTaylorExtrapolationRow: this is the row that the
    // rank owning the cell sends across a processor face, because its stencil
    // does not exist on the other side.
    void buildCellGradientDotRow
    (
        labelList& cols,
        scalarField& coeffs,
        const label cellI,
        const vector& d,
        const highOrderInterp& LREInterp
    )
    {
        const globalIndex& gc = LREInterp.globalCells();
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

        cols[selfCoeffI] = gc.toGlobal(cellI);
        coeffs[selfCoeffI] = gradCoeffs[cellI][selfCoeffI] & d;
    }

    // Added for cardiacFoam: exchange one VARIABLE-LENGTH matrix row per
    // coupled face with the rank on the other side.
    //
    // The stabilisation jump on a processor face needs the far cell's
    // reconstruction row - up to ~60 (global column, coefficient) pairs for a
    // p3 stencil in 3-D - which no field swap can provide. PstreamBuffers is
    // the standard route for variable-length data.
    //
    // Each rank builds the row for its OWN cell adjacent to the face and sends
    // it; the geometric quantity the row is evaluated at (face centre, or the
    // owner-to-neighbour vector) is shared, so no cell centres are exchanged
    // and mesh.C()'s coupled values - a slicedVolVectorField, only valid after
    // an evaluation nothing here performs - are never touched.
    //
    // Results are indexed by boundary-face position, faceI - nInternalFaces();
    // entries for non-coupled faces stay empty.
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

    // Expand the Taylor reconstruction of cell cellI evaluated at evalPoint
    // into matrix coefficients, to the polynomial order of the interpolator
    // (gradient, plus Hessian at p2+, plus third derivatives at p3). Same row
    // layout as addCellGradientDotCoeffs; the self entry additionally carries
    // the 1.0 of the zeroth-order term.
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

        // Modified for cardiacFoam: matrix columns are global cell indices.
        const globalIndex& gc = LREInterp.globalCells();

        const vector d = evalPoint - C[cellI];
        const CompactListList<label>& stencils = LREInterp.globalCellStencils();

        // Modified for cardiacFoam: CompactListList::operator[] returns a
        // UList VIEW by value, not a List. Binding it to a const labelList&
        // compiles but yields an EMPTY list, silently dropping the whole
        // stencil. Bind to the returned view type instead.
        const UList<label> curStencil = stencils[cellI];
        const label selfCoeffI = curStencil.size();

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

        const scalar selfCoeff =
            1.0 + cellTaylorEntryCoeff(cellI, selfCoeffI, d, LREInterp, twoD);

        addTripletIfNeeded(triplets, row, gc.toGlobal(cellI), scale*selfCoeff);
    }

    // Stabilisation jump on an internal face, added to one row:
    //
    //     scale*(Vm[nei] - Vm[own])
    //       - 0.5*scale*(grad Vm|_own . dPN)
    //       - 0.5*scale*(grad Vm|_nei . dPN)
    //
    // i.e. the difference between the two cell values and what the two
    // reconstructions predict. It vanishes for a field the reconstruction
    // represents exactly, so it does not change the formal order. The caller
    // passes +a/V[own] for the owner row and -a/V[nei] for the neighbour's.
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
        // Modified for cardiacFoam: these two columns were LOCAL cell indices
        // while every other column written by this assembly (including
        // addCellGradientDotCoeffs just below) is global. Identical in serial,
        // where gRowStart is 0.
        const globalIndex& gc = LREInterp.globalCells();

        addTripletIfNeeded(triplets, row, gc.toGlobal(nei), scale);
        addTripletIfNeeded(triplets, row, gc.toGlobal(own), -scale);

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

    // Added for cardiacFoam: the internal jump above, written for a face whose
    // neighbouring cell lives on another rank.
    //
    // Only the local cell's row is written; the far rank writes its own, which
    // is why the caller passes the jump seen from the local cell (dPN pointing
    // outwards, scale positive) rather than an owner/neighbour convention that
    // would have to be agreed across the partition.
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
        const globalIndex& gc = LREInterp.globalCells();

        addTripletIfNeeded(triplets, row, nbrGlobalCell, scale);
        addTripletIfNeeded(triplets, row, gc.toGlobal(own), -scale);

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
        // expressed in THIS rank's orientation, so it carries the same -0.5
        // factor as the local one.
        addRemoteRowCoeffs
        (
            triplets,
            row,
            -0.5*scale,
            nbrGradCols,
            nbrGradCoeffs
        );
    }

    // The same jump against a Dirichlet boundary face, where the prescribed
    // value goes to the right-hand side instead of the matrix, so only the
    // owner's terms are added here.
    void addRhieChowBoundaryJumpCoeffs
    (
        std::vector<Triplet>& triplets,
        const label row,
        const scalar scale,
        const label own,
        const vector& dPb,
        const highOrderInterp& LREInterp
    )
    {
        // Modified for cardiacFoam: global column (see the internal jump).
        addTripletIfNeeded(triplets, row, LREInterp.globalCells().toGlobal(own), -scale);
        addCellGradientDotCoeffs
        (
            triplets,
            row,
            -scale,
            own,
            dPb,
            LREInterp
        );
    }


    // Low-order diffusion operator (useHighOrder_Vm = false): two-point
    // orthogonal fluxes over internal, Dirichlet and coupled (processor) faces,
    // divided by the cell volume so K represents div(D grad .) itself rather
    // than its integral. Optionally stabilised with the Rhie-Chow face jump.
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

        // Modified for cardiacFoam: local rows, GLOBAL columns, as in
        // assembleHighOrderStiffnessMatrix. This assembly still used local
        // columns for its orthogonal terms - and, through
        // addCellGradientDotCoeffs, global ones for its stabilisation terms -
        // which agree only in serial.
        const globalIndex& gc = LREInterp.globalCells();

        std::vector<Triplet> triplets;
        triplets.reserve(4*mesh.nFaces());

        forAll(neighbour, faceI)
        {
            const label own = owner[faceI];
            const label nei = neighbour[faceI];

            const tensor Df = 0.5*(conductivity[own] + conductivity[nei]);
            const scalar a = orthogonalDiffusionCoeff(mesh.Sf()[faceI], C[nei] - C[own], Df);

            addTripletIfNeeded(triplets, own, gc.toGlobal(own), -a/max(V[own], SMALL));
            addTripletIfNeeded(triplets, own, gc.toGlobal(nei),  a/max(V[own], SMALL));
            addTripletIfNeeded(triplets, nei, gc.toGlobal(own),  a/max(V[nei], SMALL));
            addTripletIfNeeded(triplets, nei, gc.toGlobal(nei), -a/max(V[nei], SMALL));
        }

        const surfaceVectorField& Cf = mesh.Cf();

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

        // Added for cardiacFoam: everything the coupled (processor) faces need.
        // This assembly had no coupled branch at all, so in parallel it produced
        // partition blocks with no coupling between them - neither the
        // orthogonal flux nor the stabilisation jump. See the equivalent block
        // in assembleHighOrderStiffnessMatrix for why the far cell's data has
        // to be swapped rather than read off a patch field.
        List<tensor> nbrConductivity;
        List<label> nbrGlobalCell;
        List<labelList> nbrGradCols;
        List<scalarField> nbrGradCoeffs;

        if (Pstream::parRun())
        {
            syncTools::swapBoundaryCellList
            (
                mesh, conductivity.primitiveField(), nbrConductivity
            );

            labelList myGlobalCell(mesh.nCells());
            forAll(myGlobalCell, cellI)
            {
                myGlobalCell[cellI] = gc.toGlobal(cellI);
            }
            syncTools::swapBoundaryCellList(mesh, myGlobalCell, nbrGlobalCell);

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
                        // receiver's own owner-to-neighbour vector, so that the
                        // row arrives as grad|_here . dPN_there and can be used
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

            if (bcType == "empty" || bcType == zeroGradientFvPatchScalarField::typeName)
            {
                continue;
            }

            if (patch.coupled())
            {
                // Processor patches only, as in the high-order assembly: the
                // row exchange has no route to the far side of a cyclic.
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

                    const vector Sf = mesh.Sf().boundaryField()[patchI][faceI];
                    const vector dPN = patchDelta[faceI];
                    const tensor Df =
                        0.5*(conductivity[own] + nbrConductivity[bFaceI]);

                    const scalar a = orthogonalDiffusionCoeff(Sf, dPN, Df);

                    addTripletIfNeeded
                    (
                        triplets, own, gc.toGlobal(own), -a/max(V[own], SMALL)
                    );
                    addTripletIfNeeded
                    (
                        triplets, own, nbrGlobalCell[bFaceI], a/max(V[own], SMALL)
                    );

                    if (stabilisationAlpha > SMALL)
                    {
                        const scalar area = mag(Sf) + VSMALL;
                        const vector n = Sf/area;
                        const scalar aStab =
                            stabilisationFaceCoeff
                            (
                                stabilisationAlpha,
                                area,
                                Df,
                                n,
                                dPN
                            );

                        addRhieChowCoupledJumpCoeffs
                        (
                            triplets,
                            own,
                            aStab/max(V[own], SMALL),
                            own,
                            nbrGlobalCell[bFaceI],
                            dPN,
                            nbrGradCols[bFaceI],
                            nbrGradCoeffs[bFaceI],
                            LREInterp
                        );
                    }
                }

                continue;
            }

            if (bcType == fixedValueFvPatchScalarField::typeName || bcType == "fixedVoltage")
            {
                forAll(patch, faceI)
                {
                    const label gf = patch.start() + faceI;
                    const label own = owner[gf];

                    const scalar a =
                        orthogonalDiffusionCoeff
                        (
                            mesh.Sf().boundaryField()[patchI][faceI],
                            Cf.boundaryField()[patchI][faceI] - C[own],
                            conductivity[own]
                        );

                    addTripletIfNeeded
                    (
                        triplets, own, gc.toGlobal(own), -a/max(V[own], SMALL)
                    );

                    if (stabilisationAlpha > SMALL)
                    {
                        const vector Sf =
                            mesh.Sf().boundaryField()[patchI][faceI];
                        const scalar area = mag(Sf) + VSMALL;
                        const vector n = Sf/area;
                        const vector dPb =
                            Cf.boundaryField()[patchI][faceI] - C[own];
                        const scalar aStab =
                            stabilisationFaceCoeff
                            (
                                stabilisationAlpha,
                                area,
                                conductivity[own],
                                n,
                                dPb
                            );

                        addRhieChowBoundaryJumpCoeffs
                        (
                            triplets,
                            own,
                            aStab/max(V[own], SMALL),
                            own,
                            dPb,
                            LREInterp
                        );
                    }
                }
            }
        }

        // Modified for cardiacFoam: global column count (MPIAIJ layout).
        K.resize(mesh.nCells(), gc.totalSize());
        K.setFromTriplets(triplets.begin(), triplets.end());
        K.makeCompressed();
        std::vector<Triplet>().swap(triplets);
    }

    // High-order diffusion operator: the face flux n . (D grad Vm) integrated
    // by LRE quadrature over each face, assembled into a matrix with local rows
    // and global columns.
    //
    // Assembled in chunks of faces: on a large 3-D unstructured tet mesh the
    // monolithic triplet vector would peak in the tens of GB purely through
    // capacity doubling. Each chunk is folded into K by sparse addition, which
    // is numerically transparent because addTripletIfNeeded filters individual
    // contributions rather than accumulated ones.
    //
    // The boundary loop also handles coupled faces, both for the flux (whose
    // stencil already crosses the cut) and for the stabilisation jump (whose
    // far-side reconstruction row has to be exchanged).
    void assembleHighOrderStiffnessMatrix
    (
        const fvMesh& mesh,
        const volTensorField& conductivity,
        const highOrderInterp& LREInterp,
        const scalar stabilisationAlpha,
        const label tripletsPerFaceReserve,
        SpMat& K
    )
    {
        // Chunked assembly: avoids accumulating ~434M triplets in a single
        // std::vector (which would peak at ~19 GB transient due to repeated
        // capacity-doubling reallocations on large 3D unstructured tet
        // meshes). The chunked variant bounds the triplet buffer to a fixed
        // size per chunk and folds each chunk into the global K via sparse
        // matrix addition. Numerically equivalent to the monolithic version:
        // K = sum_i K_chunk_i, and addTripletIfNeeded filters by individual
        // value magnitude (not by accumulated value), so the chunking is
        // transparent to the final result.
        const bool twoD = mesh.nGeometricD() == 2;
        const vectorField& C = mesh.C();
        const surfaceVectorField& Cf = mesh.Cf();
        const CompactListList<point>& faceQP = LREInterp.faceQuadPoints();
        const CompactListList<scalar>& faceQW = LREInterp.faceQuadWeightPhysical();
        const CompactListList<label>& faceStencils = LREInterp.globalFaceStencils();
        const List<CompactListList<vector>>& faceGradCoeffs =
            LREInterp.QRGradFaceGPCoeffs();

        const scalarField& V = mesh.V();
        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();

        const label nCells = mesh.nCells();
        // Modified for cardiacFoam: local rows, global columns (MPIAIJ layout).
        const globalIndex& gc = LREInterp.globalCells();
        const label nGlobalCells = gc.totalSize();
        const label nInternalFaces = neighbour.size();
        const label faceChunk = 50000;


        std::vector<Triplet> triplets;
        // Reserve for a full chunk without forcing the tet/stabilised worst
        // case on flux-only hex runs. The cap is selected by the caller from
        // memoryOptimization controls; it changes capacity only, not assembly.
        triplets.reserve(faceChunk * tripletsPerFaceReserve);

        SpMat Klocal(nCells, nGlobalCells);
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
                Klocal.resize(nCells, nGlobalCells);
                kInitialised = true;
            }
            else
            {
                K += Klocal;
            }
            triplets.clear();
        };

        // --- Internal faces, processed in chunks ---
        for
        (
            label faceStart = 0;
            faceStart < nInternalFaces;
            faceStart += faceChunk
        )
        {
            const label faceEnd =
                min(faceStart + faceChunk, nInternalFaces);

            // Flux contribution
            for (label faceI = faceStart; faceI < faceEnd; ++faceI)
            {
                const label own = owner[faceI];
                const label nei = neighbour[faceI];

                const vector Sf = mesh.Sf()[faceI];
                const scalar area = mag(Sf) + VSMALL;
                const vector n = Sf/area;

                const UList<label> curStencil = faceStencils[faceI];


                forAll(faceQP[faceI], qpI)
                {
                    const scalar w = faceQW[faceI][qpI];

                    forAll(curStencil, cI)
                    {
                        const label col = curStencil[cI];
                        const vector gCoeff = faceGradCoeffs[faceI][qpI][cI];

                        // Modified for cardiacFoam: the explicit area factor is
                        // gone. LRE normalised its face quadrature weights to
                        // sum to 1, so the caller supplied |Sf|; fvMeshQuadrature
                        // returns physical weights that already sum to |Sf|.
                        const scalar fluxCoeff =
                            w*(n & (conductivity[own] & gCoeff));


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

            // Stabilisation contribution (same chunk range)
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
                    const tensor Df =
                        0.5*(conductivity[own] + conductivity[nei]);
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
        // nothing else. Invisible in serial and invisible in the operator
        // fingerprint (an omitted term is not a wrong one), it showed up only
        // as an np=1 vs np=2 deviation in Vm that did NOT shrink with the
        // nonlinear method, i.e. as a spatial-discretisation defect.
        //
        // Two pieces cross the boundary:
        //   - the far cell's conductivity, so that both ranks form the same
        //     face coefficient a (syncTools gives the CELL value; a processor
        //     patch field would give the interpolated FACE value);
        //   - the far cell's Taylor reconstruction row, evaluated at the shared
        //     face centre by the rank that owns it.
        List<tensor> nbrConductivity;
        List<labelList> nbrTaylorCols;
        List<scalarField> nbrTaylorCoeffs;

        if (stabilisationAlpha > SMALL && Pstream::parRun())
        {
            syncTools::swapBoundaryCellList
            (
                mesh, conductivity.primitiveField(), nbrConductivity
            );

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

        // --- Boundary faces in a single pass (peak is small: ~50K faces
        // × ~300 triplets ≈ 250 MB, manageable without chunking) ---
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

            // Processor patches only: exchangeCoupledFaceRows knows how to
            // reach the far rank of a processorFvPatch and nothing else. A
            // cyclic would be coupled too, and taking this branch with no row
            // to add would be worse than the omission this replaces.
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

            // Owner-to-neighbour vector across the coupled face. patch.delta()
            // is used rather than C[nei] - C[own] because mesh.C() is a
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

                const UList<label> curStencil = faceStencils[globalFaceI];

                forAll(faceQP[globalFaceI], qpI)
                {
                    const scalar w = faceQW[globalFaceI][qpI];

                    forAll(curStencil, cI)
                    {
                        const label col = curStencil[cI];
                        const vector gCoeff =
                            faceGradCoeffs[globalFaceI][qpI][cI];

                        // Modified for cardiacFoam: physical quadrature weights
                        // already carry |Sf| (see the internal-face loop above).
                        const scalar fluxCoeff =
                            w*(n & (conductivity[own] & gCoeff));

                        addTripletIfNeeded
                        (
                            triplets,
                            own,
                            col,
                            fluxCoeff/max(V[own], SMALL)
                        );
                    }
                }

                if (stabiliseCoupled)
                {
                    // Same jump as the internal-face loop, written from the
                    // point of view of the LOCAL cell: it receives
                    //
                    //     a/V[own] * ( T_far(xf) - T_own(xf) )
                    //
                    // where T_c is the Taylor reconstruction of cell c. The
                    // internal-face loop writes this with a + sign into the
                    // owner row and a - sign into the neighbour row; those two
                    // signs are the same expression seen from either cell, so
                    // no owner/neighbour convention has to be agreed across the
                    // partition. The rank on the far side adds its own copy,
                    // exactly as the internal loop adds both rows at once.
                    const label bFaceI = bStart + faceI;
                    const point xf = Cf.boundaryField()[patchI][faceI];
                    const vector dPN = patchDelta[faceI];
                    const tensor Df =
                        0.5*(conductivity[own] + nbrConductivity[bFaceI]);

                    // Identical on both ranks: the area and the normal
                    // diffusivity are orientation-independent, and dPN enters
                    // only as |dPN . n|.
                    const scalar a =
                        stabilisationFaceCoeff
                        (
                            stabilisationAlpha,
                            area,
                            Df,
                            n,
                            dPN
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
                else if
                (
                    stabilisationAlpha > SMALL
                 && (
                        bcType == fixedValueFvPatchScalarField::typeName
                     || bcType == "fixedVoltage"
                    )
                )
                {
                    const vector Sf = mesh.Sf().boundaryField()[patchI][faceI];
                    const scalar area = mag(Sf) + VSMALL;
                    const vector n = Sf/area;
                    const point xf = Cf.boundaryField()[patchI][faceI];
                    const vector dPb = xf - C[own];
                    const scalar a =
                        stabilisationFaceCoeff
                        (
                            stabilisationAlpha,
                            area,
                            conductivity[own],
                            n,
                            dPb
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

        // Ensure K has the correct shape even if there were no triplets
        // anywhere (degenerate case).
        if (!kInitialised)
        {
            K.resize(nCells, nGlobalCells);
            K.setZero();
        }

        K.makeCompressed();
        std::vector<Triplet>().swap(triplets);
    }

    // Right-hand side contribution of the Dirichlet boundaries for the
    // low-order operator: the prescribed exact value times the face
    // coefficient, which is the part of the flux that cannot go into the
    // matrix.
    EigVec assembleStandardOrthogonalBoundaryVector
    (
        const fvMesh& mesh,
        const volTensorField& conductivity,
        const highOrderInterp& LREInterp,
        const scalar stabilisationAlpha,
        const scalar t,
        const label dim
    )
    {
        EigVec b = EigVec::Zero(mesh.nCells());

        const vectorField& C = mesh.C();
        const scalarField& V = mesh.V();
        const labelUList& owner = mesh.owner();
        const surfaceVectorField& Cf = mesh.Cf();

        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            const word bcType =
                patch.lookupPatchField<volScalarField, scalar>("Vm").type();

            if (bcType != fixedValueFvPatchScalarField::typeName && bcType != "fixedVoltage")
            {
                continue;
            }

            forAll(patch, faceI)
            {
                const label gf = patch.start() + faceI;
                const label own = owner[gf];
                const point p = Cf.boundaryField()[patchI][faceI];

                const scalar a =
                    orthogonalDiffusionCoeff
                    (
                        mesh.Sf().boundaryField()[patchI][faceI],
                        p - C[own],
                        conductivity[own]
                    );

                b[own] += a/max(V[own], SMALL)*exactVm(p, t, dim);

                if (stabilisationAlpha > SMALL)
                {
                    const vector Sf =
                        mesh.Sf().boundaryField()[patchI][faceI];
                    const scalar area = mag(Sf) + VSMALL;
                    const vector n = Sf/area;
                    const vector dPb = p - C[own];
                    const scalar aStab =
                        stabilisationFaceCoeff
                        (
                            stabilisationAlpha,
                            area,
                            conductivity[own],
                            n,
                            dPb
                        );

                    b[own] +=
                        aStab/max(V[own], SMALL)*exactVm(p, t, dim);
                }
            }
        }

        return b;
    }

    // Same for the high-order operator, where the boundary value enters through
    // the ghost entry of the face stencil at every quadrature point, and again
    // through the stabilisation jump.
    EigVec assembleHighOrderBoundaryVector
    (
        const fvMesh& mesh,
        const volTensorField& conductivity,
        const highOrderInterp& LREInterp,
        const scalar stabilisationAlpha,
        const scalar t,
        const label dim
    )
    {
        EigVec b = EigVec::Zero(mesh.nCells());

        const CompactListList<point>& faceQP = LREInterp.faceQuadPoints();
        const CompactListList<scalar>& faceQW = LREInterp.faceQuadWeightPhysical();
        const CompactListList<label>& faceStencils = LREInterp.globalFaceStencils();
        const List<CompactListList<vector>>& faceGradCoeffs =
            LREInterp.QRGradFaceGPCoeffs();

        const scalarField& V = mesh.V();
        const labelUList& owner = mesh.owner();

        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            const word bcType =
                patch.lookupPatchField<volScalarField, scalar>("Vm").type();

            if
            (
                bcType != "fixedValue"
             && bcType != "fixedVoltage"
            )
            {
                continue;
            }

            forAll(patch, faceI)
            {
                const label globalFaceI = patch.start() + faceI;
                const label own = owner[globalFaceI];

                const vector Sf = mesh.Sf().boundaryField()[patchI][faceI];
                const scalar area = mag(Sf) + VSMALL;
                const vector n = Sf/area;

                const label ghostID = faceStencils[globalFaceI].size();

                forAll(faceQP[globalFaceI], qpI)
                {
                    const scalar w = faceQW[globalFaceI][qpI];
                    const scalar Vbc = exactVm(faceQP[globalFaceI][qpI], t, dim);

                    const vector gGhost =
                        faceGradCoeffs[globalFaceI][qpI][ghostID];

                    // Modified for cardiacFoam: physical quadrature weights
                    // already carry |Sf|, so the explicit area factor is gone.
                    const scalar fluxCoeff =
                        w*(n & (conductivity[own] & gGhost));

                    b[own] += fluxCoeff/max(V[own], SMALL)*Vbc;
                }

                if (stabilisationAlpha > SMALL)
                {
                    const point xf =
                        mesh.Cf().boundaryField()[patchI][faceI];
                    const vector Sf =
                        mesh.Sf().boundaryField()[patchI][faceI];
                    const scalar area = mag(Sf) + VSMALL;
                    const vector n = Sf/area;
                    const vector dPb = xf - mesh.C()[own];
                    const scalar a =
                        stabilisationFaceCoeff
                        (
                            stabilisationAlpha,
                            area,
                            conductivity[own],
                            n,
                            dPb
                        );

                    b[own] +=
                        a/max(V[own], SMALL)*exactVm(xf, t, dim);
                }
            }
        }

        return b;
    }

    // Solve A x = b with the Eigen backend: SparseLU (direct) or BiCGSTAB with
    // an ILUT preconditioner. Serial by construction - Eigen has no distributed
    // matrix - so this path is available only at one rank.
    EigVec solveSparseSystemEigen
    (
        const SpMat& A,
        const EigVec& b,
        const word& linearSolver,
        const scalar tol,
        const label maxIter,
        label& linearIterations,
        scalar& linearError
    )
    {
        if (linearSolver == "SparseLU")
        {
            Eigen::SparseLU<SpMat> solver;
            solver.analyzePattern(A);
            solver.factorize(A);

            if (solver.info() != Eigen::Success)
            {
                FatalErrorInFunction
                    << "SparseLU factorization failed"
                    << exit(FatalError);
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
            Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<scalar>> solver;
            solver.setTolerance(tol);
            solver.setMaxIterations(maxIter);
            solver.compute(A);

            if (solver.info() != Eigen::Success)
            {
                FatalErrorInFunction
                    << "BiCGSTAB setup failed"
                    << exit(FatalError);
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
    // dispatches to PETSc or Eigen according to linearSolverBackend and returns
    // the iteration count and final relative residual through its
    // out-parameters. Both backends take the same convergence controls, so a
    // case produces an algorithmically equivalent run either way. In parallel
    // only the PETSc path is usable.
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
        scalar& linearError
    )
    {
        if (usesPetscBackend(linearSolverBackend))
        {
            word kspType(petscKspType);
            word pcType(petscPcType);

            if (linearSolver == "SparseLU" || linearSolver == "LU")
            {
                kspType = "preonly";
                pcType = "lu";
            }
            else if
            (
                linearSolver == "BiCGSTAB"
             && petscKspType == linearSolver
            )
            {
                kspType = "bcgs";
            }

            PetscKspMatrixSolver solver;
            solver.reset
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

            return solver.solve(b, linearIterations, linearError);
        }

        return solveSparseSystemEigen
        (
            A,
            b,
            linearSolver,
            tol,
            maxIter,
            linearIterations,
            linearError
        );
    }

    // Gauss quadrature over one tetrahedron (a,b,c,d), appended to the running
    // lists. Uses a Duffy-style collapsed mapping of the cube onto the
    // tetrahedron, hence the r^2 s Jacobian in the weights.
    void appendTetQuadraturePointsAndWeights
    (
        const point& a,
        const point& b,
        const point& c,
        const point& d,
        DynamicList<point>& qPoints,
        DynamicList<scalar>& qWeights
    )
    {
        static const scalar xi[4] =
        {
            0.06943184420297371,
            0.33000947820757187,
            0.6699905217924281,
            0.9305681557970262
        };

        static const scalar wi[4] =
        {
            0.1739274225687269*0.5,
            0.3260725774312731*0.5,
            0.3260725774312731*0.5,
            0.1739274225687269*0.5
        };

        const scalar tetVol =
            mag((b - a) & ((c - a) ^ (d - a)))/6.0;

        if (tetVol <= VSMALL)
        {
            return;
        }

        for (label i = 0; i < 4; ++i)
        {
            const scalar r = xi[i];
            for (label j = 0; j < 4; ++j)
            {
                const scalar s = xi[j];
                for (label k = 0; k < 4; ++k)
                {
                    const scalar t = xi[k];

                    const scalar l0 = 1.0 - r;
                    const scalar l1 = r*(1.0 - s);
                    const scalar l2 = r*s*(1.0 - t);
                    const scalar l3 = r*s*t;

                    qPoints.append(l0*a + l1*b + l2*c + l3*d);
                    qWeights.append
                    (
                        tetVol*6.0*wi[i]*wi[j]*wi[k]*sqr(r)*s
                    );
                }
            }
        }
    }

    // Quadrature rule for an arbitrary polyhedral cell: the cell is decomposed
    // into tetrahedra on its face triangulation and the cell centre, and the
    // per-tet rules above are concatenated. Weights are PHYSICAL, i.e. they sum
    // to the cell volume.
    void cellQuadraturePointsAndWeights
    (
        const fvMesh& mesh,
        const label cellI,
        List<point>& qPoints,
        scalarField& qWeights
    )
    {
        DynamicList<point> qPointDyn;
        DynamicList<scalar> qWeightDyn;

        const cell& curCell = mesh.cells()[cellI];
        const faceList& faces = mesh.faces();
        const pointField& pts = mesh.points();
        const vectorField& faceCentres = mesh.faceCentres();
        const point& cellCentre = mesh.C()[cellI];

        forAll(curCell, cellFaceI)
        {
            const label faceI = curCell[cellFaceI];
            const face& f = faces[faceI];

            if (f.size() < 3)
            {
                continue;
            }

            const point& faceCentre = faceCentres[faceI];

            forAll(f, fpI)
            {
                const point& p0 = pts[f[fpI]];
                const point& p1 = pts[f[(fpI + 1) % f.size()]];

                appendTetQuadraturePointsAndWeights
                (
                    cellCentre,
                    faceCentre,
                    p0,
                    p1,
                    qPointDyn,
                    qWeightDyn
                );
            }
        }

        if (qPointDyn.size() == 0)
        {
            qPoints.setSize(1);
            qWeights.setSize(1);
            qPoints[0] = cellCentre;
            qWeights[0] = 1.0;
            return;
        }

        scalar wSum = 0.0;
        forAll(qWeightDyn, qI)
        {
            wSum += qWeightDyn[qI];
        }
        wSum = max(wSum, SMALL);

        qPoints.setSize(qPointDyn.size());
        qWeights.setSize(qWeightDyn.size());

        forAll(qPointDyn, qI)
        {
            qPoints[qI] = qPointDyn[qI];
            qWeights[qI] = qWeightDyn[qI]/wSum;
        }
    }

    // Impose the exact state values on Dirichlet patches at time t. The states
    // are local ODEs with no spatial coupling, so this matters only where the
    // reconstruction stencils reach a boundary.
    void updateStateBoundaryValues
    (
        volScalarField& u1,
        volScalarField& u2,
        volScalarField& u3,
        const scalar t,
        const label dim
    )
    {
        const fvMesh& mesh = u1.mesh();
        const surfaceVectorField& Cf = mesh.Cf();

        forAll(u1.boundaryField(), patchI)
        {
            fvPatchScalarField& u1p = u1.boundaryFieldRef()[patchI];
            fvPatchScalarField& u2p = u2.boundaryFieldRef()[patchI];
            fvPatchScalarField& u3p = u3.boundaryFieldRef()[patchI];

            if (u1p.type() == "empty")
            {
                continue;
            }

            const vectorField& CfPatch = Cf.boundaryField()[patchI];

            if (u1p.type() == fixedValueFvPatchScalarField::typeName)
            {
                forAll(u1p, faceI)
                {
                    u1p[faceI] = exactU1(CfPatch[faceI], t, dim);
                    u2p[faceI] = exactU2(CfPatch[faceI], t, dim);
                    u3p[faceI] = exactU3(CfPatch[faceI], t, dim);
                }
            }
        }
    }

    // Impose the exact Vm on Dirichlet patches at time t - the manufactured
    // boundary condition. Called at every stage where Vm changes, because the
    // reconstruction reads boundary values.
    void applyExactVmBoundaryValues
    (
        volScalarField& Vm,
        const scalar t,
        const label dim
    )
    {
        const fvMesh& mesh = Vm.mesh();
        const surfaceVectorField& Cf = mesh.Cf();

        forAll(Vm.boundaryField(), patchI)
        {
            fvPatchScalarField& Vp = Vm.boundaryFieldRef()[patchI];
            const word bcType = Vp.type();

            if
            (
                bcType == fixedValueFvPatchScalarField::typeName
             || bcType == "fixedVoltage"
            )
            {
                const vectorField& CfPatch = Cf.boundaryField()[patchI];

                forAll(Vp, faceI)
                {
                    Vp[faceI] = exactVm(CfPatch[faceI], t, dim);
                }
            }
        }

        Vm.correctBoundaryConditions();
    }

    // Evaluate the exact solution at the cell centres into the *Exact fields,
    // for the error measures and for the written output.
    void fillExactFields
    (
        volScalarField& VmExact,
        volScalarField& u1Exact,
        volScalarField& u2Exact,
        const scalar t,
        const label dim
    )
    {
        const vectorField& C = VmExact.mesh().C();

        forAll(C, cellI)
        {
            VmExact[cellI] = exactVm(C[cellI], t, dim);
            u1Exact[cellI] = exactU1(C[cellI], t, dim);
            u2Exact[cellI] = exactU2(C[cellI], t, dim);
        }

        const fvMesh& mesh = VmExact.mesh();
        const surfaceVectorField& Cf = mesh.Cf();

        forAll(VmExact.boundaryField(), patchI)
        {
            fvPatchScalarField& Vp = VmExact.boundaryFieldRef()[patchI];
            fvPatchScalarField& u1p = u1Exact.boundaryFieldRef()[patchI];
            fvPatchScalarField& u2p = u2Exact.boundaryFieldRef()[patchI];

            if (Vp.type() == "empty")
            {
                continue;
            }

            const vectorField& CfPatch = Cf.boundaryField()[patchI];

            if
            (
                Vp.type() == fixedValueFvPatchScalarField::typeName
             || Vp.type() == "fixedVoltage"
            )
            {
                forAll(Vp, faceI)
                {
                    Vp[faceI] = exactVm(CfPatch[faceI], t, dim);
                }
            }

            if (u1p.type() == fixedValueFvPatchScalarField::typeName)
            {
                forAll(u1p, faceI)
                {
                    u1p[faceI] = exactU1(CfPatch[faceI], t, dim);
                    u2p[faceI] = exactU2(CfPatch[faceI], t, dim);
                }
            }
        }

        VmExact.correctBoundaryConditions();
        u1Exact.correctBoundaryConditions();
        u2Exact.correctBoundaryConditions();
    }

    // EXPLICIT high-order diffusion: evaluates the face fluxes
    // n . (D grad Vm) by LRE quadrature and takes their divergence, writing
    // both the face flux field and the cell Laplacian. The diagnostic
    // counterpart of assembleHighOrderStiffnessMatrix, which builds the same
    // operator as a matrix.
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

        // Modified for cardiacFoam: gradScalarFaceQuad now returns a
        // CompactListList (contiguous storage); it indexes identically.
        autoPtr<CompactListList<vector>> gradVmQuadPtr =
            LREInterp_Vm.gradScalarFaceQuad(Vm);

        CompactListList<vector>& gradVmQuad = gradVmQuadPtr.ref();
        const CompactListList<scalar>& faceQuadW =
            LREInterp_Vm.faceQuadWeightPhysical();

        const surfaceVectorField nHat(mesh.Sf()/mesh.magSf());
        const scalarField& magSfInternal = mesh.magSf().internalField();

        for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
        {
            const label owner = mesh.owner()[faceI];
            const vector& faceNormal = nHat[faceI];
            const scalar faceArea = magSfInternal[faceI];

            fluxVm_HO[faceI] = 0.0;

            forAll(gradVmQuad[faceI], pI)
            {
                const vector Dg = conductivity[owner] & gradVmQuad[faceI][pI];

                // Modified for cardiacFoam: physical quadrature weights
                // already carry |Sf|, so faceArea is no longer applied here.
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
            const scalarField& pMagSf = mesh.magSf().boundaryField()[patchI];
            const vectorField& pNormals = nHat.boundaryField()[patchI];

            forAll(patchFlux, faceI)
            {
                const label globalFaceI = start + faceI;
                const label owner = mesh.owner()[globalFaceI];
                const vector& faceNormal = pNormals[faceI];
                const scalar faceArea = pMagSf[faceI];

                patchFlux[faceI] = 0.0;

                forAll(gradVmQuad[globalFaceI], pI)
                {
                    const vector Dg =
                        conductivity[owner] & gradVmQuad[globalFaceI][pI];

                    // Modified for cardiacFoam: see the internal-face loop.
                    patchFlux[faceI] +=
                        (faceNormal & Dg)*faceQuadW[globalFaceI][pI];
                }
            }
        }

        lapVm = fvc::div(fluxVm_HO);
    }

    // Ionic current evaluated once per cell at the cell centre: the low-order
    // path, and the reference the quadrature path is compared against.
    void computeCellCentredIion
    (
        const volScalarField& Vm,
        const volScalarField& u1,
        const volScalarField& u2,
        const volScalarField& u3,
        const scalar beta,
        const scalar chiVal,
        const scalar CmVal,
        volScalarField& Iion,
        const bool useOpenMP = false,
        const label openMPThreshold = 256
    )
    {
        const label nCells = Vm.internalField().size();

        #pragma omp parallel for schedule(static) if(useOpenMP && nCells >= openMPThreshold)
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            Iion[cellI] =
                ionicCurrentPDE
                (
                    Vm[cellI],
                    u1[cellI],
                    u2[cellI],
                    u3[cellI],
                    beta,
                    chiVal,
                    CmVal
                );
        }

        Iion.correctBoundaryConditions();
    }

    // Vm source evaluated once per cell at the cell centre.
    void computeCellCentredVmSource
    (
        const volScalarField& Vm,
        const volScalarField& u1,
        const volScalarField& u2,
        const volScalarField& u3,
        const scalar beta,
        const scalar chiVal,
        const scalar CmVal,
        volScalarField& sourceVm,
        const bool useOpenMP = false,
        const label openMPThreshold = 256
    )
    {
        const label nCells = Vm.internalField().size();

        #pragma omp parallel for schedule(static) if(useOpenMP && nCells >= openMPThreshold)
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            sourceVm[cellI] =
                vmSourcePDE
                (
                    Vm[cellI],
                    u1[cellI],
                    u2[cellI],
                    u3[cellI],
                    beta,
                    chiVal,
                    CmVal
                );
        }

        sourceVm.correctBoundaryConditions();
    }

    // Reconstruct Vm at every ionic quadrature point, from the cell values, at
    // the order of the interpolator. This is what makes the source term
    // high-order: evaluating a nonlinear function at the cell centre and
    // multiplying by the volume is only first-order accurate however good the
    // diffusion operator is.
    void reconstructVmAtIionIntegrationPoints
    (
        const volScalarField& Vm,
        const Switch useHighOrderVm,
        const highOrderInterp& LREInterp_Vm,
        const highOrderInterp& LREInterp_Iion,
        scalarField& VmIntegrationPoints
    )
    {
        const fvMesh& mesh = Vm.mesh();
        const vectorField& C = mesh.C();
        const CompactListList<point>& cellIionQuadP =
            LREInterp_Iion.cellQuadPoints();

        label integrationPointI = 0;

        if (useHighOrderVm)
        {
            const bool twoD = mesh.nGeometricD() == 2;

            tmp<volVectorField> tGradVm = LREInterp_Vm.grad(Vm);
            const vectorField& gradVm = tGradVm->internalField();

            tmp<volSymmTensorField> tHessVm;
            const symmTensorField* hessVm = nullptr;
            if (LREInterp_Vm.order() >= 2)
            {
                tHessVm = LREInterp_Vm.hessian(Vm);
                hessVm = &(tHessVm->internalField());
            }

            autoPtr<List<highOrderInterp::symmTensor3Order>> thirdVmPtr;
            const List<highOrderInterp::symmTensor3Order>* thirdVm = nullptr;
            if (LREInterp_Vm.order() >= 3)
            {
                thirdVmPtr = LREInterp_Vm.thirdDeriv(Vm);
                thirdVm = &thirdVmPtr();
            }

            forAll(mesh.cells(), cellI)
            {
                const scalar Vc = Vm[cellI];
                const vector& gradVc = gradVm[cellI];
                const vector& xc = C[cellI];

                const symmTensor* H = hessVm ? &((*hessVm)[cellI]) : nullptr;
                const highOrderInterp::symmTensor3Order* T3 =
                    thirdVm ? &((*thirdVm)[cellI]) : nullptr;

                forAll(cellIionQuadP[cellI], qI)
                {
                    const vector d = cellIionQuadP[cellI][qI] - xc;
                    VmIntegrationPoints[integrationPointI] =
                        reconstructFromTaylor(Vc, gradVc, H, T3, d, twoD);
                    ++integrationPointI;
                }
            }
        }
        else
        {
            tmp<volVectorField> tGradVm = fvc::grad(Vm);
            const vectorField& gradVm = tGradVm->internalField();

            forAll(mesh.cells(), cellI)
            {
                const scalar Vc = Vm[cellI];
                const vector& gradVc = gradVm[cellI];
                const vector& xc = C[cellI];

                forAll(cellIionQuadP[cellI], qI)
                {
                    const vector d = cellIionQuadP[cellI][qI] - xc;
                    VmIntegrationPoints[integrationPointI] = Vc + (gradVc & d);
                    ++integrationPointI;
                }
            }
        }
    }

    // Reconstruct the cell-centred states (u1,u2,u3) at the Iion Gauss points
    // using a dedicated high-order LRE (LREInterp_states). Mirrors the
    // high-order branch of reconstructVmAtIionIntegrationPoints, but evaluates
    // the three state fields in a single pass over the cells. Used by
    // stateIntegrationMode = cellCentredReconstruct.
    void reconstructStatesAtIionIntegrationPoints
    (
        const volScalarField& u1,
        const volScalarField& u2,
        const volScalarField& u3,
        const highOrderInterp& LREInterp_states,
        const highOrderInterp& LREInterp_Iion,
        scalarField& u1IntegrationPoints,
        scalarField& u2IntegrationPoints,
        scalarField& u3IntegrationPoints
    )
    {
        const fvMesh& mesh = u1.mesh();
        const vectorField& C = mesh.C();
        const bool twoD = mesh.nGeometricD() == 2;
        const CompactListList<point>& cellIionQuadP =
            LREInterp_Iion.cellQuadPoints();

        // Gradients (always), Hessians (order >= 2), third derivs (order >= 3)
        // of each state field from the state interpolator.
        tmp<volVectorField> tGradU1 = LREInterp_states.grad(u1);
        tmp<volVectorField> tGradU2 = LREInterp_states.grad(u2);
        tmp<volVectorField> tGradU3 = LREInterp_states.grad(u3);
        const vectorField& gradU1 = tGradU1->internalField();
        const vectorField& gradU2 = tGradU2->internalField();
        const vectorField& gradU3 = tGradU3->internalField();

        tmp<volSymmTensorField> tHessU1, tHessU2, tHessU3;
        const symmTensorField* hessU1 = nullptr;
        const symmTensorField* hessU2 = nullptr;
        const symmTensorField* hessU3 = nullptr;
        if (LREInterp_states.order() >= 2)
        {
            tHessU1 = LREInterp_states.hessian(u1);
            tHessU2 = LREInterp_states.hessian(u2);
            tHessU3 = LREInterp_states.hessian(u3);
            hessU1 = &(tHessU1->internalField());
            hessU2 = &(tHessU2->internalField());
            hessU3 = &(tHessU3->internalField());
        }

        autoPtr<List<highOrderInterp::symmTensor3Order>> thirdU1Ptr, thirdU2Ptr, thirdU3Ptr;
        const List<highOrderInterp::symmTensor3Order>* thirdU1 = nullptr;
        const List<highOrderInterp::symmTensor3Order>* thirdU2 = nullptr;
        const List<highOrderInterp::symmTensor3Order>* thirdU3 = nullptr;
        if (LREInterp_states.order() >= 3)
        {
            thirdU1Ptr = LREInterp_states.thirdDeriv(u1);
            thirdU2Ptr = LREInterp_states.thirdDeriv(u2);
            thirdU3Ptr = LREInterp_states.thirdDeriv(u3);
            thirdU1 = &thirdU1Ptr();
            thirdU2 = &thirdU2Ptr();
            thirdU3 = &thirdU3Ptr();
        }

        label integrationPointI = 0;
        forAll(mesh.cells(), cellI)
        {
            const vector& xc = C[cellI];

            const scalar u1c = u1[cellI];
            const scalar u2c = u2[cellI];
            const scalar u3c = u3[cellI];

            const vector& gU1 = gradU1[cellI];
            const vector& gU2 = gradU2[cellI];
            const vector& gU3 = gradU3[cellI];

            const symmTensor* HU1 = hessU1 ? &((*hessU1)[cellI]) : nullptr;
            const symmTensor* HU2 = hessU2 ? &((*hessU2)[cellI]) : nullptr;
            const symmTensor* HU3 = hessU3 ? &((*hessU3)[cellI]) : nullptr;

            const highOrderInterp::symmTensor3Order* TU1 =
                thirdU1 ? &((*thirdU1)[cellI]) : nullptr;
            const highOrderInterp::symmTensor3Order* TU2 =
                thirdU2 ? &((*thirdU2)[cellI]) : nullptr;
            const highOrderInterp::symmTensor3Order* TU3 =
                thirdU3 ? &((*thirdU3)[cellI]) : nullptr;

            forAll(cellIionQuadP[cellI], qI)
            {
                const vector d = cellIionQuadP[cellI][qI] - xc;
                u1IntegrationPoints[integrationPointI] =
                    reconstructFromTaylor(u1c, gU1, HU1, TU1, d, twoD);
                u2IntegrationPoints[integrationPointI] =
                    reconstructFromTaylor(u2c, gU2, HU2, TU2, d, twoD);
                u3IntegrationPoints[integrationPointI] =
                    reconstructFromTaylor(u3c, gU3, HU3, TU3, d, twoD);
                ++integrationPointI;
            }
        }
    }

    // Collapse a quantity defined at the quadrature points back to a cell
    // average, weighting by the physical quadrature weights. The inverse
    // direction of the reconstruction above.
    void averageIntegrationPointFieldToCells
    (
        const scalarField& integrationPointValues,
        const highOrderInterp& LREInterp_Iion,
        volScalarField& field
    )
    {
        const fvMesh& mesh = field.mesh();
        const CompactListList<scalar>& cellIionQuadW =
            LREInterp_Iion.cellQuadWeightPhysical();

        scalarField& cellValues = field.primitiveFieldRef();
        label integrationPointI = 0;

        forAll(mesh.cells(), cellI)
        {
            scalar valueBar = 0.0;
            scalar wSum = 0.0;

            forAll(cellIionQuadW[cellI], qI)
            {
                const scalar w = cellIionQuadW[cellI][qI];
                valueBar += w*integrationPointValues[integrationPointI];
                wSum += w;
                ++integrationPointI;
            }

            cellValues[cellI] = valueBar/max(wSum, SMALL);
        }

        field.correctBoundaryConditions();
    }

    inline void clampPhysical
    (
        scalarField& values,
        const scalar vMin,
        const scalar vMax
    )
    {
        forAll(values, i)
        {
            if (values[i] < vMin) values[i] = vMin;
            else if (values[i] > vMax) values[i] = vMax;
        }
    }

    // Seed the per-quadrature-point states from the exact solution at t = 0.
    // Necessary because those states are independent unknowns in the
    // gaussPointODE mode - they are not derived from the cell values.
    void initialiseIionIntegrationPointStates
    (
        const highOrderInterp& LREInterp_Iion,
        const scalar t,
        const label dim,
        scalarField& u1IntegrationPoints,
        scalarField& u2IntegrationPoints,
        scalarField& u3IntegrationPoints
    )
    {
        const CompactListList<point>& cellIionQuadP =
            LREInterp_Iion.cellQuadPoints();

        label integrationPointI = 0;
        forAll(cellIionQuadP, cellI)
        {
            forAll(cellIionQuadP[cellI], qI)
            {
                const point& p = cellIionQuadP[cellI][qI];
                u1IntegrationPoints[integrationPointI] = exactU1(p, t, dim);
                u2IntegrationPoints[integrationPointI] = exactU2(p, t, dim);
                u3IntegrationPoints[integrationPointI] = exactU3(p, t, dim);
                ++integrationPointI;
            }
        }
    }

    // Cell-averaged ionic current obtained by evaluating Iion at every
    // quadrature point and integrating, rather than evaluating it once at the
    // cell centre.
    void computeIionFromIntegrationPoints
    (
        const scalarField& VmIntegrationPoints,
        const scalarField& u1IntegrationPoints,
        const scalarField& u2IntegrationPoints,
        const scalarField& u3IntegrationPoints,
        const scalar beta,
        const scalar chiVal,
        const scalar CmVal,
        scalarField& IionIntegrationPoints,
        const bool useOpenMP = false,
        const label openMPThreshold = 256
    )
    {
        const label nIntegrationPoints = IionIntegrationPoints.size();

        #pragma omp parallel for schedule(static) if(useOpenMP && nIntegrationPoints >= openMPThreshold)
        for (label integrationPointI = 0; integrationPointI < nIntegrationPoints; ++integrationPointI)
        {
            IionIntegrationPoints[integrationPointI] =
                ionicCurrentPDE
                (
                    VmIntegrationPoints[integrationPointI],
                    u1IntegrationPoints[integrationPointI],
                    u2IntegrationPoints[integrationPointI],
                    u3IntegrationPoints[integrationPointI],
                    beta,
                    chiVal,
                    CmVal
                );
        }
    }

    // Same for the Vm source term.
    void computeVmSourceFromIntegrationPoints
    (
        const scalarField& VmIntegrationPoints,
        const scalarField& u1IntegrationPoints,
        const scalarField& u2IntegrationPoints,
        const scalarField& u3IntegrationPoints,
        const scalar beta,
        const scalar chiVal,
        const scalar CmVal,
        scalarField& sourceIntegrationPoints,
        const bool useOpenMP = false,
        const label openMPThreshold = 256
    )
    {
        const label nIntegrationPoints = sourceIntegrationPoints.size();

        #pragma omp parallel for schedule(static) if(useOpenMP && nIntegrationPoints >= openMPThreshold)
        for (label integrationPointI = 0; integrationPointI < nIntegrationPoints; ++integrationPointI)
        {
            sourceIntegrationPoints[integrationPointI] =
                vmSourcePDE
                (
                    VmIntegrationPoints[integrationPointI],
                    u1IntegrationPoints[integrationPointI],
                    u2IntegrationPoints[integrationPointI],
                    u3IntegrationPoints[integrationPointI],
                    beta,
                    chiVal,
                    CmVal
                );
        }
    }

    // Reaction rates with Vm linearly interpolated in time across the step,
    // which is what the state ODEs see while the PDE is advanced from t^n to
    // t^{n+1}. Without it the states would be integrated against a frozen Vm
    // and the coupling would drop to first order.
    void reactionRatesLinearVm
    (
        const scalar VmOld,
        const scalar VmNew,
        const scalar tau,
        const scalar dt,
        const scalar u1,
        const scalar u2,
        const scalar u3,
        scalar& du1dt,
        scalar& du2dt,
        scalar& du3dt
    )
    {
        const scalar alpha =
            dt > SMALL
          ? min(max(tau/dt, scalar(0.0)), scalar(1.0))
          : scalar(1.0);

        const scalar VmTau = (1.0 - alpha)*VmOld + alpha*VmNew;

        reactionRates(VmTau, u1, u2, u3, du1dt, du2dt, du3dt);
    }

    // One classical fourth-order Runge-Kutta step of the state ODEs at one
    // point. Fixed step: used when the adaptive integrator is not requested.
    void rk4StateStep
    (
        const scalar VmOld,
        const scalar VmNew,
        const scalar tau,
        const scalar h,
        const scalar dt,
        const scalar u1,
        const scalar u2,
        const scalar u3,
        scalar& u1New,
        scalar& u2New,
        scalar& u3New
    )
    {
        scalar k11 = 0.0, k12 = 0.0, k13 = 0.0;
        scalar k21 = 0.0, k22 = 0.0, k23 = 0.0;
        scalar k31 = 0.0, k32 = 0.0, k33 = 0.0;
        scalar k41 = 0.0, k42 = 0.0, k43 = 0.0;

        reactionRatesLinearVm
        (
            VmOld, VmNew, tau, dt,
            u1, u2, u3,
            k11, k12, k13
        );
        reactionRatesLinearVm
        (
            VmOld, VmNew, tau + 0.5*h, dt,
            u1 + 0.5*h*k11,
            u2 + 0.5*h*k12,
            u3 + 0.5*h*k13,
            k21, k22, k23
        );
        reactionRatesLinearVm
        (
            VmOld, VmNew, tau + 0.5*h, dt,
            u1 + 0.5*h*k21,
            u2 + 0.5*h*k22,
            u3 + 0.5*h*k23,
            k31, k32, k33
        );
        reactionRatesLinearVm
        (
            VmOld, VmNew, tau + h, dt,
            u1 + h*k31,
            u2 + h*k32,
            u3 + h*k33,
            k41, k42, k43
        );

        u1New = u1 + (h/6.0)*(k11 + 2.0*k21 + 2.0*k31 + k41);
        u2New = u2 + (h/6.0)*(k12 + 2.0*k22 + 2.0*k32 + k42);
        u3New = u3 + (h/6.0)*(k13 + 2.0*k23 + 2.0*k33 + k43);
    }

    // One Runge-Kutta-Fehlberg 4(5) step with an embedded error estimate, from
    // which the next step size is chosen. The state ODEs are stiff at the
    // upstroke, so the adaptive step is what keeps the local error below the
    // tolerance without paying for it everywhere.
    void rkf45StateStep
    (
        const scalar VmOld,
        const scalar VmNew,
        const scalar tau,
        const scalar h,
        const scalar dt,
        const scalar u1,
        const scalar u2,
        const scalar u3,
        scalar& u1Fifth,
        scalar& u2Fifth,
        scalar& u3Fifth,
        const scalar stateODEAbsTol,
        const scalar stateODERelTol,
        scalar& err
    )
    {
        scalar k11 = 0.0, k12 = 0.0, k13 = 0.0;
        scalar k21 = 0.0, k22 = 0.0, k23 = 0.0;
        scalar k31 = 0.0, k32 = 0.0, k33 = 0.0;
        scalar k41 = 0.0, k42 = 0.0, k43 = 0.0;
        scalar k51 = 0.0, k52 = 0.0, k53 = 0.0;
        scalar k61 = 0.0, k62 = 0.0, k63 = 0.0;

        reactionRatesLinearVm
        (
            VmOld, VmNew, tau, dt,
            u1, u2, u3,
            k11, k12, k13
        );
        reactionRatesLinearVm
        (
            VmOld, VmNew, tau + h/5.0, dt,
            u1 + h*(1.0/5.0)*k11,
            u2 + h*(1.0/5.0)*k12,
            u3 + h*(1.0/5.0)*k13,
            k21, k22, k23
        );
        reactionRatesLinearVm
        (
            VmOld, VmNew, tau + 3.0*h/10.0, dt,
            u1 + h*((3.0/40.0)*k11 + (9.0/40.0)*k21),
            u2 + h*((3.0/40.0)*k12 + (9.0/40.0)*k22),
            u3 + h*((3.0/40.0)*k13 + (9.0/40.0)*k23),
            k31, k32, k33
        );
        reactionRatesLinearVm
        (
            VmOld, VmNew, tau + 3.0*h/5.0, dt,
            u1 + h*((3.0/10.0)*k11 - (9.0/10.0)*k21 + (6.0/5.0)*k31),
            u2 + h*((3.0/10.0)*k12 - (9.0/10.0)*k22 + (6.0/5.0)*k32),
            u3 + h*((3.0/10.0)*k13 - (9.0/10.0)*k23 + (6.0/5.0)*k33),
            k41, k42, k43
        );
        reactionRatesLinearVm
        (
            VmOld, VmNew, tau + h, dt,
            u1 + h*((-11.0/54.0)*k11 + (5.0/2.0)*k21 - (70.0/27.0)*k31 + (35.0/27.0)*k41),
            u2 + h*((-11.0/54.0)*k12 + (5.0/2.0)*k22 - (70.0/27.0)*k32 + (35.0/27.0)*k42),
            u3 + h*((-11.0/54.0)*k13 + (5.0/2.0)*k23 - (70.0/27.0)*k33 + (35.0/27.0)*k43),
            k51, k52, k53
        );
        reactionRatesLinearVm
        (
            VmOld, VmNew, tau + 7.0*h/8.0, dt,
            u1 + h*((1631.0/55296.0)*k11 + (175.0/512.0)*k21 + (575.0/13824.0)*k31 + (44275.0/110592.0)*k41 + (253.0/4096.0)*k51),
            u2 + h*((1631.0/55296.0)*k12 + (175.0/512.0)*k22 + (575.0/13824.0)*k32 + (44275.0/110592.0)*k42 + (253.0/4096.0)*k52),
            u3 + h*((1631.0/55296.0)*k13 + (175.0/512.0)*k23 + (575.0/13824.0)*k33 + (44275.0/110592.0)*k43 + (253.0/4096.0)*k53),
            k61, k62, k63
        );

        u1Fifth = u1 + h*((37.0/378.0)*k11 + (250.0/621.0)*k31 + (125.0/594.0)*k41 + (512.0/1771.0)*k61);
        u2Fifth = u2 + h*((37.0/378.0)*k12 + (250.0/621.0)*k32 + (125.0/594.0)*k42 + (512.0/1771.0)*k62);
        u3Fifth = u3 + h*((37.0/378.0)*k13 + (250.0/621.0)*k33 + (125.0/594.0)*k43 + (512.0/1771.0)*k63);

        const scalar u1Fourth = u1 + h*((2825.0/27648.0)*k11 + (18575.0/48384.0)*k31 + (13525.0/55296.0)*k41 + (277.0/14336.0)*k51 + 0.25*k61);
        const scalar u2Fourth = u2 + h*((2825.0/27648.0)*k12 + (18575.0/48384.0)*k32 + (13525.0/55296.0)*k42 + (277.0/14336.0)*k52 + 0.25*k62);
        const scalar u3Fourth = u3 + h*((2825.0/27648.0)*k13 + (18575.0/48384.0)*k33 + (13525.0/55296.0)*k43 + (277.0/14336.0)*k53 + 0.25*k63);

        err = max
        (
            mag(u1Fifth - u1Fourth)
           /(stateODEAbsTol + stateODERelTol*max(mag(u1Fifth), mag(u1Fourth))),
            max
            (
                mag(u2Fifth - u2Fourth)
               /(stateODEAbsTol + stateODERelTol*max(mag(u2Fifth), mag(u2Fourth))),
                mag(u3Fifth - u3Fourth)
               /(stateODEAbsTol + stateODERelTol*max(mag(u3Fifth), mag(u3Fourth)))
            )
        );
    }

    // Integrate the state ODEs at one point across a full PDE time step,
    // sub-stepping adaptively until the interval is covered. Vm is interpolated
    // linearly across the step (reactionRatesLinearVm), so the sub-steps see a
    // consistent voltage.
    void advanceStateODE
    (
        const scalar VmOld,
        const scalar VmNew,
        const scalar u1Old,
        const scalar u2Old,
        const scalar u3Old,
        const scalar dt,
        const word& stateODESolver,
        const scalar stateODEInitialStep,
        const scalar stateODEAbsTol,
        const scalar stateODERelTol,
        const label stateODEMaxSteps,
        scalar& u1New,
        scalar& u2New,
        scalar& u3New
    )
    {
        if (dt <= SMALL)
        {
            u1New = u1Old;
            u2New = u2Old;
            u3New = u3Old;
            return;
        }

        if (stateODESolver == "Euler" || stateODESolver == "forwardEuler")
        {
            scalar du1 = 0.0, du2 = 0.0, du3 = 0.0;
            reactionRatesLinearVm
            (
                VmOld, VmNew, 0.0, dt,
                u1Old, u2Old, u3Old,
                du1, du2, du3
            );

            u1New = u1Old + dt*du1;
            u2New = u2Old + dt*du2;
            u3New = u3Old + dt*du3;
            return;
        }

        if (stateODESolver == "RK4")
        {
            rk4StateStep
            (
                VmOld, VmNew, 0.0, dt, dt,
                u1Old, u2Old, u3Old,
                u1New, u2New, u3New
            );
            return;
        }

        if
        (
            stateODESolver != "RKF45"
         && stateODESolver != "RKT45"
         && stateODESolver != "rkf45"
        )
        {
            FatalErrorInFunction
                << "Unknown state ODE solver: " << stateODESolver << nl
                << "Valid options are RKF45, RKT45, RK4, Euler"
                << exit(FatalError);
        }

        scalar tau = 0.0;
        scalar h =
            stateODEInitialStep > SMALL
          ? min(dt, stateODEInitialStep)
          : dt;
        scalar u1Cur = u1Old;
        scalar u2Cur = u2Old;
        scalar u3Cur = u3Old;
        const scalar absTol = max(stateODEAbsTol, scalar(SMALL));
        const scalar relTol = max(stateODERelTol, scalar(SMALL));
        const scalar hMin = max(dt*1.0e-12, SMALL);

        for (label subStep = 0; subStep < max(stateODEMaxSteps, label(1)); ++subStep)
        {
            if (tau >= dt - SMALL)
            {
                break;
            }

            h = min(h, dt - tau);

            scalar u1Trial = u1Cur;
            scalar u2Trial = u2Cur;
            scalar u3Trial = u3Cur;
            scalar err = GREAT;

            rkf45StateStep
            (
                VmOld, VmNew, tau, h, dt,
                u1Cur, u2Cur, u3Cur,
                u1Trial, u2Trial, u3Trial,
                absTol,
                relTol,
                err
            );

            if (err <= 1.0 || h <= hMin)
            {
                u1Cur = u1Trial;
                u2Cur = u2Trial;
                u3Cur = u3Trial;
                tau += h;
            }

            const scalar factor = min
            (
                scalar(4.0),
                max(scalar(0.1), scalar(0.84)*std::pow(1.0/(err + SMALL), 0.25))
            );
            h = max(hMin, min(dt - tau, h*factor));
        }

        if (tau < dt - 10.0*SMALL)
        {
            WarningInFunction
                << "State RKF45 reached maxSteps before completing dt. "
                << "tau = " << tau << ", dt = " << dt << nl;
        }

        u1New = u1Cur;
        u2New = u2Cur;
        u3New = u3Cur;
    }

    // Advance the states at every ionic quadrature point over one PDE step.
    // Each point is an independent ODE, hence the optional OpenMP loop.
    void updateIntegrationPointStatesODE
    (
        const scalarField& VmOldIntegrationPoints,
        const scalarField& VmNewIntegrationPoints,
        const scalarField& u1OldIntegrationPoints,
        const scalarField& u2OldIntegrationPoints,
        const scalarField& u3OldIntegrationPoints,
        scalarField& u1NewIntegrationPoints,
        scalarField& u2NewIntegrationPoints,
        scalarField& u3NewIntegrationPoints,
        const scalar dt,
        const word& stateODESolver,
        const scalar stateODEInitialStep,
        const scalar stateODEAbsTol,
        const scalar stateODERelTol,
        const label stateODEMaxSteps,
        const bool useOpenMP = false,
        const label openMPThreshold = 256
    )
    {
        const label nIntegrationPoints = u1NewIntegrationPoints.size();

        #pragma omp parallel for schedule(static) if(useOpenMP && nIntegrationPoints >= openMPThreshold)
        for (label integrationPointI = 0; integrationPointI < nIntegrationPoints; ++integrationPointI)
        {
            advanceStateODE
            (
                VmOldIntegrationPoints[integrationPointI],
                VmNewIntegrationPoints[integrationPointI],
                u1OldIntegrationPoints[integrationPointI],
                u2OldIntegrationPoints[integrationPointI],
                u3OldIntegrationPoints[integrationPointI],
                dt,
                stateODESolver,
                stateODEInitialStep,
                stateODEAbsTol,
                stateODERelTol,
                stateODEMaxSteps,
                u1NewIntegrationPoints[integrationPointI],
                u2NewIntegrationPoints[integrationPointI],
                u3NewIntegrationPoints[integrationPointI]
            );
        }
    }

    // Advance the cell-centred states over one PDE step - the counterpart of
    // the above for the cellCentredReconstruct mode, where the ODE is solved
    // once per cell and the states are reconstructed at the quadrature points
    // when Iion needs them.
    void updateStateFieldsODE
    (
        const volScalarField& VmOld,
        const volScalarField& VmNew,
        const volScalarField& u1Old,
        const volScalarField& u2Old,
        const volScalarField& u3Old,
        volScalarField& u1New,
        volScalarField& u2New,
        volScalarField& u3New,
        const scalar dt,
        const word& stateODESolver,
        const scalar stateODEInitialStep,
        const scalar stateODEAbsTol,
        const scalar stateODERelTol,
        const label stateODEMaxSteps,
        const bool useOpenMP = false,
        const label openMPThreshold = 256
    )
    {
        scalarField& u1NI = u1New.primitiveFieldRef();
        scalarField& u2NI = u2New.primitiveFieldRef();
        scalarField& u3NI = u3New.primitiveFieldRef();

        const scalarField& VO = VmOld.primitiveField();
        const scalarField& VN = VmNew.primitiveField();
        const scalarField& u1O = u1Old.primitiveField();
        const scalarField& u2O = u2Old.primitiveField();
        const scalarField& u3O = u3Old.primitiveField();

        const label nCells = u1NI.size();

        #pragma omp parallel for schedule(static) if(useOpenMP && nCells >= openMPThreshold)
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            advanceStateODE
            (
                VO[cellI],
                VN[cellI],
                u1O[cellI],
                u2O[cellI],
                u3O[cellI],
                dt,
                stateODESolver,
                stateODEInitialStep,
                stateODEAbsTol,
                stateODERelTol,
                stateODEMaxSteps,
                u1NI[cellI],
                u2NI[cellI],
                u3NI[cellI]
            );
        }
    }

    ReconstructionErrorSummary computeGaussReconstructedError
    (
        const volScalarField& fld,
        const scalar t,
        const ManufacturedFieldID fieldID,
        const label dim,
        const highOrderInterp& LREInterp
    )
    {
        const fvMesh& mesh = fld.mesh();
        const bool useTaylor = mesh.nGeometricD() > 1;
        const bool twoD = mesh.nGeometricD() == 2;

        tmp<volVectorField> tGrad;
        const volVectorField* gradPtr = nullptr;
        if (useTaylor)
        {
            tGrad = LREInterp.grad(fld);
            gradPtr = &tGrad();
        }

        tmp<volSymmTensorField> tHess;
        const volSymmTensorField* hessPtr = nullptr;
        if (useTaylor && LREInterp.order() >= 2)
        {
            tHess = LREInterp.hessian(fld);
            hessPtr = &tHess();
        }

        autoPtr<List<highOrderInterp::symmTensor3Order>> thirdPtr;
        const List<highOrderInterp::symmTensor3Order>* thirdList = nullptr;
        if (useTaylor && LREInterp.order() >= 3)
        {
            thirdPtr = LREInterp.thirdDeriv(fld);
            thirdList = &thirdPtr();
        }

        const vectorField& C = mesh.C();
        const scalarField& V = mesh.V();

        scalar e1 = 0.0;
        scalar e2 = 0.0;
        scalar eInf = 0.0;
        scalar n2 = 0.0;
        scalar volTot = 0.0;

        List<point> qPoints;
        scalarField qWeights;

        forAll(fld.internalField(), cellI)
        {
            cellQuadraturePointsAndWeights(mesh, cellI, qPoints, qWeights);

            const scalar cellV = V[cellI];

            forAll(qPoints, qpI)
            {
                const point& gp = qPoints[qpI];
                const scalar w = qWeights[qpI];
                const scalar meas = cellV*w;
                const vector d = gp - C[cellI];

                scalar num = fld[cellI];

                if (gradPtr)
                {
                    const symmTensor* H =
                        hessPtr ? &((*hessPtr)[cellI]) : nullptr;
                    const highOrderInterp::symmTensor3Order* T3 =
                        thirdList ? &((*thirdList)[cellI]) : nullptr;

                    num =
                        reconstructFromTaylor
                        (
                            fld[cellI],
                            (*gradPtr)[cellI],
                            H,
                            T3,
                            d,
                            twoD
                        );
                }

                const scalar ex = exactFieldValue(fieldID, gp, t, dim);
                const scalar err = num - ex;

                e1 += meas*mag(err);
                e2 += meas*sqr(err);
                eInf = max(eInf, mag(err));
                n2 += meas*sqr(ex);
                volTot += meas;
            }
        }

        e1 = returnReduce(e1, sumOp<scalar>());
        e2 = returnReduce(e2, sumOp<scalar>());
        eInf = returnReduce(eInf, maxOp<scalar>());
        n2 = returnReduce(n2, sumOp<scalar>());
        volTot = returnReduce(volTot, sumOp<scalar>());

        ReconstructionErrorSummary s;
        s.L1 = e1/max(volTot, SMALL);
        s.L2 = std::sqrt(e2/max(volTot, SMALL));
        s.Linf = eInf;
        s.normL2 = std::sqrt(n2/max(volTot, SMALL));
        s.relL2 = 100.0*s.L2/max(s.normL2, SMALL);

        return s;
    }

    ReconstructionErrorSummary cellAsReconstructionError
    (
        const FieldErrorSummary& cellErr,
        const volScalarField& exact
    )
    {
        ReconstructionErrorSummary s;
        s.L1 = cellErr.L1;
        s.L2 = cellErr.L2;
        s.Linf = cellErr.Linf;
        s.normL2 = volumeWeightedL2(exact);
        s.relL2 = 100.0*s.L2/max(s.normL2, SMALL);
        return s;
    }

    ReconstructionErrorSummary computeIntegrationPointStateError
    (
        const scalarField& stateIntegrationPoints,
        const scalar t,
        const ManufacturedFieldID fieldID,
        const label dim,
        const highOrderInterp& LREInterp_Iion,
        const fvMesh& mesh
    )
    {
        const CompactListList<point>& cellIionQuadP =
            LREInterp_Iion.cellQuadPoints();
        const CompactListList<scalar>& cellIionQuadW =
            LREInterp_Iion.cellQuadWeightPhysical();
        const scalarField& V = mesh.V();

        scalar e1 = 0.0;
        scalar e2 = 0.0;
        scalar eInf = 0.0;
        scalar n2 = 0.0;
        scalar volTot = 0.0;

        label integrationPointI = 0;
        forAll(mesh.cells(), cellI)
        {
            const scalar cellV = V[cellI];
            scalar wSum = 0.0;

            forAll(cellIionQuadW[cellI], qI)
            {
                wSum += cellIionQuadW[cellI][qI];
            }

            forAll(cellIionQuadP[cellI], qI)
            {
                const point& gp = cellIionQuadP[cellI][qI];
                const scalar meas =
                    cellV*cellIionQuadW[cellI][qI]/max(wSum, SMALL);
                const scalar ex = exactFieldValue(fieldID, gp, t, dim);
                const scalar err = stateIntegrationPoints[integrationPointI] - ex;

                e1 += meas*mag(err);
                e2 += meas*sqr(err);
                eInf = max(eInf, mag(err));
                n2 += meas*sqr(ex);
                volTot += meas;
                ++integrationPointI;
            }
        }

        e1 = returnReduce(e1, sumOp<scalar>());
        e2 = returnReduce(e2, sumOp<scalar>());
        eInf = returnReduce(eInf, maxOp<scalar>());
        n2 = returnReduce(n2, sumOp<scalar>());
        volTot = returnReduce(volTot, sumOp<scalar>());

        ReconstructionErrorSummary s;
        s.L1 = e1/max(volTot, SMALL);
        s.L2 = std::sqrt(e2/max(volTot, SMALL));
        s.Linf = eInf;
        s.normL2 = std::sqrt(n2/max(volTot, SMALL));
        s.relL2 = 100.0*s.L2/max(s.normL2, SMALL);

        return s;
    }

    // Write the run summary consumed by run_convergence.py: the error norms of
    // every field at the cell centres and at the quadrature points, the fitted
    // quantities, and the timing/memory figures. This file is the interface
    // between the solver and the convergence tables, so its layout is part of
    // the driver's contract.
    void writeSummary
    (
        const Time& runTime,
        const label nSteps,
        const scalar dt,
        const FieldErrorSummary& VmCell,
        const FieldErrorSummary& u1Cell,
        const FieldErrorSummary& u2Cell,
        const ReconstructionErrorSummary& VmHO,
        const ReconstructionErrorSummary& u1HO,
        const ReconstructionErrorSummary& u2HO,
        const volScalarField& rhsVm,
        const volScalarField& Iion,
        const fvMesh& mesh,
        const word& implicitScheme,
        const word& massMatrixType,
        const bool useHighOrderVm,
        const bool useHighOrderIion,
        const scalar stabilisationAlpha,
        const word& memoryOptimization,
        const bool memoryOptimizationEffective,
        const label memoryOptimizationCellThreshold,
        const Switch memoryOptimizationAdaptiveTripletReserve,
        const label memoryOptimizationFluxTripletsPerFace,
        const label memoryOptimizationStabilisedTripletsPerFace,
        const Switch memoryOptimizationTrimHeap,
        const Switch memoryOptimizationCompactMassAssembly,
        const bool compactMassAssembly,
        const label stiffnessTripletsPerFaceReserve,
        const label allocatedIionIntegrationPoints,
        const label lreN,
        const label lreNn,
        const scalar lreK,
        const label lreMaxStencilSize,
        const label lreIionN,
        const label lreIionNn,
        const scalar lreIionK,
        const label lreIionMaxStencilSize,
        const word& linearSolverBackend,
        const word& implicitLinearSolver,
        const word& petscLinearKspType,
        const word& petscLinearPcType,
        const scalar implicitTolerance,
        const label implicitMaxIterations,
        const label maxNonlinearIterations,
        const label minNonlinearIterations,
        const scalar nonlinearTolerance,
        const scalar nonlinearVmTolerance,
        const scalar nonlinearStatesTolerance,
        const Switch nonlinearRequireStatesConvergence,
        const Switch nonlinearAcceptUnconverged,
        const scalar nonlinearIionTolerance,
        const scalar nonlinearRelaxation,
        const label jfnkMaxKrylovIterations,
        const label jfnkMaxRestarts,
        const scalar jfnkLinearTolerance,
        const scalar jfnkEpsilon,
        const word& jfnkLinearSolverBackend,
        const word& jfnkPetscKspType,
        const word& jfnkPetscPcType,
        const label jfnkInitGuessOrder,
        const scalar jfnkInitGuessVmMin,
        const scalar jfnkInitGuessVmMax,
        const Switch jfnkClampODEInput,
        const Switch jfnkLineSearch,
        const label jfnkLineSearchMaxIter,
        const scalar jfnkLineSearchAlphaMin,
        const scalar jfnkArmijoC,
        const word& jfnkPreconditioner,
        const scalar jfnkPreconditionerDropTolerance,
        const label jfnkPreconditionerFillFactor,
        const label jfnkPreconditionerUpdateFrequency,
        const scalar diagonalIionEpsilon,
        const word& stateODESolver,
        const scalar stateODEInitialStep,
        const scalar stateODEAbsTol,
        const scalar stateODERelTol,
        const label stateODEMaxSteps,
        const Switch stateODEUseOpenMP,
        const label stateODENumThreads,
        const std::vector<NonlinearConvergenceRecord>& nonlinearHistory,
        const scalar peakMemoryKB,
        const scalar setupWallTime,
        const scalar timeLoopWallTime,
        const scalar postProcessWallTime,
        const scalar totalWallTime,
        const word& nonlinearMethod
    )
    {
        const label N = estimatedN(mesh);
        const scalar dx = characteristicDx(mesh);
        const word dimName = name(mesh.nGeometricD()) + "D";
        const bool isPicard =
            nonlinearMethod == "Picard" || nonlinearMethod == "picard";
        const bool isJFNK =
            nonlinearMethod == "JFNK" || nonlinearMethod == "jfnk";
        const bool isDiagonalIion =
            nonlinearMethod == "diagonalIion"
         || nonlinearMethod == "diagonal"
         || nonlinearMethod == "localDiagonal";

        const NonlinearConvergenceRecord* finalNonlinear =
            nonlinearHistory.empty() ? nullptr : &nonlinearHistory.back();

        label maxNonlinearIterationsUsed = 0;
        label nonConvergedNonlinearSteps = 0;
        label rolledBackNonlinearSteps = 0;
        scalar maxCoupledResidual = 0.0;
        scalar maxVmResidual = 0.0;
        scalar maxU1Residual = 0.0;
        scalar maxU2Residual = 0.0;
        scalar maxU3Residual = 0.0;
        scalar maxStateResidual = 0.0;
        scalar maxIionResidual = 0.0;

        for (const NonlinearConvergenceRecord& rec : nonlinearHistory)
        {
            maxNonlinearIterationsUsed =
                max(maxNonlinearIterationsUsed, rec.iterations);
            if (!rec.converged)
            {
                ++nonConvergedNonlinearSteps;
            }
            if (rec.rolledBack)
            {
                ++rolledBackNonlinearSteps;
            }
            maxCoupledResidual = max(maxCoupledResidual, rec.coupledResidual);
            maxVmResidual = max(maxVmResidual, rec.VmResidual);
            maxU1Residual = max(maxU1Residual, rec.u1Residual);
            maxU2Residual = max(maxU2Residual, rec.u2Residual);
            maxU3Residual = max(maxU3Residual, rec.u3Residual);
            maxStateResidual = max(maxStateResidual, rec.maxStateResidual);
            maxIionResidual = max(maxIionResidual, rec.IionResidual);
        }

        // Modified for cardiacFoam: the canonical summary goes to the global
        // case directory, written by the master; the other ranks write a
        // throwaway copy into their own processorN directory.
        //
        // Note the shape of this: every rank must still construct an OFstream.
        // OpenFOAM's file handler can be collective, so a naive
        // "if (!Pstream::master()) return;" here deadlocks the run - the master
        // enters the file operation and the others never do. That cost a hang
        // that looked like a solver problem: the time loop completed, every
        // nonlinear step converged, and then the run simply stopped.
        //
        // runTime.path() points at processorN/ on a decomposed run, which is
        // why the master needs rootPath()/globalCaseName() instead: otherwise
        // the post-processing driver finds no file in the case directory.
        const fileName summaryName
        (
            dimName + "_" + name(N) + "_cells_transient.dat"
        );

        const fileName outFile =
            Pstream::master()
          ? runTime.rootPath()/runTime.globalCaseName()/summaryName
          : runTime.path()/summaryName;

        OFstream os(outFile);

        os  << "Manufactured-solution error summary (cell-centred):" << nl
            << "Field     L1-error       L2-error       Linf-error" << nl
            << "Vm      " << VmCell.L1 << "   " << VmCell.L2 << "   "
            << VmCell.Linf << nl
            << "u1      " << u1Cell.L1 << "   " << u1Cell.L2 << "   "
            << u1Cell.Linf << nl
            << "u2      " << u2Cell.L1 << "   " << u2Cell.L2 << "   "
            << u2Cell.Linf << nl
            << "-------------------------------------------------" << nl << nl
            << "MATLAB-like cell-centred error:" << nl
            << "Field     error_cell" << nl
            << "VmC      " << VmCell.errorCell << nl
            << "u1C      " << u1Cell.errorCell << nl
            << "u2C      " << u2Cell.errorCell << nl
            << "-------------------------------------------------" << nl << nl
            << "Gauss-reconstructed error summary:" << nl
            << "Field     L1-error       L2-error       Linf-error       normL2         relL2(%)" << nl
            << "VmG      " << VmHO.L1 << "   " << VmHO.L2 << "   "
            << VmHO.Linf << "   " << VmHO.normL2 << "   "
            << VmHO.relL2 << nl
            << "u1G      " << u1HO.L1 << "   " << u1HO.L2 << "   "
            << u1HO.Linf << "   " << u1HO.normL2 << "   "
            << u1HO.relL2 << nl
            << "u2G      " << u2HO.L1 << "   " << u2HO.L2 << "   "
            << u2HO.Linf << "   " << u2HO.normL2 << "   "
            << u2HO.relL2 << nl
            << "-------------------------------------------------" << nl << nl
            << "RHS summary at final time:" << nl
            << "Field     Linf-RHS" << nl
            << "VmR      " << linfNorm(rhsVm) << nl
            << "IionR    " << linfNorm(Iion) << nl
            << "-------------------------------------------------" << nl << nl
            << "Simulation summary:" << nl
            << "-------------------" << nl
            << "Final time            = " << runTime.value() << nl
            << "Number of cells (N)   = " << N << nl
            << "Dimension             = " << dimName << nl
            << "Grid spacing (dx)     = " << dx << nl
            << "Time step (dt)        = " << dt << nl
            << "Number of steps       = " << nSteps << nl
            << "Implicit scheme       = " << implicitScheme << nl
            << "Mass matrix           = " << massMatrixType << nl
            << "useHighOrder_Vm       = " << (useHighOrderVm ? "true" : "false") << nl
            << "Stabilisation alpha   = " << stabilisationAlpha << nl
            << "memoryOptimization    = " << memoryOptimization << nl
            << "memoryOptimization effective = "
            << (memoryOptimizationEffective ? "true" : "false") << nl
            << "memoryOptimization cellThreshold = "
            << memoryOptimizationCellThreshold << nl
            << "adaptive triplet reserve = "
            << (memoryOptimizationAdaptiveTripletReserve ? "true" : "false") << nl
            << "flux triplets/face reserve = "
            << memoryOptimizationFluxTripletsPerFace << nl
            << "stabilised triplets/face reserve = "
            << memoryOptimizationStabilisedTripletsPerFace << nl
            << "memoryOptimization trimHeap = "
            << (memoryOptimizationTrimHeap ? "true" : "false") << nl
            << "memoryOptimization compactMassAssembly = "
            << (memoryOptimizationCompactMassAssembly ? "true" : "false") << nl
            << "compact mass assembly = "
            << (compactMassAssembly ? "true" : "false") << nl
            << "stiffness triplets/face reserve = "
            << stiffnessTripletsPerFaceReserve << nl
            << "Vm LRE N              = " << lreN << nl
            << "Vm LRE Nn             = " << lreNn << nl
            << "Vm LRE k              = " << lreK << nl
            << "Vm LRE maxStencilSize = " << lreMaxStencilSize << nl
            << "useHighOrder_Iion     = " << (useHighOrderIion ? "true" : "false") << nl
            << "allocated Iion IPs    = " << allocatedIionIntegrationPoints << nl
            << "Iion LRE N            = " << lreIionN << nl
            << "Iion LRE Nn           = " << lreIionNn << nl
            << "Iion LRE k            = " << lreIionK << nl
            << "Iion LRE maxStencilSize = " << lreIionMaxStencilSize << nl
            << "Nonlinear method      = " << nonlinearMethod << nl
            << "Max nonlinear iterations = " << maxNonlinearIterations << nl
            << "Min nonlinear iterations = " << minNonlinearIterations << nl
            << "Linear backend        = " << linearSolverBackend << nl
            << "Linear PDE solver     = " << implicitLinearSolver << nl
            << "PETSc linear KSP      = " << petscLinearKspType << nl
            << "PETSc linear PC       = " << petscLinearPcType << nl
            << "Linear PDE tolerance  = " << implicitTolerance << nl
            << "Linear PDE maxIter    = " << implicitMaxIterations << nl
            << "State ODE solver      = " << stateODESolver << nl
            << "State ODE initial step= " << stateODEInitialStep << nl
            << "State ODE absTol      = " << stateODEAbsTol << nl
            << "State ODE relTol      = " << stateODERelTol << nl
            << "State ODE maxSteps    = " << stateODEMaxSteps << nl
            << "State ODE OpenMP      = "
            << (stateODEUseOpenMP ? "true" : "false") << nl
            << "State ODE threads     = " << stateODENumThreads << nl;

        if (isPicard)
        {
            os  << "Vm Picard tolerance   = " << nonlinearVmTolerance << nl
                << "State Picard tolerance= " << nonlinearStatesTolerance << nl
                << "Picard relaxation     = " << nonlinearRelaxation << nl;
        }
        else if (isJFNK)
        {
            os  << "JFNK coupled tolerance= " << nonlinearTolerance << nl
                << "JFNK Vm tolerance     = " << nonlinearVmTolerance << nl
                << "JFNK Iion tolerance   = " << nonlinearIionTolerance << nl
                << "JFNK require states   = "
                << (nonlinearRequireStatesConvergence ? "true" : "false") << nl;
            if (nonlinearRequireStatesConvergence)
            {
                os  << "JFNK state tolerance  = "
                    << nonlinearStatesTolerance << nl;
            }
            os  << "JFNK accept unconverged = "
                << (nonlinearAcceptUnconverged ? "true" : "false") << nl
                << "JFNK Newton relaxation= " << nonlinearRelaxation << nl
                << "JFNK GMRES m          = " << jfnkMaxKrylovIterations << nl
                << "JFNK GMRES restarts   = " << jfnkMaxRestarts << nl
                << "JFNK GMRES tolerance  = " << jfnkLinearTolerance << nl
                << "JFNK epsilon          = " << jfnkEpsilon << nl
                << "JFNK backend          = " << jfnkLinearSolverBackend << nl
                << "JFNK PETSc KSP        = " << jfnkPetscKspType << nl
                << "JFNK PETSc PC         = " << jfnkPetscPcType << nl
                << "JFNK initGuessOrder   = " << jfnkInitGuessOrder << nl
                << "JFNK initGuessVmMin   = " << jfnkInitGuessVmMin << nl
                << "JFNK initGuessVmMax   = " << jfnkInitGuessVmMax << nl
                << "JFNK clampODEInput    = "
                << (jfnkClampODEInput ? "true" : "false") << nl
                << "JFNK lineSearch       = "
                << (jfnkLineSearch ? "true" : "false") << nl;
            if (jfnkLineSearch)
            {
                os  << "JFNK lineSearch maxIt = " << jfnkLineSearchMaxIter << nl
                    << "JFNK lineSearch alphaMin = "
                    << jfnkLineSearchAlphaMin << nl
                    << "JFNK Armijo c         = " << jfnkArmijoC << nl;
            }
            os  << "JFNK preconditioner   = " << jfnkPreconditioner << nl
                << "JFNK PC dropTol       = "
                << jfnkPreconditionerDropTolerance << nl
                << "JFNK PC fillFactor    = "
                << jfnkPreconditionerFillFactor << nl
                << "JFNK PC updateFreq    = "
                << jfnkPreconditionerUpdateFrequency << nl;
        }
        else if (isDiagonalIion)
        {
            os  << "Diagonal coupled tolerance = " << nonlinearTolerance << nl
                << "Diagonal Vm tolerance = " << nonlinearVmTolerance << nl
                << "Diagonal Iion tolerance = " << nonlinearIionTolerance << nl
                << "Diagonal require states = "
                << (nonlinearRequireStatesConvergence ? "true" : "false") << nl;
            if (nonlinearRequireStatesConvergence)
            {
                os  << "Diagonal state tolerance = "
                    << nonlinearStatesTolerance << nl;
            }
            os  << "Diagonal accept unconverged = "
                << (nonlinearAcceptUnconverged ? "true" : "false") << nl
                << "Diagonal relaxation  = " << nonlinearRelaxation << nl
                << "Diagonal Iion epsilon= " << diagonalIionEpsilon << nl;
        }

        os  << "-------------------" << nl << nl
            << nonlinearMethod << " nonlinear convergence summary:" << nl
            << "final it/max          = "
            << (finalNonlinear ? finalNonlinear->iterations : 0)
            << "/" << maxNonlinearIterations << nl
            << "final converged       = "
            << (finalNonlinear && finalNonlinear->converged ? "true" : "false") << nl
            << "final rolledBack      = "
            << (finalNonlinear && finalNonlinear->rolledBack ? "true" : "false") << nl;

        if (!isPicard)
        {
            os  << "final coupled_relL2   = "
                << (finalNonlinear ? finalNonlinear->coupledResidual : GREAT) << nl;
        }

        os  << "final Vm_relL2        = "
            << (finalNonlinear ? finalNonlinear->VmResidual : GREAT) << nl
            << "final u1_relL2        = "
            << (finalNonlinear ? finalNonlinear->u1Residual : GREAT) << nl
            << "final u2_relL2        = "
            << (finalNonlinear ? finalNonlinear->u2Residual : GREAT) << nl
            << "final u3_relL2        = "
            << (finalNonlinear ? finalNonlinear->u3Residual : GREAT) << nl
            << "final maxState_relL2  = "
            << (finalNonlinear ? finalNonlinear->maxStateResidual : GREAT) << nl;

        if (!isPicard)
        {
            os  << "final Iion_relL2      = "
                << (finalNonlinear ? finalNonlinear->IionResidual : GREAT) << nl;
        }

        os  << "max it used           = " << maxNonlinearIterationsUsed << nl
            << "non-converged steps   = " << nonConvergedNonlinearSteps << nl
            << "rolled-back steps     = " << rolledBackNonlinearSteps << nl;

        if (!isPicard)
        {
            os  << "max coupled_relL2     = " << maxCoupledResidual << nl;
        }

        os  << "max Vm_relL2          = " << maxVmResidual << nl
            << "max u1_relL2          = " << maxU1Residual << nl
            << "max u2_relL2          = " << maxU2Residual << nl
            << "max u3_relL2          = " << maxU3Residual << nl
            << "maxState_relL2        = " << maxStateResidual << nl;

        if (!isPicard)
        {
            os  << "max Iion_relL2        = " << maxIionResidual << nl;
        }

        os  << "-------------------------------------------------" << nl << nl
            << "Computational resources:" << nl
            << "peakRSS_kB            = " << peakMemoryKB << nl
            << "peakRSS_MB            = " << peakMemoryKB/1024.0 << nl
            << "-------------------" << nl << nl
            << "Wall-clock times [s]:" << nl
            << "setup                 = " << setupWallTime << nl
            << "timeLoop              = " << timeLoopWallTime << nl
            << "postProcess           = " << postProcessWallTime << nl
            << "total                 = " << totalWallTime << nl
            << "-------------------" << nl;

        Info<< "Wrote summary to " << outFile << nl
            << "Vm error: L1 = " << VmCell.L1
            << ", L2 = " << VmCell.L2
            << ", Linf = " << VmCell.Linf << nl
            << "Vm Gauss error: L1 = " << VmHO.L1
            << ", L2 = " << VmHO.L2
            << ", Linf = " << VmHO.Linf
            << ", Relative = " << VmHO.relL2 << "%" << nl
            << "Final-time RHS Linf = " << linfNorm(rhsVm) << nl
            << "Timing [s]: setup = " << setupWallTime
            << ", loop = " << timeLoopWallTime
            << ", post = " << postProcessWallTime
            << ", total = " << totalWallTime << endl;
    }
}

int main(int argc, char* argv[])
{
#ifdef __GLIBC__
    // Tighten the glibc allocator before any large allocation. The LRE
    // construction runs ~2.3M QR factorisations on a 3D unstructured tet
    // mesh with N=40 / p3; without these tunables glibc grows the heap by
    // several GB of unreclaimed fragments and the run OOMs on a 16 GB box.
    //   - M_ARENA_MAX=2: prevents glibc from creating ~64 arenas (one set
    //     per ODE OpenMP thread). LRE setup is serial; multiple arenas
    //     just inflate the RSS.
    //   - M_MMAP_THRESHOLD=64 KB: routes mid-sized allocations (Eigen
    //     temporaries during QR) through mmap so they are returned to
    //     the kernel individually when freed.
    //   - M_TRIM_THRESHOLD=64 KB: lets glibc release sbrk-heap padding
    //     back to the OS more aggressively.
    mallopt(M_ARENA_MAX, 2);
    mallopt(M_MMAP_THRESHOLD, 64*1024);
    mallopt(M_TRIM_THRESHOLD, 64*1024);
#endif

    #include "setRootCaseLists.H"
    PetscSession petscSession(argc, argv);
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    const auto tStartTotal = std::chrono::steady_clock::now();
    const auto tStartSetup = tStartTotal;

    const label dim = mesh.nGeometricD();
    const scalar dt = runTime.deltaTValue();
    const scalar beta = computeBeta(conductivity, dim);
    const scalar chiVal = chi.value();
    const scalar CmVal = Cm.value();
    const scalar chiCmVal = chiVal*CmVal;
    const scalar lapScale = 1.0/max(chiCmVal, SMALL);
    const scalar theta = thetaFromScheme(implicitScheme);
    const word massMatrixMode = normalizedMassMatrixType(massMatrixType);
    const std::string memoryOptimizationMode = lowerWord(memoryOptimization);

    if
    (
        memoryOptimizationMode != "auto"
     && memoryOptimizationMode != "on"
     && memoryOptimizationMode != "off"
    )
    {
        FatalErrorInFunction
            << "Unknown memoryOptimization '" << memoryOptimization << "'. "
            << "Valid options are auto, on and off."
            << abort(FatalError);
    }

    const bool memoryOptimizationEffective =
        memoryOptimizationMode == "on"
     || (
            memoryOptimizationMode == "auto"
         && useHighOrder_Vm
         && dim == 3
         && mesh.nCells() >= memoryOptimizationCellThreshold
        );

    const bool useAdaptiveTripletReserve =
        memoryOptimizationEffective
     && memoryOptimizationAdaptiveTripletReserve;

    const label stiffnessTripletsPerFaceReserve =
            useAdaptiveTripletReserve
        ? (
                stabilisationAlpha > SMALL
              ? memoryOptimizationStabilisedTripletsPerFace
              : memoryOptimizationFluxTripletsPerFace
          )
        : label(800);

    const bool trimHeapAfterLargeSetup =
        memoryOptimizationEffective && memoryOptimizationTrimHeap;

    const bool compactMassAssembly =
        memoryOptimizationEffective && memoryOptimizationCompactMassAssembly;

    const label allocatedIionIntegrationPoints =
            useHighOrder_Iion && dim > 1
        ? totalIionIntegrationPoints
        : label(0);

#ifdef _OPENMP
    if (stateODENumThreads > 0)
    {
        omp_set_num_threads(stateODENumThreads);
    }
#endif

    scalarField VmOldIntegrationPoints(allocatedIionIntegrationPoints, 0.0);
    scalarField VmGuessIntegrationPoints(allocatedIionIntegrationPoints, 0.0);
    scalarField sourceIntegrationPoints(allocatedIionIntegrationPoints, 0.0);
    scalarField IionIntegrationPoints(allocatedIionIntegrationPoints, 0.0);
    scalarField u1IntegrationPoints(allocatedIionIntegrationPoints, 0.0);
    scalarField u2IntegrationPoints(allocatedIionIntegrationPoints, 0.0);
    scalarField u3IntegrationPoints(allocatedIionIntegrationPoints, 0.0);

    Info<< "Running high-order manufactured FDA implicit PETSc solver" << nl
        << "Dimension = " << dim << nl
        << "beta = " << beta << nl
        << "dt = " << dt << nl
        << "implicitScheme = " << implicitScheme << nl
        << "massMatrix = " << massMatrixMode << nl
        << "linearSolverBackend = " << linearSolverBackend << nl
        << "PETSc linear KSP = " << petscLinearKspType
        << ", PC = " << petscLinearPcType << nl
        << "nonlinearMethod = " << nonlinearMethod << nl
        << "JFNK linear backend = " << jfnkLinearSolverBackend << nl
        << "JFNK PETSc KSP = " << jfnkPetscKspType
        << ", PC = " << jfnkPetscPcType << nl
        << "stateODESolver = " << stateODESolver << nl
        << "stateODEUseOpenMP = " << stateODEUseOpenMP
        << ", stateODENumThreads = " << stateODENumThreads << nl
        << "nonlinearRelaxation = " << nonlinearRelaxation << nl
        << "useHighOrder_Vm = " << useHighOrder_Vm << nl
        << "useHighOrder_Iion = " << useHighOrder_Iion << nl
        << "stateIntegrationMode = " << stateIntegrationMode << nl
        << "stabilisationAlpha = " << stabilisationAlpha << nl
        << "memoryOptimization = " << memoryOptimization
        << " (effective: "
        << (memoryOptimizationEffective ? "true" : "false") << ")" << nl
        << "stiffnessTripletsPerFaceReserve = "
        << stiffnessTripletsPerFaceReserve << nl
        << "compactMassAssembly = "
        << (compactMassAssembly ? "true" : "false") << nl
        << "jfnkPreconditioner = " << jfnkPreconditioner << nl
        << "jfnkPreconditionerUpdateFrequency = "
        << jfnkPreconditionerUpdateFrequency << nl
        << "Iion integration points = " << totalIionIntegrationPoints
        << " (allocated: " << allocatedIionIntegrationPoints << ")"
        << endl;

    if (reconstructStatesFromCellCentres && useHighOrder_Iion && dim > 1)
    {
        Info<< "LREInterp_states order = " << LREInterp_statesPtr().order()
            << " (states evolved at cell centres, reconstructed at Gauss points)"
            << endl;
    }

    const bool usePicard =
        nonlinearMethod == "Picard"
     || nonlinearMethod == "picard";

    const bool useDiagonalIion =
        nonlinearMethod == "diagonalIion"
     || nonlinearMethod == "diagonal"
     || nonlinearMethod == "localDiagonal";

    const bool useJFNK =
        nonlinearMethod == "JFNK"
     || nonlinearMethod == "jfnk";

    if (!usePicard && !useDiagonalIion && !useJFNK)
    {
        FatalErrorInFunction
            << "Unknown nonlinearMethod '" << nonlinearMethod << "'. "
            << "Valid options are Picard, JFNK and diagonalIion."
            << abort(FatalError);
    }

    const bool useJfnkPreconditioner =
        jfnkPreconditioner != "none"
     && jfnkPreconditioner != "None"
     && jfnkPreconditioner != "off"
     && jfnkPreconditioner != "false";

    const bool useJfnkDiagonalIionPreconditioner =
        jfnkPreconditioner == "diagonalIion"
     || jfnkPreconditioner == "diagonal"
     || jfnkPreconditioner == "localDiagonal";

    const bool useJfnkDiffusionPreconditioner =
        jfnkPreconditioner == "diffusion"
     || jfnkPreconditioner == "linear"
     || jfnkPreconditioner == "AImplicit";

    if
    (
        useJfnkPreconditioner
     && !useJfnkDiagonalIionPreconditioner
     && !useJfnkDiffusionPreconditioner
    )
    {
        FatalErrorInFunction
            << "Unknown jfnkPreconditioner '" << jfnkPreconditioner << "'. "
            << "Valid options are none, diffusion and diagonalIion."
            << abort(FatalError);
    }

    fillExactFields(VmExact, u1Exact, u2Exact, runTime.value(), dim);

    Vm.primitiveFieldRef() = VmExact.primitiveField();
    u1.primitiveFieldRef() = u1Exact.primitiveField();
    u2.primitiveFieldRef() = u2Exact.primitiveField();
    u3 = dimensionedScalar("zero", dimless, 0.0);

    applyExactVmBoundaryValues(Vm, runTime.value(), dim);
    updateStateBoundaryValues(u1, u2, u3, runTime.value(), dim);
    u1.correctBoundaryConditions();
    u2.correctBoundaryConditions();
    u3.correctBoundaryConditions();

    if (useHighOrder_Iion && dim > 1)
    {
        initialiseIionIntegrationPointStates
        (
            LREInterp_Iion,
            runTime.value(),
            dim,
            u1IntegrationPoints,
            u2IntegrationPoints,
            u3IntegrationPoints
        );
    }

    // Added for cardiacFoam: publish this rank's global row offset for the
    // linear-solver wrappers. The matrices carry local rows and global columns,
    // so every PETSc insertion has to shift the row index by this amount. Zero
    // in serial, which is why the serial path is unaffected.
    gRowStart = LREInterp_Vm.globalCells().localStart();

    Info<< "Distributed linear algebra: rows ["
        << gRowStart << ", " << gRowStart + mesh.nCells() << ") of "
        << LREInterp_Vm.globalCells().totalSize() << " global cells" << endl;

    SpMat M;
    SpMat K;
    SpMat AImplicit;
    SpMat BImplicit;

    if (profileTimings)
    {
        logMemoryCheckpoint("before mass assembly");
    }

    if (!useHighOrder_Vm || massMatrixMode == "lumped")
    {
        assembleDiagonalMassMatrix(mesh, LREInterp_Vm.globalCells(), 1.0, M);
    }
    else
    {
        assembleConsistentMassMatrixHO
        (
            mesh,
            LREInterp_Vm,
            1.0,
            compactMassAssembly,
            M
        );
    }

    if (profileTimings)
    {
        logMemoryCheckpoint("after mass assembly");
    }

#ifdef __GLIBC__
    // Reclaim the heap fragmentation left by the ~290k per-cell QR
    // factorisations triggered when LREInterp_Vm populated its cell
    // gradient / Hessian / 3rd-derivative coefficient tables.
    if (trimHeapAfterLargeSetup)
    {
        malloc_trim(0);
    }
#endif

    if (profileTimings)
    {
        logMemoryCheckpoint("before stiffness assembly");
    }

    if (useHighOrder_Vm)
    {
        assembleHighOrderStiffnessMatrix
        (
            mesh,
            conductivity,
            LREInterp_Vm,
            stabilisationAlpha,
            stiffnessTripletsPerFaceReserve,
            K
        );
    }
    else
    {
        assembleStandardOrthogonalStiffnessMatrix
        (
            mesh,
            conductivity,
            LREInterp_Vm,
            stabilisationAlpha,
            K
        );
    }

    if (profileTimings)
    {
        logMemoryCheckpoint("after stiffness assembly");
    }

    // Added for cardiacFoam: consistency check on the assembled diffusion
    // operator. A constant field has zero gradient, so with no-flux boundaries
    // every row of K must sum to zero. This caught a silent failure during the
    // LRE -> movingLeastSquares port: CompactListList::operator[] returns a
    // UList VIEW by value, and binding it to a const labelList& compiled but
    // produced an empty stencil, so K came out with no entries at all. The
    // symptom was a fixed error that did not reduce under refinement, which is
    // easy to mistake for a discretisation problem. Cheap (one pass over the
    // non-zeros) and it doubles as the conservation probe in parallel.
    {
        scalar maxRowSum = 0.0;
        for (label r = 0; r < K.outerSize(); ++r)
        {
            scalar rowSum = 0.0;
            for (SpMat::InnerIterator it(K, r); it; ++it)
            {
                rowSum += it.value();
            }
            maxRowSum = max(maxRowSum, mag(rowSum));
        }
        reduce(maxRowSum, maxOp<scalar>());

        const label nnzGlobal = returnReduce(label(K.nonZeros()), sumOp<label>());

        Info<< "Stiffness matrix: " << nnzGlobal << " non-zeros, "
            << "max |row sum| = " << maxRowSum
            << " (must be ~0: K must annihilate constants)" << endl;

        if (nnzGlobal == 0)
        {
            FatalErrorInFunction
                << "The stiffness matrix is empty. The diffusion operator would"
                << " be absent from the system." << abort(FatalError);
        }
    }

    // Added for cardiacFoam: with CF_DUMP_STENCILS set, write the MLS stencil
    // MEMBERSHIP - one line per cell, "globalCell : global columns" - so that a
    // serial and a parallel build can be diffed cell by cell. This is one level
    // below CF_DUMP_K: it says whether a differing operator row comes from
    // different stencil members or from different coefficients over the same
    // members.
    if (getenv("CF_DUMP_STENCILS"))
    {
        const globalIndex& gcD = LREInterp_Vm.globalCells();
        const CompactListList<label>& cellSt = LREInterp_Vm.globalCellStencils();

        OFstream os
        (
            "stencil_dump_rank" + Foam::name(Pstream::myProcNo()) + ".dat"
        );
        // Cell centres alongside, at full precision: the selection is a
        // nearest-N by distance, so a centroid that differs in the last bits
        // between a serial and a decomposed mesh is enough to reorder a tie.
        OFstream osc
        (
            "cellCentre_dump_rank" + Foam::name(Pstream::myProcNo()) + ".dat"
        );
        osc.precision(17);
        forAll(mesh.C(), cellI)
        {
            osc << gcD.toGlobal(cellI) << ' ' << mesh.C()[cellI].x() << ' '
                << mesh.C()[cellI].y() << ' ' << mesh.C()[cellI].z() << nl;
        }

        forAll(mesh.C(), cellI)
        {
            os  << gcD.toGlobal(cellI) << " :";
            const UList<label> st = cellSt[cellI];
            forAll(st, i)
            {
                os  << ' ' << st[i];
            }
            os  << nl;
        }

        // Face stencils, keyed by LOCAL face index: the reader maps them with
        // processorN/constant/polyMesh/faceProcAddressing.
        const CompactListList<label>& faceSt = LREInterp_Vm.globalFaceStencils();

        OFstream osf
        (
            "faceStencil_dump_rank" + Foam::name(Pstream::myProcNo()) + ".dat"
        );
        forAll(faceSt, faceI)
        {
            osf << faceI << " :";
            const UList<label> st = faceSt[faceI];
            forAll(st, i)
            {
                osf << ' ' << st[i];
            }
            osf << nl;
        }
    }

    // Added for cardiacFoam: with CF_DUMP_K set in the environment, write K as
    // (global row, global column, value) triples, one file per rank.
    //
    // This is the check that a row-sum or nnz fingerprint cannot make: mapping
    // the parallel indices back through processorN/constant/polyMesh/
    // cellProcAddressing gives an entry-by-entry comparison against the serial
    // assembly, which is what separates a genuine discretisation difference
    // from partition-dependent preconditioning. Off by default, no cost.
    if (getenv("CF_DUMP_K"))
    {
        OFstream os
        (
            "K_dump_rank" + Foam::name(Pstream::myProcNo()) + ".dat"
        );
        os.precision(17);
        for (label r = 0; r < K.outerSize(); ++r)
        {
            for (SpMat::InnerIterator it(K, r); it; ++it)
            {
                os  << (gRowStart + r) << ' ' << it.col() << ' '
                    << it.value() << nl;
            }
        }
    }

#ifdef __GLIBC__
    // Reclaim the heap fragmentation from the ~2.3M per-(face, GP) QR
    // factorisations during LREInterp_Vm.QRGradFaceGPCoeffs() population.
    if (trimHeapAfterLargeSetup)
    {
        malloc_trim(0);
    }
#endif

    AImplicit = M;
    AImplicit *= (1.0/dt);
    AImplicit -= (theta*lapScale)*K;

    // Move M into BImplicit instead of copying: M is not referenced again
    // after this block, so we can hand its storage over and free the
    // ~280 MB it occupied as a separate matrix.
    BImplicit = std::move(M);
    BImplicit *= (1.0/dt);

    if (theta < 1.0 - SMALL)
    {
        BImplicit += ((1.0 - theta)*lapScale)*K;
    }

    // Added for cardiacFoam: PETSc mirrors of the two implicit operators, used
    // for every matrix-vector product. Eigen cannot do A*x any more, because A
    // has local rows and global columns while x only holds this rank's block;
    // MatMult performs the halo exchange internally. Built once here, since
    // AImplicit and BImplicit are constant while dt is (dt only changes on a
    // final clipped step, which this solver does not take).
    DistributedMatVec applyA;
    DistributedMatVec applyB;
    DistributedMatVec applyK;
    applyA.reset(AImplicit);
    applyB.reset(BImplicit);
    applyK.reset(K);

    if (profileTimings)
    {
        logMemoryCheckpoint("after implicit matrix assembly");
    }

#ifdef __GLIBC__
    if (trimHeapAfterLargeSetup)
    {
        malloc_trim(0);
    }
#endif

    PetscKspMatrixSolver cachedPicardPetscSolver;
    const bool useCachedPicardPetscSolver =
        usePicard && usesPetscBackend(linearSolverBackend);

    if (useCachedPicardPetscSolver)
    {
        cachedPicardPetscSolver.reset
        (
            AImplicit,
            petscLinearKspType,
            petscLinearPcType,
            implicitTolerance,
            implicitMaxIterations,
            petscLinearRestart,
            petscLinearOptionsPrefix,
            petscUseOptions
        );
    }

    // Persistent JFNK preconditioner state, kept across the time loop so that
    // the PETSc Mat, KSP and ILUT symbolic factorisation are reused. Memory
    // peak per PC update is dominated by the freshly allocated ILUT factor
    // (~fillFactor * nnz(AImplicit)); reusing the cached factor avoids
    // repeated alloc/free that fragments the glibc heap on long runs.
    Eigen::IncompleteLUT<scalar> jfnkIlu;
    PetscKspMatrixSolver jfnkPetscPcSolver;
    SpMat jfnkP;
    bool jfnkPMatInitialised = false;

    // Persistent JFNK shell-Krylov solver (PETSc MatShell + KSP). The Krylov
    // subspace and KSP structures are allocated once and reused across every
    // Newton iteration, avoiding the ~25-40 MB alloc/free cycle that was
    // happening once per Newton step in the old solvePetscShellSystem call.
    // The mat-vec / apply-PC callback payloads are refreshed per solve().
    PetscShellKspSolver jfnkPetscShellSolver;

    autoPtr<OFstream> nonlinearResidualFilePtr;
    // Modified for cardiacFoam: master only, global case directory. See the
    // transient summary above.
    if (writeNonlinearResiduals && Pstream::master())
    {
        const fileName residualDir
        (
            runTime.rootPath()/runTime.globalCaseName()
          / "postProcessing"/"highOrderManufacturedFDAImplicitPETSc"
        );
        mkDir(residualDir);
        nonlinearResidualFilePtr.reset
        (
            new OFstream(residualDir/"nonlinearResiduals.dat")
        );
        nonlinearResidualFilePtr()
            << "# time_s step nonlinearMethod stateODESolver iter linearIterations linearError "
            << "newtonResidual Vm_relL2 u1_relL2 u2_relL2 u3_relL2 Iion_relL2 "
            << "lineSearchIters converged"
            << nl;
    }

    auto updateNumericalLaplacian =
    [&](const scalar evalTime, const bool updateFluxField)
    {
        if (useHighOrder_Vm && updateFluxField)
        {
            computeHighOrderLaplacian
            (
                Vm,
                conductivity,
                LREInterp_Vm,
                fluxVm_HO,
                lapVm
            );
        }

        EigVec bcNow =
                useHighOrder_Vm
            ? assembleHighOrderBoundaryVector
              (
                  mesh,
                  conductivity,
                  LREInterp_Vm,
                  stabilisationAlpha,
                  evalTime,
                  dim
              )
            : assembleStandardOrthogonalBoundaryVector
              (
                  mesh,
                  conductivity,
                  LREInterp_Vm,
                  stabilisationAlpha,
                  evalTime,
                  dim
              );

        // Modified for cardiacFoam: distributed product, see DistributedMatVec.
        EigVec lapNow = applyK(fieldToEigVec(Vm)) + bcNow;
        eigVecToField(lapNow, lapVm);
        lapVm.correctBoundaryConditions();
    };

    const auto tEndSetup = std::chrono::steady_clock::now();
    const auto tStartLoop = tEndSetup;

    label nSteps = 0;
    scalarField VmTwoStepsAgo(mesh.nCells(), 0.0);
    scalar dtPrev = 0.0;
    bool hasPrevStep = false;
    std::vector<NonlinearConvergenceRecord> nonlinearHistory;

    scalar nonlinearEvalWallTime = 0.0;
    scalar gmresWallTime = 0.0;
    scalar sparseLinearSolveWallTime = 0.0;
    scalar preconditionerSetupWallTime = 0.0;
    scalar preconditionerApplyWallTime = 0.0;
    label nonlinearEvalCalls = 0;
    label gmresCalls = 0;
    label sparseLinearSolveCalls = 0;
    label preconditionerSetups = 0;
    label preconditionerApplications = 0;

    while (runTime.value() < runTime.endTime().value() - SMALL)
    {
        const scalar t = runTime.value();

        volScalarField VmOld
        (
            IOobject
            (
                "VmOld",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Vm
        );

        volScalarField u1Old
        (
            IOobject
            (
                "u1Old",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            u1
        );

        volScalarField u2Old
        (
            IOobject
            (
                "u2Old",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            u2
        );

        volScalarField u3Old
        (
            IOobject
            (
                "u3Old",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            u3
        );

        scalarField u1OldIntegrationPoints(u1IntegrationPoints);
        scalarField u2OldIntegrationPoints(u2IntegrationPoints);
        scalarField u3OldIntegrationPoints(u3IntegrationPoints);
        scalarField u1GuessIntegrationPoints(u1OldIntegrationPoints);
        scalarField u2GuessIntegrationPoints(u2OldIntegrationPoints);
        scalarField u3GuessIntegrationPoints(u3OldIntegrationPoints);

        if (useHighOrder_Iion && dim > 1)
        {
            reconstructVmAtIionIntegrationPoints
            (
                VmOld,
                useHighOrder_Vm,
                LREInterp_Vm,
                LREInterp_Iion,
                VmOldIntegrationPoints
            );

            computeVmSourceFromIntegrationPoints
            (
                VmOldIntegrationPoints,
                u1OldIntegrationPoints,
                u2OldIntegrationPoints,
                u3OldIntegrationPoints,
                beta,
                chiVal,
                CmVal,
                sourceIntegrationPoints,
                stateODEUseOpenMP,
                stateODEOpenMPThreshold
            );

            averageIntegrationPointFieldToCells
            (
                sourceIntegrationPoints,
                LREInterp_Iion,
                sourceVm
            );
        }
        else
        {
            computeCellCentredVmSource
            (
                VmOld,
                u1Old,
                u2Old,
                u3Old,
                beta,
                chiVal,
                CmVal,
                sourceVm,
                stateODEUseOpenMP,
                stateODEOpenMPThreshold
            );
        }

        const EigVec Vn = fieldToEigVec(VmOld);
        const EigVec sourceN = sourceToEigVec(sourceVm);

        EigVec bcN =
                useHighOrder_Vm
            ? assembleHighOrderBoundaryVector
              (
                  mesh,
                  conductivity,
                  LREInterp_Vm,
                  stabilisationAlpha,
                  t,
                  dim
              )
            : assembleStandardOrthogonalBoundaryVector
              (
                  mesh,
                  conductivity,
                  LREInterp_Vm,
                  stabilisationAlpha,
                  t,
                  dim
              );

            EigVec bcNp1 =
                useHighOrder_Vm
            ? assembleHighOrderBoundaryVector
              (
                  mesh,
                  conductivity,
                  LREInterp_Vm,
                  stabilisationAlpha,
                  t + dt,
                  dim
              )
            : assembleStandardOrthogonalBoundaryVector
              (
                  mesh,
                  conductivity,
                  LREInterp_Vm,
                  stabilisationAlpha,
                  t + dt,
                  dim
              );


        bcN *= lapScale;
        bcNp1 *= lapScale;

        volScalarField VmGuess
        (
            IOobject
            (
                "VmGuess",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            VmOld
        );

        volScalarField u1Guess
        (
            IOobject
            (
                "u1Guess",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            u1Old
        );

        volScalarField u2Guess
        (
            IOobject
            (
                "u2Guess",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            u2Old
        );

        volScalarField u3Guess
        (
            IOobject
            (
                "u3Guess",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            u3Old
        );

        volScalarField IionGuess
        (
            IOobject
            (
                "IionGuess",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Iion
        );

        if (useJFNK && jfnkInitGuessOrder >= 1 && hasPrevStep && dtPrev > SMALL)
        {
            scalarField& Vg = VmGuess.primitiveFieldRef();
            const scalarField& VnInternal = VmOld.primitiveField();
            const scalar ratio = dt/dtPrev;

            forAll(Vg, cellI)
            {
                scalar v =
                    VnInternal[cellI]
                  + ratio*(VnInternal[cellI] - VmTwoStepsAgo[cellI]);

                v = min(max(v, jfnkInitGuessVmMin), jfnkInitGuessVmMax);
                Vg[cellI] = v;
            }
        }
        applyExactVmBoundaryValues(VmGuess, t + dt, dim);

        VmTwoStepsAgo = VmOld.primitiveField();
        dtPrev = dt;
        hasPrevStep = true;

        volScalarField sourceDerivative
        (
            IOobject
            (
                "sourceDerivative",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("sourceDerivative", dimless/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );

        auto evaluateNonlinearFields =
        [&]
        (
            const volScalarField& VmCandidate,
            volScalarField& u1Candidate,
            volScalarField& u2Candidate,
            volScalarField& u3Candidate,
            scalarField& VmCandidateIntegrationPoints,
            scalarField& u1CandidateIntegrationPoints,
            scalarField& u2CandidateIntegrationPoints,
            scalarField& u3CandidateIntegrationPoints,
            volScalarField& sourceCandidate,
            volScalarField& IionCandidate
        )
        {
            const auto evalStart = std::chrono::steady_clock::now();

            if (useHighOrder_Iion && dim > 1)
            {
                reconstructVmAtIionIntegrationPoints
                (
                    VmCandidate,
                    useHighOrder_Vm,
                    LREInterp_Vm,
                    LREInterp_Iion,
                    VmCandidateIntegrationPoints
                );

                if (jfnkClampODEInput)
                {
                    clampPhysical
                    (
                        VmCandidateIntegrationPoints,
                        jfnkInitGuessVmMin,
                        jfnkInitGuessVmMax
                    );
                }

                if (reconstructStatesFromCellCentres)
                {
                    // ODE once per cell (cell-centred), then reconstruct the
                    // states at the Iion Gauss points with LREInterp_states.
                    updateStateFieldsODE
                    (
                        VmOld,
                        VmCandidate,
                        u1Old,
                        u2Old,
                        u3Old,
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        dt,
                        stateODESolver,
                        stateODEInitialStep,
                        stateODEAbsTol,
                        stateODERelTol,
                        stateODEMaxSteps,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );

                    updateStateBoundaryValues
                    (
                        u1Candidate, u2Candidate, u3Candidate, t + dt, dim
                    );
                    u1Candidate.correctBoundaryConditions();
                    u2Candidate.correctBoundaryConditions();
                    u3Candidate.correctBoundaryConditions();

                    reconstructStatesAtIionIntegrationPoints
                    (
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        LREInterp_statesPtr(),
                        LREInterp_Iion,
                        u1CandidateIntegrationPoints,
                        u2CandidateIntegrationPoints,
                        u3CandidateIntegrationPoints
                    );
                }
                else
                {
                    // Legacy: integrate the ODE independently at every Gauss
                    // point, then average the states back to the cell centres.
                    updateIntegrationPointStatesODE
                    (
                        VmOldIntegrationPoints,
                        VmCandidateIntegrationPoints,
                        u1OldIntegrationPoints,
                        u2OldIntegrationPoints,
                        u3OldIntegrationPoints,
                        u1CandidateIntegrationPoints,
                        u2CandidateIntegrationPoints,
                        u3CandidateIntegrationPoints,
                        dt,
                        stateODESolver,
                        stateODEInitialStep,
                        stateODEAbsTol,
                        stateODERelTol,
                        stateODEMaxSteps,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );

                    averageIntegrationPointFieldToCells
                    (
                        u1CandidateIntegrationPoints,
                        LREInterp_Iion,
                        u1Candidate
                    );
                    averageIntegrationPointFieldToCells
                    (
                        u2CandidateIntegrationPoints,
                        LREInterp_Iion,
                        u2Candidate
                    );
                    averageIntegrationPointFieldToCells
                    (
                        u3CandidateIntegrationPoints,
                        LREInterp_Iion,
                        u3Candidate
                    );
                }

                computeVmSourceFromIntegrationPoints
                (
                    VmCandidateIntegrationPoints,
                    u1CandidateIntegrationPoints,
                    u2CandidateIntegrationPoints,
                    u3CandidateIntegrationPoints,
                    beta,
                    chiVal,
                    CmVal,
                    sourceIntegrationPoints,
                    stateODEUseOpenMP,
                    stateODEOpenMPThreshold
                );
                averageIntegrationPointFieldToCells
                (
                    sourceIntegrationPoints,
                    LREInterp_Iion,
                    sourceCandidate
                );

                computeIionFromIntegrationPoints
                (
                    VmCandidateIntegrationPoints,
                    u1CandidateIntegrationPoints,
                    u2CandidateIntegrationPoints,
                    u3CandidateIntegrationPoints,
                    beta,
                    chiVal,
                    CmVal,
                    IionIntegrationPoints,
                    stateODEUseOpenMP,
                    stateODEOpenMPThreshold
                );
                averageIntegrationPointFieldToCells
                (
                    IionIntegrationPoints,
                    LREInterp_Iion,
                    IionCandidate
                );
            }
            else
            {
                if (jfnkClampODEInput)
                {
                    scalarField VmClampedInternal(VmCandidate.internalField());
                    clampPhysical
                    (
                        VmClampedInternal,
                        jfnkInitGuessVmMin,
                        jfnkInitGuessVmMax
                    );

                    volScalarField VmClamped
                    (
                        IOobject
                        (
                            "VmClamped",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        VmCandidate
                    );
                    VmClamped.primitiveFieldRef() = VmClampedInternal;
                    VmClamped.correctBoundaryConditions();

                    updateStateFieldsODE
                    (
                        VmOld,
                        VmClamped,
                        u1Old,
                        u2Old,
                        u3Old,
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        dt,
                        stateODESolver,
                        stateODEInitialStep,
                        stateODEAbsTol,
                        stateODERelTol,
                        stateODEMaxSteps,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );

                    computeCellCentredVmSource
                    (
                        VmClamped,
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        beta,
                        chiVal,
                        CmVal,
                        sourceCandidate,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );

                    computeCellCentredIion
                    (
                        VmClamped,
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        beta,
                        chiVal,
                        CmVal,
                        IionCandidate,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );
                }
                else
                {
                    updateStateFieldsODE
                    (
                        VmOld,
                        VmCandidate,
                        u1Old,
                        u2Old,
                        u3Old,
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        dt,
                        stateODESolver,
                        stateODEInitialStep,
                        stateODEAbsTol,
                        stateODERelTol,
                        stateODEMaxSteps,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );

                    computeCellCentredVmSource
                    (
                        VmCandidate,
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        beta,
                        chiVal,
                        CmVal,
                        sourceCandidate,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );

                    computeCellCentredIion
                    (
                        VmCandidate,
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        beta,
                        chiVal,
                        CmVal,
                        IionCandidate,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );
                }
            }

            updateStateBoundaryValues(u1Candidate, u2Candidate, u3Candidate, t + dt, dim);
            u1Candidate.correctBoundaryConditions();
            u2Candidate.correctBoundaryConditions();
            u3Candidate.correctBoundaryConditions();

            if (profileTimings)
            {
                const auto evalEnd = std::chrono::steady_clock::now();
                nonlinearEvalWallTime +=
                    std::chrono::duration<scalar>(evalEnd - evalStart).count();
                ++nonlinearEvalCalls;
            }
        };

        auto computeSourceDerivative = [&]()
        {
            volScalarField VmPlus
            (
                IOobject
                (
                    "VmPlus",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                VmGuess
            );

            VmPlus.primitiveFieldRef() += diagonalIionEpsilon;
            applyExactVmBoundaryValues(VmPlus, t + dt, dim);

            volScalarField u1Plus
            (
                IOobject
                (
                    "u1Plus",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                u1Old
            );
            volScalarField u2Plus
            (
                IOobject
                (
                    "u2Plus",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                u2Old
            );
            volScalarField u3Plus
            (
                IOobject
                (
                    "u3Plus",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                u3Old
            );
            volScalarField sourcePlus
            (
                IOobject
                (
                    "sourcePlus",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                sourceVm
            );
            volScalarField IionPlus
            (
                IOobject
                (
                    "IionPlus",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                IionGuess
            );

            scalarField VmPlusIntegrationPoints
            (
                allocatedIionIntegrationPoints,
                0.0
            );
            scalarField u1PlusIntegrationPoints(u1OldIntegrationPoints);
            scalarField u2PlusIntegrationPoints(u2OldIntegrationPoints);
            scalarField u3PlusIntegrationPoints(u3OldIntegrationPoints);

            evaluateNonlinearFields
            (
                VmPlus,
                u1Plus,
                u2Plus,
                u3Plus,
                VmPlusIntegrationPoints,
                u1PlusIntegrationPoints,
                u2PlusIntegrationPoints,
                u3PlusIntegrationPoints,
                sourcePlus,
                IionPlus
            );

            scalarField& deriv = sourceDerivative.primitiveFieldRef();
            forAll(deriv, cellI)
            {
                deriv[cellI] =
                    (sourcePlus[cellI] - sourceVm[cellI])
                   /max(diagonalIionEpsilon, SMALL);
            }
            sourceDerivative.correctBoundaryConditions();
        };

        auto rhsFromSource = [&](const volScalarField& sourceCandidate)
        {
            const EigVec sourceNp1 = sourceToEigVec(sourceCandidate);

            EigVec rhs =
                applyB(Vn)
              + theta*(sourceNp1 + bcNp1)
              + (1.0 - theta)*(sourceN + bcN);

            return rhs;
        };

        if (useJFNK)
        {
            evaluateNonlinearFields
            (
                VmGuess,
                u1Guess,
                u2Guess,
                u3Guess,
                VmGuessIntegrationPoints,
                u1GuessIntegrationPoints,
                u2GuessIntegrationPoints,
                u3GuessIntegrationPoints,
                sourceVm,
                IionGuess
            );
        }

        bool nonlinearConverged = false;
        label nonlinearIters = 0;
        scalar finalCoupledResidual = GREAT;
        scalar finalVmResidual = GREAT;
        scalar finalU1Residual = GREAT;
        scalar finalU2Residual = GREAT;
        scalar finalU3Residual = GREAT;
        scalar finalMaxStateResidual = GREAT;
        scalar finalIionResidual = GREAT;
        label finalLinearIterations = 0;
        scalar finalLinearError = GREAT;
        label finalLineSearchIterations = -1;

        auto writeAndPrintResiduals =
        [&]
        (
            const label corr,
            const label linearIterations,
            const scalar linearError,
            const scalar newtonResidual,
            const scalar VmResidual,
            const scalar u1Residual,
            const scalar u2Residual,
            const scalar u3Residual,
            const scalar IionResidual,
            const bool converged,
            const label lineSearchIters = -1
        )
        {
            Info<< "Time = " << t << "," << nl
                << "       Nonlinear method:          " << nonlinearMethod << nl
                << "       Implicit FDA solver:       " << implicitLinearSolver
                << "; iterations = " << linearIterations
                << ", estimated error = " << linearError << nl
                << "       NonLinSolver:              iter = " << corr
                << ", Vm residual = " << VmResidual
                << ", u1 residual = " << u1Residual
                << ", u2 residual = " << u2Residual
                << ", u3 residual = " << u3Residual
                << ", Iion residual = " << IionResidual
                << ", coupled residual = " << newtonResidual;
            if (lineSearchIters >= 0)
            {
                Info<< ", lineSearchIters = " << lineSearchIters;
            }
            Info<< ", converged = " << converged << nl;

            if (writeNonlinearResiduals && nonlinearResidualFilePtr.valid())
            {
                nonlinearResidualFilePtr()
                    << t << ' ' << (nSteps + 1) << ' ' << nonlinearMethod << ' '
                    << stateODESolver << ' '
                    << corr << ' ' << linearIterations << ' ' << linearError << ' '
                    << newtonResidual << ' '
                    << VmResidual << ' ' << u1Residual << ' '
                    << u2Residual << ' ' << u3Residual << ' '
                    << IionResidual << ' '
                    << lineSearchIters << ' ' << converged << nl;
            }

            finalCoupledResidual = newtonResidual;
            finalVmResidual = VmResidual;
            finalU1Residual = u1Residual;
            finalU2Residual = u2Residual;
            finalU3Residual = u3Residual;
            finalMaxStateResidual =
                max(u1Residual, max(u2Residual, u3Residual));
            finalIionResidual = IionResidual;
            finalLinearIterations = linearIterations;
            finalLinearError = linearError;
            finalLineSearchIterations = lineSearchIters;
        };

        if (useJFNK)
        {
            EigVec x = fieldToEigVec(VmGuess);

            volScalarField VmTmp
            (
                IOobject
                (
                    "VmTmp",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                VmGuess
            );
            volScalarField u1Tmp
            (
                IOobject
                (
                    "u1Tmp",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                u1Guess
            );
            volScalarField u2Tmp
            (
                IOobject
                (
                    "u2Tmp",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                u2Guess
            );
            volScalarField u3Tmp
            (
                IOobject
                (
                    "u3Tmp",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                u3Guess
            );
            volScalarField sourceTmp
            (
                IOobject
                (
                    "sourceTmp",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                sourceVm
            );
            volScalarField IionTmp
            (
                IOobject
                (
                    "IionTmp",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                IionGuess
            );

            scalarField VmTmpIntegrationPoints
            (
                allocatedIionIntegrationPoints,
                0.0
            );
            scalarField u1TmpIntegrationPoints(u1OldIntegrationPoints);
            scalarField u2TmpIntegrationPoints(u2OldIntegrationPoints);
            scalarField u3TmpIntegrationPoints(u3OldIntegrationPoints);

            auto residualFor =
            [&](const EigVec& xCandidate, const bool updateGuess) -> EigVec
            {
                if (updateGuess)
                {
                    eigVecToField(xCandidate, VmGuess);
                    applyExactVmBoundaryValues(VmGuess, t + dt, dim);
                    evaluateNonlinearFields
                    (
                        VmGuess,
                        u1Guess,
                        u2Guess,
                        u3Guess,
                        VmGuessIntegrationPoints,
                        u1GuessIntegrationPoints,
                        u2GuessIntegrationPoints,
                        u3GuessIntegrationPoints,
                        sourceVm,
                        IionGuess
                    );

                    const EigVec xApplied = fieldToEigVec(VmGuess);
                    return applyA(xApplied) - rhsFromSource(sourceVm);
                }

                VmTmp = VmGuess;
                u1Tmp = u1Guess;
                u2Tmp = u2Guess;
                u3Tmp = u3Guess;
                sourceTmp = sourceVm;
                IionTmp = IionGuess;

                VmTmpIntegrationPoints = 0.0;
                u1TmpIntegrationPoints = u1OldIntegrationPoints;
                u2TmpIntegrationPoints = u2OldIntegrationPoints;
                u3TmpIntegrationPoints = u3OldIntegrationPoints;

                eigVecToField(xCandidate, VmTmp);
                applyExactVmBoundaryValues(VmTmp, t + dt, dim);

                evaluateNonlinearFields
                (
                    VmTmp,
                    u1Tmp,
                    u2Tmp,
                    u3Tmp,
                    VmTmpIntegrationPoints,
                    u1TmpIntegrationPoints,
                    u2TmpIntegrationPoints,
                    u3TmpIntegrationPoints,
                    sourceTmp,
                    IionTmp
                );

                const EigVec xApplied = fieldToEigVec(VmTmp);
                return applyA(xApplied) - rhsFromSource(sourceTmp);
            };

            const bool usePetscJfnkBackend =
                usesPetscBackend(jfnkLinearSolverBackend);

            // Track when the persistent JFNK PC was last refreshed within
            // the current timestep so that jfnkPreconditionerUpdateFrequency
            // throttling still works on a per-timestep basis. The Eigen ILUT
            // (fallback path) is reused across timesteps only via PETSc;
            // Eigen has no in-place refactor, so the Eigen branch still pays
            // the full compute() cost on each update.
            bool jfnkPcReadyThisStep = false;
            label jfnkPcSetupCorr = -1;

            for
            (
                label corr = 0;
                corr < max(implicitNonlinearIterations, label(1));
                ++corr
            )
            {
                const scalarField VmPrevious(VmGuess.primitiveField());
                const scalarField u1Previous(u1Guess.primitiveField());
                const scalarField u2Previous(u2Guess.primitiveField());
                const scalarField u3Previous(u3Guess.primitiveField());
                const scalarField u1IPPrevious(u1GuessIntegrationPoints);
                const scalarField u2IPPrevious(u2GuessIntegrationPoints);
                const scalarField u3IPPrevious(u3GuessIntegrationPoints);
                const scalarField IionPrevious(IionGuess.primitiveField());

                const EigVec R = residualFor(x, true);
                const EigVec rhsCurrent = rhsFromSource(sourceVm);
                const scalar newtonResidual = relativeL2Norm(R, rhsCurrent);
                const scalar VmResidualInitial =
                    relativeL2Difference(VmGuess.primitiveField(), VmPrevious);
                scalar u1ResidualInitial = GREAT;
                scalar u2ResidualInitial = GREAT;
                scalar u3ResidualInitial = GREAT;
                if (useHighOrder_Iion && dim > 1 && !reconstructStatesFromCellCentres)
                {
                    u1ResidualInitial =
                        relativeL2Difference(u1GuessIntegrationPoints, u1IPPrevious);
                    u2ResidualInitial =
                        relativeL2Difference(u2GuessIntegrationPoints, u2IPPrevious);
                    u3ResidualInitial =
                        relativeL2Difference(u3GuessIntegrationPoints, u3IPPrevious);
                }
                else
                {
                    // cellCentredReconstruct (and the low-order path) carry the
                    // authoritative states at cell centres -> measure there.
                    u1ResidualInitial =
                        relativeL2Difference(u1Guess.primitiveField(), u1Previous);
                    u2ResidualInitial =
                        relativeL2Difference(u2Guess.primitiveField(), u2Previous);
                    u3ResidualInitial =
                        relativeL2Difference(u3Guess.primitiveField(), u3Previous);
                }
                const scalar IionResidualInitial =
                    relativeL2Difference(IionGuess.primitiveField(), IionPrevious);

                const bool residualConverged =
                    corr + 1 >= implicitMinNonlinearIterations
                 && newtonResidual <= nonlinearTolerance
                 && VmResidualInitial <= nonlinearVmTolerance
                 && IionResidualInitial <= nonlinearIionTolerance
                 && (
                        !nonlinearRequireStatesConvergence
                     || (
                            u1ResidualInitial <= nonlinearStatesTolerance
                         && u2ResidualInitial <= nonlinearStatesTolerance
                         && u3ResidualInitial <= nonlinearStatesTolerance
                        )
                    );

                if (residualConverged)
                {
                    writeAndPrintResiduals
                    (
                        corr + 1,
                        0,
                        0.0,
                        newtonResidual,
                        VmResidualInitial,
                        u1ResidualInitial,
                        u2ResidualInitial,
                        u3ResidualInitial,
                        IionResidualInitial,
                        true
                    );
                    nonlinearConverged = true;
                    nonlinearIters = corr + 1;
                    break;
                }

                // Modified for cardiacFoam: globally reduced. An unreduced
                // norm here makes the finite-difference step, and therefore the
                // Jacobian approximation, depend on the rank count.
                const scalar epsScale = jfnkEpsilon*(1.0 + gNorm(x));

                auto matVec = [&](const EigVec& v) -> EigVec
                {
                    const scalar eps = epsScale/max(gNorm(v), SMALL);
                    return (residualFor(x + eps*v, false) - R)/eps;
                };

                label gmresIterations = 0;
                scalar gmresError = GREAT;
                EigVec delta(R.size());

                const auto gmresStart = std::chrono::steady_clock::now();

                if (useJfnkPreconditioner)
                {
                    const bool updatePreconditioner =
                        !jfnkPcReadyThisStep
                     || (
                            corr - jfnkPcSetupCorr
                         >= jfnkPreconditionerUpdateFrequency
                        );

                    if (updatePreconditioner)
                    {
                        const auto pcSetupStart =
                            std::chrono::steady_clock::now();

                        // Reuse the persistent jfnkP buffer. On the first
                        // update we materialise it as a full copy of
                        // AImplicit; afterwards we overwrite just the value
                        // array (the sparsity pattern is invariant because
                        // M, K and AImplicit are time-constant and the
                        // diagonalIion correction only touches existing
                        // diagonal entries).
                        if (!jfnkPMatInitialised)
                        {
                            jfnkP = AImplicit;
                            jfnkPMatInitialised = true;
                        }
                        else
                        {
                            // Same sparsity pattern: copy only the value
                            // buffer. Avoids the ~50 MB allocation that
                            // "SpMat P = AImplicit" used to incur per update.
                            std::copy
                            (
                                AImplicit.valuePtr(),
                                AImplicit.valuePtr() + AImplicit.nonZeros(),
                                jfnkP.valuePtr()
                            );
                        }

                        if (useJfnkDiagonalIionPreconditioner)
                        {
                            computeSourceDerivative();
                            addDiagonalToMatrix
                            (
                                sourceDerivative.primitiveField(),
                                LREInterp_Vm.globalCells(),
                               -theta,
                                jfnkP
                            );
                        }

                        if (usePetscJfnkBackend)
                        {
                            if (!jfnkPetscPcSolver.isInitialised())
                            {
                                // First-ever setup: full reset() builds the
                                // PETSc Mat, KSP, PC and ILUT factor.
                                jfnkPetscPcSolver.reset
                                (
                                    jfnkP,
                                    jfnkPreconditionerKspType,
                                    jfnkPreconditionerPcType,
                                    jfnkPreconditionerTolerance,
                                    jfnkPreconditionerMaxIterations,
                                    jfnkMaxKrylovIterations,
                                    jfnkPreconditionerOptionsPrefix,
                                    petscUseOptions,
                                    scalar(jfnkPreconditionerFillFactor),
                                    jfnkPreconditionerDropTolerance
                                );
                            }
                            else
                            {
                                // Subsequent updates: refresh values in the
                                // cached PETSc Mat and force the PC to
                                // recompute. The ILUT symbolic factorisation
                                // is reused — only the numerical factor is
                                // recomputed. Eliminates the per-update
                                // alloc/free of ~500-700 MB that otherwise
                                // fragments the heap on long runs.
                                jfnkPetscPcSolver.updateValues(jfnkP);
                            }
                        }
                        else
                        {
                            jfnkIlu.setDroptol(jfnkPreconditionerDropTolerance);
                            jfnkIlu.setFillfactor(jfnkPreconditionerFillFactor);
                            jfnkIlu.compute(jfnkP);

                            if (jfnkIlu.info() != Eigen::Success)
                            {
                                FatalErrorInFunction
                                    << "JFNK ILUT preconditioner setup failed"
                                    << exit(FatalError);
                            }
                        }

                        jfnkPcReadyThisStep = true;
                        jfnkPcSetupCorr = corr;

                        if (profileTimings)
                        {
                            const auto pcSetupEnd =
                                std::chrono::steady_clock::now();
                            preconditionerSetupWallTime +=
                                std::chrono::duration<scalar>
                                (
                                    pcSetupEnd - pcSetupStart
                                ).count();
                            ++preconditionerSetups;
                        }
                    }

                    std::function<EigVec(const EigVec&)> applyPreconditioner =
                    [&](const EigVec& r) -> EigVec
                    {
                        const auto pcApplyStart =
                            std::chrono::steady_clock::now();

                        EigVec z(r.size());
                        if (usePetscJfnkBackend)
                        {
                            label pcIterations = 0;
                            scalar pcError = GREAT;
                            z = jfnkPetscPcSolver.solve
                            (
                                r,
                                pcIterations,
                                pcError
                            );
                        }
                        else
                        {
                            z = jfnkIlu.solve(r);
                        }

                        if (profileTimings)
                        {
                            const auto pcApplyEnd =
                                std::chrono::steady_clock::now();
                            preconditionerApplyWallTime +=
                                std::chrono::duration<scalar>
                                (
                                    pcApplyEnd - pcApplyStart
                                ).count();
                            ++preconditionerApplications;
                        }

                        return z;
                    };

                    if (usePetscJfnkBackend)
                    {
                        if (!jfnkPetscShellSolver.isInitialised())
                        {
                            jfnkPetscShellSolver.initialise
                            (
                                static_cast<label>(R.size()),
                                jfnkPetscKspType,
                                jfnkPetscPcType,
                                jfnkPetscRestart,
                                max
                                (
                                    jfnkMaxKrylovIterations
                                   *(jfnkMaxRestarts + 1),
                                    label(1)
                                ),
                                jfnkLinearTolerance,
                                jfnkPetscOptionsPrefix,
                                petscUseOptions,
                                true /* withShellPC */
                            );
                        }

                        delta = jfnkPetscShellSolver.solve
                        (
                            matVec,
                            &applyPreconditioner,
                            -R,
                            gmresIterations,
                            gmresError
                        );
                    }
                    else
                    {
                        delta = solveLeftPreconditionedGMRES
                        (
                            matVec,
                            applyPreconditioner,
                            -R,
                            jfnkMaxKrylovIterations,
                            jfnkMaxRestarts,
                            jfnkLinearTolerance,
                            gmresIterations,
                            gmresError
                        );
                    }
                }
                else
                {
                    if (usePetscJfnkBackend)
                    {
                        if (!jfnkPetscShellSolver.isInitialised())
                        {
                            jfnkPetscShellSolver.initialise
                            (
                                static_cast<label>(R.size()),
                                jfnkPetscKspType,
                                jfnkPetscPcType,
                                jfnkPetscRestart,
                                max
                                (
                                    jfnkMaxKrylovIterations
                                   *(jfnkMaxRestarts + 1),
                                    label(1)
                                ),
                                jfnkLinearTolerance,
                                jfnkPetscOptionsPrefix,
                                petscUseOptions,
                                false /* withShellPC */
                            );
                        }

                        delta = jfnkPetscShellSolver.solve
                        (
                            matVec,
                            nullptr,
                            -R,
                            gmresIterations,
                            gmresError
                        );
                    }
                    else
                    {
                        delta = solveGMRES
                        (
                            matVec,
                            -R,
                            jfnkMaxKrylovIterations,
                            jfnkMaxRestarts,
                            jfnkLinearTolerance,
                            gmresIterations,
                            gmresError
                        );
                    }
                }

                if (profileTimings)
                {
                    const auto gmresEnd = std::chrono::steady_clock::now();
                    gmresWallTime +=
                        std::chrono::duration<scalar>
                        (
                            gmresEnd - gmresStart
                        ).count();
                    ++gmresCalls;
                }

                label lineSearchIters = 0;
                EigVec Rnew(R.size());
                scalar newtonResidualNew = GREAT;

                if (jfnkLineSearch)
                {
                    const EigVec xBackup = x;
                    const scalar Rold2 = gSquaredNorm(R);   // Modified for cardiacFoam
                    scalar alpha = nonlinearRelaxation;

                    while (true)
                    {
                        x = xBackup + alpha*delta;
                        Rnew = residualFor(x, true);
                        const scalar Rnew2 = gSquaredNorm(Rnew);   // Modified for cardiacFoam

                        if (Rnew2 <= (1.0 - 2.0*jfnkArmijoC*alpha)*Rold2)
                        {
                            break;
                        }

                        ++lineSearchIters;
                        if (lineSearchIters >= jfnkLineSearchMaxIter)
                        {
                            break;
                        }

                        alpha *= 0.5;
                        if (alpha < jfnkLineSearchAlphaMin)
                        {
                            alpha = jfnkLineSearchAlphaMin;
                            x = xBackup + alpha*delta;
                            Rnew = residualFor(x, true);
                            break;
                        }
                    }
                }
                else
                {
                    x += nonlinearRelaxation*delta;
                    Rnew = residualFor(x, true);
                }

                const EigVec rhsNew = rhsFromSource(sourceVm);
                newtonResidualNew = relativeL2Norm(Rnew, rhsNew);

                const scalar VmResidual =
                    relativeL2Difference(VmGuess.primitiveField(), VmPrevious);
                scalar u1Residual = GREAT;
                scalar u2Residual = GREAT;
                scalar u3Residual = GREAT;
                if (useHighOrder_Iion && dim > 1 && !reconstructStatesFromCellCentres)
                {
                    u1Residual =
                        relativeL2Difference(u1GuessIntegrationPoints, u1IPPrevious);
                    u2Residual =
                        relativeL2Difference(u2GuessIntegrationPoints, u2IPPrevious);
                    u3Residual =
                        relativeL2Difference(u3GuessIntegrationPoints, u3IPPrevious);
                }
                else
                {
                    // cellCentredReconstruct (and the low-order path) carry the
                    // authoritative states at cell centres -> measure there.
                    u1Residual =
                        relativeL2Difference(u1Guess.primitiveField(), u1Previous);
                    u2Residual =
                        relativeL2Difference(u2Guess.primitiveField(), u2Previous);
                    u3Residual =
                        relativeL2Difference(u3Guess.primitiveField(), u3Previous);
                }
                const scalar IionResidual =
                    relativeL2Difference(IionGuess.primitiveField(), IionPrevious);

                const bool converged =
                    corr + 1 >= implicitMinNonlinearIterations
                 && newtonResidualNew <= nonlinearTolerance
                 && VmResidual <= nonlinearVmTolerance
                 && IionResidual <= nonlinearIionTolerance
                 && (
                        !nonlinearRequireStatesConvergence
                     || (
                            u1Residual <= nonlinearStatesTolerance
                         && u2Residual <= nonlinearStatesTolerance
                         && u3Residual <= nonlinearStatesTolerance
                        )
                    );

                writeAndPrintResiduals
                (
                    corr + 1,
                    gmresIterations,
                    gmresError,
                    newtonResidualNew,
                    VmResidual,
                    u1Residual,
                    u2Residual,
                    u3Residual,
                    IionResidual,
                    converged,
                    lineSearchIters
                );

                nonlinearIters = corr + 1;
                if (converged)
                {
                    nonlinearConverged = true;
                    break;
                }
            }
        }
        else if (usePicard || useDiagonalIion)
        {
            // Persistent ACurrent buffer reused across the corr iterations
            // (and across timesteps). Only allocated when actually needed —
            // pure Picard with a cached PETSc solver never touches it. For
            // useDiagonalIion the value buffer is refreshed each iter from
            // AImplicit and then the diagonal source derivative correction
            // is applied in-place; sparsity is invariant.
            SpMat ACurrent;
            bool ACurrentInitialised = false;
            const bool needsACurrent =
                useDiagonalIion
             || (usePicard && !useCachedPicardPetscSolver);

            for
            (
                label corr = 0;
                corr < max(implicitNonlinearIterations, label(1));
                ++corr
            )
            {
                const scalarField VmPrevious(VmGuess.primitiveField());
                const scalarField u1Previous(u1Guess.primitiveField());
                const scalarField u2Previous(u2Guess.primitiveField());
                const scalarField u3Previous(u3Guess.primitiveField());
                const scalarField u1IPPrevious(u1GuessIntegrationPoints);
                const scalarField u2IPPrevious(u2GuessIntegrationPoints);
                const scalarField u3IPPrevious(u3GuessIntegrationPoints);
                const scalarField IionPrevious(IionGuess.primitiveField());

                evaluateNonlinearFields
                (
                    VmGuess,
                    u1Guess,
                    u2Guess,
                    u3Guess,
                    VmGuessIntegrationPoints,
                    u1GuessIntegrationPoints,
                    u2GuessIntegrationPoints,
                    u3GuessIntegrationPoints,
                    sourceVm,
                    IionGuess
                );

                EigVec rhs = rhsFromSource(sourceVm);

                if (needsACurrent)
                {
                    if (!ACurrentInitialised)
                    {
                        ACurrent = AImplicit;
                        ACurrentInitialised = true;
                    }
                    else
                    {
                        // Refresh values only; sparsity pattern is invariant
                        // (M, K and AImplicit are time-constant and only the
                        // diagonal source term gets added below).
                        std::copy
                        (
                            AImplicit.valuePtr(),
                            AImplicit.valuePtr() + AImplicit.nonZeros(),
                            ACurrent.valuePtr()
                        );
                    }
                }

                if (useDiagonalIion)
                {
                    computeSourceDerivative();
                    const EigVec sourceDerivativeVec =
                        fieldToEigVec(sourceDerivative);
                    const EigVec VmLinearisationPoint =
                        fieldToEigVec(VmGuess);

                    addDiagonalToMatrix
                    (
                        sourceDerivative.primitiveField(),
                        LREInterp_Vm.globalCells(),
                       -theta,
                        ACurrent
                    );

                    rhs -= theta
                       * sourceDerivativeVec.cwiseProduct(VmLinearisationPoint);
                }

                label linearIterations = 0;
                scalar linearError = GREAT;

                const auto sparseSolveStart =
                    std::chrono::steady_clock::now();
                EigVec Vsol(rhs.size());
                if (useCachedPicardPetscSolver && usePicard)
                {
                    Vsol = cachedPicardPetscSolver.solve
                    (
                        rhs,
                        linearIterations,
                        linearError
                    );
                }
                else
                {
                    Vsol =
                        solveSparseSystem
                        (
                            ACurrent,
                            rhs,
                            linearSolverBackend,
                            implicitLinearSolver,
                            petscLinearKspType,
                            petscLinearPcType,
                            petscLinearRestart,
                            petscLinearOptionsPrefix,
                            petscUseOptions,
                            implicitTolerance,
                            implicitMaxIterations,
                            linearIterations,
                            linearError
                        );
                }
                if (profileTimings)
                {
                    const auto sparseSolveEnd =
                        std::chrono::steady_clock::now();
                    sparseLinearSolveWallTime +=
                        std::chrono::duration<scalar>
                        (
                            sparseSolveEnd - sparseSolveStart
                        ).count();
                    ++sparseLinearSolveCalls;
                }

                const EigVec Vprev = fieldToEigVec(VmGuess);
                const EigVec Vrelaxed =
                    Vprev + nonlinearRelaxation*(Vsol - Vprev);

                eigVecToField(Vrelaxed, VmGuess);
                applyExactVmBoundaryValues(VmGuess, t + dt, dim);

                if (useDiagonalIion)
                {
                    evaluateNonlinearFields
                    (
                        VmGuess,
                        u1Guess,
                        u2Guess,
                        u3Guess,
                        VmGuessIntegrationPoints,
                        u1GuessIntegrationPoints,
                        u2GuessIntegrationPoints,
                        u3GuessIntegrationPoints,
                        sourceVm,
                        IionGuess
                    );
                }
                else
                {
                    updateStateBoundaryValues(u1Guess, u2Guess, u3Guess, t + dt, dim);
                    u1Guess.correctBoundaryConditions();
                    u2Guess.correctBoundaryConditions();
                    u3Guess.correctBoundaryConditions();
                }

                const EigVec nonlinearResidualVec =
                    applyA(fieldToEigVec(VmGuess)) - rhsFromSource(sourceVm);
                const scalar newtonResidual =
                    relativeL2Norm(nonlinearResidualVec, rhsFromSource(sourceVm));

                const scalar VmResidual =
                    relativeL2Difference(VmGuess.primitiveField(), VmPrevious);
                scalar u1Residual = GREAT;
                scalar u2Residual = GREAT;
                scalar u3Residual = GREAT;
                if (useHighOrder_Iion && dim > 1 && !reconstructStatesFromCellCentres)
                {
                    u1Residual =
                        relativeL2Difference(u1GuessIntegrationPoints, u1IPPrevious);
                    u2Residual =
                        relativeL2Difference(u2GuessIntegrationPoints, u2IPPrevious);
                    u3Residual =
                        relativeL2Difference(u3GuessIntegrationPoints, u3IPPrevious);
                }
                else
                {
                    // cellCentredReconstruct (and the low-order path) carry the
                    // authoritative states at cell centres -> measure there.
                    u1Residual =
                        relativeL2Difference(u1Guess.primitiveField(), u1Previous);
                    u2Residual =
                        relativeL2Difference(u2Guess.primitiveField(), u2Previous);
                    u3Residual =
                        relativeL2Difference(u3Guess.primitiveField(), u3Previous);
                }
                const scalar IionResidual =
                    relativeL2Difference(IionGuess.primitiveField(), IionPrevious);

                const bool statesConverged =
                    u1Residual <= nonlinearStatesTolerance
                 && u2Residual <= nonlinearStatesTolerance
                 && u3Residual <= nonlinearStatesTolerance;

                const bool converged =
                    corr + 1 >= implicitMinNonlinearIterations
                 && VmResidual <= nonlinearVmTolerance
                 && (
                        usePicard
                      ? statesConverged
                      : (
                            newtonResidual <= nonlinearTolerance
                         && IionResidual <= nonlinearIionTolerance
                         && (
                                !nonlinearRequireStatesConvergence
                             || statesConverged
                            )
                        )
                    );

                writeAndPrintResiduals
                (
                    corr + 1,
                    linearIterations,
                    linearError,
                    newtonResidual,
                    VmResidual,
                    u1Residual,
                    u2Residual,
                    u3Residual,
                    IionResidual,
                    converged
                );

                nonlinearIters = corr + 1;
                if (converged)
                {
                    nonlinearConverged = true;
                    break;
                }
            }
        }

        bool nonlinearRolledBack = false;

        if (!nonlinearConverged && usePicard)
        {
            WarningInFunction
                << "Picard nonlinear solver did not converge at t = " << t
                << " after " << nonlinearIters << " iterations."
                << " Accepting the last Picard iterate, matching the"
                << " highOrderManufacturedFDAImplicit behaviour." << endl;
        }
        else if (!nonlinearConverged && !nonlinearAcceptUnconverged)
        {
            WarningInFunction
                << "Nonlinear solver did not converge at t = " << t
                << " after " << nonlinearIters << " iterations."
                << " Rolling back Vm, Iion and states to the values"
                << " at the beginning of the time step." << nl
                << "Set 'nonlinearAcceptUnconverged true;' in"
                << " spatialIntegrationProperties to accept unconverged"
                << " iterates instead." << endl;

            VmGuess.primitiveFieldRef() = VmOld.primitiveField();
            applyExactVmBoundaryValues(VmGuess, t, dim);
            u1Guess.primitiveFieldRef() = u1Old.primitiveField();
            u2Guess.primitiveFieldRef() = u2Old.primitiveField();
            u3Guess.primitiveFieldRef() = u3Old.primitiveField();
            updateStateBoundaryValues(u1Guess, u2Guess, u3Guess, t, dim);
            u1Guess.correctBoundaryConditions();
            u2Guess.correctBoundaryConditions();
            u3Guess.correctBoundaryConditions();
            IionGuess.primitiveFieldRef() = Iion.primitiveField();
            IionGuess.correctBoundaryConditions();
            nonlinearRolledBack = true;

            if (useHighOrder_Iion && dim > 1)
            {
                u1GuessIntegrationPoints = u1OldIntegrationPoints;
                u2GuessIntegrationPoints = u2OldIntegrationPoints;
                u3GuessIntegrationPoints = u3OldIntegrationPoints;
            }
        }
        else if (!nonlinearConverged && nonlinearAcceptUnconverged)
        {
            WarningInFunction
                << "Nonlinear solver did not converge at t = " << t
                << " after " << nonlinearIters << " iterations."
                << " Accepting the unconverged iterate"
                << " (nonlinearAcceptUnconverged=true)." << endl;
        }

        NonlinearConvergenceRecord nonlinearRecord;
        nonlinearRecord.time = t + dt;
        nonlinearRecord.step = nSteps + 1;
        nonlinearRecord.iterations = nonlinearIters;
        nonlinearRecord.maxIterations =
            max(implicitNonlinearIterations, label(1));
        nonlinearRecord.converged = nonlinearConverged;
        nonlinearRecord.rolledBack = nonlinearRolledBack;
        nonlinearRecord.coupledResidual = finalCoupledResidual;
        nonlinearRecord.VmResidual = finalVmResidual;
        nonlinearRecord.u1Residual = finalU1Residual;
        nonlinearRecord.u2Residual = finalU2Residual;
        nonlinearRecord.u3Residual = finalU3Residual;
        nonlinearRecord.maxStateResidual = finalMaxStateResidual;
        nonlinearRecord.IionResidual = finalIionResidual;
        nonlinearRecord.linearIterations = finalLinearIterations;
        nonlinearRecord.linearError = finalLinearError;
        nonlinearRecord.lineSearchIterations = finalLineSearchIterations;
        nonlinearHistory.push_back(nonlinearRecord);

        Vm.primitiveFieldRef() = VmGuess.primitiveField();
        u1.primitiveFieldRef() = u1Guess.primitiveField();
        u2.primitiveFieldRef() = u2Guess.primitiveField();
        u3.primitiveFieldRef() = u3Guess.primitiveField();

        if (useHighOrder_Iion && dim > 1)
        {
            u1IntegrationPoints = u1GuessIntegrationPoints;
            u2IntegrationPoints = u2GuessIntegrationPoints;
            u3IntegrationPoints = u3GuessIntegrationPoints;
        }

        const scalar committedTime = nonlinearRolledBack ? t : t + dt;
        applyExactVmBoundaryValues(Vm, committedTime, dim);
        updateStateBoundaryValues(u1, u2, u3, committedTime, dim);
        u1.correctBoundaryConditions();
        u2.correctBoundaryConditions();
        u3.correctBoundaryConditions();

        if (useHighOrder_Iion && dim > 1)
        {
            reconstructVmAtIionIntegrationPoints
            (
                Vm,
                useHighOrder_Vm,
                LREInterp_Vm,
                LREInterp_Iion,
                VmGuessIntegrationPoints
            );

            if (reconstructStatesFromCellCentres)
            {
                // Keep the Gauss-point states consistent with the committed
                // cell-centred states (matters after a rollback too).
                reconstructStatesAtIionIntegrationPoints
                (
                    u1,
                    u2,
                    u3,
                    LREInterp_statesPtr(),
                    LREInterp_Iion,
                    u1IntegrationPoints,
                    u2IntegrationPoints,
                    u3IntegrationPoints
                );
            }

            computeIionFromIntegrationPoints
            (
                VmGuessIntegrationPoints,
                u1IntegrationPoints,
                u2IntegrationPoints,
                u3IntegrationPoints,
                beta,
                chiVal,
                CmVal,
                IionIntegrationPoints,
                stateODEUseOpenMP,
                stateODEOpenMPThreshold
            );

            averageIntegrationPointFieldToCells
            (
                IionIntegrationPoints,
                LREInterp_Iion,
                Iion
            );
        }
        else
        {
            computeCellCentredIion
            (
                Vm,
                u1,
                u2,
                u3,
                beta,
                chiVal,
                CmVal,
                Iion,
                stateODEUseOpenMP,
                stateODEOpenMPThreshold
            );
        }

        ++nSteps;
        ++runTime;

        const bool needsRhsVm =
            nSteps % 50 == 0 || nSteps <= 5 || runTime.outputTime();

        if (needsRhsVm)
        {
            updateNumericalLaplacian(t + dt, runTime.outputTime());
            rhsVm = lapVm/(chi*Cm) - Iion;
        }

        if (nSteps % 50 == 0 || nSteps <= 5)
        {
            Info<< "Time step " << nSteps
                << " : t = " << runTime.value()
                << " , RHS Linf = " << linfNorm(rhsVm)
                << endl;
        }

        if (runTime.outputTime())
        {
            fillExactFields(VmExact, u1Exact, u2Exact, runTime.value(), dim);
            VmError = Vm - VmExact;
            u1Error = u1 - u1Exact;
            u2Error = u2 - u2Exact;
            runTime.write();
        }
    }

    const auto tEndLoop = std::chrono::steady_clock::now();
    const auto tStartPost = tEndLoop;

    fillExactFields(VmExact, u1Exact, u2Exact, runTime.value(), dim);

    if (useHighOrder_Iion && dim > 1)
    {
        reconstructVmAtIionIntegrationPoints
        (
            Vm,
            useHighOrder_Vm,
            LREInterp_Vm,
            LREInterp_Iion,
            VmGuessIntegrationPoints
        );

        computeIionFromIntegrationPoints
        (
            VmGuessIntegrationPoints,
            u1IntegrationPoints,
            u2IntegrationPoints,
            u3IntegrationPoints,
            beta,
            chiVal,
            CmVal,
            IionIntegrationPoints,
            stateODEUseOpenMP,
            stateODEOpenMPThreshold
        );

        averageIntegrationPointFieldToCells
        (
            IionIntegrationPoints,
            LREInterp_Iion,
            Iion
        );
    }
    else
    {
        computeCellCentredIion
        (
            Vm,
            u1,
            u2,
            u3,
            beta,
            chiVal,
            CmVal,
            Iion,
            stateODEUseOpenMP,
            stateODEOpenMPThreshold
        );
    }

    updateNumericalLaplacian(runTime.value(), true);
    rhsVm = lapVm/(chi*Cm) - Iion;

    VmError = Vm - VmExact;
    u1Error = u1 - u1Exact;
    u2Error = u2 - u2Exact;

    const FieldErrorSummary VmCell = computeFieldErrorSummary(VmError);
    const FieldErrorSummary u1Cell = computeFieldErrorSummary(u1Error);
    const FieldErrorSummary u2Cell = computeFieldErrorSummary(u2Error);

    const ReconstructionErrorSummary VmHO =
            useHighOrder_Vm && dim > 1
        ? computeGaussReconstructedError(Vm, runTime.value(), mfVm, dim, LREInterp_Vm)
        : cellAsReconstructionError(VmCell, VmExact);

    const ReconstructionErrorSummary u1HO =
            useHighOrder_Iion && dim > 1
        ? computeIntegrationPointStateError
          (
              u1IntegrationPoints,
              runTime.value(),
              mfU1,
              dim,
              LREInterp_Iion,
              mesh
          )
        : cellAsReconstructionError(u1Cell, u1Exact);

    const ReconstructionErrorSummary u2HO =
            useHighOrder_Iion && dim > 1
        ? computeIntegrationPointStateError
          (
              u2IntegrationPoints,
              runTime.value(),
              mfU2,
              dim,
              LREInterp_Iion,
              mesh
          )
        : cellAsReconstructionError(u2Cell, u2Exact);

    const auto tEndPost = std::chrono::steady_clock::now();

    const scalar setupWallTime =
        std::chrono::duration_cast<std::chrono::duration<scalar>>
        (
            tEndSetup - tStartSetup
        ).count();

    const scalar timeLoopWallTime =
        std::chrono::duration_cast<std::chrono::duration<scalar>>
        (
            tEndLoop - tStartLoop
        ).count();

    const scalar postProcessWallTime =
        std::chrono::duration_cast<std::chrono::duration<scalar>>
        (
            tEndPost - tStartPost
        ).count();

    const scalar totalWallTime =
        std::chrono::duration_cast<std::chrono::duration<scalar>>
        (
            tEndPost - tStartTotal
        ).count();

    if (profileTimings)
    {
        Info<< "Fine-grained timing [s]:" << nl
            << "  nonlinear evaluations: calls = " << nonlinearEvalCalls
            << ", wall = " << nonlinearEvalWallTime << nl
            << "  JFNK KSP/GMRES: calls = " << gmresCalls
            << ", wall = " << gmresWallTime << nl
            << "  JFNK preconditioner setup: calls = " << preconditionerSetups
            << ", wall = " << preconditionerSetupWallTime << nl
            << "  JFNK preconditioner apply: calls = "
            << preconditionerApplications
            << ", wall = " << preconditionerApplyWallTime << nl
            << "  sparse PDE solves: calls = " << sparseLinearSolveCalls
            << ", wall = " << sparseLinearSolveWallTime << endl;
    }

    runTime.write();

    writeSummary
    (
        runTime,
        nSteps,
        dt,
        VmCell,
        u1Cell,
        u2Cell,
        VmHO,
        u1HO,
        u2HO,
        rhsVm,
        Iion,
        mesh,
        implicitScheme,
        massMatrixMode,
        useHighOrder_Vm,
        useHighOrder_Iion,
        stabilisationAlpha,
        memoryOptimization,
        memoryOptimizationEffective,
        memoryOptimizationCellThreshold,
        memoryOptimizationAdaptiveTripletReserve,
        memoryOptimizationFluxTripletsPerFace,
        memoryOptimizationStabilisedTripletsPerFace,
        memoryOptimizationTrimHeap,
        memoryOptimizationCompactMassAssembly,
        compactMassAssembly,
        stiffnessTripletsPerFaceReserve,
        allocatedIionIntegrationPoints,
        lreN,
        lreNn,
        lreK,
        lreMaxStencil,
        lreIionN,
        lreIionNn,
        lreIionK,
        lreIionMaxStencil,
        linearSolverBackend,
        implicitLinearSolver,
        petscLinearKspType,
        petscLinearPcType,
        implicitTolerance,
        implicitMaxIterations,
        max(implicitNonlinearIterations, label(1)),
        implicitMinNonlinearIterations,
        nonlinearTolerance,
        nonlinearVmTolerance,
        nonlinearStatesTolerance,
        nonlinearRequireStatesConvergence,
        nonlinearAcceptUnconverged,
        nonlinearIionTolerance,
        nonlinearRelaxation,
        jfnkMaxKrylovIterations,
        jfnkMaxRestarts,
        jfnkLinearTolerance,
        jfnkEpsilon,
        jfnkLinearSolverBackend,
        jfnkPetscKspType,
        jfnkPetscPcType,
        jfnkInitGuessOrder,
        jfnkInitGuessVmMin,
        jfnkInitGuessVmMax,
        jfnkClampODEInput,
        jfnkLineSearch,
        jfnkLineSearchMaxIter,
        jfnkLineSearchAlphaMin,
        jfnkArmijoC,
        jfnkPreconditioner,
        jfnkPreconditionerDropTolerance,
        jfnkPreconditionerFillFactor,
        jfnkPreconditionerUpdateFrequency,
        diagonalIionEpsilon,
        stateODESolver,
        stateODEInitialStep,
        stateODEAbsTol,
        stateODERelTol,
        stateODEMaxSteps,
        stateODEUseOpenMP,
        stateODENumThreads,
        nonlinearHistory,
        peakResidentSetSizeKB(),
        setupWallTime,
        timeLoopWallTime,
        postProcessWallTime,
        totalWallTime,
        nonlinearMethod
    );

    Info<< "End" << nl << endl;
    return 0;
}
