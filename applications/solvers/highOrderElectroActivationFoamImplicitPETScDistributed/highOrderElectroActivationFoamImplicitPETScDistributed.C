/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

Solver
    highOrderElectroActivationFoamImplicitPETSc

Description
    Implicit monodomain solver for cardiac electrophysiology with Jacobian-Free
    Newton-Krylov (JFNK) nonlinear coupling.

    PDE solved (monodomain):

        chi * Cm * dVm/dt  =  div( D . grad Vm )  -  Iion(Vm, s)  +  Iext
                  ds/dt    =  f_ionic(Vm, s)

    where Vm is the transmembrane potential, D the conductivity tensor, chi
    the surface-to-volume ratio, Cm the membrane capacitance, Iion the
    (nonlinear) ionic current, s the vector of ionic states, and Iext an
    external stimulus current.

    Numerical features:
      * Time integration: theta-method (Backward Euler or Crank-Nicolson).
      * Mass matrix: lumped (diagonal) or LRE-consistent high-order.
      * Diffusion operator: standard FV orthogonal stencil or LRE high-order
        face-quadrature reconstruction (Taylor expansion up to 3rd order).
      * Ionic source: optionally evaluated at LRE cell quadrature points to
        match the spatial accuracy of the diffusion operator.
      * Nonlinear solve: JFNK (matrix-free Newton + restarted GMRES with
        modified Gram-Schmidt) or Picard / diagonal-linearised variants.
      * Linear solve (non-JFNK branches): Eigen SparseLU or BiCGSTAB.
      * Robustness: if the nonlinear iteration fails to converge, Vm, Iion
        and ionic states are rolled back to the values at the beginning of
        the time step (configurable via nonlinearAcceptUnconverged).

    Diagnostics:
      * activationTime: first crossing of activationThreshold (default 0 V),
        recorded via linear interpolation across the time step.
      * P8 point and P1-P8 diagonal activation samples (Niederer 2012
        benchmark) written to postProcessing.
      * Per-time-step nonlinear residuals written to nonlinearResiduals.dat.

\*---------------------------------------------------------------------------*/

#include <petscksp.h>
#include "fvCFD.H"
// Modified for cardiacFoam: LRE (serial) -> highOrderInterp, a thin
// adapter over solids4foam's parallel movingLeastSquares. The API
// mapping, and the two places where semantics changed (quadrature
// weights, stencil storage), are documented in
// src/highOrderAdapter/highOrderInterp.H.
#include "highOrderInterp.H"
#include "globalIndex.H"
#include "syncTools.H"
// Added for cardiacFoam: needed by the coupled-patch (processor face)
// stabilisation exchange in the stiffness assemblies.
#include "PstreamBuffers.H"
#include "processorFvPatch.H"
#include "Field.H"
#include "volFields.H"
#include <cmath>
#include <chrono>
#include <algorithm>
#include <cctype>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
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

    // Added for cardiacFoam: global row offset of this rank's block of cells,
    // i.e. globalIndex::localStart(). Constant for the whole run (one mesh, one
    // decomposition), so it is set once in main() rather than threaded through
    // every linear-solver signature. Zero in serial.
    label gRowStart = 0;

    // Added for cardiacFoam: global cell addressing used for matrix COLUMNS.
    // Owned by main(); the assembly routines are free functions, so a
    // file-scope handle avoids threading an extra argument through all of them.
    autoPtr<globalIndex> gCellsPtr;

    // Added for cardiacFoam: local cell index -> global column index.
    inline label gCol(const label localCellI)
    {
        return gCellsPtr->toGlobal(localCellI);
    }

    // Added for cardiacFoam: total number of columns (global cell count).
    inline label gNCols()
    {
        return gCellsPtr->totalSize();
    }

    using SpMat = Eigen::SparseMatrix<scalar, Eigen::RowMajor>;
    using Triplet = Eigen::Triplet<scalar>;
    using EigVec = Eigen::Matrix<scalar, Eigen::Dynamic, 1>;

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
        Eigen::IncompleteLUT<scalar> jfnkPrec;
        bool jfnkPrecComputed = false;
    };

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

    // Lower-case copy of a dictionary word, so that user input is matched
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

    // Map the solver names used by the dictionaries onto PETSc KSP type names,
    // so that the same case file drives either backend. Names PETSc already
    // knows pass through unchanged. Note that the direct solvers (sparseLU/lu)
    // become "preonly", i.e. apply the preconditioner once - it is PCLU that
    // does the factorisation.
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
    // spellings the dictionaries have accumulated (off/false/none/nopc, ilut,
    // jakobi) and the AMG aliases. Anything else passes through, so any PC
    // PETSc supports can be named directly.
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
        // Phase C3: algebraic multigrid convenience aliases. For the assembled
        // (Picard / diagonalIion) branch the monodomain operator is SPD, so
        // "kspType cg" + "pcType gamg" (PETSc native AMG) or "hypre" (BoomerAMG)
        // scales far better than ILU/BiCGSTAB/LU on large 3D meshes. Any other
        // PETSc PC name (gamg, hypre, gasm, bjacobi, ...) already passes through
        // unchanged below, so no explicit mapping is required for them.
        if (s == "amg")
        {
            return "gamg";
        }
        if (s == "boomeramg")
        {
            return "hypre";
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

    // True when the requested KSP is restarted GMRES, which is the only one
    // that needs a restart length to be configured.
    bool isGmresType(const word& kspType)
    {
        return petscKspTypeName(kspType) == "gmres";
    }

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

    // Persistent PETSc KSP for an ASSEMBLED matrix, kept alive across time
    // steps and nonlinear iterations.
    //
    // The point is the MatrixUpdate policy: rebuilding the matrix and its
    // preconditioner on every solve dominates the cost when, as in Picard with
    // a constant operator, the matrix does not change at all. Holding the Mat,
    // the KSP and the two work vectors lets a solve be reduced to Reuse (solve
    // only) or ValuesOnly (refresh values, re-factorise, keep the pattern).
    class PetscKspMatrixSolver
    {
        Mat A_;
        KSP ksp_;
        Vec b_;
        Vec x_;
        label n_;

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
            const PetscInt nRows = static_cast<PetscInt>(A.rows());

            // Modified for cardiacFoam: was MatCreateSeqAIJ, which PETSc
            // rejects on a communicator of size > 1 ("Comm must be of size 1").
            // buildDistributedMat produces MATAIJ - sequential on one rank,
            // MPIAIJ on several - from the same local-rows/global-columns Eigen
            // matrix, with the diagonal/off-diagonal preallocation split that
            // stops assembly degenerating into repeated reallocation.
            buildDistributedMat(A, A_);

            // Modified for cardiacFoam: distributed vector with nRows LOCAL entries.
            checkPetscError
            (
                VecCreateMPI(PETSC_COMM_WORLD, nRows, PETSC_DETERMINE, &b_),
                "VecCreateMPI(b)"
            );
            checkPetscError(VecDuplicate(b_, &x_), "VecDuplicate(x)");

            checkPetscError(KSPCreate(PETSC_COMM_WORLD, &ksp_), "KSPCreate");
            checkPetscError(KSPSetOperators(ksp_, A_, A_), "KSPSetOperators");

            const std::string kspName = petscKspTypeName(kspType);
            // Modified for cardiacFoam: ilu/lu have no MPIAIJ implementation
            // and KSPSetUp fails with PETSC_ERR_SUP. In parallel they are
            // wrapped in block-Jacobi, which applies the requested
            // factorisation on each rank's diagonal block. The driver's
            // default pcType is ilu, so this is needed from the first run.
            std::string subPcName;
            const std::string pcName = petscParallelPcTypeName(pcType, subPcName);
            checkPetscError(KSPSetType(ksp_, kspName.c_str()), "KSPSetType");

            PC pc = nullptr;
            checkPetscError(KSPGetPC(ksp_, &pc), "KSPGetPC");
            checkPetscError(PCSetType(pc, pcName.c_str()), "PCSetType");

            // Added for cardiacFoam: configure the per-block factorisation.
            if (!subPcName.empty())
            {
                checkPetscError(KSPSetUp(ksp_), "KSPSetUp(bjacobi)");
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
                    checkPetscError(KSPGetPC(subKsp[b], &subPc), "KSPGetPC(sub)");
                    checkPetscError
                    (
                        PCSetType(subPc, subPcName.c_str()), "PCSetType(sub)"
                    );
                }
            }
            if (pcName == "ilu" || pcName == "lu")
            {
                checkPetscError
                (
                    PCFactorSetFill(pc, max(factorFill, scalar(1.0))),
                    "PCFactorSetFill"
                );
            }
            if (pcName == "ilu" && dropTolerance > 0.0)
            {
                checkPetscError
                (
                    PCFactorSetDropTolerance(pc, dropTolerance, dropTolerance, 1000),
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

        void updateValues(const SpMat& A, const bool reusePreconditioner = false)
        {
            if (!isInitialised())
            {
                FatalErrorInFunction
                    << "PetscKspMatrixSolver::updateValues requires reset()"
                    << exit(FatalError);
            }

            // Modified for cardiacFoam: n_ is the LOCAL row count, while
            // A.cols() is the GLOBAL column count - they are only equal in
            // serial. Comparing them aborted every diagonalIion run in
            // parallel ("320x640 but cached PETSc Mat is 320x320").
            if (A.rows() != n_ || A.cols() != gNCols())
            {
                FatalErrorInFunction
                    << "PetscKspMatrixSolver::updateValues called with "
                    << "matrix of size " << A.rows() << "x" << A.cols()
                    << " but expected " << n_ << "x" << gNCols()
                    << exit(FatalError);
            }

            checkPetscError(MatZeroEntries(A_), "MatZeroEntries(updateValues)");

            const PetscInt nRows = static_cast<PetscInt>(A.rows());
            for (PetscInt row = 0; row < nRows; ++row)
            {
                for (SpMat::InnerIterator it(A, row); it; ++it)
                {
                    checkPetscError
                    (
                        MatSetValue
                        (
                            A_,
                            // Modified for cardiacFoam: PETSc addresses rows
                            // GLOBALLY. it.row() is this rank's local index, so
                            // it needs the same gRowStart offset that
                            // buildDistributedMat() applies; without it every
                            // rank but the first wrote into the wrong rows.
                            static_cast<PetscInt>(gRowStart + it.row()),
                            static_cast<PetscInt>(it.col()),
                            static_cast<PetscScalar>(it.value()),
                            INSERT_VALUES
                        ),
                        "MatSetValue(updateValues)"
                    );
                }
            }

            checkPetscError(MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin(updateValues)");
            checkPetscError(MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd(updateValues)");
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

            checkPetscError(KSPGetIterationNumber(ksp_, &petscIterations), "KSPGetIterationNumber");
            checkPetscError(KSPGetResidualNorm(ksp_, &residualNorm), "KSPGetResidualNorm");
            checkPetscError(KSPGetConvergedReason(ksp_, &reason), "KSPGetConvergedReason");

            if (reason < 0)
            {
                FatalErrorInFunction
                    << "PETSc KSP diverged with reason " << reason
                    << exit(FatalError);
            }

            EigVec result(rhs.size());
            copyPetscVecToEigVec(x_, result);

            iterations = static_cast<label>(petscIterations);
            estimatedError =
                static_cast<scalar>(residualNorm)/max(gNorm(rhs), scalar(SMALL));

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

    // Payload for a PCShell: the preconditioner application and the local row
    // count. Lets an arbitrary C++ operator (here an ILU of the diffusion
    // matrix) act as a PETSc preconditioner.
    struct PetscShellPCContext
    {
        std::function<EigVec(const EigVec&)> apply;
        label n;
    };

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

    // Persistent PETSc KSP over a MATRIX-FREE operator, used for the JFNK inner
    // Krylov solve.
    //
    // Same lifetime argument as PetscKspMatrixSolver, but here there is no
    // matrix to reuse: what is worth keeping is the KSP, the work vectors and
    // the shell contexts. The preconditioner is optional and is itself a shell,
    // so JFNK can run unpreconditioned or against an ILU of AImplicit.
    class PetscShellKspSolver
    {
        Mat A_;
        KSP ksp_;
        Vec b_;
        Vec x_;
        PetscShellMatVecContext matContext_;
        PetscShellPCContext pcContext_;
        label n_;
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

            checkPetscError
            (
                // Modified for cardiacFoam: nP is the LOCAL row/col count; the global
                // sizes are left to PETSc to sum across ranks.
                MatCreateShell
                (
                    PETSC_COMM_WORLD, nP, nP,
                    PETSC_DETERMINE, PETSC_DETERMINE,
                    &matContext_, &A_
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

            // Modified for cardiacFoam: distributed vector, nP LOCAL entries.
            checkPetscError
            (
                VecCreateMPI(PETSC_COMM_WORLD, nP, PETSC_DETERMINE, &b_),
                "VecCreateMPI(cached shell b)"
            );
            checkPetscError(VecDuplicate(b_, &x_), "VecDuplicate(cached shell x)");

            checkPetscError(KSPCreate(PETSC_COMM_WORLD, &ksp_), "KSPCreate(cached shell)");
            checkPetscError(KSPSetOperators(ksp_, A_, A_), "KSPSetOperators(cached shell)");

            const std::string kspName = petscKspTypeName(kspType);
            checkPetscError(KSPSetType(ksp_, kspName.c_str()), "KSPSetType(cached shell)");
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
                checkPetscError(KSPSetFromOptions(ksp_), "KSPSetFromOptions(cached shell)");
                checkPetscError(KSPGetPC(ksp_, &pc), "KSPGetPC(cached shell, after options)");
            }

            if (withShellPC)
            {
                checkPetscError(PCSetType(pc, PCSHELL), "PCSetType(cached shell PCSHELL)");
                checkPetscError(PCShellSetContext(pc, &pcContext_), "PCShellSetContext(cached)");
                checkPetscError(PCShellSetApply(pc, petscShellPCApply), "PCShellSetApply(cached)");
                // Right preconditioning so the KSP monitors the true residual
                // ||b - A x|| (matching the Eigen solveGMRES path), rather than
                // the left-preconditioned residual ||M^{-1}(b - A x)||.
                checkPetscError(KSPSetPCSide(ksp_, PC_RIGHT), "KSPSetPCSide(cached shell)");
                hasShellPC_ = true;
            }
            else
            {
                // Modified for cardiacFoam: see the note on bjacobi above.
                std::string subPcNameShell;
                const std::string pcName =
                    petscParallelPcTypeName(pcType, subPcNameShell);
                checkPetscError(PCSetType(pc, pcName.c_str()), "PCSetType(cached shell scalar PC)");
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
                    << "PetscShellKspSolver::solve shell-PC configuration changed"
                    << exit(FatalError);
            }

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

            checkPetscError(KSPGetIterationNumber(ksp_, &petscIterations), "KSPGetIterationNumber(cached shell)");
            checkPetscError(KSPGetResidualNorm(ksp_, &residualNorm), "KSPGetResidualNorm(cached shell)");
            checkPetscError(KSPGetConvergedReason(ksp_, &reason), "KSPGetConvergedReason(cached shell)");

            if (reason < 0)
            {
                FatalErrorInFunction
                    << "PETSc cached shell KSP diverged with reason " << reason
                    << exit(FatalError);
            }

            EigVec result(rhs.size());
            copyPetscVecToEigVec(x_, result);

            iterations = static_cast<label>(petscIterations);
            estimatedError =
                static_cast<scalar>(residualNorm)/max(gNorm(rhs), scalar(SMALL));

            return result;
        }
    };


    // ----------------------------------------------------------------- //
    // Embedded ten Tusscher--Noble--Noble--Panfilov 2004 ionic model.
    //
    // The original electro solver used cardiacFoam's runtime-selected
    // ionicModel/TNNP classes from src/.  This solver keeps the same TNNP
    // equations locally so the numerical PDE/source coupling can be kept
    // aligned with the manufactured JFNK solver without linking against
    // cardiacFoam-pc/src.  Units follow the generated CellML code:
    // Vm is stored internally in mV and time in ms.  Since 1 mV/ms = 1 V/s,
    // Iion_cm can be inserted directly into the monodomain source used by
    // the PDE, whose Vm field is in V.
    // ----------------------------------------------------------------- //

enum CONSTANTS_INDEX{
    R,              // 0  : gas constant
    T,              // 1  : temperature
    F,              // 2  : Faraday constant
    Cm,             // 3  : membrane capacitance (pF)
    V_c,            // 4  : cell volume (um^3)

    // Stimulus (S1)
    stim_start,         // 5
    stim_period_S1,     // 6
    stim_duration,      // 7
    stim_amplitude,     // 8

    // Reversal potentials & ion concentrations
    P_kna,          // 9
    K_o,            // 10
    Na_o,           // 11
    Ca_o,           // 12

    // Conductances
    g_K1,           // 13
    g_Kr,           // 14
    g_Ks,           // 15
    g_Na,           // 16
    g_bNa,          // 17
    g_CaL,          // 18
    g_bCa,          // 19
    g_to,           // 20

    // Na/K Pump
    P_NaK,          // 21
    K_mk,           // 22
    K_mNa,          // 23

    // NCX exchanger
    K_NaCa,         // 24
    K_sat,          // 25
    alpha_NCX,      // 26
    gamma_NCX,      // 27
    Km_Ca,          // 28
    Km_Nai,         // 29

    // Calcium pump + potassium pump
    g_pCa,          // 30
    K_pCa,          // 31
    g_pK,           // 32

    // Calcium dynamics (SR release / uptake)
    tau_g,          // 33
    a_rel,          // 34
    b_rel,          // 35
    c_rel,          // 36
    K_up,           // 37
    V_leak,         // 38
    Vmax_up,        // 39

    // Buffers
    Buf_c,          // 40
    K_buf_c,        // 41
    Buf_sr,         // 42
    K_buf_sr,       // 43
    V_sr,           // 44
    tau_fCa,        // 45

    // Added multi-stimulus parameters (S1/S2 protocol)
    nstim1,         // 46
    stim_period_S2, // 47
    nstim2,         // 48

    NUM_CONSTANTS
};


enum STATES_INDEX {
    V,      // membrane voltage
    K_i,    // intracellular potassium
    Na_i,   // intracellular sodium
    Ca_i,   // intracellular calcium
    Xr1,
    Xr2,
    Xs,
    m,
    h,
    j,
    d,
    f,
    fCa,
    s,
    r,
    Ca_SR,
    g,
    NUM_STATES
};


enum ALGEBRAIC_INDEX {
    Istim,
    xr1_inf,
    xr2_inf,
    xs_inf,
    m_inf,
    h_inf,
    j_inf,
    d_inf,
    f_inf,
    alpha_fCa,
    s_inf,
    r_inf,
    g_inf,
    E_Na,
    alpha_xr1,
    alpha_xr2,
    alpha_xs,
    alpha_m,
    alpha_h,
    alpha_j,
    alpha_d,
    tau_f,
    beta_fCa,
    tau_s,
    tau_r,
    d_g,
    E_K,
    beta_xr1,
    beta_xr2,
    beta_xs,
    beta_m,
    beta_h,
    beta_j,
    gama_fCa,
    beta_d,
    E_Ks,
    tau_xr1,
    tau_xr2,
    tau_xs,
    tau_m,
    tau_h,
    tau_j,
    gamma_d,
    fCa_inf,
    E_Ca,
    tau_d,
    d_fCa,
    alpha_K1,
    beta_K1,
    xK1_inf,
    i_K1,
    i_Kr,
    i_Ks,
    i_Na,
    i_b_Na,
    i_CaL,
    i_b_Ca,
    i_to,
    i_NaK,
    i_NaCa,
    i_p_Ca,
    i_p_K,
    i_rel,
    i_up,
    i_leak,
    ddt_Ca_i_total,
    ddt_Ca_sr_total,
    f_JCa_i_free,
    f_JCa_sr_free,
    Iion_cm,

    NUM_ALGEBRAIC
};

static const char* TNNP_STATES_NAMES[17] = {
    "V",        // 0: membrane voltage (millivolt)
    "K_i",      // 1: intracellular potassium (millimolar)
    "Na_i",     // 2: intracellular sodium (millimolar)
    "Ca_i",     // 3: intracellular calcium (millimolar)
    "Xr1",      // 4: Xr1 gate (dimensionless)
    "Xr2",      // 5: Xr2 gate (dimensionless)
    "Xs",       // 6: Xs gate (dimensionless)
    "m",        // 7: m gate (dimensionless)
    "h",        // 8: h gate (dimensionless)
    "j",        // 9: j gate (dimensionless)
    "d",        // 10: d gate (dimensionless)
    "f",        // 11: f gate (dimensionless)
    "fCa",      // 12: fCa gate (dimensionless)
    "s",        // 13: s gate (dimensionless)
    "r",        // 14: r gate (dimensionless)
    "Ca_SR",    // 15: sarcoplasmic reticulum calcium (millimolar)
    "g"         // 16: g gate (dimensionless)
};
static const char* TNNP_ALGEBRAIC_NAMES[70] = {
    "Istim",                    // 0
    "xr1_inf",                  // 1
    "xr2_inf",                  // 2
    "xs_inf",                   // 3
    "m_inf",                    // 4
    "h_inf",                    // 5
    "j_inf",                    // 6
    "d_inf",                    // 7
    "f_inf",                    // 8
    "alpha_fCa",               // 9
    "s_inf",                    // 10
    "r_inf",                    // 11
    "g_inf",                    // 12
    "E_Na",                     // 13
    "alpha_xr1",               // 14
    "alpha_xr2",               // 15
    "alpha_xs",                // 16
    "alpha_m",                 // 17
    "alpha_h",                 // 18
    "alpha_j",                 // 19
    "alpha_d",                 // 20
    "tau_f",                    // 21
    "beta_fCa",                // 22
    "tau_s",                    // 23
    "tau_r",                    // 24
    "d_g",                      // 25
    "E_K",                      // 26
    "beta_xr1",                // 27
    "beta_xr2",                // 28
    "beta_xs",                 // 29
    "beta_m",                  // 30
    "beta_h",                  // 31
    "beta_j",                  // 32
    "gama_fCa",                // 33
    "beta_d",                  // 34
    "E_Ks",                     // 35
    "tau_xr1",                 // 36
    "tau_xr2",                 // 37
    "tau_xs",                  // 38
    "tau_m",                   // 39
    "tau_h",                   // 40
    "tau_j",                   // 41
    "gamma_d",                 // 42
    "fCa_inf",                 // 43
    "E_Ca",                     // 44
    "tau_d",                   // 45
    "d_fCa",                   // 46
    "alpha_K1",                // 47
    "beta_K1",                 // 48
    "xK1_inf",                 // 49
    "i_K1",                     // 50
    "i_Kr",                     // 51
    "i_Ks",                     // 52
    "i_Na",                     // 53
    "i_b_Na",                   // 54
    "i_CaL",                    // 55
    "i_b_Ca",                   // 56
    "i_to",                     // 57
    "i_NaK",                    // 58
    "i_NaCa",                   // 59
    "i_p_Ca",                   // 60
    "i_p_K",                    // 61
    "i_rel",                    // 62
    "i_up",                     // 63
    "i_leak",                   // 64
    "ddt_Ca_i_total",           // 65
    "ddt_Ca_sr_total",          // 66
    "f_JCa_i_free",             // 67
    "f_JCa_sr_free",            // 68
    "Iion_cm"                   // 69
};

    // S1/S2 pacing waveform of the embedded TNNP kernel: returns the stimulus
    // current at time VOI as a train of nStim1 pulses of period stimPeriodS1
    // followed by nStim2 of period stimPeriodS2, each of stimDuration and
    // amplitude stimAmplitude. This is the CELL-level stimulus of the ionic
    // model, distinct from applyStimulus() below, which is the TISSUE-level one
    // applied to a region of the mesh.
    scalar embeddedComputeStimulus
    (
        const scalar VOI,
        const scalar stimStart,
        const scalar stimPeriodS1,
        const scalar stimDuration,
        const scalar stimAmplitude,
        const label nStim1,
        const scalar stimPeriodS2,
        const label nStim2
    )
    {
        scalar Istim = 0.0;

        if (stimPeriodS1 > 0.0 && nStim1 > 0)
        {
            const scalar tp = VOI - stimStart;
            if (tp >= 0.0 && tp <= stimPeriodS1*nStim1)
            {
                const scalar phase = tp - std::floor(tp/stimPeriodS1)*stimPeriodS1;
                if (phase >= 0.0 && phase <= stimDuration)
                {
                    Istim = -stimAmplitude;
                }
            }
        }

        if (Istim == 0.0 && stimPeriodS2 > 0.0 && nStim2 > 0)
        {
            const scalar tS1End = stimStart + stimPeriodS1*nStim1;
            const scalar tp = VOI - tS1End;
            if (tp >= 0.0 && tp <= stimPeriodS2*nStim2)
            {
                const scalar phase = tp - std::floor(tp/stimPeriodS2)*stimPeriodS2;
                if (phase >= 0.0 && phase <= stimDuration)
                {
                    Istim = -stimAmplitude;
                }
            }
        }

        return Istim;
    }

    inline Foam::scalar computeIstim(Foam::scalar t, const double* C)
    {
        return embeddedComputeStimulus
        (
            t,
            C[stim_start],
            C[stim_period_S1],
            C[stim_duration],
            C[stim_amplitude],
            Foam::label(C[nstim1]),
            C[stim_period_S2],
            Foam::label(C[nstim2])
        );
    }

void
TNNPinitConsts(double* CONSTANTS, double* RATES, double *STATES, int tissueFlag, const Foam::dictionary& stimulus)

{
STATES[V] = -86.2;
CONSTANTS[R] = 8.314;
CONSTANTS[T] = 310;
CONSTANTS[F] = 96.485;
CONSTANTS[Cm] = 185;
CONSTANTS[V_c] = 16404;


CONSTANTS[P_kna] = 0.03;
CONSTANTS[K_o] = 5.4;
CONSTANTS[Na_o] = 140;
STATES[K_i] = 138.3;
STATES[Na_i] = 11.6;
CONSTANTS[Ca_o] = 2;
STATES[Ca_i] = 0.0002;
CONSTANTS[g_K1] = 5.405;
CONSTANTS[g_Kr] = 0.096;
STATES[Xr1] = 0;
STATES[Xr2] = 1;
//Conditional tissue conductance Gks
CONSTANTS[g_Ks] = (tissueFlag == 2)
          ? 0.062
          : 0.245;

STATES[Xs] = 0;
CONSTANTS[g_Na] = 14.838;
STATES[m] = 0;
STATES[h] = 0.75;
STATES[j] = 0.75;
CONSTANTS[g_bNa] = 0.00029;
CONSTANTS[g_CaL] = 0.175;
STATES[d] = 0;
STATES[f] = 1;
STATES[fCa] = 1;
CONSTANTS[g_bCa] = 0.000592;
//Conditional tissue conductance Gto
CONSTANTS[g_to] = (tissueFlag == 3)
          ? 0.073
          : 0.294;

STATES[s] = 1;
STATES[r] = 0;
CONSTANTS[P_NaK] = 1.362;
CONSTANTS[K_mk] = 1;
CONSTANTS[K_mNa] = 40;
CONSTANTS[K_NaCa] = 1000;
CONSTANTS[K_sat] = 0.1;
CONSTANTS[alpha_NCX] = 2.5;
CONSTANTS[gamma_NCX] = 0.35;
CONSTANTS[Km_Ca] = 1.38;
CONSTANTS[Km_Nai] = 87.5;
CONSTANTS[g_pCa] = 0.825;
CONSTANTS[K_pCa] = 0.0005;
CONSTANTS[g_pK] = 0.0146;
STATES[Ca_SR] = 0.2;
STATES[g] = 1;
CONSTANTS[tau_g] = 2;
CONSTANTS[a_rel] = 0.016464;
CONSTANTS[b_rel] = 0.25;
CONSTANTS[c_rel] = 0.008232;
CONSTANTS[K_up] = 0.00025;
CONSTANTS[V_leak] = 8e-5;
CONSTANTS[Vmax_up] = 0.000425;
CONSTANTS[Buf_c] = 0.15;
CONSTANTS[K_buf_c] = 0.001;
CONSTANTS[Buf_sr] = 10;
CONSTANTS[K_buf_sr] = 0.3;
CONSTANTS[V_sr] = 1094;
CONSTANTS[tau_fCa] = 2.00000;

}

void
TNNPcomputeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag, bool solveVmWithinODESolver)
{

//INa
ALGEBRAIC[m_inf] = 1.00000/std::pow(1.00000+std::exp((- 56.8600 - STATES[V])/9.03000), 2.00000);
ALGEBRAIC[alpha_m] = 1.00000/(1.00000+std::exp((- 60.0000 - STATES[V])/5.00000));
ALGEBRAIC[beta_m] = 0.100000/(1.00000+std::exp((STATES[V]+35.0000)/5.00000))+0.100000/(1.00000+std::exp((STATES[V] - 50.0000)/200.000));
ALGEBRAIC[tau_m] =  1.00000*ALGEBRAIC[alpha_m]*ALGEBRAIC[beta_m];
RATES[m] = (ALGEBRAIC[m_inf] - STATES[m])/ALGEBRAIC[tau_m];
ALGEBRAIC[h_inf] = 1.00000/std::pow(1.00000+std::exp((STATES[V]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[alpha_h] = (STATES[V]<- 40.0000 ?  0.0570000*std::exp(- (STATES[V]+80.0000)/6.80000) : 0.00000);
ALGEBRAIC[beta_h] = (STATES[V]<- 40.0000 ?  2.70000*std::exp( 0.0790000*STATES[V])+ 310000.*std::exp( 0.348500*STATES[V]) : 0.770000/( 0.130000*(1.00000+std::exp((STATES[V]+10.6600)/- 11.1000))));
ALGEBRAIC[tau_h] = 1.00000/(ALGEBRAIC[alpha_h]+ALGEBRAIC[beta_h]);
RATES[h] = (ALGEBRAIC[h_inf] - STATES[h])/ALGEBRAIC[tau_h];
ALGEBRAIC[j_inf] = 1.00000/std::pow(1.00000+std::exp((STATES[V]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[alpha_j] = (STATES[V]<- 40.0000 ? (( ( - 25428.0*std::exp( 0.244400*STATES[V]) -  6.94800e-06*std::exp( - 0.0439100*STATES[V]))*(STATES[V]+37.7800))/1.00000)/(1.00000+std::exp( 0.311000*(STATES[V]+79.2300))) : 0.00000);
ALGEBRAIC[beta_j] = (STATES[V]<- 40.0000 ? ( 0.0242400*std::exp( - 0.0105200*STATES[V]))/(1.00000+std::exp( - 0.137800*(STATES[V]+40.1400))) : ( 0.600000*std::exp( 0.0570000*STATES[V]))/(1.00000+std::exp( - 0.100000*(STATES[V]+32.0000))));
ALGEBRAIC[tau_j] = 1.00000/(ALGEBRAIC[alpha_j]+ALGEBRAIC[beta_j]);
RATES[j] = (ALGEBRAIC[j_inf] - STATES[j])/ALGEBRAIC[tau_j];

//ICaL
ALGEBRAIC[d_inf] = 1.00000/(1.00000+std::exp((- 5.00000 - STATES[V])/7.50000));
ALGEBRAIC[alpha_d] = 1.40000/(1.00000+std::exp((- 35.0000 - STATES[V])/13.0000))+0.250000;
ALGEBRAIC[gama_fCa] = 1.40000/(1.00000+std::exp((STATES[V]+5.00000)/5.00000));
ALGEBRAIC[gamma_d] = 1.00000/(1.00000+std::exp((50.0000 - STATES[V])/20.0000));
ALGEBRAIC[tau_d] =  1.00000*ALGEBRAIC[alpha_d]*ALGEBRAIC[gama_fCa]+ALGEBRAIC[gamma_d];
RATES[d] = (ALGEBRAIC[d_inf] - STATES[d])/ALGEBRAIC[tau_d];
ALGEBRAIC[f_inf] = 1.00000/(1.00000+std::exp((STATES[V]+20.0000)/7.00000));
ALGEBRAIC[tau_f] =  1125.00*std::exp(- std::pow(STATES[V]+27.0000, 2.00000)/240.000)+80.0000+165.000/(1.00000+std::exp((25.0000 - STATES[V])/10.0000));
RATES[f] = (ALGEBRAIC[f_inf] - STATES[f])/ALGEBRAIC[tau_f];

ALGEBRAIC[alpha_fCa] = 1.00000/(1.00000+std::pow(STATES[Ca_i]/0.000325000, 8.00000));
ALGEBRAIC[beta_fCa] = 0.100000/(1.00000+std::exp((STATES[Ca_i] - 0.000500000)/0.000100000));
ALGEBRAIC[beta_d] = 0.200000/(1.00000+std::exp((STATES[Ca_i] - 0.000750000)/0.000800000));
ALGEBRAIC[fCa_inf] = (ALGEBRAIC[alpha_fCa]+ALGEBRAIC[beta_fCa]+ALGEBRAIC[beta_d]+0.230000)/1.46000;
ALGEBRAIC[d_fCa] = (ALGEBRAIC[fCa_inf] - STATES[fCa])/CONSTANTS[tau_fCa];
RATES[fCa] = (ALGEBRAIC[fCa_inf]>STATES[fCa]&&STATES[V]>- 60.0000 ? 0.00000 : ALGEBRAIC[d_fCa]);

//Ito
//Conditional tissue flag for the Ito current
ALGEBRAIC[s_inf] = (tissueFlag == 1)
    ? 1.00000 / (1.00000 + std::exp((STATES[V] + 20.0000) / 5.00000))
    : 1.10000 / (1.00000 + std::exp((STATES[V] + 28.0000) / 6.00000));
ALGEBRAIC[tau_s] = (tissueFlag == 1)
     ? 85.0000*std::exp(- std::pow(STATES[V]+45.0000, 2.00000)/320.000)+5.00000/(1.00000+std::exp((STATES[V] - 20.0000)/5.00000))+3.00000
     : 1000.0000*std::exp(- std::pow(STATES[V]+67.0000, 2.00000)/1000.000)+8.00000;
RATES[s] = (ALGEBRAIC[s_inf] - STATES[s])/ALGEBRAIC[tau_s];
ALGEBRAIC[r_inf] = 1.00000/(1.00000+std::exp((20.0000 - STATES[V])/6.00000));
ALGEBRAIC[tau_r] =  9.50000*std::exp(- std::pow(STATES[V]+40.0000, 2.00000)/1800.00)+0.800000;
RATES[r] = (ALGEBRAIC[r_inf] - STATES[r])/ALGEBRAIC[tau_r];

//IKr
ALGEBRAIC[xr1_inf] = 1.00000/(1.00000+std::exp((- 26.0000 - STATES[V])/7.00000));
ALGEBRAIC[alpha_xr1] = 450.000/(1.00000+std::exp((- 45.0000 - STATES[V])/10.0000));
ALGEBRAIC[beta_xr1] = 6.00000/(1.00000+std::exp((STATES[V]+30.0000)/11.5000));
ALGEBRAIC[tau_xr1] =  1.00000*ALGEBRAIC[alpha_xr1]*ALGEBRAIC[beta_xr1];
RATES[Xr1] = (ALGEBRAIC[xr1_inf] - STATES[Xr1])/ALGEBRAIC[tau_xr1];
ALGEBRAIC[xr2_inf] = 1.00000/(1.00000+std::exp((STATES[V]+88.0000)/24.0000));
ALGEBRAIC[alpha_xr2] = 3.00000/(1.00000+std::exp((- 60.0000 - STATES[V])/20.0000));
ALGEBRAIC[beta_xr2] = 1.12000/(1.00000+std::exp((STATES[V] - 60.0000)/20.0000));
ALGEBRAIC[tau_xr2] =  1.00000*ALGEBRAIC[alpha_xr2]*ALGEBRAIC[beta_xr2];
RATES[Xr2] = (ALGEBRAIC[xr2_inf] - STATES[Xr2])/ALGEBRAIC[tau_xr2];

//IKs
ALGEBRAIC[xs_inf] = 1.00000/(1.00000+std::exp((- 5.00000 - STATES[V])/14.0000));
ALGEBRAIC[alpha_xs] = 1100.00/ std::pow((1.00000+std::exp((- 10.0000 - STATES[V])/6.00000)), 1.0 / 2);
ALGEBRAIC[beta_xs] = 1.00000/(1.00000+std::exp((STATES[V] - 60.0000)/20.0000));
ALGEBRAIC[tau_xs] =  1.00000*ALGEBRAIC[alpha_xs]*ALGEBRAIC[beta_xs];
RATES[Xs] = (ALGEBRAIC[xs_inf] - STATES[Xs])/ALGEBRAIC[tau_xs];

//SR calcium RyR gate g
ALGEBRAIC[g_inf] = (STATES[Ca_i]<0.000350000 ? 1.00000/(1.00000+std::pow(STATES[Ca_i]/0.000350000, 6.00000)) : 1.00000/(1.00000+std::pow(STATES[Ca_i]/0.000350000, 16.0000)));
ALGEBRAIC[d_g] = (ALGEBRAIC[g_inf] - STATES[g])/CONSTANTS[tau_g];
RATES[g] = (ALGEBRAIC[g_inf]>STATES[g]&&STATES[V]>- 60.0000 ? 0.00000 : ALGEBRAIC[d_g]);


ALGEBRAIC[i_NaK] = (( (( CONSTANTS[P_NaK]*CONSTANTS[K_o])/(CONSTANTS[K_o]+CONSTANTS[K_mk]))*STATES[Na_i])/(STATES[Na_i]+CONSTANTS[K_mNa]))/(1.00000+ 0.124500*std::exp(( - 0.100000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))+ 0.0353000*std::exp(( - STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])));
ALGEBRAIC[E_Na] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log(CONSTANTS[Na_o]/STATES[Na_i]);
ALGEBRAIC[i_Na] =  CONSTANTS[g_Na]*std::pow(STATES[m], 3.00000)*STATES[h]*STATES[j]*(STATES[V] - ALGEBRAIC[E_Na]);
ALGEBRAIC[i_b_Na] =  CONSTANTS[g_bNa]*(STATES[V] - ALGEBRAIC[E_Na]);
ALGEBRAIC[i_NaCa] = ( CONSTANTS[K_NaCa]*( std::exp(( CONSTANTS[gamma_NCX]*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))*std::pow(STATES[Na_i], 3.00000)*CONSTANTS[Ca_o] -  std::exp(( (CONSTANTS[gamma_NCX] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))*std::pow(CONSTANTS[Na_o], 3.00000)*STATES[Ca_i]*CONSTANTS[alpha_NCX]))/( (std::pow(CONSTANTS[Km_Nai], 3.00000)+std::pow(CONSTANTS[Na_o], 3.00000))*(CONSTANTS[Km_Ca]+CONSTANTS[Ca_o])*(1.00000+ CONSTANTS[K_sat]*std::exp(( (CONSTANTS[gamma_NCX] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))));
RATES[Na_i] = ( - (ALGEBRAIC[i_Na]+ALGEBRAIC[i_b_Na]+ 3.00000*ALGEBRAIC[i_NaK]+ 3.00000*ALGEBRAIC[i_NaCa])*CONSTANTS[Cm])/( CONSTANTS[V_c]*CONSTANTS[F]);

ALGEBRAIC[E_K] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log(CONSTANTS[K_o]/STATES[K_i]);
ALGEBRAIC[alpha_K1] = 0.100000/(1.00000+std::exp( 0.0600000*((STATES[V] - ALGEBRAIC[E_K]) - 200.000)));
ALGEBRAIC[beta_K1] = ( 3.00000*std::exp( 0.000200000*((STATES[V] - ALGEBRAIC[E_K])+100.000))+ 1.00000*std::exp( 0.100000*((STATES[V] - ALGEBRAIC[E_K]) - 10.0000)))/(1.00000+std::exp( - 0.500000*(STATES[V] - ALGEBRAIC[E_K])));
ALGEBRAIC[xK1_inf] = ALGEBRAIC[alpha_K1]/(ALGEBRAIC[alpha_K1]+ALGEBRAIC[beta_K1]);
ALGEBRAIC[i_K1] =  CONSTANTS[g_K1]*ALGEBRAIC[xK1_inf]* std::pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*(STATES[V] - ALGEBRAIC[E_K]);

ALGEBRAIC[i_to] =  CONSTANTS[g_to]*STATES[r]*STATES[s]*(STATES[V] - ALGEBRAIC[E_K]);

ALGEBRAIC[i_Kr] =  CONSTANTS[g_Kr]* std::pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*STATES[Xr1]*STATES[Xr2]*(STATES[V] - ALGEBRAIC[E_K]);
ALGEBRAIC[E_Ks] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log((CONSTANTS[K_o]+ CONSTANTS[P_kna]*CONSTANTS[Na_o])/(STATES[K_i]+ CONSTANTS[P_kna]*STATES[Na_i]));
ALGEBRAIC[i_Ks] =  CONSTANTS[g_Ks]*std::pow(STATES[Xs], 2.00000)*(STATES[V] - ALGEBRAIC[E_Ks]);
ALGEBRAIC[i_CaL] = ( (( CONSTANTS[g_CaL]*STATES[d]*STATES[f]*STATES[fCa]*4.00000*STATES[V]*std::pow(CONSTANTS[F], 2.00000))/( CONSTANTS[R]*CONSTANTS[T]))*( STATES[Ca_i]*std::exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])) -  0.341000*CONSTANTS[Ca_o]))/(std::exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])) - 1.00000);
ALGEBRAIC[E_Ca] =  (( 0.500000*CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log(CONSTANTS[Ca_o]/STATES[Ca_i]);
ALGEBRAIC[i_b_Ca] =  CONSTANTS[g_bCa]*(STATES[V] - ALGEBRAIC[E_Ca]);
ALGEBRAIC[i_p_K] = ( CONSTANTS[g_pK]*(STATES[V] - ALGEBRAIC[E_K]))/(1.00000+std::exp((25.0000 - STATES[V])/5.98000));
ALGEBRAIC[i_p_Ca] = ( CONSTANTS[g_pCa]*STATES[Ca_i])/(STATES[Ca_i]+CONSTANTS[K_pCa]);

ALGEBRAIC[Iion_cm] = ALGEBRAIC[i_K1]+ALGEBRAIC[i_to]+ALGEBRAIC[i_Kr]+
                ALGEBRAIC[i_Ks]+ALGEBRAIC[i_CaL]+ALGEBRAIC[i_NaK]+
                ALGEBRAIC[i_Na]+ALGEBRAIC[i_b_Na]+ALGEBRAIC[i_NaCa]+
                ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_K]+ALGEBRAIC[i_p_Ca];

ALGEBRAIC[Istim] = computeIstim(VOI, CONSTANTS);

RATES[V] = 0.0;
if (solveVmWithinODESolver)
    {
        RATES[V] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }

RATES[K_i] = ( - ((ALGEBRAIC[i_K1]+ALGEBRAIC[i_to]+ALGEBRAIC[i_Kr]+ALGEBRAIC[i_Ks]+ALGEBRAIC[i_p_K]+ALGEBRAIC[Istim]) -  2.00000*ALGEBRAIC[i_NaK])*CONSTANTS[Cm])/( CONSTANTS[V_c]*CONSTANTS[F]);
ALGEBRAIC[i_rel] =  (( CONSTANTS[a_rel]*std::pow(STATES[Ca_SR], 2.00000))/(std::pow(CONSTANTS[b_rel], 2.00000)+std::pow(STATES[Ca_SR], 2.00000))+CONSTANTS[c_rel])*STATES[d]*STATES[g];
ALGEBRAIC[i_up] = CONSTANTS[Vmax_up]/(1.00000+std::pow(CONSTANTS[K_up], 2.00000)/std::pow(STATES[Ca_i], 2.00000));
ALGEBRAIC[i_leak] =  CONSTANTS[V_leak]*(STATES[Ca_SR] - STATES[Ca_i]);
ALGEBRAIC[ddt_Ca_i_total] = (( (- ((ALGEBRAIC[i_CaL]+ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_Ca]) -  2.00000*ALGEBRAIC[i_NaCa])/( 2.00000*CONSTANTS[V_c]*CONSTANTS[F]))*CONSTANTS[Cm]+ALGEBRAIC[i_leak]) - ALGEBRAIC[i_up])+ALGEBRAIC[i_rel];
ALGEBRAIC[f_JCa_i_free] = 1.00000/(1.00000+( CONSTANTS[Buf_c]*CONSTANTS[K_buf_c])/std::pow(STATES[Ca_i]+CONSTANTS[K_buf_c], 2.00000));
RATES[Ca_i] =  ALGEBRAIC[ddt_Ca_i_total]*ALGEBRAIC[f_JCa_i_free];
ALGEBRAIC[ddt_Ca_sr_total] =  (CONSTANTS[V_c]/CONSTANTS[V_sr])*(ALGEBRAIC[i_up] - (ALGEBRAIC[i_rel]+ALGEBRAIC[i_leak]));
ALGEBRAIC[f_JCa_sr_free] = 1.00000/(1.00000+( CONSTANTS[Buf_sr]*CONSTANTS[K_buf_sr])/std::pow(STATES[Ca_SR]+CONSTANTS[K_buf_sr], 2.00000));
RATES[Ca_SR] =  ALGEBRAIC[ddt_Ca_sr_total]*ALGEBRAIC[f_JCa_sr_free];
}





void
TNNPcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver)
{

ALGEBRAIC[f_inf] = 1.00000/(1.00000+std::exp((STATES[V]+20.0000)/7.00000));
ALGEBRAIC[tau_f] =  1125.00*std::exp(- std::pow(STATES[V]+27.0000, 2.00000)/240.000)+80.0000+165.000/(1.00000+std::exp((25.0000 - STATES[V])/10.0000));

//Conditional tissue flag for the Ito current
ALGEBRAIC[s_inf] = (tissueFlag == 1)
    ? 1.00000 / (1.00000 + std::exp((STATES[V] + 20.0000) / 5.00000))
    : 1.10000 / (1.00000 + std::exp((STATES[V] + 28.0000) / 6.00000));

ALGEBRAIC[tau_s] = (tissueFlag == 1)
     ? 85.0000*std::exp(- std::pow(STATES[V]+45.0000, 2.00000)/320.000)+5.00000/(1.00000+std::exp((STATES[V] - 20.0000)/5.00000))+3.00000
     : 1000.0000*std::exp(- std::pow(STATES[V]+67.0000, 2.00000)/1000.000)+8.00000;

ALGEBRAIC[r_inf] = 1.00000/(1.00000+std::exp((20.0000 - STATES[V])/6.00000));
ALGEBRAIC[tau_r] =  9.50000*std::exp(- std::pow(STATES[V]+40.0000, 2.00000)/1800.00)+0.800000;
ALGEBRAIC[g_inf] = (STATES[Ca_i]<0.000350000 ? 1.00000/(1.00000+std::pow(STATES[Ca_i]/0.000350000, 6.00000)) : 1.00000/(1.00000+std::pow(STATES[Ca_i]/0.000350000, 16.0000)));
ALGEBRAIC[d_g] = (ALGEBRAIC[g_inf] - STATES[g])/CONSTANTS[tau_g];
ALGEBRAIC[xr1_inf] = 1.00000/(1.00000+std::exp((- 26.0000 - STATES[V])/7.00000));
ALGEBRAIC[alpha_xr1] = 450.000/(1.00000+std::exp((- 45.0000 - STATES[V])/10.0000));
ALGEBRAIC[beta_xr1] = 6.00000/(1.00000+std::exp((STATES[V]+30.0000)/11.5000));
ALGEBRAIC[tau_xr1] =  1.00000*ALGEBRAIC[alpha_xr1]*ALGEBRAIC[beta_xr1];
ALGEBRAIC[xr2_inf] = 1.00000/(1.00000+std::exp((STATES[V]+88.0000)/24.0000));
ALGEBRAIC[alpha_xr2] = 3.00000/(1.00000+std::exp((- 60.0000 - STATES[V])/20.0000));
ALGEBRAIC[beta_xr2] = 1.12000/(1.00000+std::exp((STATES[V] - 60.0000)/20.0000));
ALGEBRAIC[tau_xr2] =  1.00000*ALGEBRAIC[alpha_xr2]*ALGEBRAIC[beta_xr2];
ALGEBRAIC[xs_inf] = 1.00000/(1.00000+std::exp((- 5.00000 - STATES[V])/14.0000));
ALGEBRAIC[alpha_xs] = 1100.00/ std::pow((1.00000+std::exp((- 10.0000 - STATES[V])/6.00000)), 1.0 / 2);
ALGEBRAIC[beta_xs] = 1.00000/(1.00000+std::exp((STATES[V] - 60.0000)/20.0000));
ALGEBRAIC[tau_xs] =  1.00000*ALGEBRAIC[alpha_xs]*ALGEBRAIC[beta_xs];
ALGEBRAIC[m_inf] = 1.00000/std::pow(1.00000+std::exp((- 56.8600 - STATES[V])/9.03000), 2.00000);
ALGEBRAIC[alpha_m] = 1.00000/(1.00000+std::exp((- 60.0000 - STATES[V])/5.00000));
ALGEBRAIC[beta_m] = 0.100000/(1.00000+std::exp((STATES[V]+35.0000)/5.00000))+0.100000/(1.00000+std::exp((STATES[V] - 50.0000)/200.000));
ALGEBRAIC[tau_m] =  1.00000*ALGEBRAIC[alpha_m]*ALGEBRAIC[beta_m];
ALGEBRAIC[h_inf] = 1.00000/std::pow(1.00000+std::exp((STATES[V]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[alpha_h] = (STATES[V]<- 40.0000 ?  0.0570000*std::exp(- (STATES[V]+80.0000)/6.80000) : 0.00000);
ALGEBRAIC[beta_h] = (STATES[V]<- 40.0000 ?  2.70000*std::exp( 0.0790000*STATES[V])+ 310000.*std::exp( 0.348500*STATES[V]) : 0.770000/( 0.130000*(1.00000+std::exp((STATES[V]+10.6600)/- 11.1000))));
ALGEBRAIC[tau_h] = 1.00000/(ALGEBRAIC[alpha_h]+ALGEBRAIC[beta_h]);
ALGEBRAIC[j_inf] = 1.00000/std::pow(1.00000+std::exp((STATES[V]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[alpha_j] = (STATES[V]<- 40.0000 ? (( ( - 25428.0*std::exp( 0.244400*STATES[V]) -  6.94800e-06*std::exp( - 0.0439100*STATES[V]))*(STATES[V]+37.7800))/1.00000)/(1.00000+std::exp( 0.311000*(STATES[V]+79.2300))) : 0.00000);
ALGEBRAIC[beta_j] = (STATES[V]<- 40.0000 ? ( 0.0242400*std::exp( - 0.0105200*STATES[V]))/(1.00000+std::exp( - 0.137800*(STATES[V]+40.1400))) : ( 0.600000*std::exp( 0.0570000*STATES[V]))/(1.00000+std::exp( - 0.100000*(STATES[V]+32.0000))));
ALGEBRAIC[tau_j] = 1.00000/(ALGEBRAIC[alpha_j]+ALGEBRAIC[beta_j]);
ALGEBRAIC[d_inf] = 1.00000/(1.00000+std::exp((- 5.00000 - STATES[V])/7.50000));
ALGEBRAIC[alpha_d] = 1.40000/(1.00000+std::exp((- 35.0000 - STATES[V])/13.0000))+0.250000;
ALGEBRAIC[gama_fCa] = 1.40000/(1.00000+std::exp((STATES[V]+5.00000)/5.00000));
ALGEBRAIC[gamma_d] = 1.00000/(1.00000+std::exp((50.0000 - STATES[V])/20.0000));
ALGEBRAIC[tau_d] =  1.00000*ALGEBRAIC[alpha_d]*ALGEBRAIC[gama_fCa]+ALGEBRAIC[gamma_d];
ALGEBRAIC[alpha_fCa] = 1.00000/(1.00000+std::pow(STATES[Ca_i]/0.000325000, 8.00000));
ALGEBRAIC[beta_fCa] = 0.100000/(1.00000+std::exp((STATES[Ca_i] - 0.000500000)/0.000100000));
ALGEBRAIC[beta_d] = 0.200000/(1.00000+std::exp((STATES[Ca_i] - 0.000750000)/0.000800000));
ALGEBRAIC[fCa_inf] = (ALGEBRAIC[alpha_fCa]+ALGEBRAIC[beta_fCa]+ALGEBRAIC[beta_d]+0.230000)/1.46000;
ALGEBRAIC[d_fCa] = (ALGEBRAIC[fCa_inf] - STATES[fCa])/CONSTANTS[tau_fCa];
ALGEBRAIC[i_NaK] = (( (( CONSTANTS[P_NaK]*CONSTANTS[K_o])/(CONSTANTS[K_o]+CONSTANTS[K_mk]))*STATES[Na_i])/(STATES[Na_i]+CONSTANTS[K_mNa]))/(1.00000+ 0.124500*std::exp(( - 0.100000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))+ 0.0353000*std::exp(( - STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])));
ALGEBRAIC[E_Na] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log(CONSTANTS[Na_o]/STATES[Na_i]);
ALGEBRAIC[i_Na] =  CONSTANTS[g_Na]*std::pow(STATES[m], 3.00000)*STATES[h]*STATES[j]*(STATES[V] - ALGEBRAIC[E_Na]);
ALGEBRAIC[i_b_Na] =  CONSTANTS[g_bNa]*(STATES[V] - ALGEBRAIC[E_Na]);
ALGEBRAIC[i_NaCa] = ( CONSTANTS[K_NaCa]*( std::exp(( CONSTANTS[gamma_NCX]*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))*std::pow(STATES[Na_i], 3.00000)*CONSTANTS[Ca_o] -  std::exp(( (CONSTANTS[gamma_NCX] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))*std::pow(CONSTANTS[Na_o], 3.00000)*STATES[Ca_i]*CONSTANTS[alpha_NCX]))/( (std::pow(CONSTANTS[Km_Nai], 3.00000)+std::pow(CONSTANTS[Na_o], 3.00000))*(CONSTANTS[Km_Ca]+CONSTANTS[Ca_o])*(1.00000+ CONSTANTS[K_sat]*std::exp(( (CONSTANTS[gamma_NCX] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))));
ALGEBRAIC[E_K] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log(CONSTANTS[K_o]/STATES[K_i]);
ALGEBRAIC[alpha_K1] = 0.100000/(1.00000+std::exp( 0.0600000*((STATES[V] - ALGEBRAIC[E_K]) - 200.000)));
ALGEBRAIC[beta_K1] = ( 3.00000*std::exp( 0.000200000*((STATES[V] - ALGEBRAIC[E_K])+100.000))+ 1.00000*std::exp( 0.100000*((STATES[V] - ALGEBRAIC[E_K]) - 10.0000)))/(1.00000+std::exp( - 0.500000*(STATES[V] - ALGEBRAIC[E_K])));
ALGEBRAIC[xK1_inf] = ALGEBRAIC[alpha_K1]/(ALGEBRAIC[alpha_K1]+ALGEBRAIC[beta_K1]);
ALGEBRAIC[i_K1] =  CONSTANTS[g_K1]*ALGEBRAIC[xK1_inf]* std::pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*(STATES[V] - ALGEBRAIC[E_K]);
ALGEBRAIC[i_to] =  CONSTANTS[g_to]*STATES[r]*STATES[s]*(STATES[V] - ALGEBRAIC[E_K]);
ALGEBRAIC[i_Kr] =  CONSTANTS[g_Kr]* std::pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*STATES[Xr1]*STATES[Xr2]*(STATES[V] - ALGEBRAIC[E_K]);
ALGEBRAIC[E_Ks] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log((CONSTANTS[K_o]+ CONSTANTS[P_kna]*CONSTANTS[Na_o])/(STATES[K_i]+ CONSTANTS[P_kna]*STATES[Na_i]));
ALGEBRAIC[i_Ks] =  CONSTANTS[g_Ks]*std::pow(STATES[Xs], 2.00000)*(STATES[V] - ALGEBRAIC[E_Ks]);
ALGEBRAIC[i_CaL] = ( (( CONSTANTS[g_CaL]*STATES[d]*STATES[f]*STATES[fCa]*4.00000*STATES[V]*std::pow(CONSTANTS[F], 2.00000))/( CONSTANTS[R]*CONSTANTS[T]))*( STATES[Ca_i]*std::exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])) -  0.341000*CONSTANTS[Ca_o]))/(std::exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])) - 1.00000);
ALGEBRAIC[E_Ca] =  (( 0.500000*CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log(CONSTANTS[Ca_o]/STATES[Ca_i]);
ALGEBRAIC[i_b_Ca] =  CONSTANTS[g_bCa]*(STATES[V] - ALGEBRAIC[E_Ca]);
ALGEBRAIC[i_p_K] = ( CONSTANTS[g_pK]*(STATES[V] - ALGEBRAIC[E_K]))/(1.00000+std::exp((25.0000 - STATES[V])/5.98000));
ALGEBRAIC[i_p_Ca] = ( CONSTANTS[g_pCa]*STATES[Ca_i])/(STATES[Ca_i]+CONSTANTS[K_pCa]);

ALGEBRAIC[i_rel] =  (( CONSTANTS[a_rel]*std::pow(STATES[Ca_SR], 2.00000))/(std::pow(CONSTANTS[b_rel], 2.00000)+std::pow(STATES[Ca_SR], 2.00000))+CONSTANTS[c_rel])*STATES[d]*STATES[g];
ALGEBRAIC[i_up] = CONSTANTS[Vmax_up]/(1.00000+std::pow(CONSTANTS[K_up], 2.00000)/std::pow(STATES[Ca_i], 2.00000));
ALGEBRAIC[i_leak] =  CONSTANTS[V_leak]*(STATES[Ca_SR] - STATES[Ca_i]);
ALGEBRAIC[ddt_Ca_i_total] = (( (- ((ALGEBRAIC[i_CaL]+ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_Ca]) -  2.00000*ALGEBRAIC[i_NaCa])/( 2.00000*CONSTANTS[V_c]*CONSTANTS[F]))*CONSTANTS[Cm]+ALGEBRAIC[i_leak]) - ALGEBRAIC[i_up])+ALGEBRAIC[i_rel];
ALGEBRAIC[f_JCa_i_free] = 1.00000/(1.00000+( CONSTANTS[Buf_c]*CONSTANTS[K_buf_c])/std::pow(STATES[Ca_i]+CONSTANTS[K_buf_c], 2.00000));
ALGEBRAIC[ddt_Ca_sr_total] =  (CONSTANTS[V_c]/CONSTANTS[V_sr])*(ALGEBRAIC[i_up] - (ALGEBRAIC[i_rel]+ALGEBRAIC[i_leak]));
ALGEBRAIC[f_JCa_sr_free] = 1.00000/(1.00000+( CONSTANTS[Buf_sr]*CONSTANTS[K_buf_sr])/std::pow(STATES[Ca_SR]+CONSTANTS[K_buf_sr], 2.00000));

ALGEBRAIC[Iion_cm] = ALGEBRAIC[i_K1]+ALGEBRAIC[i_to]+ALGEBRAIC[i_Kr]+
                ALGEBRAIC[i_Ks]+ALGEBRAIC[i_CaL]+ALGEBRAIC[i_NaK]+
                ALGEBRAIC[i_Na]+ALGEBRAIC[i_b_Na]+ALGEBRAIC[i_NaCa]+
                ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_K]+ALGEBRAIC[i_p_Ca];

ALGEBRAIC[Istim] = computeIstim(VOI, CONSTANTS);
}


    // ten Tusscher-Noble-Noble-Panfilov (2004) ionic model, embedded in this
    // translation unit rather than taken from the ionicModels library.
    //
    // It owns one state vector per INTEGRATION POINT, not per cell: the
    // high-order path evaluates Iion at the cell quadrature points, so the ODE
    // is integrated there. The scratch arrays cover two disjoint banks - Gauss
    // points and cell centres - because the front-aware hybrid needs both alive
    // at once; see solveODE() and the note on cellSlotOffset_ above.
    class EmbeddedTNNPModel
    {
        dictionary dict_;
        label nPoints_;

        // Added for cardiacFoam: number of scratch slots reserved for
        // cell-centred integration, and the index at which that second bank
        // starts. The per-point scratch below is a single array covering two
        // disjoint point sets: Gauss-point slots occupy [0, nPoints_) and
        // cell-centre slots occupy [cellSlotOffset_, cellSlotOffset_+nCells_).
        // See solveODE() for why the two banks must not overlap.
        label nCellSlots_;
        label cellSlotOffset_;

        scalarField constants_;
        Field<Field<scalar>> algebraic_;
        Field<Field<scalar>> rates_;
        Field<Field<scalar>> oldStates_;
        scalarList stepMs_;
        wordList exportedNames_;
        wordList debugNames_;
        label tissue_;
        word odeSolverName_;
        scalar odeInitialStepMs_;
        scalar odeAbsTol_;
        scalar odeRelTol_;
        label odeMaxSteps_;
        Switch numericalProtection_;
        Switch clampGates_;
        Switch clampVmForRates_;
        Switch logNumericalProtection_;
        scalar concentrationFloor_;
        scalar vmMinForRates_;
        scalar vmMaxForRates_;
        scalar vmZeroEps_;
        label maxProtectionLogEntries_;
        word protectionLogName_;
        label protectionCorrections_;

        // Added for cardiacFoam: OpenMP over the ionic ODEs. Every point is an
        // independent initial-value problem sharing no state, so the loop is
        // embarrassingly parallel - this is the same treatment the MMS solver
        // already applies in updateIntegrationPointStatesODE.
        //
        // Off by default: it interacts with the MPI rank count, and oversubscribing
        // a node is slower than either alone.
        Switch parallelODE_;
        label parallelODEThreshold_;

        static label tissueFlag(const word& tissueName)
        {
            if (tissueName == "epicardialCells") return 1;
            if (tissueName == "mCells") return 2;
            if (tissueName == "endocardialCells") return 3;
            if (tissueName == "myocyte") return 4;

            FatalErrorInFunction
                << "Unsupported TNNP tissue '" << tissueName << "'. Valid options are "
                << "epicardialCells, mCells, endocardialCells, myocyte."
                << exit(FatalError);

            return 1;
        }

        static label stateIndex(const word& name)
        {
            for (label i = 0; i < NUM_STATES; ++i)
            {
                if (name == TNNP_STATES_NAMES[i]) return i;
            }
            return -1;
        }

        static label algebraicIndex(const word& name)
        {
            for (label i = 0; i < NUM_ALGEBRAIC; ++i)
            {
                if (name == TNNP_ALGEBRAIC_NAMES[i]) return i;
            }
            return -1;
        }

        void loadStimulusProtocol()
        {
            constants_[stim_start] = dict_.lookupOrDefault<scalar>("stim_start", 0.0);
            constants_[stim_period_S1] = dict_.lookupOrDefault<scalar>("stim_period_S1", 0.0);
            constants_[stim_duration] = dict_.lookupOrDefault<scalar>("stim_duration", 0.0);
            constants_[stim_amplitude] = dict_.lookupOrDefault<scalar>("stim_amplitude", 0.0);
            constants_[nstim1] = dict_.lookupOrDefault<scalar>("nstim1", 0.0);
            constants_[stim_period_S2] = dict_.lookupOrDefault<scalar>("stim_period_S2", 0.0);
            constants_[nstim2] = dict_.lookupOrDefault<scalar>("nstim2", 0.0);
        }

        // Modified for cardiacFoam: the only shared mutable state on the ODE
        // path. The bound test and the increment must be one indivisible step,
        // so this is "critical" rather than "atomic" - an atomic increment would
        // still race with the test. It only executes when a correction actually
        // fires, which is rare, so it does not serialise the loop.
        void recordCorrection()
        {
            #ifdef _OPENMP
            #pragma omp critical(cardiacFoamProtectionCount)
            #endif
            {
                if (protectionCorrections_ < maxProtectionLogEntries_)
                {
                    ++protectionCorrections_;
                }
            }
        }

        void protectState
        (
            scalarField& state,
            const scalar,
            const label,
            const char*
        )
        {
            if (!numericalProtection_) return;

            auto correct = [&](const label i, const scalar value)
            {
                if (state[i] != value)
                {
                    state[i] = value;
                    recordCorrection();
                }
            };

            if (!std::isfinite(state[V]))
            {
                correct(V, -86.2);
            }
            if (clampVmForRates_)
            {
                correct(V, min(max(state[V], vmMinForRates_), vmMaxForRates_));
            }
            if (mag(state[V]) < vmZeroEps_)
            {
                correct(V, state[V] < 0.0 ? -vmZeroEps_ : vmZeroEps_);
            }

            const label concentrationIDs[] = {K_i, Na_i, Ca_i, Ca_SR};
            for (const label id : concentrationIDs)
            {
                if (!std::isfinite(state[id]) || state[id] < concentrationFloor_)
                {
                    correct(id, concentrationFloor_);
                }
            }

            if (clampGates_)
            {
                const label gateIDs[] = {Xr1, Xr2, Xs, m, h, j, d, f, fCa, s, r, g};
                for (const label id : gateIDs)
                {
                    if (!std::isfinite(state[id]))
                    {
                        correct(id, 0.0);
                    }
                    else
                    {
                        correct(id, min(max(state[id], scalar(0.0)), scalar(1.0)));
                    }
                }
            }
        }

        void computeRates
        (
            const scalar tMs,
            const scalarField& y,
            scalarField& dydt,
            scalarField& algebraic,
            const label pointI,
            const char* phase
        )
        {
            if (numericalProtection_)
            {
                scalarField yProtected(y);
                protectState(yProtected, tMs, pointI, phase);
                TNNPcomputeRates
                (
                    tMs,
                    constants_.data(),
                    dydt.data(),
                    yProtected.data(),
                    algebraic.data(),
                    tissue_,
                    false
                );
            }
            else
            {
                scalarField yCopy(y);
                TNNPcomputeRates
                (
                    tMs,
                    constants_.data(),
                    dydt.data(),
                    yCopy.data(),
                    algebraic.data(),
                    tissue_,
                    false
                );
            }
        }

        void computeVariables
        (
            const scalar tMs,
            scalarField& y,
            scalarField& rates,
            scalarField& algebraic
        )
        {
            TNNPcomputeVariables
            (
                tMs,
                constants_.data(),
                rates.data(),
                y.data(),
                algebraic.data(),
                tissue_,
                false
            );
        }

        void eulerStep
        (
            const scalar tMs,
            const scalar hMs,
            scalarField& y,
            scalarField& rates,
            scalarField& algebraic,
            const label pointI
        )
        {
            computeRates(tMs, y, rates, algebraic, pointI, "Euler");
            forAll(y, i)
            {
                y[i] += hMs*rates[i];
            }
        }

        void rk4Step
        (
            const scalar tMs,
            const scalar hMs,
            scalarField& y,
            scalarField& rates,
            scalarField& algebraic,
            const label pointI
        )
        {
            scalarField k1(NUM_STATES, 0.0), k2(NUM_STATES, 0.0);
            scalarField k3(NUM_STATES, 0.0), k4(NUM_STATES, 0.0);
            scalarField yt(NUM_STATES, 0.0);

            computeRates(tMs, y, k1, algebraic, pointI, "RK4_k1");
            forAll(y, i) yt[i] = y[i] + 0.5*hMs*k1[i];
            computeRates(tMs + 0.5*hMs, yt, k2, algebraic, pointI, "RK4_k2");
            forAll(y, i) yt[i] = y[i] + 0.5*hMs*k2[i];
            computeRates(tMs + 0.5*hMs, yt, k3, algebraic, pointI, "RK4_k3");
            forAll(y, i) yt[i] = y[i] + hMs*k3[i];
            computeRates(tMs + hMs, yt, k4, algebraic, pointI, "RK4_k4");

            forAll(y, i)
            {
                y[i] += (hMs/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
            }
            rates = k4;
        }

        void rkf45Trial
        (
            const scalar tMs,
            const scalar hMs,
            const scalarField& y,
            scalarField& yFifth,
            scalarField& rates,
            scalarField& algebraic,
            const label pointI,
            scalar& err
        )
        {
            scalarField k1(NUM_STATES, 0.0), k2(NUM_STATES, 0.0), k3(NUM_STATES, 0.0);
            scalarField k4(NUM_STATES, 0.0), k5(NUM_STATES, 0.0), k6(NUM_STATES, 0.0);
            scalarField yt(NUM_STATES, 0.0), yFourth(NUM_STATES, 0.0);

            computeRates(tMs, y, k1, algebraic, pointI, "RKF45_k1");
            forAll(y, i) yt[i] = y[i] + hMs*(1.0/5.0)*k1[i];
            computeRates(tMs + hMs/5.0, yt, k2, algebraic, pointI, "RKF45_k2");
            forAll(y, i) yt[i] = y[i] + hMs*((3.0/40.0)*k1[i] + (9.0/40.0)*k2[i]);
            computeRates(tMs + 3.0*hMs/10.0, yt, k3, algebraic, pointI, "RKF45_k3");
            forAll(y, i) yt[i] = y[i] + hMs*((3.0/10.0)*k1[i] - (9.0/10.0)*k2[i] + (6.0/5.0)*k3[i]);
            computeRates(tMs + 3.0*hMs/5.0, yt, k4, algebraic, pointI, "RKF45_k4");
            forAll(y, i) yt[i] = y[i] + hMs*((-11.0/54.0)*k1[i] + (5.0/2.0)*k2[i] - (70.0/27.0)*k3[i] + (35.0/27.0)*k4[i]);
            computeRates(tMs + hMs, yt, k5, algebraic, pointI, "RKF45_k5");
            forAll(y, i) yt[i] = y[i] + hMs*((1631.0/55296.0)*k1[i] + (175.0/512.0)*k2[i] + (575.0/13824.0)*k3[i] + (44275.0/110592.0)*k4[i] + (253.0/4096.0)*k5[i]);
            computeRates(tMs + 7.0*hMs/8.0, yt, k6, algebraic, pointI, "RKF45_k6");

            err = 0.0;
            forAll(y, i)
            {
                yFifth[i] = y[i] + hMs*((37.0/378.0)*k1[i] + (250.0/621.0)*k3[i] + (125.0/594.0)*k4[i] + (512.0/1771.0)*k6[i]);
                yFourth[i] = y[i] + hMs*((2825.0/27648.0)*k1[i] + (18575.0/48384.0)*k3[i] + (13525.0/55296.0)*k4[i] + (277.0/14336.0)*k5[i] + 0.25*k6[i]);
                const scalar scale = odeAbsTol_ + odeRelTol_*max(mag(yFifth[i]), mag(yFourth[i]));
                err = max(err, mag(yFifth[i] - yFourth[i])/max(scale, SMALL));
            }
            rates = k6;
        }

        void advancePoint
        (
            const scalar tStartMs,
            const scalar dtMs,
            const scalar VmVolt,
            scalarField& state,
            scalarField& algebraic,
            scalarField& rates,
            scalar& stepMs,
            const label pointI
        )
        {
            if (state.size() != NUM_STATES)
            {
                state.setSize(NUM_STATES, 0.0);
                TNNPinitConsts(constants_.data(), rates.data(), state.data(), tissue_, dict_);
                loadStimulusProtocol();
            }

            state[V] = 1000.0*VmVolt;
            protectState(state, tStartMs, pointI, "solveODE_start");

            if (dtMs <= SMALL)
            {
                computeVariables(tStartMs, state, rates, algebraic);
                return;
            }

            if (odeSolverName_ == "Euler" || odeSolverName_ == "forwardEuler")
            {
                eulerStep(tStartMs, dtMs, state, rates, algebraic, pointI);
                protectState(state, tStartMs + dtMs, pointI, "solveODE_end");
                computeVariables(tStartMs + dtMs, state, rates, algebraic);
                stepMs = dtMs;
                return;
            }

            if (odeSolverName_ == "RK4")
            {
                rk4Step(tStartMs, dtMs, state, rates, algebraic, pointI);
                protectState(state, tStartMs + dtMs, pointI, "solveODE_end");
                computeVariables(tStartMs + dtMs, state, rates, algebraic);
                stepMs = dtMs;
                return;
            }

            if
            (
                odeSolverName_ != "RKF45"
             && odeSolverName_ != "RKT45"
             && odeSolverName_ != "rkf45"
            )
            {
                FatalErrorInFunction
                    << "Unknown embedded TNNP ODE solver: " << odeSolverName_ << nl
                    << "Valid options are RKF45, RKT45, RK4, Euler"
                    << exit(FatalError);
            }

            scalar tau = 0.0;
            scalar h = stepMs > SMALL ? min(dtMs, stepMs) : min(dtMs, odeInitialStepMs_);
            if (h <= SMALL) h = dtMs;
            const scalar hMin = max(dtMs*1.0e-12, SMALL);
            scalarField yTrial(NUM_STATES, 0.0);

            for (label subStep = 0; subStep < max(odeMaxSteps_, label(1)); ++subStep)
            {
                if (tau >= dtMs - SMALL) break;
                h = min(h, dtMs - tau);

                scalar err = GREAT;
                rkf45Trial
                (
                    tStartMs + tau,
                    h,
                    state,
                    yTrial,
                    rates,
                    algebraic,
                    pointI,
                    err
                );

                if (err <= 1.0 || h <= hMin)
                {
                    state = yTrial;
                    tau += h;
                }

                const scalar factor = min
                (
                    scalar(4.0),
                    max(scalar(0.1), scalar(0.84)*std::pow(1.0/(err + SMALL), 0.25))
                );
                h = max(hMin, min(dtMs - tau, h*factor));
                stepMs = h;
            }

            if (tau < dtMs - 10.0*SMALL)
            {
                WarningInFunction
                    << "Embedded TNNP RKF45 reached maxSteps before completing dt. "
                    << "tau = " << tau << " ms, dt = " << dtMs << " ms" << nl;
            }

            protectState(state, tStartMs + dtMs, pointI, "solveODE_end");
            computeVariables(tStartMs + dtMs, state, rates, algebraic);
        }

    public:

        // Added for cardiacFoam: the "no slot map" default for solveODE().
        // A function-local static is used rather than labelUList::null() so
        // that .empty() is always well defined.
        static const labelList& emptyScratchIDs()
        {
            static const labelList empty;
            return empty;
        }

        // Added for cardiacFoam: first scratch slot of the cell-centre bank.
        // Callers integrating cell-centred points alongside Gauss points build
        // their slot map as cellSlotOffset() + cellI.
        label cellSlotOffset() const
        {
            return cellSlotOffset_;
        }

        label nCellSlots() const
        {
            return nCellSlots_;
        }

        // Modified for cardiacFoam: takes nCellSlots (= mesh.nCells()) and sizes
        // the per-point scratch to nIntegrationPoints + nCellSlots, so that the
        // front-aware hybrid's Gauss-point pass and cell-centre pass get
        // disjoint scratch banks instead of both starting at slot 0.
        // oldStates_ is deliberately NOT extended: it is only ever assigned
        // whole-array (oldStates_ = states), so it resizes itself to the state
        // array in use and is not indexed by scratch slot.
        EmbeddedTNNPModel
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const label nCellSlots,
            const scalar initialDeltaT
        )
        :
            dict_(dict),
            nPoints_(nIntegrationPoints),
            nCellSlots_(nCellSlots),
            cellSlotOffset_(nIntegrationPoints),
            constants_(NUM_CONSTANTS, 0.0),
            algebraic_
            (
                nIntegrationPoints + nCellSlots,
                Field<scalar>(NUM_ALGEBRAIC, 0.0)
            ),
            rates_
            (
                nIntegrationPoints + nCellSlots,
                Field<scalar>(NUM_STATES, 0.0)
            ),
            oldStates_(nIntegrationPoints, Field<scalar>(NUM_STATES, 0.0)),
            stepMs_
            (
                nIntegrationPoints + nCellSlots,
                max(1000.0*initialDeltaT, SMALL)
            ),
            exportedNames_(),
            debugNames_(),
            tissue_(tissueFlag(dict.lookupOrDefault<word>("tissue", "epicardialCells"))),
            odeSolverName_(dict.lookupOrDefault<word>("solver", dict.lookupOrDefault<word>("stateODESolver", "RKF45"))),
            odeInitialStepMs_(dict.lookupOrDefault<scalar>("initialODEStep", dict.lookupOrDefault<scalar>("stateODEInitialStep", 1.0e-5))),
            odeAbsTol_(max(dict.lookupOrDefault<scalar>("absTol", dict.lookupOrDefault<scalar>("stateODEAbsTol", 1.0e-9)), scalar(SMALL))),
            odeRelTol_(max(dict.lookupOrDefault<scalar>("relTol", dict.lookupOrDefault<scalar>("stateODERelTol", 1.0e-6)), scalar(SMALL))),
            odeMaxSteps_(dict.lookupOrDefault<label>("maxSteps", dict.lookupOrDefault<label>("stateODEMaxSteps", 10000))),
            numericalProtection_(dict.lookupOrDefault<Switch>("TNNPNumericalProtection", false)),
            clampGates_(dict.lookupOrDefault<Switch>("TNNPClampGates", true)),
            clampVmForRates_(dict.lookupOrDefault<Switch>("TNNPClampVmForRates", true)),
            logNumericalProtection_(dict.lookupOrDefault<Switch>("TNNPLogNumericalProtection", true)),
            concentrationFloor_(max(dict.lookupOrDefault<scalar>("TNNPConcentrationFloor", 1.0e-12), scalar(SMALL))),
            vmMinForRates_(dict.lookupOrDefault<scalar>("TNNPVMinForRates", -200.0)),
            vmMaxForRates_(dict.lookupOrDefault<scalar>("TNNPVMaxForRates", 100.0)),
            vmZeroEps_(max(dict.lookupOrDefault<scalar>("TNNPVoltageZeroEps", 1.0e-8), scalar(SMALL))),
            maxProtectionLogEntries_(dict.lookupOrDefault<label>("TNNPMaxProtectionLogEntries", 1000000)),
            protectionLogName_(dict.lookupOrDefault<word>("TNNPProtectionLog", "TNNP_numericalProtection_summary.dat")),
            protectionCorrections_(0),
            // Added for cardiacFoam. The key name matches the one the solver
            // report documents and the one the tutorial driver already writes
            // into this dictionary - which until now nothing read.
            parallelODE_(dict.lookupOrDefault<Switch>("parallelODE", false)),
            parallelODEThreshold_
            (
                dict.lookupOrDefault<label>("parallelODEThreshold", 256)
            )
        {
            if (dict_.found("exportedVariables"))
            {
                dict_.lookup("exportedVariables") >> exportedNames_;
            }
            if (dict_.found("debugPrintVariables"))
            {
                dict_.lookup("debugPrintVariables") >> debugNames_;
            }

            Info<< "Using embedded TNNP ionic model" << nl
                << "    tissue flag = " << tissue_ << nl
                << "    ODE solver = " << odeSolverName_ << nl
                << "    absTol = " << odeAbsTol_ << nl
                << "    relTol = " << odeRelTol_ << nl
                << "    maxSteps = " << odeMaxSteps_ << nl;

            if (numericalProtection_)
            {
                Info<< "    numerical protection enabled" << nl;
            }
        }

        ~EmbeddedTNNPModel()
        {
            if (numericalProtection_ && logNumericalProtection_)
            {
                // Modified for cardiacFoam: protectionLogName_ is a bare
                // file name, and mpirun starts every rank in the case
                // directory, so all ranks used to open and truncate the SAME
                // file. Each non-master rank now writes inside its own
                // processorN. Every rank still constructs the stream - see the
                // deadlock note on postProcessingDir().
                const fileName protectionLogPath
                (
                    Pstream::parRun() && !Pstream::master()
                  ? fileName("processor" + Foam::name(Pstream::myProcNo()))
                        /protectionLogName_
                  : fileName(protectionLogName_)
                );

                OFstream os(protectionLogPath);
                os  << "# Embedded TNNP numerical protection summary" << nl
                    << "corrections " << protectionCorrections_ << nl;
            }
        }

        label nEqns() const
        {
            return NUM_STATES;
        }

        wordList exportedFieldNames() const
        {
            return exportedNames_;
        }

        // Initialise the ionic states on the point set the caller owns.
        //
        // The caller sizes the state array to the point set the ODE is
        // integrated on, and that is NOT always nPoints_: with
        // stateIntegrationMode = cellCentredReconstruct the ODE is driven by
        // the cell-centred Vm, so the array holds mesh.nCells() entries while
        // nPoints_ counts the Iion Gauss points. Resizing to nPoints_ here
        // used to make the size guard in solveODE()/calculateCurrent()
        // permanently unsatisfiable in that mode, so the guard re-fired on
        // every call and silently reset the whole cell model to its resting
        // state at every time step (no state history could ever accumulate).
        //
        // The per-point scratch (rates_, algebraic_, stepMs_) is sized in the
        // constructor to nPoints_, the largest point count the model is ever
        // driven with, and is left alone here.
        void initialiseStates(Field<Field<scalar>>& states)
        {
            if (states.empty())
            {
                states.setSize(nPoints_);
            }

            // The scratch is indexed by the point set being integrated, which
            // may be larger than states.size(); keep it covering all nPoints_
            // entries regardless of how many states the caller holds.
            forAll(rates_, pointI)
            {
                rates_[pointI].setSize(NUM_STATES, 0.0);
                algebraic_[pointI].setSize(NUM_ALGEBRAIC, 0.0);
            }

            forAll(states, pointI)
            {
                states[pointI].setSize(NUM_STATES, 0.0);
                TNNPinitConsts
                (
                    constants_.data(),
                    rates_[pointI].data(),
                    states[pointI].data(),
                    tissue_,
                    dict_
                );
                loadStimulusProtocol();
                computeVariables(0.0, states[pointI], rates_[pointI], algebraic_[pointI]);
            }
            oldStates_ = states;
        }

        void updateStatesOld(const Field<Field<scalar>>& states)
        {
            oldStates_ = states;
        }

        void resetStatesToStatesOld(Field<Field<scalar>>& states)
        {
            states = oldStates_;
        }

        void calculateCurrent
        (
            const scalar stepStartTime,
            const scalar,
            const scalarField& Vm,
            scalarField& Im,
            Field<Field<scalar>>& states
        )
        {
            const scalar tStartMs = 1000.0*stepStartTime;
            if (Im.size() != Vm.size())
            {
                FatalErrorInFunction
                    << "Im.size() != Vm.size()" << abort(FatalError);
            }

            // A mismatch here is a programming error, not something to repair
            // on the fly: re-initialising would discard the accumulated state
            // history. The state array must already be sized to the point set
            // it is integrated on (mesh.nCells() for cellCentredReconstruct,
            // the Iion Gauss points otherwise).
            if (states.size() != Vm.size())
            {
                FatalErrorInFunction
                    << "states.size() = " << states.size()
                    << " != Vm.size() = " << Vm.size() << nl
                    << "The ionic state array must be sized to the point set "
                    << "it is integrated on." << abort(FatalError);
            }

            forAll(Vm, pointI)
            {
                scalarField& y = states[pointI];
                y[V] = 1000.0*Vm[pointI];
                protectState(y, tStartMs, pointI, "calculateCurrent");
                computeVariables(tStartMs, y, rates_[pointI], algebraic_[pointI]);
                Im[pointI] = algebraic_[pointI][Iion_cm];
            }
        }

        void solveODE
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im,
            Field<Field<scalar>>& states,
            const labelUList& scratchIDs = emptyScratchIDs()
        )
        {
            const scalar tStartMs = 1000.0*stepStartTime;
            const scalar dtMs = 1000.0*deltaT;
            if (Im.size() != Vm.size())
            {
                FatalErrorInFunction
                    << "Im.size() != Vm.size()" << abort(FatalError);
            }

            // Added for cardiacFoam: the caller may pass an explicit scratch
            // slot for each point. This is what keeps the front-aware hybrid
            // correct. That path calls solveODE twice per iteration, once with
            // a compacted array of nF front Gauss points and once with a
            // compacted array of nS smooth cells, both indexed from 0. Without
            // a slot map both passes write algebraic_/rates_/stepMs_[0..n),
            // so the second silently overwrites the first's scratch even
            // though the two arrays describe different physical points.
            if (!scratchIDs.empty())
            {
                if (scratchIDs.size() != Vm.size())
                {
                    FatalErrorInFunction
                        << "scratchIDs.size() = " << scratchIDs.size()
                        << " != Vm.size() = " << Vm.size()
                        << abort(FatalError);
                }

                forAll(scratchIDs, i)
                {
                    if (scratchIDs[i] < 0 || scratchIDs[i] >= rates_.size())
                    {
                        FatalErrorInFunction
                            << "scratchIDs[" << i << "] = " << scratchIDs[i]
                            << " out of range [0, " << rates_.size() << ")"
                            << abort(FatalError);
                    }
                }
            }

            // A mismatch here is a programming error, not something to repair
            // on the fly: re-initialising would discard the accumulated state
            // history. The state array must already be sized to the point set
            // it is integrated on (mesh.nCells() for cellCentredReconstruct,
            // the Iion Gauss points otherwise).
            if (states.size() != Vm.size())
            {
                FatalErrorInFunction
                    << "states.size() = " << states.size()
                    << " != Vm.size() = " << Vm.size() << nl
                    << "The ionic state array must be sized to the point set "
                    << "it is integrated on." << abort(FatalError);
            }

            // Modified for cardiacFoam: OpenMP across the points.
            //
            // Safe because every iteration touches only data of its own point:
            // states[pointI] and Im[pointI] are distinct entries, and the
            // scratch is addressed by SLOT, with the slot map injective (the
            // identity when absent, and the front/smooth maps of the hybrid are
            // disjoint by construction - see the two scratch banks above). The
            // only shared mutable state on this path is the protection counter,
            // which recordCorrection() now serialises.
            //
            // Written as a canonical for-loop rather than forAll because OpenMP
            // cannot parallelise the macro's form. Iteration order does not
            // affect the result: the points do not communicate.
            //
            // Caveat worth stating: a FatalError raised inside the loop escapes
            // a parallel region, which OpenMP does not define. It aborts the
            // process either way, but the message may interleave. The MMS
            // solver's OpenMP loops carry the same exposure.
            const label nPointsLocal = Vm.size();

            #ifdef _OPENMP
            #pragma omp parallel for schedule(static) \
                if (parallelODE_ && nPointsLocal >= parallelODEThreshold_)
            #endif
            for (label pointI = 0; pointI < nPointsLocal; ++pointI)
            {
                // The scratch is addressed by slot, not by position in the
                // passed array. With no slot map the two coincide, so every
                // other call site is unaffected.
                const label slotI =
                    scratchIDs.empty() ? pointI : scratchIDs[pointI];

                advancePoint
                (
                    tStartMs,
                    dtMs,
                    Vm[pointI],
                    states[pointI],
                    algebraic_[slotI],
                    rates_[slotI],
                    stepMs_[slotI],
                    slotI
                );
                Im[pointI] = algebraic_[slotI][Iion_cm];
            }
        }

        void exportStates
        (
            const Field<Field<scalar>>& states,
            PtrList<volScalarField>& outFields
        )
        {
            forAll(outFields, outI)
            {
                const word& name = exportedNames_[outI];
                const label sI = stateIndex(name);
                const label aI = algebraicIndex(name);
                scalarField& fld = outFields[outI].primitiveFieldRef();

                forAll(fld, cellI)
                {
                    if (sI >= 0 && cellI < states.size())
                    {
                        fld[cellI] = states[cellI][sI];
                    }
                    else if (aI >= 0 && cellI < algebraic_.size())
                    {
                        fld[cellI] = algebraic_[cellI][aI];
                    }
                    else
                    {
                        fld[cellI] = 0.0;
                    }
                }
                outFields[outI].correctBoundaryConditions();
            }
        }

        void exportStatesIntegrationPoints
        (
            const Field<Field<scalar>>& states,
            PtrList<volScalarField>& outFields,
            const CompactListList<scalar>& cellIionQuadW
        )
        {
            forAll(outFields, outI)
            {
                const word& name = exportedNames_[outI];
                const label sI = stateIndex(name);
                const label aI = algebraicIndex(name);
                scalarField& fld = outFields[outI].primitiveFieldRef();

                label integrationPointI = 0;
                forAll(cellIionQuadW, cellI)
                {
                    scalar wSum = 0.0;
                    scalar value = 0.0;
                    forAll(cellIionQuadW[cellI], qI)
                    {
                        const scalar w = cellIionQuadW[cellI][qI];
                        wSum += w;
                        if (sI >= 0 && integrationPointI < states.size())
                        {
                            value += w*states[integrationPointI][sI];
                        }
                        else if (aI >= 0 && integrationPointI < algebraic_.size())
                        {
                            value += w*algebraic_[integrationPointI][aI];
                        }
                        ++integrationPointI;
                    }
                    fld[cellI] = value/max(wSum, SMALL);
                }
                outFields[outI].correctBoundaryConditions();
            }
        }
    };

    // Representative cell size: the geometric mean edge of a cell if the domain
    // were divided uniformly, i.e. (active volume / total cells)^(1/dim). Used
    // only to size the stable time step, so a global average is the right
    // notion. The cell count is reduced, so the value does not depend on the
    // decomposition.
    scalar characteristicDx(const fvMesh& mesh)
    {
        const boundBox bb(mesh.points());
        const vector span = bb.max() - bb.min();
        const label dim = max(mesh.nGeometricD(), label(1));

        scalar activeMeasure = 1.0;
        if (dim >= 1)
        {
            activeMeasure *= max(span.x(), SMALL);
        }
        if (dim >= 2)
        {
            activeMeasure *= max(span.y(), SMALL);
        }
        if (dim >= 3)
        {
            activeMeasure *= max(span.z(), SMALL);
        }

        const scalar nCellsGlobal = scalar(returnReduce(mesh.nCells(), sumOp<label>()));

        return std::pow(activeMeasure/max(nCellsGlobal, scalar(1.0)), 1.0/scalar(dim));
    }

    // Explicit-diffusion stability limit dt = CFL*dx^2/(dim*D_eff), with
    // D_eff = max(D)/(chi*Cm). Reported (and optionally used) as a reference
    // even though the solver is implicit, because it is the scale at which the
    // nonlinear coupling stops being stiff. Reduced with min so all ranks agree.
    scalar computeStableDeltaT
    (
        const volTensorField& conductivity,
        const scalar chiVal,
        const scalar CmVal,
        const scalar CFL,
        const scalar dx,
        const label dim
    )
    {
        const scalar Dxx = max(conductivity.component(tensor::XX)).value();
        const scalar Dyy = max(conductivity.component(tensor::YY)).value();
        const scalar Dzz = max(conductivity.component(tensor::ZZ)).value();
        const scalar Dmax = max(Dxx, max(Dyy, Dzz));
        const scalar DefMax = Dmax/(chiVal*CmVal);

        scalar dtStable = GREAT;
        if (dim > 0 && DefMax > SMALL)
        {
            dtStable = CFL*sqr(dx)/(scalar(dim)*DefMax);
        }

        reduce(dtStable, minOp<scalar>());
        return dtStable;
    }

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

    // Tissue-level external stimulus: sets the current to stimulusIntensity in
    // the cells of every stimulus region whose window [tStart, tStart+duration]
    // contains t0, and to zero everywhere else. Regions and their start times
    // come from constant/stimulusProtocol.
    void applyStimulus
    (
        const scalar t0,
        volScalarField& externalStimulusCurrent,
        const List<labelList>& stimulusCellIDsList,
        const List<scalar>& stimulusStartTimes,
        const scalar stimulusIntensity,
        const scalar stimulusDuration
    )
    {
        scalarField& extI = externalStimulusCurrent.primitiveFieldRef();
        extI = 0.0;

        forAll(stimulusCellIDsList, bI)
        {
            const scalar tStart = stimulusStartTimes[bI];
            if (t0 < tStart || t0 > (tStart + stimulusDuration))
            {
                continue;
            }

            const labelList& stimulusCellIDs = stimulusCellIDsList[bI];
            forAll(stimulusCellIDs, cI)
            {
                extI[stimulusCellIDs[cI]] = stimulusIntensity;
            }
        }

        externalStimulusCurrent.correctBoundaryConditions();
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
    // the code actually branches on, lumped and consistent.
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
            << "Valid options are lumped/diagonal and consistent/consistentHO"
            << exit(FatalError);

        return "lumped";
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

    // ||current - previous|| / ||current||, globally reduced. The nonlinear
    // stopping test for the state fields.
    scalar relativeL2Difference
    (
        const scalarField& current,
        const scalarField& previous
    )
    {
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

        return std::sqrt(num)/(std::sqrt(den) + SMALL);
    }

    // ||r|| / ||reference||, both global. The nonlinear stopping test for Vm
    // and the linear one for the Krylov solvers.
    scalar relativeL2Norm(const EigVec& r, const EigVec& reference)
    {
        // Modified for cardiacFoam: both norms are over LOCAL blocks only.
        return gNorm(r)/(gNorm(reference) + SMALL);
    }

    // Worst relative L2 change over all state variables between two nonlinear
    // iterations, also filling the per-state residuals for reporting. The
    // states converge at different rates, so the stopping test uses the maximum
    // rather than an average.
    scalar maxStateRelativeL2Difference
    (
        const Field<Field<scalar>>& current,
        const Field<Field<scalar>>& previous,
        scalarField& stateResiduals
    )
    {
        stateResiduals = 0.0;

        if (current.empty() || previous.empty())
        {
            return 0.0;
        }

        const label nStates = min(current[0].size(), stateResiduals.size());
        scalar maxResidual = 0.0;

        for (label stateI = 0; stateI < nStates; ++stateI)
        {
            scalar num = 0.0;
            scalar den = 0.0;

            forAll(current, pointI)
            {
                if (stateI >= current[pointI].size() || stateI >= previous[pointI].size())
                {
                    continue;
                }

                const scalar c = current[pointI][stateI];
                const scalar diff = c - previous[pointI][stateI];
                num += diff*diff;
                den += c*c;
            }

            reduce(num, sumOp<scalar>());
            reduce(den, sumOp<scalar>());

            stateResiduals[stateI] = std::sqrt(num)/(std::sqrt(den) + SMALL);
            maxResidual = max(maxResidual, stateResiduals[stateI]);
        }

        return maxResidual;
    }

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
                const label col = it.col() - gRowStart;
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

    // Added for cardiacFoam: directory for this run's postProcessing output.
    //
    // In a parallel run runTime.path() is <case>/processorN, so every one of
    // these files used to be written inside a processor directory where the
    // tutorial driver does not look for them.
    //
    // EVERY rank must still construct the stream. Guarding with
    // "if (!Pstream::master()) return;" before opening DEADLOCKS: OpenFOAM's
    // file handler can be collective, so the master blocks inside the
    // collective open while the other ranks have already moved on. That failure
    // was diagnosed in the MMS solver as a run which completed every time step,
    // printed its timings and then hung. Only the destination differs here.
    fileName postProcessingDir(const Time& runTime)
    {
        const fileName base =
        (
            Pstream::master()
          ? runTime.rootPath()/runTime.globalCaseName()
          : runTime.path()
        );

        return base/"postProcessing"/"highOrderElectroActivationFoamImplicitPETSc";
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

    // High-order consistent mass matrix assembled with LRE cell quadrature.
    //
    //   For every cell i and every cell j in its LRE stencil we accumulate
    //
    //       M(i, j) = coefficient * sum_{qp in cell_i} w_qp * phi_j(x_qp)
    //
    //   where phi_j(x_qp) is the LRE Taylor-basis function of stencil cell j
    //   evaluated at quadrature point x_qp using the cell-centred Taylor
    //   expansion (value + grad + Hessian + 3rd-order derivative). The
    //   resulting M is effectively a high-order weighted-averaging operator:
    //   (M * Vm)[i] approximates the cell-averaged Vm at cell i to the
    //   order of the LRE reconstruction.
    void assembleConsistentMassMatrixHO
    (
        const fvMesh& mesh,
        const highOrderInterp& LREInterp,
        const scalar coefficient,
        const bool compactRows,
        SpMat& M
    )
    {
        const bool twoD = mesh.nGeometricD() == 2;
        const vectorField& C = mesh.C();

        // Stencils and quadrature data from the LRE object.
        const CompactListList<label>& stencils =
            LREInterp.globalCellStencils();
        const CompactListList<point>& cellQP = LREInterp.cellQuadPoints();
        const CompactListList<scalar>& cellQW = LREInterp.cellQuadWeightPhysical();
        // Linear, quadratic and cubic Taylor-reconstruction coefficients.
        // For each stencil cell cI of cell cellI, gradCoeffs[cellI][cI] gives
        // the contribution of Vm_{stencil_cI} to grad Vm at cell cellI's
        // centre; similarly for the Hessian and 3rd-order derivatives.
        const CompactListList<vector>& gradCoeffs = LREInterp.QRGradCoeffs();
        const CompactListList<symmTensor>& hessCoeffs =
            LREInterp.cellHessianCoeffs();
        const CompactListList<highOrderInterp::symmTensor3Order>& thirdCoeffs =
            LREInterp.cellThirdDerivCoeffs();

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
            triplets.reserve(mesh.nCells()*40);
        }

        forAll(stencils, cellI)
        {
            // Modified for cardiacFoam: must be UList<label> BY VALUE, not
                // a labelList reference - see the note on
                // CompactListList in the header comment.
                const UList<label> curStencil = stencils[cellI];
            // LRE stores stencil-cell coefficients in [0..nStencil-1] and
            // the self (centre) coefficients at index nStencil.
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

            // Normalise quadrature weights so they sum to one. The mass
            // matrix then approximates the cell-averaged operator
            //   (M Vm)[i] ~ (1/|Omega_i|) integral_{Omega_i} Vm dV
            scalar wSum = 0.0;
            forAll(cellQW[cellI], qpI)
            {
                wSum += cellQW[cellI][qpI];
            }
            wSum = max(wSum, SMALL);

            forAll(cellQP[cellI], qpI)
            {
                const scalar w = cellQW[cellI][qpI]/wSum;
                // Displacement from cell centre to quadrature point.
                const vector d = cellQP[cellI][qpI] - C[cellI];

                // Off-diagonal entries: basis function of every stencil
                // neighbour evaluated at this quadrature point.
                forAll(curStencil, cI)
                {
                    // phi_j(x_qp) = grad . d + (1/2) H : (d x d) + (1/6) T3 :: (d x d x d)
                    scalar coeff = (gradCoeffs[cellI][cI] & d);

                    if (LREInterp.order() >= 2)
                    {
                        coeff += 0.5*quadraticForm(hessCoeffs[cellI][cI], d);
                    }

                    if (LREInterp.order() >= 3)
                    {
                        coeff +=
                            (1.0/6.0)
                           *highOrderInterp::cubicForm(thirdCoeffs[cellI][cI], d, twoD);
                    }

                    addTripletIfNeeded
                    (
                        triplets,
                        cellI,
                        curStencil[cI],
                        coefficient*w*coeff
                    );
                }

                // Diagonal (self) entry: starts at 1.0 because phi_i is the
                // partition-of-unity reconstruction centred at cell i.
                scalar selfCoeff = 1.0 + (gradCoeffs[cellI][selfCoeffI] & d);

                if (LREInterp.order() >= 2)
                {
                    selfCoeff +=
                        0.5*quadraticForm(hessCoeffs[cellI][selfCoeffI], d);
                }

                if (LREInterp.order() >= 3)
                {
                    selfCoeff +=
                        (1.0/6.0)
                       *highOrderInterp::cubicForm(thirdCoeffs[cellI][selfCoeffI], d, twoD);
                }

                addTripletIfNeeded
                (
                    triplets, cellI, gCol(cellI), coefficient*w*selfCoeff
                );
            }
        }

        M.resize(mesh.nCells(), gNCols());
        M.setFromTriplets(triplets.begin(), triplets.end());
        M.makeCompressed();

        reportOperatorFingerprint("M (consistent)", M);
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


    // Same jump against a Dirichlet boundary face, where the "neighbour" value
    // is prescribed and therefore goes to the right-hand side rather than into
    // the matrix; only the owner's terms are added here.
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
        addTripletIfNeeded(triplets, row, gCol(own), -scale);
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
    void assembleStandardOrthogonalStiffnessMatrix
    (
        const fvMesh& mesh,
        const volTensorField& conductivity,
        const highOrderInterp* LREInterp,
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

        if (stabilisationAlpha > SMALL && LREInterp == nullptr)
        {
            FatalErrorInFunction
                << "stabilisationAlpha > 0 requires LRECoeffs_Vm so the "
                << "Rhie-Chow stabilisation jump can be reconstructed."
                << abort(FatalError);
        }

        forAll(neighbour, faceI)
        {
            const label own = owner[faceI];
            const label nei = neighbour[faceI];

            const tensor Df = 0.5*(conductivity[own] + conductivity[nei]);
            const scalar a =
                orthogonalDiffusionCoeff(mesh.Sf()[faceI], C[nei] - C[own], Df);

            addTripletIfNeeded(triplets, own, gCol(own), -a/max(V[own], SMALL));
            addTripletIfNeeded(triplets, own, gCol(nei),  a/max(V[own], SMALL));
            addTripletIfNeeded(triplets, nei, gCol(own),  a/max(V[nei], SMALL));
            addTripletIfNeeded(triplets, nei, gCol(nei), -a/max(V[nei], SMALL));
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
                    *LREInterp
                );

                addRhieChowInternalJumpCoeffs
                (
                    triplets,
                    nei,
                    -a/max(V[nei], SMALL),
                    own,
                    nei,
                    dPN,
                    *LREInterp
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

            // Added for cardiacFoam: the far cell's gradient row, for the
            // stabilisation jump. Only under alpha > 0 - the flux terms above
            // need no reconstruction, and LREInterp may legitimately be null
            // when alpha is zero (the guard at the top of this function only
            // fatals when alpha > 0), so the lambda must not be built at all.
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
                            *LREInterp
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
                            *LREInterp
                        );
                    }
                }

                continue;
            }

            if
            (
                bcType == fixedValueFvPatchScalarField::typeName
             || bcType == "fixedVoltage"
            )
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

                    addTripletIfNeeded(triplets, own, gCol(own), -a/max(V[own], SMALL));

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
                            *LREInterp
                        );
                    }
                }
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
        const label faceChunk = 50000;

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
                        const scalar fluxCoeff =
                            // Modified for cardiacFoam: no area
                            // factor, w is a physical weight.
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

        // -------- Boundary faces -------------------------------------------
        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            const word bcType =
                patch.lookupPatchField<volScalarField, scalar>("Vm").type();

            // Insulating (Neumann) and empty patches contribute zero flux.
            // For Dirichlet-style patches we still need to assemble the
            // stencil contributions; the boundary value enters through the
            // RHS elsewhere.
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
                            w*(n & (conductivity[own] & gCoeff));

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
                else if
                (
                    stabilisationAlpha > SMALL
                 && (
                        bcType == fixedValueFvPatchScalarField::typeName
                     || bcType == "fixedVoltage"
                    )
                )
                {
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

        if (!kInitialised)
        {
            K.resize(nCells, gNCols());
            K.setZero();
        }

        K.makeCompressed();

        reportOperatorFingerprint("K (high order)", K);
        std::vector<Triplet>().swap(triplets);
    }

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

    // Evaluate Vm at every Iion integration point using a Taylor expansion
    // around the cell centre:
    //
    //     Vm(x_qp) = Vm_c + grad . d + (1/2) H : (d x d) + (1/6) T3 :: d^3
    //
    // The expansion order is set by the LRE object for Vm. When high-order
    // is disabled we fall back to a plain linear extrapolation from the
    // standard FVM cell gradient (still 2nd-order accurate on smooth fields
    // but without LRE's extended stencil).
    void reconstructVmAtIionIntegrationPoints
    (
        const volScalarField& Vm,
        const Switch useHighOrderVm,
        const highOrderInterp* LREInterp_Vm,
        const highOrderInterp& LREInterp_Iion,
        scalarField& VmIntegrationPoints
    )
    {
        const fvMesh& mesh = Vm.mesh();
        const vectorField& C = mesh.C();
        const CompactListList<point>& cellIionQuadP =
            LREInterp_Iion.cellQuadPoints();

        // Flat counter walks through every cell's quadrature points in the
        // same order they were laid out by the LRE constructor.
        label integrationPointI = 0;

        if (useHighOrderVm && LREInterp_Vm)
        {
            const bool twoD = mesh.nGeometricD() == 2;

            tmp<volVectorField> tGradVm = LREInterp_Vm->grad(Vm);
            const vectorField& gradVm = tGradVm->internalField();

            tmp<volSymmTensorField> tHessVm;
            const symmTensorField* hessVm = nullptr;
            if (LREInterp_Vm->order() >= 2)
            {
                tHessVm = LREInterp_Vm->hessian(Vm);
                hessVm = &(tHessVm->internalField());
            }

            autoPtr<List<highOrderInterp::symmTensor3Order>> thirdVmPtr;
            const List<highOrderInterp::symmTensor3Order>* thirdVm = nullptr;
            if (LREInterp_Vm->order() >= 3)
            {
                thirdVmPtr = LREInterp_Vm->thirdDeriv(Vm);
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


    // Reconstruct the cell-centred ionic states at the Iion Gauss points using
    // a dedicated high-order LRE (LREInterp_states). Generic companion to
    // reconstructVmAtIionIntegrationPoints() for the stateIntegrationMode =
    // cellCentredReconstruct path: each of the nStates state columns is loaded
    // into a reusable scratch volScalarField, its LRE Taylor coefficients are
    // computed, and the field is evaluated at every quadrature point.
    //
    //   statesCells [cellI][stateI]  ->  statesIP [ipI][stateI]
    void reconstructStatesAtIionIntegrationPoints
    (
        const Field<Field<scalar>>& statesCells,
        const label nStates,
        const highOrderInterp& LREInterp_states,
        const highOrderInterp& LREInterp_Iion,
        volScalarField& scratch,
        Field<Field<scalar>>& statesIP,
        const boolList* writeCellMask = nullptr,
        const bool limitToNeighbourRange = true
    )
    {
        const fvMesh& mesh = scratch.mesh();
        const vectorField& C = mesh.C();
        const bool twoD = mesh.nGeometricD() == 2;
        const CompactListList<point>& cellIionQuadP =
            LREInterp_Iion.cellQuadPoints();
        // Barth--Jespersen-style limiter: bound each reconstructed Gauss value
        // to the [min,max] of the cell's own value and its face neighbours.
        // High-order (p3) reconstruction of the stiff TNNP states across an
        // under-resolved front overshoots (Gibbs) and can drive the ionic
        // algebra to overflow (SIGFPE); clamping to the local data range
        // removes the overshoot while leaving smooth regions untouched.
        const labelListList& cellCells = mesh.cellCells();

        scalarField& s = scratch.primitiveFieldRef();

        for (label stateI = 0; stateI < nStates; ++stateI)
        {
            forAll(s, cellI)
            {
                s[cellI] = statesCells[cellI][stateI];
            }
            scratch.correctBoundaryConditions();

            // Added for cardiacFoam: the limiter's admissible range must also
            // include neighbours on the far side of a processor boundary.
            // mesh.cellCells() stops at the cut, so cells next to one saw a
            // narrower [min,max] and were clamped harder than the same cell in
            // serial. correctBoundaryConditions() has just run, so the coupled
            // patch fields already hold the neighbouring cell values.
            scalarField coupledLo(mesh.nCells(), GREAT);
            scalarField coupledHi(mesh.nCells(), -GREAT);
            if (Pstream::parRun())
            {
                forAll(mesh.boundary(), patchI)
                {
                    const fvPatch& p = mesh.boundary()[patchI];
                    if (!p.coupled())
                    {
                        continue;
                    }

                    const scalarField pnf
                    (
                        scratch.boundaryField()[patchI].patchNeighbourField()
                    );
                    const labelUList& faceCells = p.faceCells();

                    forAll(faceCells, i)
                    {
                        const label cellI = faceCells[i];
                        coupledLo[cellI] = min(coupledLo[cellI], pnf[i]);
                        coupledHi[cellI] = max(coupledHi[cellI], pnf[i]);
                    }
                }
            }

            tmp<volVectorField> tGrad = LREInterp_states.grad(scratch);
            const vectorField& grad = tGrad->internalField();

            tmp<volSymmTensorField> tHess;
            const symmTensorField* hess = nullptr;
            if (LREInterp_states.order() >= 2)
            {
                tHess = LREInterp_states.hessian(scratch);
                hess = &(tHess->internalField());
            }

            autoPtr<List<highOrderInterp::symmTensor3Order>> thirdPtr;
            const List<highOrderInterp::symmTensor3Order>* third = nullptr;
            if (LREInterp_states.order() >= 3)
            {
                thirdPtr = LREInterp_states.thirdDeriv(scratch);
                third = &thirdPtr();
            }

            label integrationPointI = 0;
            forAll(mesh.cells(), cellI)
            {
                const bool writeThisCell =
                    (writeCellMask == nullptr) || (*writeCellMask)[cellI];

                const scalar sc = s[cellI];
                const vector& gradc = grad[cellI];
                const vector& xc = C[cellI];

                const symmTensor* H = hess ? &((*hess)[cellI]) : nullptr;
                const highOrderInterp::symmTensor3Order* T3 =
                    third ? &((*third)[cellI]) : nullptr;

                // Local admissible range for the limiter.
                scalar sLo = sc, sHi = sc;
                if (writeThisCell && limitToNeighbourRange)
                {
                    const labelList& nbr = cellCells[cellI];
                    forAll(nbr, j)
                    {
                        const scalar sn = s[nbr[j]];
                        sLo = min(sLo, sn);
                        sHi = max(sHi, sn);
                    }

                    // Added for cardiacFoam: neighbours across a processor
                    // face. The initialisers are GREAT / -GREAT, so a cell with
                    // no coupled face is left exactly as it was.
                    sLo = min(sLo, coupledLo[cellI]);
                    sHi = max(sHi, coupledHi[cellI]);
                }

                forAll(cellIionQuadP[cellI], qI)
                {
                    if (writeThisCell)
                    {
                        const vector d = cellIionQuadP[cellI][qI] - xc;
                        scalar sq =
                            reconstructFromTaylor(sc, gradc, H, T3, d, twoD);
                        if (limitToNeighbourRange)
                        {
                            sq = min(max(sq, sLo), sHi);
                        }
                        statesIP[integrationPointI][stateI] = sq;
                    }
                    ++integrationPointI;
                }
            }
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

    // This is the discrete cell-average operator and it is the natural
    // companion to reconstructVmAtIionIntegrationPoints(): together they
    // form a Gauss-quadrature evaluation of  (1/|Omega|) int_Omega Iion(Vm).
    void averageIionIntegrationPointsToCells
    (
        const scalarField& IionIntegrationPoints,
        const highOrderInterp& LREInterp_Iion,
        volScalarField& Iion
    )
    {
        const fvMesh& mesh = Iion.mesh();
        const CompactListList<scalar>& cellIionQuadW =
            LREInterp_Iion.cellQuadWeightPhysical();

        scalarField& IionCells = Iion.primitiveFieldRef();
        label integrationPointI = 0;

        forAll(mesh.cells(), cellI)
        {
            scalar iBar = 0.0;
            scalar wSum = 0.0;

            forAll(cellIionQuadW[cellI], qI)
            {
                const scalar w = cellIionQuadW[cellI][qI];
                iBar += w*IionIntegrationPoints[integrationPointI];
                wSum += w;
                ++integrationPointI;
            }

            // Normalise by sum of weights so the result is an unbiased
            // average regardless of whether quadrature weights sum to 1.
            IionCells[cellI] = iBar/max(wSum, SMALL);
        }

        Iion.correctBoundaryConditions();
    }

    // Quadrature-average the per-Gauss-point state vectors to cell centres,
    // component by component. Companion to averageIionIntegrationPointsToCells
    // for the Field<Field<scalar>> state carrier; used by the front-aware
    // hybrid to build cell-centred states from the Gauss-point carrier.
    void averageStatesIntegrationPointsToCells
    (
        const Field<Field<scalar>>& statesIP,
        const label nStates,
        const highOrderInterp& LREInterp_Iion,
        Field<Field<scalar>>& statesCells,
        const boolList* cellMask = nullptr
    )
    {
        const CompactListList<scalar>& cellIionQuadW =
            LREInterp_Iion.cellQuadWeightPhysical();

        label integrationPointI = 0;
        forAll(cellIionQuadW, cellI)
        {
            const bool updateThisCell =
                (cellMask == nullptr) || (*cellMask)[cellI];

            if (!updateThisCell)
            {
                integrationPointI += cellIionQuadW[cellI].size();
                continue;
            }

            Field<scalar>& sc = statesCells[cellI];
            sc = 0.0;
            scalar wSum = 0.0;

            forAll(cellIionQuadW[cellI], qI)
            {
                const scalar w = cellIionQuadW[cellI][qI];
                const Field<scalar>& sq = statesIP[integrationPointI];
                for (label k = 0; k < nStates; ++k)
                {
                    sc[k] += w*sq[k];
                }
                wSum += w;
                ++integrationPointI;
            }

            const scalar invW = 1.0/max(wSum, SMALL);
            for (label k = 0; k < nStates; ++k)
            {
                sc[k] *= invW;
            }
        }
    }

    // Flag "front" cells for the front-aware hybrid: those whose sub-cell Vm
    // range (over their Iion Gauss points) exceeds a threshold, or which are
    // stimulated. This classification is computed ONCE per time step (from the
    // committed Vm) and frozen for the whole nonlinear solve -- recomputing it
    // per nonlinear iteration makes the front/smooth switch flip between
    // iterations and turns the Picard/Newton fixed point into a period-2 limit
    // cycle that never converges.
    void flagFrontCells
    (
        const fvMesh& mesh,
        const scalarField& VmIntegrationPoints,
        const highOrderInterp& LREInterp_Iion,
        const volScalarField& externalStimulusCurrent,
        const scalar vmRangeThreshold,
        boolList& frontCell
    )
    {
        const CompactListList<point>& cellQP = LREInterp_Iion.cellQuadPoints();
        label ip = 0;
        forAll(mesh.cells(), cellI)
        {
            scalar vMin = GREAT, vMax = -GREAT;
            forAll(cellQP[cellI], qI)
            {
                const scalar v = VmIntegrationPoints[ip];
                vMin = min(vMin, v);
                vMax = max(vMax, v);
                ++ip;
            }
            frontCell[cellI] =
                (vMax - vMin > vmRangeThreshold)
             || (mag(externalStimulusCurrent[cellI]) > SMALL);
        }
    }

    // Grow the set of front cells by nLayers cell layers (face neighbours),
    // so the moving front stays inside the Gauss-resolved band for the whole
    // time step and the front/smooth transition is padded.
    // Added for cardiacFoam: opt-in front-cell count diagnostic.
    bool debugFrontCellCount = false;

    // Grow the front flag by nLayers layers of face neighbours, so that the
    // expensive per-Gauss-point ODE integration is applied to a band around the
    // propagating front rather than exactly on it.
    void dilateFrontCells
    (
        const fvMesh& mesh,
        const label nLayers,
        boolList& frontCell
    )
    {
        if (nLayers <= 0)
        {
            return;
        }

        const labelListList& cellCells = mesh.cellCells();

        // Added for cardiacFoam: marker field used to carry the front flag
        // across processor boundaries. mesh.cellCells() stops at a processor
        // face, so without this the front band simply failed to grow across a
        // partition. Unlike the stabilisation gap this is an O(1) error, not a
        // shrinking surface effect: a cell classified front integrates a
        // different set of ODEs from one classified smooth.
        //
        // Constructing the field with zeroGradient does NOT put zeroGradient on
        // the processor patches: fvPatchField::New overrides the requested type
        // with the patch's own constraint type, so processor patches receive a
        // real processorFvPatchScalarField and correctBoundaryConditions()
        // performs the halo exchange. One exchange per dilation layer.
        autoPtr<volScalarField> markerPtr;
        if (Pstream::parRun())
        {
            markerPtr.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "frontCellMarker",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("zero", dimless, 0.0),
                    zeroGradientFvPatchScalarField::typeName
                )
            );
        }

        for (label layer = 0; layer < nLayers; ++layer)
        {
            boolList grown(frontCell);
            forAll(frontCell, cellI)
            {
                if (!frontCell[cellI])
                {
                    continue;
                }
                const labelList& nbr = cellCells[cellI];
                forAll(nbr, j)
                {
                    grown[nbr[j]] = true;
                }
            }

            if (Pstream::parRun())
            {
                volScalarField& marker = markerPtr();
                scalarField& m = marker.primitiveFieldRef();
                forAll(frontCell, cellI)
                {
                    m[cellI] = frontCell[cellI] ? 1.0 : 0.0;
                }
                marker.correctBoundaryConditions();

                forAll(mesh.boundary(), patchI)
                {
                    const fvPatch& p = mesh.boundary()[patchI];
                    if (!p.coupled())
                    {
                        continue;
                    }

                    const scalarField nbrMarker
                    (
                        marker.boundaryField()[patchI].patchNeighbourField()
                    );
                    const labelUList& faceCells = p.faceCells();

                    forAll(faceCells, i)
                    {
                        if (nbrMarker[i] > 0.5)
                        {
                            grown[faceCells[i]] = true;
                        }
                    }
                }
            }

            frontCell = grown;
        }

        // Added for cardiacFoam: the sharpest available check that the front
        // classification is partition-independent. It is an integer, so it must
        // match EXACTLY between rank counts; before the exchange above it
        // differed by roughly the number of cells adjacent to a cut.
        if (debugFrontCellCount)
        {
            label nFront = 0;
            forAll(frontCell, cellI)
            {
                if (frontCell[cellI]) ++nFront;
            }
            reduce(nFront, sumOp<label>());
            Info<< "    [front] cells = " << nFront << endl;
        }
    }

    // Detect the first time Vm crosses the activation threshold in each
    // cell and record it via linear interpolation across the time step:
    //
    //     w = (Vthr - Vm^n) / (Vm^{n+1} - Vm^n)   in [0, 1]
    //     activationTime = t_n + w * dt
    //
    // Once a cell is marked activated its flag in calculateActivationTime
    // is cleared so subsequent threshold re-crossings (e.g. recovery and
    // re-excitation) are ignored.
    void updateActivationTimes
    (
        const scalarField& VmOld,
        const volScalarField& Vm,
        volScalarField& activationTime,
        boolList& calculateActivationTime,
        const scalar t0,
        const scalar dt,
        const scalar activationThreshold
    )
    {
        const scalarField& VmI = Vm.primitiveField();
        scalarField& activationTimeI = activationTime.primitiveFieldRef();

        forAll(activationTimeI, cellI)
        {
            if (!calculateActivationTime[cellI])
            {
                continue;
            }

            if (VmOld[cellI] < activationThreshold && VmI[cellI] >= activationThreshold)
            {
                const scalar denom = VmI[cellI] - VmOld[cellI];
                scalar w = 0.0;
                if (mag(denom) > VSMALL)
                {
                    w = (activationThreshold - VmOld[cellI])/denom;
                }

                w = min(max(w, scalar(0.0)), scalar(1.0));
                activationTimeI[cellI] = t0 + w*dt;
                calculateActivationTime[cellI] = false;
            }
        }

        activationTime.correctBoundaryConditions();
    }

    // Index of the cell whose centre is nearest a sample point, or -1 on the
    // ranks that do not own it. Used to locate benchmark probe points.
    label nearestCellToPoint
    (
        const fvMesh& mesh,
        const point& samplePoint
    )
    {
        const vectorField& C = mesh.C();
        label bestCellI = -1;
        scalar bestD2 = GREAT;

        forAll(C, cellI)
        {
            const scalar d2 = magSqr(C[cellI] - samplePoint);
            if (d2 < bestD2)
            {
                bestD2 = d2;
                bestCellI = cellI;
            }
        }

        // Modified for cardiacFoam: the search above is over LOCAL cells only,
        // so every rank used to believe it owned the sample point. Downstream
        // that meant each rank triggered the stop condition on ITS OWN nearest
        // cell, clipped effectiveEndTime differently, and left the time loop at
        // a different step - some ranks dropping out of the collective solve
        // while the others blocked in it. A hang, preceded by corrupted physics,
        // because effectiveEndTime also sets the final step's dt.
        //
        // Exactly one rank now returns a cell; the others return -1. Ties are
        // broken by lowest rank, so the choice is deterministic and independent
        // of the decomposition. reduce() is preferred over Pstream::broadcast,
        // which broadcasts from masterNo() rather than from an arbitrary root.
        if (Pstream::parRun())
        {
            scalar globalBestD2 = bestD2;
            reduce(globalBestD2, minOp<scalar>());

            label winner =
                (bestCellI >= 0 && bestD2 <= globalBestD2)
              ? Pstream::myProcNo()
              : Pstream::nProcs();
            reduce(winner, minOp<label>());

            if (winner >= Pstream::nProcs())
            {
                FatalErrorInFunction
                    << "Could not find a cell for point " << samplePoint
                    << " on any rank" << exit(FatalError);
            }

            return (winner == Pstream::myProcNo()) ? bestCellI : -1;
        }

        if (bestCellI < 0)
        {
            FatalErrorInFunction
                << "Could not find a cell for point " << samplePoint
                << exit(FatalError);
        }

        return bestCellI;
    }


    // Sample a field at an arbitrary point by inverse-distance weighting over
    // its k nearest cell centres, with an exact-hit shortcut. This is how the
    // Niederer benchmark reads activation times at points that are not cell
    // centres.
    scalar sampleIDW
    (
        const volScalarField& field,
        const point& samplePoint,
        const label requestedNeighbours,
        const scalar requestedPower
    )
    {
        const vectorField& C = field.mesh().C();
        const scalarField& values = field.primitiveField();

        // Modified for cardiacFoam: the k nearest cells must be the k nearest
        // GLOBALLY. The search below is over local cells, so on N ranks this
        // used to return "the k nearest on my rank" - N different answers, none
        // of them right unless the sample point happened to sit well inside one
        // partition.
        //
        // Each rank builds its own top-k, the candidates are gathered, and every
        // rank merges the same set with the same ordering, so all ranks agree by
        // construction. Ties are broken by global cell ID, which reproduces the
        // serial result exactly: the serial loop tests d2 < bestD2 strictly, so
        // among equal distances it keeps the LOWEST cell index, and in serial
        // the global ID is the cell index.
        //
        // Note the early return on an exact hit has moved AFTER the merge. Left
        // inside the search loop it would make one rank return immediately while
        // the others entered the gather, which deadlocks.
        const label nGlobalCells = returnReduce(C.size(), sumOp<label>());
        if (nGlobalCells == 0)
        {
            return 0.0;
        }

        const label nNbrs =
            min(max(requestedNeighbours, label(1)), nGlobalCells);
        const scalar power = max(requestedPower, scalar(SMALL));

        List<scalar> bestD2(nNbrs, GREAT);
        List<label> bestId(nNbrs, -1);
        List<scalar> bestVal(nNbrs, 0.0);

        forAll(C, cellI)
        {
            const scalar d2 = magSqr(C[cellI] - samplePoint);

            for (label i = 0; i < nNbrs; ++i)
            {
                if (d2 < bestD2[i])
                {
                    for (label j = nNbrs - 1; j > i; --j)
                    {
                        bestD2[j] = bestD2[j - 1];
                        bestId[j] = bestId[j - 1];
                        bestVal[j] = bestVal[j - 1];
                    }

                    bestD2[i] = d2;
                    bestId[i] = gCol(cellI);
                    bestVal[i] = values[cellI];
                    break;
                }
            }
        }

        // Added for cardiacFoam: merge every rank's candidates identically.
        if (Pstream::parRun())
        {
            List<List<scalar>> allD2(Pstream::nProcs());
            List<List<label>> allId(Pstream::nProcs());
            List<List<scalar>> allVal(Pstream::nProcs());
            allD2[Pstream::myProcNo()] = bestD2;
            allId[Pstream::myProcNo()] = bestId;
            allVal[Pstream::myProcNo()] = bestVal;
            Pstream::allGatherList(allD2);
            Pstream::allGatherList(allId);
            Pstream::allGatherList(allVal);

            DynamicList<label> order;
            DynamicList<scalar> d2Flat;
            DynamicList<label> idFlat;
            DynamicList<scalar> valFlat;
            forAll(allD2, procI)
            {
                forAll(allD2[procI], i)
                {
                    if (allId[procI][i] < 0)
                    {
                        continue;
                    }
                    order.append(d2Flat.size());
                    d2Flat.append(allD2[procI][i]);
                    idFlat.append(allId[procI][i]);
                    valFlat.append(allVal[procI][i]);
                }
            }

            // Sort by (distance, global cell ID). The ID makes the order total,
            // so every rank produces the same list.
            labelList idx(order);
            std::stable_sort
            (
                idx.begin(),
                idx.end(),
                [&](const label a, const label b)
                {
                    if (d2Flat[a] != d2Flat[b]) return d2Flat[a] < d2Flat[b];
                    return idFlat[a] < idFlat[b];
                }
            );

            bestD2 = List<scalar>(nNbrs, GREAT);
            bestId = List<label>(nNbrs, -1);
            bestVal = List<scalar>(nNbrs, 0.0);
            const label nKeep = min(nNbrs, label(idx.size()));
            for (label i = 0; i < nKeep; ++i)
            {
                bestD2[i] = d2Flat[idx[i]];
                bestId[i] = idFlat[idx[i]];
                bestVal[i] = valFlat[idx[i]];
            }
        }

        // Exact hit on a cell centre: return that cell's value. Deferred to
        // here so that every rank has already taken part in the gather.
        if (bestId[0] >= 0 && bestD2[0] < VSMALL)
        {
            return bestVal[0];
        }

        scalar wSum = 0.0;
        scalar vSum = 0.0;
        forAll(bestId, i)
        {
            if (bestId[i] < 0)
            {
                continue;
            }

            const scalar d = std::sqrt(max(bestD2[i], scalar(VSMALL)));
            const scalar w = 1.0/std::pow(d, power);
            wSum += w;
            vSum += w*bestVal[i];
        }

        return vSum/max(wSum, scalar(SMALL));
    }

    // Write the benchmark activation-time outputs: the value at the P8 corner
    // point and a profile sampled along the diagonal of the slab, both by
    // inverse-distance weighting, into the solver's postProcessing directory.
    void writeActivationSamples
    (
        const Time& runTime,
        const volScalarField& activationTime,
        const point& diagonalStart,
        const point& diagonalEnd,
        const point& p8Point,
        const label nDiagonalSamplesRequested,
        const label sampleNeighbours,
        const scalar samplePower,
        const bool verbose = true
    )
    {
        // Modified for cardiacFoam: see postProcessingDir().
        const fileName sampleDir(postProcessingDir(runTime));
        mkDir(sampleDir);

        const label nDiagonalSamples = max(nDiagonalSamplesRequested, label(2));

        const scalar p8Time = sampleIDW
        (
            activationTime,
            p8Point,
            sampleNeighbours,
            samplePower
        );

        OFstream p8File(sampleDir/"P8_activationTime.dat");
        p8File<< "# label x y z activationTime_s activationTime_ms" << nl;
        p8File<< "P8 "
              << p8Point.x() << ' ' << p8Point.y() << ' ' << p8Point.z() << ' '
              << p8Time << ' ' << 1000.0*p8Time << nl;

        OFstream diagFile(sampleDir/"diagonalActivationTime.dat");
        diagFile<< "# index x y z distance_m distance_mm activationTime_s activationTime_ms" << nl;

        const vector diag = diagonalEnd - diagonalStart;
        for (label i = 0; i < nDiagonalSamples; ++i)
        {
            const scalar s = scalar(i)/scalar(nDiagonalSamples - 1);
            const point p = diagonalStart + s*diag;
            const scalar distance = mag(p - diagonalStart);
            const scalar activation = sampleIDW
            (
                activationTime,
                p,
                sampleNeighbours,
                samplePower
            );

            diagFile<< i << ' '
                    << p.x() << ' ' << p.y() << ' ' << p.z() << ' '
                    << distance << ' ' << 1000.0*distance << ' '
                    << activation << ' ' << 1000.0*activation << nl;
        }

        if (verbose)
        {
            Info<< "P8 activation time = " << p8Time << " s ("
                << 1000.0*p8Time << " ms)" << nl
                << "Activation samples written to " << sampleDir << nl;
        }
    }
}

int main(int argc, char* argv[])
{
#ifdef __GLIBC__
    mallopt(M_ARENA_MAX, 2);
    mallopt(M_MMAP_THRESHOLD, 64*1024);
    mallopt(M_TRIM_THRESHOLD, 64*1024);
#endif

    #include "setRootCaseLists.H"
    PetscSession petscSession(argc, argv);
    #include "createTime.H"
    #include "createMesh.H"
    // Added for cardiacFoam: build the global cell addressing before anything
    // assembles a matrix. In serial toGlobal() is the identity and totalSize()
    // is mesh.nCells(), so every column and dimension below is unchanged.
    gCellsPtr.reset(new globalIndex(mesh.nCells()));

    // Added for cardiacFoam: this rank's first global row. Zero in serial.
    gRowStart = gCellsPtr->localStart();

    if (Pstream::parRun())
    {
        Info<< "Distributed: " << gCellsPtr->totalSize() << " cells over "
            << Pstream::nProcs() << " ranks" << endl;
    }

    #include "createFields.H"

    const auto tStartTotal = std::chrono::steady_clock::now();
    const auto tStartSetup = tStartTotal;

    const label nStates = ionicModel->nEqns();
    EmbeddedTNNPModel* tnnpModel = &ionicModel();

    // Front-aware hybrid: Gauss-point ODEs in "front" cells, cell-centre ODE +
    // reconstruction elsewhere. It carries the states at the Gauss points (like
    // gaussPointODE), so it is excluded from statesAtCells below.
    const bool useFrontHybrid =
        frontAwareHybrid
     && useHighOrder_Iion
     && mesh.nGeometricD() > 1;

    // stateIntegrationMode = cellCentredReconstruct keeps the persistent states
    // at cell centres (1 ODE/cell) and reconstructs them at the Iion Gauss
    // points on demand; gaussPointODE (and the hybrid) keep them at the Gauss
    // points.
    const bool statesAtCells =
        reconstructStatesFromCellCentres
     && !useFrontHybrid
     && useHighOrder_Iion
     && mesh.nGeometricD() > 1;

    Field<Field<scalar>> states
    (
        statesAtCells ? mesh.nCells() : totalIionIntegrationPoints,
        Field<scalar>(nStates, 0.0)
    );
    ionicModel->initialiseStates(states);

    // Transient IP-sized buffer + scratch field used when reconstructing the
    // cell-centred states at the Iion Gauss points (cellCentredReconstruct).
    Field<Field<scalar>> statesIP
    (
        statesAtCells ? totalIionIntegrationPoints : 0,
        Field<scalar>(nStates, 0.0)
    );

    // Front-aware hybrid working set (per-cell front flag and a reusable
    // cell-centred state buffer for the smooth-cell ODE + reconstruction).
    boolList frontCell(useFrontHybrid ? mesh.nCells() : 0, false);
    Field<Field<scalar>> statesCellBuf
    (
        useFrontHybrid ? mesh.nCells() : 0,
        Field<scalar>(nStates, 0.0)
    );

    autoPtr<volScalarField> stateScratchPtr;
    if (statesAtCells || useFrontHybrid)
    {
        stateScratchPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "stateScratch",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("stateScratch", dimless, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }

    scalarField VmIntegrationPoints(totalIionIntegrationPoints, 0.0);
    scalarField IionIntegrationPoints(totalIionIntegrationPoints, 0.0);

    const label dim = max(mesh.nGeometricD(), label(1));
    const scalar dx = characteristicDx(mesh);
    const scalar CFL = timeIntegrationProperties.lookupOrDefault<scalar>("CFL", 0.1);
    const scalar dtExplicitReference = computeStableDeltaT
    (
        conductivity,
        chi.value(),
        Cm.value(),
        CFL,
        dx,
        dim
    );

    scalar dt = runTime.deltaTValue();
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

    Info<< "Running high-order electro activation implicit PETSc solver" << nl
        << "Dimension = " << dim << nl
        << "dx = " << dx << nl
        << "requested dt = " << dt << nl
        << "explicit stable dt reference = " << dtExplicitReference << nl
        << "implicitScheme = " << implicitScheme << nl
        << "massMatrix = " << massMatrixMode << nl
        << "linearSolverBackend = " << linearSolverBackend << nl
        << "PETSc linear KSP = " << petscLinearKspType
        << ", PC = " << petscLinearPcType << nl
        << "useHighOrder_Vm = " << useHighOrder_Vm << nl
        << "useHighOrder_Iion = " << useHighOrder_Iion << nl
        << "stateIntegrationMode = " << stateIntegrationMode << nl
        << "frontAwareHybrid = " << frontAwareHybrid
        << (frontAwareHybrid ? " (hybridFrontVmRange = " : "")
        << (frontAwareHybrid ? Foam::name(hybridFrontVmRange) : word(""))
        << (frontAwareHybrid ? " V, hybridFrontDilation = " : "")
        << (frontAwareHybrid ? Foam::name(hybridFrontDilation) : word(""))
        << (frontAwareHybrid ? ")" : "") << nl
        << "stateReconstructLimiter = " << stateReconstructLimiter << nl
        << "stabilisationAlpha = " << stabilisationAlpha << nl
        << "memoryOptimization = " << memoryOptimization
        << " (effective: "
        << (memoryOptimizationEffective ? "true" : "false") << ")" << nl
        << "stiffnessTripletsPerFaceReserve = "
        << stiffnessTripletsPerFaceReserve << nl
        << "compactMassAssembly = "
        << (compactMassAssembly ? "true" : "false") << nl
        << "Iion integration points = " << totalIionIntegrationPoints << nl
        << "nonlinearMethod = " << nonlinearMethod << nl
        << "JFNK linear backend = " << jfnkLinearSolverBackend << nl
        << "JFNK PETSc KSP = " << jfnkPetscKspType
        << ", PC = " << jfnkPetscPcType << nl
        << "state ODE solver = " << timeIntegrationProperties.lookupOrDefault<word>("solver", "default") << nl
        << "nonlinearRelaxation = " << nonlinearRelaxation << nl
        << "nonlinearStatesTolerance = " << nonlinearStatesTolerance << nl
        << "nonlinearRequireStatesConvergence = "
        << nonlinearRequireStatesConvergence << nl
        << "activationThreshold = " << activationThreshold << " V" << nl
        << "stopAfterPointActivation = " << stopAfterPointActivation << nl
        << "stopActivationPoint = " << stopActivationPoint << nl
        << "stopDelayAfterActivation = " << stopDelayAfterActivation << " s" << nl
        << "implicitNonlinearIterations = " << implicitNonlinearIterations << nl
        << "implicitLinearSolver = " << implicitLinearSolver << nl << endl;

    if (statesAtCells)
    {
        Info<< "LREInterp_states order = " << LREInterp_statesPtr().order()
            << " (states evolved at cell centres, reconstructed at Gauss points)"
            << endl;
    }
    if (useFrontHybrid)
    {
        Info<< "LREInterp_states order = " << LREInterp_statesPtr().order()
            << " (front-aware hybrid: Gauss-point ODEs on the front,"
            << " cell-centre + reconstruction elsewhere)" << endl;
    }

    const label stopActivationCellI =
        stopAfterPointActivation
      ? nearestCellToPoint(mesh, stopActivationPoint)
      : -1;

    scalar effectiveEndTime = runTime.endTime().value();
    bool stopPointTriggered = false;
    scalar stopPointActivationTime = -GREAT;

    if (stopAfterPointActivation)
    {
        Info<< "Stop point nearest cell = " << stopActivationCellI
            << ", centre = " << mesh.C()[stopActivationCellI]
            << ", distance = "
            << mag(mesh.C()[stopActivationCellI] - stopActivationPoint)
            << " m" << nl << endl;
    }

    const wordList exportNames = ionicModel->exportedFieldNames();
    if (!exportNames.empty())
    {
        Info<< "Exporting ionic fields: " << exportNames << nl;
    }

    SpMat M;
    SpMat K;
    SpMat L;

    if (profileTimings)
    {
        logMemoryCheckpoint("before mass assembly");
    }

    if (!useHighOrder_Vm || dim <= 1 || massMatrixMode == "lumped")
    {
        if (massMatrixMode == "consistent" && (!useHighOrder_Vm || dim <= 1))
        {
            Info<< "Consistent mass requested without high-order Vm; "
                << "using the diagonal/lumped mass matrix." << nl;
        }
        assembleDiagonalMassMatrix(mesh, 1.0, M);
    }
    else
    {
        assembleConsistentMassMatrixHO
        (
            mesh,
            LREInterp_VmPtr(),
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
    if (trimHeapAfterLargeSetup)
    {
        malloc_trim(0);
    }
#endif

    if (profileTimings)
    {
        logMemoryCheckpoint("before stiffness assembly");
    }

    if (useHighOrder_Vm && dim > 1)
    {
        assembleHighOrderStiffnessMatrix
        (
            mesh,
            conductivity,
            LREInterp_VmPtr(),
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
            LREInterp_VmPtr.valid() ? &LREInterp_VmPtr() : nullptr,
            stabilisationAlpha,
            K
        );
    }

    if (profileTimings)
    {
        logMemoryCheckpoint("after stiffness assembly");
    }

#ifdef __GLIBC__
    if (trimHeapAfterLargeSetup)
    {
        malloc_trim(0);
    }
#endif

    L = K;
    L *= 1.0/(chi.value()*Cm.value());

    autoPtr<OFstream> nonlinearResidualFilePtr;
    if (writeNonlinearResiduals)
    {
        // Modified for cardiacFoam: see postProcessingDir().
        const fileName residualDir(postProcessingDir(runTime));
        mkDir(residualDir);
        nonlinearResidualFilePtr.reset
        (
            new OFstream(residualDir/"nonlinearResiduals.dat")
        );
        nonlinearResidualFilePtr()
            << "# time_s step nonlinearMethod stateODESolver iter linearIterations linearError "
            << "coupledResidual Vm_relL2 Iion_relL2 maxState_relL2 "
            << "lineSearchIters converged";
        for (label stateI = 0; stateI < nStates; ++stateI)
        {
            nonlinearResidualFilePtr() << " state" << stateI << "_relL2";
        }
        nonlinearResidualFilePtr() << nl;
    }

    PetscShellKspSolver jfnkPetscShellSolver;

    const auto tEndSetup = std::chrono::steady_clock::now();
    const auto tStartLoop = tEndSetup;

    label nSteps = 0;

    // ----------------------------------------------------------------- //
    // State for the linear-extrapolation Newton initial guess.
    //
    //   Vm_init = Vm_n + (Vm_n - Vm_{n-1}) * (dt / dt_prev)
    //
    // Requires Vm from the previous time step (VmTwoStepsAgo) and the
    // previous time-step length (dtPrev). On the first step we have no
    // history so we fall back to x_0 = Vm_n (legacy behaviour).
    // ----------------------------------------------------------------- //
    scalarField VmTwoStepsAgo(mesh.nCells(), 0.0);
    scalar dtPrev = 0.0;
    bool hasPrevStep = false;

    // Live activation output: write the diagonal/P8 samples once up front (so
    // the file exists from t=0 with every point still pending) and then again
    // whenever a new point crosses the activation threshold, so the run can be
    // monitored as it progresses. writeActivationSamples() overwrites the files.
    label nActivatedPrev = 0;
    forAll(calculateActivationTime, cellI)
    {
        if (!calculateActivationTime[cellI])
        {
            ++nActivatedPrev;
        }
    }
    reduce(nActivatedPrev, sumOp<label>());

    // Activation-failure tracking (see abortOnActivationFailure). The status is
    // written to postProcessing/.../activationStatus.dat at the end.
    scalar tLastActivationChange = runTime.value();
    word activationStatus("RUNNING");
    bool activationAborted = false;

    writeActivationSamples
    (
        runTime,
        activationTime,
        diagonalStart,
        diagonalEnd,
        p8Point,
        nDiagonalSamples,
        sampleNeighbours,
        samplePower,
        false
    );

    // ----------------------------------------------------------------------- //
    // Phase B: persistent linear-solver state, reused across time steps and
    // nonlinear iterations. AImplicit/BImplicit only depend on the (constant)
    // M, L and theta and on dt, which is fixed except on the final clipped
    // step; they are therefore assembled once and rebuilt only when dt changes.
    // The factorisations (PETSc KSP or Eigen) live for the whole run and are
    // refreshed according to the MatrixUpdate policy computed each solve.
    // ----------------------------------------------------------------------- //
    SpMat AImplicit;
    SpMat BImplicit;

    // Added for cardiacFoam: PETSc-backed operators for the two products that
    // used to be plain Eigen expressions. They are rebuilt only when the
    // matrices are (on the first step and whenever dt changes), so the extra
    // PETSc assembly is not on the per-iteration path.
    DistributedMatVec applyAImplicitOp;
    DistributedMatVec applyBImplicitOp;

    auto applyAImplicit = [&](const EigVec& x) { return applyAImplicitOp(x); };
    auto applyBImplicit = [&](const EigVec& x) { return applyBImplicitOp(x); };
    scalar assembledDt = -1.0;
    bool linearSolverNeedsRebuild = true;
    PetscKspMatrixSolver persistentPetscSolver;
    PersistentEigenSolvers persistentEigenSolvers;

    while (runTime.value() < effectiveEndTime - SMALL)
    {
        const scalar t0 = runTime.value();
        const scalar remaining = effectiveEndTime - t0;

        if (remaining <= SMALL)
        {
            break;
        }

        dt = min(runTime.deltaTValue(), remaining);
        runTime.setDeltaT(dt);

        // ----------------------------------------------------------------- //
        // Build the theta-scheme implicit operators.
        //
        //   The monodomain PDE in semi-discrete form is
        //
        //       dVm/dt = L * Vm + s(Vm, states)
        //
        //   where L = K / (chi * Cm) is the discrete diffusion operator and
        //   s = -Iion/(chi*Cm) + Iext/(chi*Cm) lumps the reactive sources.
        //   Applying the theta rule (theta=1: BE, theta=0.5: CN) and dividing
        //   by dt yields
        //
        //       AImplicit * Vm^{n+1} = BImplicit * Vm^n + source-blend
        //
        //   with  AImplicit = M/dt - theta*L
        //         BImplicit = M/dt + (1-theta)*L
        //
        //   M is either the lumped (diagonal) or consistent high-order mass
        //   matrix; L already carries the 1/(chi*Cm) scaling.
        //
        //   Phase B: AImplicit/BImplicit are constant while dt is unchanged, so
        //   they are only reassembled on the first step and whenever dt changes
        //   (the final clipped step). A dt change also forces the persistent
        //   factorisation to be rebuilt on the next solve.
        // ----------------------------------------------------------------- //
        if (mag(dt - assembledDt) > SMALL)
        {
            AImplicit = M;
            AImplicit *= (1.0/dt);
            AImplicit -= theta*L;

            BImplicit = M;
            BImplicit *= (1.0/dt);
            if (theta < 1.0 - SMALL)
            {
                // Crank-Nicolson (or any theta < 1) needs the explicit half of
                // the diffusion operator on the right-hand side. For pure
                // Backward Euler (theta=1) this term vanishes and is skipped.
                BImplicit += (1.0 - theta)*L;
            }

            // Added for cardiacFoam: mirror the freshly assembled matrices
            // into PETSc. Must follow every reassembly, which is exactly here.
            reportOperatorFingerprint("AImplicit", AImplicit);
            reportOperatorFingerprint("BImplicit", BImplicit);

            applyAImplicitOp.reset(AImplicit);
            applyBImplicitOp.reset(BImplicit);

            assembledDt = dt;
            linearSolverNeedsRebuild = true;
        }

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

        const scalarField VmOldValues(Vm.primitiveField());
        Field<Field<scalar>> statesOld(states);

        volScalarField IionOld
        (
            IOobject
            (
                "IionOld",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Iion
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

        applyStimulus
        (
            t0,
            externalStimulusCurrent,
            stimulusCellIDsList,
            stimulusStartTimes,
            stimulusIntensity.value(),
            stimulusDuration.value()
        );

        if (useHighOrder_Iion)
        {
            reconstructVmAtIionIntegrationPoints
            (
                VmOld,
                useHighOrder_Vm,
                useHighOrder_Vm ? &LREInterp_VmPtr() : nullptr,
                LREInterp_IionPtr(),
                VmIntegrationPoints
            );

            if (statesAtCells)
            {
                // statesOld lives at cell centres -> reconstruct to the Gauss
                // points before evaluating Iion there.
                reconstructStatesAtIionIntegrationPoints
                (
                    statesOld,
                    nStates,
                    LREInterp_statesPtr(),
                    LREInterp_IionPtr(),
                    stateScratchPtr(),
                    statesIP,
                    nullptr,
                    stateReconstructLimiter
                );
            }

            ionicModel->calculateCurrent
            (
                t0,
                dt,
                VmIntegrationPoints,
                IionIntegrationPoints,
                statesAtCells ? statesIP : statesOld
            );

            averageIionIntegrationPointsToCells
            (
                IionIntegrationPoints,
                LREInterp_IionPtr(),
                IionOld
            );
        }
        else
        {
            ionicModel->calculateCurrent
            (
                t0,
                dt,
                VmOld.internalField(),
                IionOld,
                statesOld
            );
            IionOld.correctBoundaryConditions();
        }

        // Front-aware hybrid: classify the front cells ONCE per time step (from
        // the committed Vm, reconstructed above into VmIntegrationPoints) and
        // freeze the flag for the whole nonlinear solve. Dilating a couple of
        // cell layers covers the front's motion within the step and keeps the
        // fixed-point map continuous (no front/smooth flip between iterations).
        if (useFrontHybrid)
        {
            flagFrontCells
            (
                mesh,
                VmIntegrationPoints,
                LREInterp_IionPtr(),
                externalStimulusCurrent,
                hybridFrontVmRange,
                frontCell
            );
            dilateFrontCells(mesh, hybridFrontDilation, frontCell);
        }

        scalarField sourceOld(mesh.nCells(), 0.0);
        scalarField sourceGuess(mesh.nCells(), 0.0);
        forAll(sourceOld, cellI)
        {
            sourceOld[cellI] =
               -IionOld[cellI]
              + externalStimulusCurrent[cellI]/(chi.value()*Cm.value());
        }

        const EigVec Vn = fieldToEigVec(VmOld);

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

        // ------------------------------------------------------------- //
        // Newton initial guess: linear extrapolation Vm_n + (Vm_n -
        // Vm_{n-1}) * (dt / dt_prev), clamped to a physical bracket.
        // First time step has no history -> identity (Vm_n).
        // ------------------------------------------------------------- //
        if (jfnkInitGuessOrder >= 1 && hasPrevStep && dtPrev > SMALL)
        {
            scalarField& Vg = VmGuess.primitiveFieldRef();
            const scalarField& Vn_ = VmOld.primitiveField();
            const scalar ratio = dt / dtPrev;
            forAll(Vg, cellI)
            {
                scalar v = Vn_[cellI] + ratio*(Vn_[cellI] - VmTwoStepsAgo[cellI]);
                if (v < jfnkInitGuessVmMin) v = jfnkInitGuessVmMin;
                if (v > jfnkInitGuessVmMax) v = jfnkInitGuessVmMax;
                Vg[cellI] = v;
            }
            VmGuess.correctBoundaryConditions();
        }

        // Snapshot Vm_n into VmTwoStepsAgo *before* Newton mutates the
        // field, so the next iteration of the time loop sees the correct
        // "previous step" Vm even after rollback.
        VmTwoStepsAgo = VmOld.primitiveField();
        dtPrev = dt;
        hasPrevStep = true;

        if (tnnpModel)
        {
            tnnpModel->updateStatesOld(statesOld);
        }

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
            Vm.boundaryField().types()
        );

        // Evaluate the coupled (Vm, states, Iion) triple for a given
        // membrane-potential candidate.
        //
        // When resetStatesToTimeStart is true the ionic state variables are
        // restored to their values at the beginning of the current time step
        // BEFORE integrating the ODE. This is required during nonlinear
        // iterations: each Newton/Picard step must integrate from t_n to
        // t_{n+1} starting from the same state, otherwise the iterates
        // drift along the ODE trajectory instead of converging.
        auto evaluateNonlinearFields =
        [&](
            const volScalarField& VmCandidate,
            Field<Field<scalar>>& statesCandidate,
            volScalarField& IionCandidate,
            const bool resetStatesToTimeStart
        )
        {
            if (resetStatesToTimeStart)
            {
                // Generic reset: works for every ionic model since `states`
                // is the externally-visible state vector field.
                statesCandidate = statesOld;

                if (tnnpModel)
                {
                    // TNNP additionally keeps its own internal state snapshot
                    // (STATES_OLD_) used for stiff-gating protection; refresh
                    // it here so the next ODE step starts from a consistent
                    // pair of (external, internal) states.
                    tnnpModel->resetStatesToStatesOld(statesCandidate);
                }
            }

            if (useHighOrder_Iion)
            {
                reconstructVmAtIionIntegrationPoints
                (
                    VmCandidate,
                    useHighOrder_Vm,
                    useHighOrder_Vm ? &LREInterp_VmPtr() : nullptr,
                    LREInterp_IionPtr(),
                    VmIntegrationPoints
                );

                // Protect TNNP from unphysical probes coming from the
                // matrix-free JFNK Jacobian. VmIntegrationPoints is a
                // working buffer for the Iion evaluation so we clamp it in
                // place; this does not affect VmCandidate or any persistent
                // field.
                if (jfnkClampODEInput)
                {
                    clampPhysical
                    (
                        VmIntegrationPoints,
                        jfnkInitGuessVmMin,
                        jfnkInitGuessVmMax
                    );
                }

                if (useFrontHybrid)
                {
                    // ------------------------------------------------------- //
                    // Front-aware hybrid. statesCandidate is Gauss-sized.
                    //   front cells  : one ODE per Gauss point (local Vm),
                    //   smooth cells : one ODE at the cell centre + LRE
                    //                  reconstruction of the states to Gauss.
                    // VmIntegrationPoints is already (optionally) clamped above.
                    // ------------------------------------------------------- //
                    const CompactListList<point>& cellQP =
                        LREInterp_IionPtr().cellQuadPoints();

                    // (a) frontCell is FROZEN for the whole time step (computed
                    //     once before the nonlinear loop from the committed Vm,
                    //     then dilated). It must NOT be recomputed here: doing so
                    //     per iteration flips the front/smooth switch and makes
                    //     the fixed point a non-convergent period-2 limit cycle.

                    // (b) statesCellBuf <- cell-average of the reset states.
                    averageStatesIntegrationPointsToCells
                    (
                        statesCandidate, nStates, LREInterp_IionPtr(),
                        statesCellBuf
                    );

                    // (c) Front cells: gather their Gauss points, advance one
                    //     ODE per point with the local Vm, scatter back.
                    DynamicList<label> frontIPs;
                    {
                        label ip = 0;
                        forAll(mesh.cells(), cellI)
                        {
                            const label nq = cellQP[cellI].size();
                            if (frontCell[cellI])
                            {
                                for (label k = 0; k < nq; ++k)
                                {
                                    frontIPs.append(ip + k);
                                }
                            }
                            ip += nq;
                        }
                    }
                    if (frontIPs.size() > 0)
                    {
                        const label nF = frontIPs.size();
                        scalarField VmF(nF), ImF(nF, 0.0);
                        Field<Field<scalar>> sF(nF);
                        forAll(frontIPs, i)
                        {
                            VmF[i] = VmIntegrationPoints[frontIPs[i]];
                            sF[i] = statesCandidate[frontIPs[i]];
                        }
                        // Modified for cardiacFoam: frontIPs are the true Gauss
                        // point indices of these front points, so passing them
                        // as the scratch map puts each point's algebraic_,
                        // rates_ and stepMs_ in the slot that belongs to it.
                        // Previously this pass used slots 0..nF-1, which the
                        // smooth-cell pass below then overwrote.
                        ionicModel->solveODE(t0, dt, VmF, ImF, sF, frontIPs);
                        forAll(frontIPs, i)
                        {
                            statesCandidate[frontIPs[i]] = sF[i];
                        }
                    }

                    // (d) Smooth cells: one ODE at the cell-centre Vm, from the
                    //     cell-average old state; result overwrites statesCellBuf.
                    boolList smoothCell(frontCell.size());
                    forAll(frontCell, cellI)
                    {
                        smoothCell[cellI] = !frontCell[cellI];
                    }
                    DynamicList<label> smoothCells;
                    forAll(smoothCell, cellI)
                    {
                        if (smoothCell[cellI]) smoothCells.append(cellI);
                    }
                    if (smoothCells.size() > 0)
                    {
                        const label nS = smoothCells.size();
                        const scalarField& VmC = VmCandidate.internalField();
                        scalarField VmS(nS), ImS(nS, 0.0);
                        Field<Field<scalar>> sS(nS);
                        forAll(smoothCells, i)
                        {
                            scalar v = VmC[smoothCells[i]];
                            if (jfnkClampODEInput)
                            {
                                if (v < jfnkInitGuessVmMin) v = jfnkInitGuessVmMin;
                                else if (v > jfnkInitGuessVmMax) v = jfnkInitGuessVmMax;
                            }
                            VmS[i] = v;
                            sS[i] = statesCellBuf[smoothCells[i]];
                        }

                        // Added for cardiacFoam: these are cell-centred points,
                        // so they take slots from the second scratch bank. The
                        // map is injective and disjoint from the front pass's
                        // Gauss slots, so the two passes no longer alias.
                        labelList smoothSlots(nS);
                        forAll(smoothCells, i)
                        {
                            smoothSlots[i] =
                                ionicModel->cellSlotOffset() + smoothCells[i];
                        }

                        ionicModel->solveODE(t0, dt, VmS, ImS, sS, smoothSlots);
                        forAll(smoothCells, i)
                        {
                            statesCellBuf[smoothCells[i]] = sS[i];
                        }
                    }

                    // (e) Front cells' new cell state = average of their advanced
                    //     Gauss states (needed as a smooth reconstruction stencil).
                    averageStatesIntegrationPointsToCells
                    (
                        statesCandidate, nStates, LREInterp_IionPtr(),
                        statesCellBuf, &frontCell
                    );

                    // (f) Reconstruct the new cell states to the Gauss points in
                    //     smooth cells only (front Gauss states already advanced).
                    reconstructStatesAtIionIntegrationPoints
                    (
                        statesCellBuf,
                        nStates,
                        LREInterp_statesPtr(),
                        LREInterp_IionPtr(),
                        stateScratchPtr(),
                        statesCandidate,
                        &smoothCell,
                        stateReconstructLimiter
                    );

                    // (g) Iion at all Gauss points from the assembled states.
                    ionicModel->calculateCurrent
                    (
                        t0,
                        dt,
                        VmIntegrationPoints,
                        IionIntegrationPoints,
                        statesCandidate
                    );
                }
                else if (statesAtCells)
                {
                    // cellCentredReconstruct: one ODE per cell (driven by the
                    // cell-centred Vm), then reconstruct the states at the
                    // Gauss points and evaluate Iion there (compute-only).
                    if (jfnkClampODEInput)
                    {
                        scalarField VmClamped(VmCandidate.internalField());
                        clampPhysical
                        (
                            VmClamped,
                            jfnkInitGuessVmMin,
                            jfnkInitGuessVmMax
                        );
                        ionicModel->solveODE
                        (
                            t0, dt, VmClamped, IionCandidate, statesCandidate
                        );
                    }
                    else
                    {
                        ionicModel->solveODE
                        (
                            t0,
                            dt,
                            VmCandidate.internalField(),
                            IionCandidate,
                            statesCandidate
                        );
                    }

                    reconstructStatesAtIionIntegrationPoints
                    (
                        statesCandidate,
                        nStates,
                        LREInterp_statesPtr(),
                        LREInterp_IionPtr(),
                        stateScratchPtr(),
                        statesIP,
                        nullptr,
                        stateReconstructLimiter
                    );

                    ionicModel->calculateCurrent
                    (
                        t0,
                        dt,
                        VmIntegrationPoints,
                        IionIntegrationPoints,
                        statesIP
                    );
                }
                else
                {
                    // gaussPointODE: integrate the ODE independently at every
                    // Gauss point.
                    ionicModel->solveODE
                    (
                        t0,
                        dt,
                        VmIntegrationPoints,
                        IionIntegrationPoints,
                        statesCandidate
                    );
                }

                averageIionIntegrationPointsToCells
                (
                    IionIntegrationPoints,
                    LREInterp_IionPtr(),
                    IionCandidate
                );
            }
            else
            {
                if (jfnkClampODEInput)
                {
                    // VmCandidate is const, so we copy into a clamped
                    // scratch buffer before invoking the ODE. The copy is
                    // O(nCells) and negligible compared to the integrator
                    // sub-steps inside solveODE.
                    scalarField VmClamped(VmCandidate.internalField());
                    clampPhysical
                    (
                        VmClamped,
                        jfnkInitGuessVmMin,
                        jfnkInitGuessVmMax
                    );
                    ionicModel->solveODE
                    (
                        t0,
                        dt,
                        VmClamped,
                        IionCandidate,
                        statesCandidate
                    );
                }
                else
                {
                    ionicModel->solveODE
                    (
                        t0,
                        dt,
                        VmCandidate.internalField(),
                        IionCandidate,
                        statesCandidate
                    );
                }
                IionCandidate.correctBoundaryConditions();
            }
        };

        // Build the right-hand side of the implicit system for a given Iion
        // candidate.
        //
        //   rhs = BImplicit * Vm^n  +  theta*src^{n+1}  +  (1-theta)*src^n
        //
        // No explicit dt factor is needed on the source terms: after writing
        // (Vm^{n+1} - Vm^n)/dt = ... and dividing both sides by dt, the
        // sources naturally appear with weights theta and (1-theta) only.
        // The reaction source is  src = -Iion/(chi*Cm) + Iext/(chi*Cm).
        auto rhsFromIion = [&](const volScalarField& IionCandidate)
        {
            // Modified for cardiacFoam: see applyAImplicit.
            EigVec rhs = applyBImplicit(Vn);
            forAll(sourceGuess, cellI)
            {
                const scalar sourceNp1 =
                   -IionCandidate[cellI]
                  + externalStimulusCurrent[cellI]/(chi.value()*Cm.value());

                rhs[cellI] +=
                    theta*sourceNp1
                  + (1.0 - theta)*sourceOld[cellI];
            }
            return rhs;
        };

        auto computeSourceDerivative =
        [&](
            const volScalarField& VmLinearisationPoint,
            const volScalarField& IionLinearisationPoint
        )
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
                VmLinearisationPoint
            );

            VmPlus.primitiveFieldRef() += diagonalIionEpsilon;
            VmPlus.correctBoundaryConditions();

            Field<Field<scalar>> statesPlus(statesOld);
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
                IionLinearisationPoint
            );

            evaluateNonlinearFields(VmPlus, statesPlus, IionPlus, true);

            scalarField& deriv = sourceDerivative.primitiveFieldRef();
            forAll(deriv, cellI)
            {
                deriv[cellI] =
                  -(IionPlus[cellI] - IionLinearisationPoint[cellI])
                   /max(diagonalIionEpsilon, SMALL);
            }
            sourceDerivative.correctBoundaryConditions();
        };

        const bool useDiagonalIion =
            nonlinearMethod == "diagonalIion"
         || nonlinearMethod == "diagonal"
         || nonlinearMethod == "localDiagonal";

        const bool useJFNK =
            nonlinearMethod == "JFNK"
         || nonlinearMethod == "jfnk";

        Field<Field<scalar>> statesGuess(statesOld);
        evaluateNonlinearFields(VmGuess, statesGuess, IionGuess, true);

        scalarField stateResiduals(nStates, 0.0);

        auto writeAndPrintResiduals =
        [&](
            const label corr,
            const label linearIterations,
            const scalar linearError,
            const scalar coupledResidual,
            const scalar VmResidual,
            const scalar IionResidual,
            const scalar maxStateResidual,
            const scalarField& stateResidualValues,
            const bool converged,
            const label lineSearchIters = -1  // -1 = N/A (Picard / early exit)
        )
        {
            Info<< "Time = " << t0 << "," << nl
                << "       Nonlinear method:          " << nonlinearMethod << nl
                << "       Implicit electro solver:   " << implicitLinearSolver
                << "; iterations = " << linearIterations
                << ", estimated error = " << linearError << nl
                << "       NonLinSolver:              iter = " << corr
                << ", Vm residual = " << VmResidual
                << ", Iion residual = " << IionResidual
                << ", max state residual = " << maxStateResidual
                << ", coupled residual = " << coupledResidual;
            if (lineSearchIters >= 0)
            {
                Info<< ", lineSearchIters = " << lineSearchIters;
            }
            Info<< ", converged = " << converged << nl
                << "       State residuals:";

            forAll(stateResidualValues, stateI)
            {
                Info<< " s" << stateI << '=' << stateResidualValues[stateI];
            }
            Info<< nl;

            if (writeNonlinearResiduals && nonlinearResidualFilePtr.valid())
            {
                nonlinearResidualFilePtr()
                    << t0 << ' ' << (nSteps + 1) << ' ' << nonlinearMethod << ' '
                    << timeIntegrationProperties.lookupOrDefault<word>("solver", "default") << ' '
                    << corr << ' ' << linearIterations << ' ' << linearError << ' '
                    << coupledResidual << ' '
                    << VmResidual << ' ' << IionResidual << ' '
                    << maxStateResidual << ' ' << lineSearchIters << ' ' << converged;
                forAll(stateResidualValues, stateI)
                {
                    nonlinearResidualFilePtr() << ' ' << stateResidualValues[stateI];
                }
                nonlinearResidualFilePtr() << nl;
            }
        };

        // Track whether the nonlinear solve converged in either branch
        // (JFNK or Picard/diagonal). Used after both loops to decide whether
        // to commit the iterate or roll back to the time-step-start state.
        bool nonlinearConverged = false;
        label nonlinearIters = 0;

        if (useJFNK)
        {
            EigVec x = fieldToEigVec(VmGuess);

            // ------------------------------------------------------------- //
            // Phase A3: diffusion preconditioner for the matrix-free Jacobian.
            //
            // M = AImplicit is the constant, SPD linear part of the Jacobian.
            // We LU-factorise it once (reused across time steps; refreshed only
            // when dt changes, tracked by linearSolverNeedsRebuild) and expose
            // r -> M^{-1} r as a closure. Each Krylov step then costs one cheap
            // triangular solve instead of an extra ionic-ODE integration, so the
            // matvec count (the dominant cost) drops sharply. Disabled by default
            // (jfnkPreconditioner == "none"), preserving the legacy solve path.
            // ------------------------------------------------------------- //
            const bool useJfnkPreconditioner = (jfnkPreconditioner == "diffusion");

            if (useJfnkPreconditioner && linearSolverNeedsRebuild)
            {
                // ILUT setup is cheap and has bounded fill; recomputed only when
                // dt changes (linearSolverNeedsRebuild).
                // Modified for cardiacFoam: factorise the square diagonal
                // block, not AImplicit itself - see diagonalBlock().
                persistentEigenSolvers.jfnkPrec.compute(diagonalBlock(AImplicit));
                if (persistentEigenSolvers.jfnkPrec.info() != Eigen::Success)
                {
                    FatalErrorInFunction
                        << "JFNK diffusion preconditioner: AImplicit ILUT "
                        << "factorization failed" << exit(FatalError);
                }
                persistentEigenSolvers.jfnkPrecComputed = true;
                linearSolverNeedsRebuild = false;
            }

            const std::function<EigVec(const EigVec&)> jfnkPrecApply =
                [&](const EigVec& r) -> EigVec
                {
                    return persistentEigenSolvers.jfnkPrec.solve(r);
                };
            const std::function<EigVec(const EigVec&)>* jfnkPrecApplyPtr =
                useJfnkPreconditioner ? &jfnkPrecApply : nullptr;

            // IMPORTANT: residualFor and standardResidual must use an explicit
            // `-> EigVec` return type. Without it, the deduced return type is
            // Eigen's expression-template (CwiseBinaryOp<..., Product<...>>),
            // which holds *references* to the temporaries produced by
            // fieldToEigVec(...) and rhsFromIion(...). Those temporaries die
            // when the lambda returns, leaving the expression with dangling
            // pointers; later materialisation reads garbage sizes and
            // triggers std::bad_alloc (140 TB allocations have been observed).
            // Returning by value as EigVec forces the expression to be
            // evaluated *inside* the lambda while the temporaries are alive.
            auto residualFor =
            [&](const EigVec& xCandidate, const bool updateGuess) -> EigVec
            {
                // Standard JfNK residual: F(Vm) = A*Vm - rhs(Iion(Vm))
                // This avoids an inner linear solve inside every matVec evaluation.
                auto standardResidual =
                [&](volScalarField& VmCandidate,
                    Field<Field<scalar>>& statesCandidate,
                    volScalarField& IionCandidate) -> EigVec
                {
                    evaluateNonlinearFields(VmCandidate, statesCandidate, IionCandidate, true);
                    // Modified for cardiacFoam: AImplicit has LOCAL rows and GLOBAL
                    // columns, so Eigen cannot form this product - it would
                    // need the whole global vector. PETSc's MatMult does the
                    // halo exchange internally.
                    return applyAImplicit(fieldToEigVec(VmCandidate))
                         - rhsFromIion(IionCandidate);
                };

                if (updateGuess)
                {
                    eigVecToField(xCandidate, VmGuess);
                    VmGuess.correctBoundaryConditions();
                    return standardResidual(VmGuess, statesGuess, IionGuess);
                }

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
                Field<Field<scalar>> statesTmp(statesOld);

                eigVecToField(xCandidate, VmTmp);
                VmTmp.correctBoundaryConditions();
                return standardResidual(VmTmp, statesTmp, IionTmp);
            };

            for
            (
                label corr = 0;
                corr < max(implicitNonlinearIterations, label(1));
                ++corr
            )
            {
                const scalarField VmPrevious(VmGuess.primitiveField());
                const scalarField IionPrevious(IionGuess.primitiveField());
                Field<Field<scalar>> statesPrevious(statesGuess);

                const EigVec R = residualFor(x, true);
                const scalar coupledResidual = relativeL2Norm(R, x);
                const scalar VmResidualInitial =
                    relativeL2Difference(VmGuess.primitiveField(), VmPrevious);
                const scalar IionResidualInitial =
                    relativeL2Difference(IionGuess.primitiveField(), IionPrevious);
                const scalar maxStateResidualInitial = maxStateRelativeL2Difference
                (
                    statesGuess,
                    statesPrevious,
                    stateResiduals
                );

                // Early-exit check: if the *current* iterate already
                // satisfies all tolerances, skip the GMRES solve entirely.
                // This typically fires at corr == 0 when the previous time
                // step's solution is a good initial guess.
                const bool residualConverged =
                    corr + 1 >= implicitMinNonlinearIterations
                 && coupledResidual <= nonlinearTolerance
                 && VmResidualInitial <= nonlinearVmTolerance
                 && IionResidualInitial <= nonlinearIionTolerance
                 && (
                        !nonlinearRequireStatesConvergence
                     || maxStateResidualInitial <= nonlinearStatesTolerance
                    );

                if (residualConverged)
                {
                    writeAndPrintResiduals
                    (
                        corr + 1,
                        0,
                        0.0,
                        coupledResidual,
                        VmResidualInitial,
                        IionResidualInitial,
                        maxStateResidualInitial,
                        stateResiduals,
                        true
                    );
                    nonlinearConverged = true;
                    nonlinearIters = corr + 1;
                    break;
                }

                // Finite-difference approximation to the Jacobian-vector
                // product:    J*v ~= (F(x + eps*v) - F(x)) / eps
                // The forward-difference perturbation follows the standard
                // JFNK formula  eps = sqrt(eps_mach) * (1 + ||x||) / ||v|| ,
                // which balances truncation and round-off error and is
                // invariant under scaling of v.
                //
                // The explicit `-> EigVec` return type is critical: see the
                // long comment on residualFor() above. Without it, the
                // returned CwiseBinaryOp captures references to R (safe,
                // it lives outside) AND to the temporary returned by
                // residualFor(..., false), causing dangling-pointer
                // bad_alloc inside Eigen when the result is finally evaluated.
                const scalar epsScale = jfnkEpsilon*(1.0 + gNorm(x));
                auto matVec = [&](const EigVec& v) -> EigVec
                {
                    const scalar eps = epsScale/max(gNorm(v), SMALL);
                    return (residualFor(x + eps*v, false) - R)/eps;
                };

                // Solve  J * delta = -R  with restarted GMRES.
                // The matVec closure above re-evaluates the nonlinear
                // residual every Krylov step, so each iteration triggers a
                // full ionic-ODE integration; this dominates the cost.
                label gmresIterations = 0;
                scalar gmresError = GREAT;
                EigVec delta(R.size());

                if (usesPetscBackend(jfnkLinearSolverBackend))
                {
                    if (!jfnkPetscShellSolver.isInitialised())
                    {
                        jfnkPetscShellSolver.initialise
                        (
                            R.size(),
                            jfnkPetscKspType,
                            jfnkPetscPcType,
                            jfnkPetscRestart,
                            jfnkMaxKrylovIterations*(jfnkMaxRestarts + 1),
                            jfnkLinearTolerance,
                            jfnkPetscOptionsPrefix,
                            petscUseOptions,
                            useJfnkPreconditioner
                        );
                    }

                    delta = jfnkPetscShellSolver.solve
                    (
                        std::function<EigVec(const EigVec&)>(matVec),
                        jfnkPrecApplyPtr,
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
                        gmresError,
                        jfnkPrecApplyPtr
                    );
                }

                // ----------------------------------------------------- //
                // Newton update with optional Armijo backtracking.
                //
                // The unconditional damped step  x <- x + omega * delta
                // is fragile when the linearisation of F at the current
                // iterate fails to capture the true nonlinearity (typical
                // case: cardiac AP upstroke with large dt). We optionally
                // wrap the step in a backtracking line search that halves
                // alpha until the sufficient-decrease Armijo condition
                //
                //     ||F(x + alpha*delta)||^2
                //         <= (1 - 2*c*alpha) * ||F(x)||^2
                //
                // is satisfied. This guarantees monotone descent in the
                // residual norm and globalises Newton.
                // ----------------------------------------------------- //
                label lineSearchIters = 0;
                EigVec Rnew(R.size());
                scalar coupledResidualNew = GREAT;

                if (jfnkLineSearch)
                {
                    // Snapshot of the current iterate; restored on each
                    // backtracking trial.
                    const EigVec xBackup = x;
                    const scalar Rold2 = gSquaredNorm(R);
                    scalar alpha = nonlinearRelaxation;

                    while (true)
                    {
                        // x = xBackup + alpha * delta
                        x = xBackup + alpha*delta;

                        // F(x_new); commits VmGuess/statesGuess/IionGuess
                        Rnew = residualFor(x, true);
                        const scalar Rnew2 = gSquaredNorm(Rnew);

                        // Armijo sufficient-decrease test
                        if (Rnew2 <= (1.0 - 2.0*jfnkArmijoC*alpha)*Rold2)
                        {
                            break;  // accept trial
                        }

                        ++lineSearchIters;
                        if (lineSearchIters >= jfnkLineSearchMaxIter)
                        {
                            break;  // budget exhausted
                        }

                        alpha *= 0.5;
                        if (alpha < jfnkLineSearchAlphaMin)
                        {
                            // alpha floored; do one final trial with the
                            // floor and accept whatever comes out.
                            alpha = jfnkLineSearchAlphaMin;
                            x = xBackup + alpha*delta;
                            Rnew = residualFor(x, true);
                            break;
                        }
                    }
                    coupledResidualNew = relativeL2Norm(Rnew, x);
                }
                else
                {
                    // Legacy path: unconditional damped step.
                    x += nonlinearRelaxation*delta;
                    Rnew = residualFor(x, true);
                    coupledResidualNew = relativeL2Norm(Rnew, x);
                }

                const scalar VmResidual =
                    relativeL2Difference(VmGuess.primitiveField(), VmPrevious);
                const scalar IionResidual =
                    relativeL2Difference(IionGuess.primitiveField(), IionPrevious);
                const scalar maxStateResidual = maxStateRelativeL2Difference
                (
                    statesGuess,
                    statesPrevious,
                    stateResiduals
                );

                const bool converged =
                    corr + 1 >= implicitMinNonlinearIterations
                 && coupledResidualNew <= nonlinearTolerance
                 && VmResidual <= nonlinearVmTolerance
                 && IionResidual <= nonlinearIionTolerance
                 && (
                        !nonlinearRequireStatesConvergence
                     || maxStateResidual <= nonlinearStatesTolerance
                    );

                writeAndPrintResiduals
                (
                    corr + 1,
                    gmresIterations,
                    gmresError,
                    coupledResidualNew,
                    VmResidual,
                    IionResidual,
                    maxStateResidual,
                    stateResiduals,
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
        else
        {
            for
            (
                label corr = 0;
                corr < max(implicitNonlinearIterations, label(1));
                ++corr
            )
            {
                const scalarField VmPrevious(VmGuess.primitiveField());
                const scalarField IionPrevious(IionGuess.primitiveField());
                Field<Field<scalar>> statesPrevious(statesGuess);

                evaluateNonlinearFields(VmGuess, statesGuess, IionGuess, true);
                EigVec rhs = rhsFromIion(IionGuess);
                SpMat ACurrent = AImplicit;

                if (useDiagonalIion)
                {
                    computeSourceDerivative(VmGuess, IionGuess);
                    const EigVec sourceDerivativeVec = fieldToEigVec(sourceDerivative);
                    const EigVec VmLinearisationPoint = fieldToEigVec(VmGuess);

                    addDiagonalToMatrix
                    (
                        sourceDerivative.primitiveField(),
                       -theta,
                        ACurrent
                    );

                    rhs -= theta
                       * sourceDerivativeVec.cwiseProduct(VmLinearisationPoint);
                }

                label linearIterations = 0;
                scalar linearError = GREAT;

                // Phase B: decide how the persistent solver treats the matrix.
                //   - first solve, or dt changed  -> Rebuild (full setup)
                //   - diagonalIion (diagonal moves every iteration) -> ValuesOnly
                //   - Picard with a constant operator -> Reuse (solve only)
                MatrixUpdate updatePolicy;
                if (linearSolverNeedsRebuild)
                {
                    updatePolicy = MatrixUpdate::Rebuild;
                    linearSolverNeedsRebuild = false;
                }
                else if (useDiagonalIion)
                {
                    updatePolicy = MatrixUpdate::ValuesOnly;
                }
                else
                {
                    updatePolicy = MatrixUpdate::Reuse;
                }

                const EigVec Vsol = solveSparseSystem
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
                    linearError,
                    persistentPetscSolver,
                    persistentEigenSolvers,
                    updatePolicy,
                    petscReusePreconditioner
                );

                const EigVec Vprev = fieldToEigVec(VmGuess);
                const EigVec Vrelaxed =
                    Vprev + nonlinearRelaxation*(Vsol - Vprev);

                eigVecToField(Vrelaxed, VmGuess);
                VmGuess.correctBoundaryConditions();
                evaluateNonlinearFields(VmGuess, statesGuess, IionGuess, true);

                const EigVec nonlinearResidualVec =
                    applyAImplicit(fieldToEigVec(VmGuess)) - rhsFromIion(IionGuess);
                const scalar coupledResidual =
                    relativeL2Norm(nonlinearResidualVec, rhsFromIion(IionGuess));
                const scalar VmResidual =
                    relativeL2Difference(VmGuess.primitiveField(), VmPrevious);
                const scalar IionResidual =
                    relativeL2Difference(IionGuess.primitiveField(), IionPrevious);
                const scalar maxStateResidual = maxStateRelativeL2Difference
                (
                    statesGuess,
                    statesPrevious,
                    stateResiduals
                );

                const bool converged =
                    corr + 1 >= implicitMinNonlinearIterations
                 && coupledResidual <= nonlinearTolerance
                 && VmResidual <= nonlinearVmTolerance
                 && IionResidual <= nonlinearIionTolerance
                 && (
                        !nonlinearRequireStatesConvergence
                     || maxStateResidual <= nonlinearStatesTolerance
                    );

                writeAndPrintResiduals
                (
                    corr + 1,
                    linearIterations,
                    linearError,
                    coupledResidual,
                    VmResidual,
                    IionResidual,
                    maxStateResidual,
                    stateResiduals,
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

        // ----------------------------------------------------------------- //
        // Commit or roll back the nonlinear iterate.
        //
        // If the Newton/Picard loop exhausted its iteration budget without
        // reaching the convergence tolerances, we have two unsafe options
        // available:
        //
        //   (a) accept the last iterate, which may be far from the implicit
        //       solution and pollute every subsequent time step;
        //   (b) discard the time step's progress and reset to t_n values.
        //
        // We choose (b) by default: print a clear warning and restore Vm,
        // Iion and the ionic states from their time-step-start snapshots.
        // The user can opt back into (a) with `nonlinearAcceptUnconverged`
        // in spatialIntegrationProperties.
        // ----------------------------------------------------------------- //
        if (!nonlinearConverged && !nonlinearAcceptUnconverged)
        {
            WarningInFunction
                << "Nonlinear solver did not converge at t = " << t0
                << " after " << nonlinearIters << " iterations."
                << " Rolling back Vm, Iion and ionic states to the values"
                << " at the beginning of the time step." << nl
                << "Set 'nonlinearAcceptUnconverged true;' in"
                << " spatialIntegrationProperties to accept unconverged"
                << " iterates instead." << endl;

            VmGuess.primitiveFieldRef() = VmOldValues;
            VmGuess.correctBoundaryConditions();
            IionGuess.primitiveFieldRef() = IionOld.primitiveField();
            IionGuess.correctBoundaryConditions();
            statesGuess = statesOld;

            if (tnnpModel)
            {
                tnnpModel->resetStatesToStatesOld(statesGuess);
            }
        }
        else if (!nonlinearConverged && nonlinearAcceptUnconverged)
        {
            WarningInFunction
                << "Nonlinear solver did not converge at t = " << t0
                << " after " << nonlinearIters << " iterations."
                << " Accepting the unconverged iterate"
                << " (nonlinearAcceptUnconverged=true)." << endl;
        }

        // Commit the (possibly rolled-back) iterate into the persistent
        // fields used by subsequent time steps and post-processing.
        Vm.primitiveFieldRef() = VmGuess.primitiveField();
        Vm.correctBoundaryConditions();
        Iion.primitiveFieldRef() = IionGuess.primitiveField();
        Iion.correctBoundaryConditions();
        states = statesGuess;

        if (outFields.size())
        {
            if (useHighOrder_Iion && !statesAtCells)
            {
                ionicModel->exportStatesIntegrationPoints
                (
                    states,
                    outFields,
                    LREInterp_IionPtr().cellQuadWeightPhysical()
                );
            }
            else
            {
                // cellCentredReconstruct (and the low-order path) keep the
                // authoritative states at cell centres.
                ionicModel->exportStates(states, outFields);
            }
        }

        if (useHighOrder_Vm && dim > 1)
        {
            computeHighOrderLaplacian
            (
                Vm,
                conductivity,
                LREInterp_VmPtr(),
                fluxVm_HO,
                lapVm
            );
        }
        else
        {
            lapVm = fvc::laplacian(conductivity, Vm);
        }
        rhsVm = lapVm/(chi*Cm) - Iion + externalStimulusCurrent/(chi*Cm);

        updateActivationTimes
        (
            VmOldValues,
            Vm,
            activationTime,
            calculateActivationTime,
            t0,
            dt,
            activationThreshold
        );

        // Refresh the live activation samples whenever new points crossed the
        // threshold this step (files are overwritten; quiet to avoid log spam).
        label nActivated = 0;
        forAll(calculateActivationTime, cellI)
        {
            if (!calculateActivationTime[cellI])
            {
                ++nActivated;
            }
        }
        reduce(nActivated, sumOp<label>());

        if (nActivated > nActivatedPrev)
        {
            nActivatedPrev = nActivated;
            tLastActivationChange = runTime.value();
            writeActivationSamples
            (
                runTime,
                activationTime,
                diagonalStart,
                diagonalEnd,
                p8Point,
                nDiagonalSamples,
                sampleNeighbours,
                samplePower,
                false
            );
        }

        // Early-abort on activation failure (recorded in activationStatus.dat).
        if (abortOnActivationFailure && !stopPointTriggered)
        {
            const scalar tNow = runTime.value();
            if (nActivated == 0)
            {
                if (tNow > activationIgnitionTimeout)
                {
                    activationStatus = "FAILED_NO_IGNITION";
                    activationAborted = true;
                }
            }
            else if
            (
                stopAfterPointActivation
             && (tNow - tLastActivationChange) > activationStallTime
            )
            {
                activationStatus = "FAILED_CONDUCTION_BLOCK";
                activationAborted = true;
            }

            if (activationAborted)
            {
                WarningInFunction
                    << "Activation failure detected (" << activationStatus
                    << ") at t = " << tNow << " s, nActivated = " << nActivated
                    << "; aborting the run early. See activationStatus.dat."
                    << nl << endl;
                break;
            }
        }

        // Modified for cardiacFoam: only ONE rank owns the stop cell now, so
        // the test can no longer be gated on stopActivationCellI >= 0 - the
        // other ranks would never trigger. The owner's value is broadcast with
        // a max-reduction (every other rank contributes -GREAT), which makes
        // stopPointTriggered, stopPointActivationTime and effectiveEndTime
        // uniform BY CONSTRUCTION, and with them the loop's exit step,
        // FAILED_CONDUCTION_BLOCK and the final activationStatus.
        if (stopAfterPointActivation && !stopPointTriggered)
        {
            scalar actTime =
                (stopActivationCellI >= 0)
              ? activationTime[stopActivationCellI]
              : -GREAT;
            reduce(actTime, maxOp<scalar>());

            if (actTime > SMALL)
            {
                stopPointTriggered = true;
                stopPointActivationTime = actTime;
                effectiveEndTime = min
                (
                    effectiveEndTime,
                    stopPointActivationTime + stopDelayAfterActivation
                );

                Info<< "Stop activation point crossed threshold at t = "
                    << stopPointActivationTime << " s; solver will stop at t = "
                    << effectiveEndTime << " s" << nl << endl;
            }
        }

        activationVelocity = fvc::grad
        (
            1.0/(activationTime + dimensionedScalar("SMALL", dimTime, SMALL))
        );

        ++nSteps;
        ++runTime;

        if (nSteps <= 5 || nSteps % 100 == 0)
        {
            Info<< "Step " << nSteps
                << " time = " << runTime.value()
                << " min(Vm) = " << gMin(Vm.primitiveField())
                << " max(Vm) = " << gMax(Vm.primitiveField()) << nl;
        }

        runTime.write();
    }

    activationVelocity = fvc::grad
    (
        1.0/(activationTime + dimensionedScalar("SMALL", dimTime, SMALL))
    );
    activationTime.write();
    activationVelocity.write();
    Vm.write();
    Iion.write();
    externalStimulusCurrent.write();
    lapVm.write();
    rhsVm.write();
    forAll(outFields, i)
    {
        outFields[i].write();
    }

    if (stopAfterPointActivation && !stopPointTriggered)
    {
        Info<< "Stop activation point did not cross the threshold before t = "
            << runTime.value() << " s" << nl << endl;
    }

    writeActivationSamples
    (
        runTime,
        activationTime,
        diagonalStart,
        diagonalEnd,
        p8Point,
        nDiagonalSamples,
        sampleNeighbours,
        samplePower
    );

    // Finalise and record the activation status so a driver can tell a healthy
    // run from a non-igniting / conduction-blocked one without parsing the log.
    if (!activationAborted)
    {
        if (stopPointTriggered)
        {
            activationStatus = "OK";
        }
        else if (nActivatedPrev == 0)
        {
            activationStatus = "FAILED_NO_IGNITION";
        }
        else
        {
            activationStatus = "INCOMPLETE";
        }
    }
    {
        // Modified for cardiacFoam: see postProcessingDir().
        const fileName statusDir(postProcessingDir(runTime));
        mkDir(statusDir);
        OFstream statusFile(statusDir/"activationStatus.dat");
        statusFile
            << "# status finalTime_s nActivated stopPointActivated" << nl
            << activationStatus << ' ' << runTime.value() << ' '
            << nActivatedPrev << ' ' << (stopPointTriggered ? 1 : 0) << nl;
        Info<< "Activation status: " << activationStatus
            << " (t = " << runTime.value() << " s, nActivated = "
            << nActivatedPrev << ")" << nl << endl;
    }

    const auto tEndTotal = std::chrono::steady_clock::now();
    const scalar setupWallTime =
        std::chrono::duration<scalar>(tEndSetup - tStartSetup).count();
    const scalar loopWallTime =
        std::chrono::duration<scalar>(tEndTotal - tStartLoop).count();
    const scalar totalWallTime =
        std::chrono::duration<scalar>(tEndTotal - tStartTotal).count();
    const scalar peakRSSMB = currentPeakRSSMB();

    const fileName performanceDir
    (
        postProcessingDir(runTime)   // Modified for cardiacFoam
    );
    mkDir(performanceDir);
    OFstream perfFile(performanceDir/"solverPerformance.dat");
    perfFile
        << "# solver finalTime_s nSteps setupWall_s loopWall_s totalWall_s peakRSS_MB "
        << "memoryOptimization memoryOptimizationEffective "
        << "linearSolverBackend petscLinearKspType petscLinearPcType "
        << "jfnkLinearSolverBackend jfnkPetscKspType jfnkPetscPcType"
        << nl;
    perfFile
        << "highOrderElectroActivationFoamImplicitPETSc "
        << runTime.value() << ' '
        << nSteps << ' '
        << setupWallTime << ' '
        << loopWallTime << ' '
        << totalWallTime << ' '
        << peakRSSMB << ' '
        << memoryOptimization << ' '
        << (memoryOptimizationEffective ? "true" : "false") << ' '
        << linearSolverBackend << ' '
        << petscLinearKspType << ' '
        << petscLinearPcType << ' '
        << jfnkLinearSolverBackend << ' '
        << jfnkPetscKspType << ' '
        << jfnkPetscPcType
        << nl;

    Info<< "Solver performance written to " << performanceDir
        << "/solverPerformance.dat" << nl
        << "Wall time: setup = " << setupWallTime
        << " s, loop = " << loopWallTime
        << " s, total = " << totalWallTime << " s" << nl
        << "Peak RSS = " << peakRSSMB << " MB" << nl << endl;

    runTime.printExecutionTime(Info);

    Info<< "End" << nl << endl;

    return 0;
}
