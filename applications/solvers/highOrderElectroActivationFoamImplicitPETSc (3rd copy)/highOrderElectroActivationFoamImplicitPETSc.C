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
#include "LRE.H"
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
    using SpMat = Eigen::SparseMatrix<scalar, Eigen::RowMajor>;
    using Triplet = Eigen::Triplet<scalar>;
    using EigVec = Eigen::Matrix<scalar, Eigen::Dynamic, 1>;

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

    bool usesPetscBackend(const word& backend)
    {
        const std::string b = lowerWord(backend);
        return b == "petsc" || b == "ksp";
    }

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

    bool isGmresType(const word& kspType)
    {
        return petscKspTypeName(kspType) == "gmres";
    }

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
            const PetscInt nCols = static_cast<PetscInt>(A.cols());

            std::vector<PetscInt> rowNnz(nRows, 0);
            for (PetscInt row = 0; row < nRows; ++row)
            {
                rowNnz[row] =
                    static_cast<PetscInt>(A.outerIndexPtr()[row + 1])
                  - static_cast<PetscInt>(A.outerIndexPtr()[row]);
            }

            checkPetscError
            (
                MatCreateSeqAIJ
                (
                    PETSC_COMM_SELF,
                    nRows,
                    nCols,
                    0,
                    rowNnz.empty() ? nullptr : rowNnz.data(),
                    &A_
                ),
                "MatCreateSeqAIJ"
            );

            checkPetscError
            (
                MatSetOption(A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE),
                "MatSetOption"
            );

            for (PetscInt row = 0; row < nRows; ++row)
            {
                for (SpMat::InnerIterator it(A, row); it; ++it)
                {
                    const PetscInt r = static_cast<PetscInt>(it.row());
                    const PetscInt c = static_cast<PetscInt>(it.col());
                    const PetscScalar v = static_cast<PetscScalar>(it.value());

                    checkPetscError
                    (
                        MatSetValue(A_, r, c, v, INSERT_VALUES),
                        "MatSetValue"
                    );
                }
            }

            checkPetscError
            (
                MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY),
                "MatAssemblyBegin"
            );
            checkPetscError
            (
                MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY),
                "MatAssemblyEnd"
            );

            checkPetscError(VecCreateSeq(PETSC_COMM_SELF, nRows, &b_), "VecCreateSeq(b)");
            checkPetscError(VecDuplicate(b_, &x_), "VecDuplicate(x)");

            checkPetscError(KSPCreate(PETSC_COMM_SELF, &ksp_), "KSPCreate");
            checkPetscError(KSPSetOperators(ksp_, A_, A_), "KSPSetOperators");

            const std::string kspName = petscKspTypeName(kspType);
            const std::string pcName = petscPcTypeName(pcType);
            checkPetscError(KSPSetType(ksp_, kspName.c_str()), "KSPSetType");

            PC pc = nullptr;
            checkPetscError(KSPGetPC(ksp_, &pc), "KSPGetPC");
            checkPetscError(PCSetType(pc, pcName.c_str()), "PCSetType");
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

        void updateValues(const SpMat& A)
        {
            if (!isInitialised())
            {
                FatalErrorInFunction
                    << "PetscKspMatrixSolver::updateValues requires reset()"
                    << exit(FatalError);
            }

            if (A.rows() != n_ || A.cols() != n_)
            {
                FatalErrorInFunction
                    << "PetscKspMatrixSolver::updateValues called with "
                    << "matrix of size " << A.rows() << "x" << A.cols()
                    << " but cached PETSc Mat is " << n_ << "x" << n_
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
                            static_cast<PetscInt>(it.row()),
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
            checkPetscError(KSPSetReusePreconditioner(ksp_, PETSC_FALSE), "KSPSetReusePreconditioner(updateValues)");
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
                static_cast<scalar>(residualNorm)/max(rhs.norm(), scalar(SMALL));

            return result;
        }
    };

    struct PetscShellMatVecContext
    {
        std::function<EigVec(const EigVec&)> matVec;
        label n;
    };

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

    struct PetscShellPCContext
    {
        std::function<EigVec(const EigVec&)> apply;
        label n;
    };

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
                MatCreateShell(PETSC_COMM_SELF, nP, nP, nP, nP, &matContext_, &A_),
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

            checkPetscError(VecCreateSeq(PETSC_COMM_SELF, nP, &b_), "VecCreateSeq(cached shell b)");
            checkPetscError(VecDuplicate(b_, &x_), "VecDuplicate(cached shell x)");

            checkPetscError(KSPCreate(PETSC_COMM_SELF, &ksp_), "KSPCreate(cached shell)");
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
                hasShellPC_ = true;
            }
            else
            {
                const std::string pcName = petscPcTypeName(pcType);
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
                static_cast<scalar>(residualNorm)/max(rhs.norm(), scalar(SMALL));

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


    class EmbeddedTNNPModel
    {
        dictionary dict_;
        label nPoints_;
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

        void recordCorrection()
        {
            if (protectionCorrections_ < maxProtectionLogEntries_)
            {
                ++protectionCorrections_;
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
        EmbeddedTNNPModel
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT
        )
        :
            dict_(dict),
            nPoints_(nIntegrationPoints),
            constants_(NUM_CONSTANTS, 0.0),
            algebraic_(nIntegrationPoints, Field<scalar>(NUM_ALGEBRAIC, 0.0)),
            rates_(nIntegrationPoints, Field<scalar>(NUM_STATES, 0.0)),
            oldStates_(nIntegrationPoints, Field<scalar>(NUM_STATES, 0.0)),
            stepMs_(nIntegrationPoints, max(1000.0*initialDeltaT, SMALL)),
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
            protectionCorrections_(0)
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
                OFstream os(protectionLogName_);
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

        void initialiseStates(Field<Field<scalar>>& states)
        {
            states.setSize(nPoints_);
            forAll(states, pointI)
            {
                states[pointI].setSize(NUM_STATES, 0.0);
                rates_[pointI].setSize(NUM_STATES, 0.0);
                algebraic_[pointI].setSize(NUM_ALGEBRAIC, 0.0);
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

            if (states.size() != Vm.size())
            {
                initialiseStates(states);
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
            Field<Field<scalar>>& states
        )
        {
            const scalar tStartMs = 1000.0*stepStartTime;
            const scalar dtMs = 1000.0*deltaT;
            if (Im.size() != Vm.size())
            {
                FatalErrorInFunction
                    << "Im.size() != Vm.size()" << abort(FatalError);
            }

            if (states.size() != Vm.size())
            {
                initialiseStates(states);
            }

            forAll(Vm, pointI)
            {
                advancePoint
                (
                    tStartMs,
                    dtMs,
                    Vm[pointI],
                    states[pointI],
                    algebraic_[pointI],
                    rates_[pointI],
                    stepMs_[pointI],
                    pointI
                );
                Im[pointI] = algebraic_[pointI][Iion_cm];
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

    void logMemoryCheckpoint(const word& stage)
    {
        Info<< "Memory checkpoint [" << stage << "]: peakRSS_MB = "
            << currentPeakRSSMB() << nl;
    }

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

    void computeHighOrderLaplacian
    (
        const volScalarField& Vm,
        const volTensorField& conductivity,
        const LRE& LREInterp_Vm,
        surfaceScalarField& fluxVm_HO,
        volScalarField& lapVm
    )
    {
        const fvMesh& mesh = Vm.mesh();

        autoPtr<List<List<vector>>> gradVmQuadPtr = LREInterp_Vm.gradScalarFaceQuad(Vm);
        List<List<vector>>& gradVmQuad = gradVmQuadPtr.ref();
        const CompactListList<scalar>& faceQuadW = LREInterp_Vm.faceQuadWeight();

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
                fluxVm_HO[faceI] += faceArea*(faceNormal & Dg)*faceQuadW[faceI][pI];
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
                    const vector Dg = conductivity[owner] & gradVmQuad[globalFaceI][pI];
                    patchFlux[faceI] += faceArea*(faceNormal & Dg)*faceQuadW[globalFaceI][pI];
                }
            }
        }

        lapVm = fvc::div(fluxVm_HO);
    }

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

    scalar reconstructFromTaylor
    (
        const scalar c,
        const vector& grad,
        const symmTensor* H,
        const LRE::symmTensor3Order* T3,
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
            val += (1.0/6.0)*LRE::cubicForm(*T3, d, twoD);
        }

        return val;
    }


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

    void eigVecToField(const EigVec& v, volScalarField& fld)
    {
        scalarField& f = fld.primitiveFieldRef();

        forAll(f, cellI)
        {
            f[cellI] = v[cellI];
        }
    }

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

    scalar relativeL2Norm(const EigVec& r, const EigVec& reference)
    {
        return r.norm()/(reference.norm() + SMALL);
    }

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

    void addDiagonalToMatrix
    (
        const scalarField& diagonal,
        const scalar coefficient,
        SpMat& A
    )
    {
        forAll(diagonal, cellI)
        {
            A.coeffRef(cellI, cellI) += coefficient*diagonal[cellI];
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
        scalar& estimatedError
    )
    {
        const label n = b.size();
        const label m = min(max(krylovDim, label(1)), label(n));
        const label maxOuter = max(maxRestarts, label(0));

        // Initial iterate x_0 = 0  =>  initial residual r_0 = b - A*0 = b
        EigVec x = EigVec::Zero(n);
        EigVec r = b;
        const scalar bNorm = b.norm();
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

                // Arnoldi step: form w = A * V[j]
                EigVec w = matVec(V[j]);

                // Modified Gram-Schmidt orthogonalisation.
                //
                //   Each H(i,j) is computed against the *currently updated* w,
                //   so loss-of-orthogonality is bounded by O(eps * kappa(A_m))
                //   instead of O(eps * kappa(A_m)^2) as in classical GS.
                //   This is important for stiff ionic systems (TNNP) where
                //   the Jacobian can be poorly conditioned.
                for (label i = 0; i <= j; ++i)
                {
                    H(i, j) = V[i].dot(w);
                    w -= H(i, j)*V[i];
                }

                // Sub-diagonal entry of the Hessenberg matrix.
                // Happy breakdown (H(j+1,j) ~= 0) means the Krylov space is
                // already A-invariant and the current iterate is exact.
                H(j + 1, j) = w.norm();
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

                // Reconstruct candidate iterate: x_j = x_restart + V_j * y
                EigVec xj = x;
                for (label i = 0; i <= j; ++i)
                {
                    xj += y[i]*V[i];
                }

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
            beta = r.norm();
            estimatedError = beta/max(bNorm, SMALL);
            if (beta <= SMALL*bNorm) break;
        }

        return x;
    }

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
            triplets.emplace_back(cellI, cellI, coefficient);
        }

        M.resize(mesh.nCells(), mesh.nCells());
        M.setFromTriplets(triplets.begin(), triplets.end());
        M.makeCompressed();
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
        const LRE& LREInterp,
        const scalar coefficient,
        const bool compactRows,
        SpMat& M
    )
    {
        const bool twoD = mesh.nGeometricD() == 2;
        const vectorField& C = mesh.C();

        // Stencils and quadrature data from the LRE object.
        const labelListList& stencils = LREInterp.globalCellStencils();
        const CompactListList<point>& cellQP = LREInterp.cellQuadPoints();
        const CompactListList<scalar>& cellQW = LREInterp.cellQuadWeight();
        // Linear, quadratic and cubic Taylor-reconstruction coefficients.
        // For each stencil cell cI of cell cellI, gradCoeffs[cellI][cI] gives
        // the contribution of Vm_{stencil_cI} to grad Vm at cell cellI's
        // centre; similarly for the Hessian and 3rd-order derivatives.
        const CompactListList<vector>& gradCoeffs = LREInterp.QRGradCoeffs();
        const CompactListList<symmTensor>& hessCoeffs =
            LREInterp.cellHessianCoeffs();
        const CompactListList<LRE::symmTensor3Order>& thirdCoeffs =
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
            const labelList& curStencil = stencils[cellI];
            // LRE stores stencil-cell coefficients in [0..nStencil-1] and
            // the self (centre) coefficients at index nStencil.
            const label selfCoeffI = curStencil.size();

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
                           *LRE::cubicForm(thirdCoeffs[cellI][cI], d, twoD);
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
                       *LRE::cubicForm(thirdCoeffs[cellI][selfCoeffI], d, twoD);
                }

                addTripletIfNeeded(triplets, cellI, cellI, coefficient*w*selfCoeff);
            }
        }

        M.resize(mesh.nCells(), mesh.nCells());
        M.setFromTriplets(triplets.begin(), triplets.end());
        M.makeCompressed();
    }

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

    scalar normalDiffusivity(const tensor& D, const vector& n)
    {
        return max(mag(n & (D & n)), SMALL);
    }

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

    void addCellGradientDotCoeffs
    (
        std::vector<Triplet>& triplets,
        const label row,
        const scalar scale,
        const label cellI,
        const vector& d,
        const LRE& LREInterp
    )
    {
        if (mag(scale) <= SMALL)
        {
            return;
        }

        const labelListList& stencils = LREInterp.globalCellStencils();
        const CompactListList<vector>& gradCoeffs = LREInterp.QRGradCoeffs();
        const labelList& curStencil = stencils[cellI];
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
            cellI,
            scale*(gradCoeffs[cellI][selfCoeffI] & d)
        );
    }

    void addCellTaylorExtrapolationCoeffs
    (
        std::vector<Triplet>& triplets,
        const label row,
        const scalar scale,
        const label cellI,
        const point& evalPoint,
        const LRE& LREInterp,
        const vectorField& C,
        const bool twoD
    )
    {
        if (mag(scale) <= SMALL)
        {
            return;
        }

        const vector d = evalPoint - C[cellI];
        const labelListList& stencils = LREInterp.globalCellStencils();
        const CompactListList<vector>& gradCoeffs = LREInterp.QRGradCoeffs();
        const CompactListList<symmTensor>& hessCoeffs =
            LREInterp.cellHessianCoeffs();
        const CompactListList<LRE::symmTensor3Order>& thirdCoeffs =
            LREInterp.cellThirdDerivCoeffs();

        const labelList& curStencil = stencils[cellI];
        const label selfCoeffI = curStencil.size();

        forAll(curStencil, cI)
        {
            scalar coeff = gradCoeffs[cellI][cI] & d;

            if (LREInterp.order() >= 2)
            {
                coeff += 0.5*quadraticForm(hessCoeffs[cellI][cI], d);
            }

            if (LREInterp.order() >= 3)
            {
                coeff +=
                    (1.0/6.0)
                   *LRE::cubicForm(thirdCoeffs[cellI][cI], d, twoD);
            }

            addTripletIfNeeded(triplets, row, curStencil[cI], scale*coeff);
        }

        scalar selfCoeff = 1.0 + (gradCoeffs[cellI][selfCoeffI] & d);

        if (LREInterp.order() >= 2)
        {
            selfCoeff += 0.5*quadraticForm(hessCoeffs[cellI][selfCoeffI], d);
        }

        if (LREInterp.order() >= 3)
        {
            selfCoeff +=
                (1.0/6.0)
               *LRE::cubicForm(thirdCoeffs[cellI][selfCoeffI], d, twoD);
        }

        addTripletIfNeeded(triplets, row, cellI, scale*selfCoeff);
    }

    void addRhieChowInternalJumpCoeffs
    (
        std::vector<Triplet>& triplets,
        const label row,
        const scalar scale,
        const label own,
        const label nei,
        const vector& dPN,
        const LRE& LREInterp
    )
    {
        addTripletIfNeeded(triplets, row, nei, scale);
        addTripletIfNeeded(triplets, row, own, -scale);

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

    void addRhieChowBoundaryJumpCoeffs
    (
        std::vector<Triplet>& triplets,
        const label row,
        const scalar scale,
        const label own,
        const vector& dPb,
        const LRE& LREInterp
    )
    {
        addTripletIfNeeded(triplets, row, own, -scale);
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

    void assembleStandardOrthogonalStiffnessMatrix
    (
        const fvMesh& mesh,
        const volTensorField& conductivity,
        const LRE* LREInterp,
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

            addTripletIfNeeded(triplets, own, own, -a/max(V[own], SMALL));
            addTripletIfNeeded(triplets, own, nei,  a/max(V[own], SMALL));
            addTripletIfNeeded(triplets, nei, own,  a/max(V[nei], SMALL));
            addTripletIfNeeded(triplets, nei, nei, -a/max(V[nei], SMALL));
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

                    addTripletIfNeeded(triplets, own, own, -a/max(V[own], SMALL));

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

        K.resize(mesh.nCells(), mesh.nCells());
        K.setFromTriplets(triplets.begin(), triplets.end());
        K.makeCompressed();
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
        const LRE& LREInterp,
        const scalar stabilisationAlpha,
        const label tripletsPerFaceReserve,
        SpMat& K
    )
    {
        const bool twoD = mesh.nGeometricD() == 2;
        const vectorField& C = mesh.C();
        const surfaceVectorField& Cf = mesh.Cf();
        const CompactListList<point>& faceQP = LREInterp.faceQuadPoints();
        const CompactListList<scalar>& faceQW = LREInterp.faceQuadWeight();
        const labelListList& faceStencils = LREInterp.globalFaceStencils();
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

        SpMat Klocal(nCells, nCells);
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
                Klocal.resize(nCells, nCells);
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
                const labelList& curStencil = faceStencils[faceI];

                forAll(faceQP[faceI], qpI)
                {
                    const scalar w = faceQW[faceI][qpI];

                    forAll(curStencil, cI)
                    {
                        const label col = curStencil[cI];
                        const vector gCoeff = faceGradCoeffs[faceI][qpI][cI];
                        const scalar fluxCoeff =
                            area*w*(n & (conductivity[own] & gCoeff));

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

            forAll(patch, faceI)
            {
                const label globalFaceI = patch.start() + faceI;
                const label own = owner[globalFaceI];
                const vector Sf = mesh.Sf().boundaryField()[patchI][faceI];
                const scalar area = mag(Sf) + VSMALL;
                const vector n = Sf/area;
                const labelList& curStencil = faceStencils[globalFaceI];

                forAll(faceQP[globalFaceI], qpI)
                {
                    const scalar w = faceQW[globalFaceI][qpI];

                    forAll(curStencil, cI)
                    {
                        const label col = curStencil[cI];
                        const vector gCoeff =
                            faceGradCoeffs[globalFaceI][qpI][cI];
                        const scalar fluxCoeff =
                            area*w*(n & (conductivity[own] & gCoeff));

                        addTripletIfNeeded(triplets, own, col, fluxCoeff/max(V[own], SMALL));
                    }
                }

                if
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
            K.resize(nCells, nCells);
            K.setZero();
        }

        K.makeCompressed();
        std::vector<Triplet>().swap(triplets);
    }

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
        const LRE* LREInterp_Vm,
        const LRE& LREInterp_Iion,
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

            autoPtr<List<LRE::symmTensor3Order>> thirdVmPtr;
            const List<LRE::symmTensor3Order>* thirdVm = nullptr;
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
                const LRE::symmTensor3Order* T3 =
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
        const LRE& LREInterp_states,
        const LRE& LREInterp_Iion,
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

            tmp<volVectorField> tGrad = LREInterp_states.grad(scratch);
            const vectorField& grad = tGrad->internalField();

            tmp<volSymmTensorField> tHess;
            const symmTensorField* hess = nullptr;
            if (LREInterp_states.order() >= 2)
            {
                tHess = LREInterp_states.hessian(scratch);
                hess = &(tHess->internalField());
            }

            autoPtr<List<LRE::symmTensor3Order>> thirdPtr;
            const List<LRE::symmTensor3Order>* third = nullptr;
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
                const LRE::symmTensor3Order* T3 =
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
        const LRE& LREInterp_Iion,
        volScalarField& Iion
    )
    {
        const fvMesh& mesh = Iion.mesh();
        const CompactListList<scalar>& cellIionQuadW =
            LREInterp_Iion.cellQuadWeight();

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
        const LRE& LREInterp_Iion,
        Field<Field<scalar>>& statesCells,
        const boolList* cellMask = nullptr
    )
    {
        const CompactListList<scalar>& cellIionQuadW =
            LREInterp_Iion.cellQuadWeight();

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
        const LRE& LREInterp_Iion,
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
            frontCell = grown;
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

        if (bestCellI < 0)
        {
            FatalErrorInFunction
                << "Could not find a cell for point " << samplePoint
                << exit(FatalError);
        }

        return bestCellI;
    }


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

        if (C.size() == 0)
        {
            return 0.0;
        }

        const label nNbrs = min(max(requestedNeighbours, label(1)), C.size());
        const scalar power = max(requestedPower, scalar(SMALL));

        List<scalar> bestD2(nNbrs, GREAT);
        List<label> bestId(nNbrs, -1);

        forAll(C, cellI)
        {
            const scalar d2 = magSqr(C[cellI] - samplePoint);
            if (d2 < VSMALL)
            {
                return values[cellI];
            }

            for (label i = 0; i < nNbrs; ++i)
            {
                if (d2 < bestD2[i])
                {
                    for (label j = nNbrs - 1; j > i; --j)
                    {
                        bestD2[j] = bestD2[j - 1];
                        bestId[j] = bestId[j - 1];
                    }

                    bestD2[i] = d2;
                    bestId[i] = cellI;
                    break;
                }
            }
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
            vSum += w*values[bestId[i]];
        }

        return vSum/max(wSum, scalar(SMALL));
    }

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
        const fileName sampleDir
        (
            runTime.path()/"postProcessing"/"highOrderElectroActivationFoamImplicitPETSc"
        );
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
        const fileName residualDir
        (
            runTime.path()/"postProcessing"/"highOrderElectroActivationFoamImplicitPETSc"
        );
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
        // ----------------------------------------------------------------- //
        SpMat AImplicit = M;
        AImplicit *= (1.0/dt);
        AImplicit -= theta*L;

        SpMat BImplicit = M;
        BImplicit *= (1.0/dt);
        if (theta < 1.0 - SMALL)
        {
            // Crank-Nicolson (or any theta < 1) needs the explicit half of
            // the diffusion operator on the right-hand side. For pure
            // Backward Euler (theta=1) this term vanishes and is skipped.
            BImplicit += (1.0 - theta)*L;
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
                        ionicModel->solveODE(t0, dt, VmF, ImF, sF);
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
                        ionicModel->solveODE(t0, dt, VmS, ImS, sS);
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
            EigVec rhs = BImplicit*Vn;
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
                    return AImplicit * fieldToEigVec(VmCandidate) - rhsFromIion(IionCandidate);
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
                const scalar epsScale = jfnkEpsilon*(1.0 + x.norm());
                auto matVec = [&](const EigVec& v) -> EigVec
                {
                    const scalar eps = epsScale/max(v.norm(), SMALL);
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
                            false
                        );
                    }

                    delta = jfnkPetscShellSolver.solve
                    (
                        std::function<EigVec(const EigVec&)>(matVec),
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
                    const scalar Rold2 = R.squaredNorm();
                    scalar alpha = nonlinearRelaxation;

                    while (true)
                    {
                        // x = xBackup + alpha * delta
                        x = xBackup + alpha*delta;

                        // F(x_new); commits VmGuess/statesGuess/IionGuess
                        Rnew = residualFor(x, true);
                        const scalar Rnew2 = Rnew.squaredNorm();

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
                    linearError
                );

                const EigVec Vprev = fieldToEigVec(VmGuess);
                const EigVec Vrelaxed =
                    Vprev + nonlinearRelaxation*(Vsol - Vprev);

                eigVecToField(Vrelaxed, VmGuess);
                VmGuess.correctBoundaryConditions();
                evaluateNonlinearFields(VmGuess, statesGuess, IionGuess, true);

                const EigVec nonlinearResidualVec =
                    AImplicit*fieldToEigVec(VmGuess) - rhsFromIion(IionGuess);
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
                    LREInterp_IionPtr().cellQuadWeight()
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
        {
            label nAct = 0;
            forAll(calculateActivationTime, cellI)
            {
                if (!calculateActivationTime[cellI])
                {
                    ++nAct;
                }
            }
            reduce(nAct, sumOp<label>());

            if (nAct > nActivatedPrev)
            {
                nActivatedPrev = nAct;
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
        }

        if
        (
            stopAfterPointActivation
         && !stopPointTriggered
         && stopActivationCellI >= 0
        )
        {
            const scalar actTime = activationTime[stopActivationCellI];
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
        runTime.path()/"postProcessing"/"highOrderElectroActivationFoamImplicitPETSc"
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
