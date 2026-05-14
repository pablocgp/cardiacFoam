/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

Solver
    highOrderElectroActivationFoamImplicitJfNK

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

// PETSc replaces Eigen entirely in this solver. We use Mat for sparse
// matrices, Vec for dense vectors, KSP for linear solves and (optionally)
// SNES for the Newton-Krylov nonlinear solve.
//
// IMPORTANT: PETSc headers MUST come before fvCFD.H / any OpenFOAM header.
// PETSc's petscmath.h defines macros (e.g. PetscPowReal) that expand to
// unqualified pow() calls; if fvCFD.H is included first it pulls in
// `using namespace Foam;` plus Foam::pow overloads, and the macro expansion
// then sees an ambiguous overload set (std::pow + Foam::pow) and refuses
// to compile.
#include <petscksp.h>
#include <petscsnes.h>
#include <petscmat.h>
#include <petscvec.h>

#include "fvCFD.H"
#include "ionicModel.H"
#include "TNNP.H"
#include "LRE.H"
#include "Field.H"
#include "volFields.H"
#include <cmath>
#include <vector>

// Convenience macro: abort on PETSc errors. PETSc functions return non-zero
// on failure; chaining checks via CHKERRABORT keeps the call sites compact.
#define PCHECK(ierr) CHKERRABORT(PETSC_COMM_SELF, (ierr))

namespace
{
    // RAII wrappers that own a PETSc Mat / Vec and call MatDestroy / VecDestroy
    // on destruction. This eliminates manual cleanup, especially in the
    // presence of FatalError exits, and avoids leaking PETSc objects.
    //
    // The underlying PETSc handle (Mat / Vec) is itself a pointer type, so
    // these wrappers store it by value and expose it via raw() and operator&.
    class PetscMatRAII
    {
        Mat m_;
    public:
        PetscMatRAII() : m_(nullptr) {}
        ~PetscMatRAII() { if (m_) MatDestroy(&m_); }
        PetscMatRAII(const PetscMatRAII&) = delete;
        PetscMatRAII& operator=(const PetscMatRAII&) = delete;
        Mat raw() const { return m_; }
        // Reference suitable for assembly functions that need Mat&.
        // Destroys any pre-existing handle so each call starts clean.
        Mat& slot() { reset(); return m_; }
        operator Mat() const { return m_; }
        void reset()
        {
            if (m_) { MatDestroy(&m_); m_ = nullptr; }
        }
        // Move-take ownership of a raw Mat.
        void adopt(Mat m) { reset(); m_ = m; }
        Mat release() { Mat t = m_; m_ = nullptr; return t; }
    };

    class PetscVecRAII
    {
        Vec v_;
    public:
        PetscVecRAII() : v_(nullptr) {}
        ~PetscVecRAII() { if (v_) VecDestroy(&v_); }
        PetscVecRAII(const PetscVecRAII&) = delete;
        PetscVecRAII& operator=(const PetscVecRAII&) = delete;
        Vec raw() const { return v_; }
        Vec& slot() { reset(); return v_; }
        operator Vec() const { return v_; }
        void reset()
        {
            if (v_) { VecDestroy(&v_); v_ = nullptr; }
        }
        void adopt(Vec v) { reset(); v_ = v; }
        Vec release() { Vec t = v_; v_ = nullptr; return t; }
    };

    // Create a sequential dense PETSc vector of size n filled with zeros.
    Vec createPetscVec(const label n)
    {
        Vec v;
        PetscErrorCode ierr;
        ierr = VecCreateSeq(PETSC_COMM_SELF, n, &v); PCHECK(ierr);
        ierr = VecSet(v, 0.0); PCHECK(ierr);
        return v;
    }

    // Create a sequential AIJ (compressed sparse row) matrix of size n x n
    // with nzPerRow non-zeros per row reserved. PETSc grows dynamically if
    // exceeded but it is far cheaper to size correctly up front.
    Mat createPetscMat(const label n, const label nzPerRow)
    {
        Mat A;
        PetscErrorCode ierr;
        ierr = MatCreate(PETSC_COMM_SELF, &A); PCHECK(ierr);
        ierr = MatSetSizes(A, n, n, n, n); PCHECK(ierr);
        ierr = MatSetType(A, MATSEQAIJ); PCHECK(ierr);
        ierr = MatSeqAIJSetPreallocation(A, nzPerRow, nullptr); PCHECK(ierr);
        ierr = MatSetOption
        (
            A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE
        ); PCHECK(ierr);
        return A;
    }

    // Copy an OpenFOAM scalarField into a PETSc Vec (sequential, same size).
    void fieldToPetscVec(const scalarField& f, Vec v)
    {
        PetscErrorCode ierr;
        PetscScalar* vArr = nullptr;
        ierr = VecGetArray(v, &vArr); PCHECK(ierr);
        forAll(f, cellI)
        {
            vArr[cellI] = f[cellI];
        }
        ierr = VecRestoreArray(v, &vArr); PCHECK(ierr);
    }

    void fieldToPetscVec(const volScalarField& fld, Vec v)
    {
        fieldToPetscVec(fld.primitiveField(), v);
    }

    // Copy a PETSc Vec back into an OpenFOAM volScalarField (internal).
    void petscVecToField(const Vec v, volScalarField& fld)
    {
        PetscErrorCode ierr;
        const PetscScalar* vArr = nullptr;
        ierr = VecGetArrayRead(v, &vArr); PCHECK(ierr);
        scalarField& f = fld.primitiveFieldRef();
        forAll(f, cellI)
        {
            f[cellI] = vArr[cellI];
        }
        ierr = VecRestoreArrayRead(v, &vArr); PCHECK(ierr);
    }

    // Convenience: copy a Vec into a freshly created Vec of the same size.
    void copyVec(Vec src, Vec dst)
    {
        PetscErrorCode ierr;
        ierr = VecCopy(src, dst); PCHECK(ierr);
    }

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

    // PETSc-native L2 norm of (current - previous) over an OpenFOAM scalarField,
    // relative to ||current||. Computed on the OpenFOAM side directly: avoids
    // creating temporary Vecs for what is just two passes over a flat array.
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

    // PETSc-Vec analogue of  ||r|| / ||reference||
    // VecNorm reduces in parallel for MPI runs; in serial it is a fast loop.
    scalar relativeL2Norm(Vec r, Vec reference)
    {
        PetscReal rn = 0.0, refn = 0.0;
        PetscErrorCode ierr;
        ierr = VecNorm(r, NORM_2, &rn); PCHECK(ierr);
        ierr = VecNorm(reference, NORM_2, &refn); PCHECK(ierr);
        return scalar(rn)/(scalar(refn) + SMALL);
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

    // Add a diagonal contribution to a PETSc matrix:  A_ii += coeff * d_i
    // We accept that the matrix has already been assembled and we are
    // re-opening it for further modifications, hence the explicit assembly
    // call at the end.
    void addDiagonalToMatrix
    (
        const scalarField& diagonal,
        const scalar coefficient,
        Mat A
    )
    {
        PetscErrorCode ierr;
        forAll(diagonal, cellI)
        {
            const PetscInt idx = static_cast<PetscInt>(cellI);
            const PetscScalar v = coefficient*diagonal[cellI];
            ierr = MatSetValues(A, 1, &idx, 1, &idx, &v, ADD_VALUES);
            PCHECK(ierr);
        }
        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); PCHECK(ierr);
        ierr = MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY); PCHECK(ierr);
    }

    // ----------------------------------------------------------------------
    // The legacy hand-rolled GMRES(m) used Eigen sparse matrices and dense
    // Hessenberg QR. It has been replaced by PETSc's KSP framework, which
    // provides a tuned, parallel-ready GMRES (and FGMRES, BiCGSTAB, ...)
    // with proper restart, multiple orthogonalisation flavours and a wide
    // selection of preconditioners. See solveLinearKSP() below.
    // ----------------------------------------------------------------------
    // Linear solve  A * x = b  with PETSc KSP.
    //
    //   - solverType    : KSP type string (e.g. KSPGMRES, KSPBCGS, KSPFGMRES,
    //                     or KSPPREONLY when combined with a direct LU PC)
    //   - pcType        : preconditioner type (PCILU, PCJACOBI, PCBJACOBI,
    //                     PCASM, PCLU, ...). For SparseLU-equivalent: PCLU.
    //   - krylovDim     : restart parameter for GMRES
    //   - tolerance     : relative residual tolerance ||r||/||b|| < tol
    //   - maxIterations : KSP iteration budget
    //
    // Inputs A, b are unchanged; x is overwritten with the solution.
    void solveLinearKSP
    (
        Mat A,
        Vec b,
        Vec x,
        const word& solverType,
        const word& pcType,
        const label krylovDim,
        const label maxRestarts,
        const scalar tolerance,
        const label maxIterations,
        label& iterations,
        scalar& estimatedError
    )
    {
        PetscErrorCode ierr;

        KSP ksp;
        ierr = KSPCreate(PETSC_COMM_SELF, &ksp); PCHECK(ierr);
        ierr = KSPSetOperators(ksp, A, A); PCHECK(ierr);

        // Translate the user-facing solver string into a KSP type token.
        const char* kspType = KSPGMRES;
        if      (solverType == "GMRES")    kspType = KSPGMRES;
        else if (solverType == "FGMRES")   kspType = KSPFGMRES;
        else if (solverType == "BiCGSTAB") kspType = KSPBCGS;
        else if (solverType == "BCGS")     kspType = KSPBCGS;
        else if (solverType == "PREONLY")  kspType = KSPPREONLY; // direct via PC
        else if (solverType == "SparseLU") kspType = KSPPREONLY; // alias
        ierr = KSPSetType(ksp, kspType); PCHECK(ierr);

        // For GMRES family, configure the restart cycle (= max Krylov dim
        // before restart). The total iteration budget is the product of
        // (krylovDim) * (maxRestarts + 1), capped by maxIterations.
        if
        (
            std::strcmp(kspType, KSPGMRES) == 0
         || std::strcmp(kspType, KSPFGMRES) == 0
        )
        {
            ierr = KSPGMRESSetRestart(ksp, krylovDim); PCHECK(ierr);
            // Modified Gram-Schmidt: more stable than classical GS for
            // ill-conditioned Jacobians (e.g. stiff ionic systems).
            ierr = KSPGMRESSetOrthogonalization
            (
                ksp, KSPGMRESModifiedGramSchmidtOrthogonalization
            ); PCHECK(ierr);
        }

        PC pc;
        ierr = KSPGetPC(ksp, &pc); PCHECK(ierr);
        const char* pType = PCILU;
        if      (pcType == "ILU"     || pcType == "ilu")     pType = PCILU;
        else if (pcType == "JACOBI"  || pcType == "jacobi")  pType = PCJACOBI;
        else if (pcType == "BJACOBI" || pcType == "bjacobi") pType = PCBJACOBI;
        else if (pcType == "LU"      || pcType == "lu"
             ||  pcType == "SparseLU")                       pType = PCLU;
        else if (pcType == "NONE"    || pcType == "none")    pType = PCNONE;
        else if (pcType == "SOR"     || pcType == "sor")     pType = PCSOR;
        else if (pcType == "ASM"     || pcType == "asm")     pType = PCASM;
        ierr = PCSetType(pc, pType); PCHECK(ierr);

        const PetscInt totalIters =
            static_cast<PetscInt>
            (
                max(maxIterations, krylovDim*(maxRestarts + 1))
            );
        ierr = KSPSetTolerances
        (
            ksp,
            tolerance,            // relative
            PETSC_DEFAULT,        // absolute
            PETSC_DEFAULT,        // divergence
            totalIters
        ); PCHECK(ierr);

        // Allow command-line overrides (-ksp_view, -ksp_monitor, ...).
        ierr = KSPSetFromOptions(ksp); PCHECK(ierr);

        ierr = KSPSolve(ksp, b, x); PCHECK(ierr);

        PetscInt its = 0;
        PetscReal rnorm = 0.0;
        ierr = KSPGetIterationNumber(ksp, &its); PCHECK(ierr);
        ierr = KSPGetResidualNorm(ksp, &rnorm); PCHECK(ierr);
        iterations = static_cast<label>(its);
        estimatedError = scalar(rnorm);

        ierr = KSPDestroy(&ksp); PCHECK(ierr);
    }

    // Buffer accumulator that batches MatSetValues calls per row.
    //
    //   PETSc's MatSetValues with ADD_VALUES is efficient when many entries
    //   in the same row are inserted at once, so we collect (col, value)
    //   pairs into a small per-row buffer and flush at row boundaries.
    //   Drops entries whose magnitude is below SMALL to keep the sparsity
    //   pattern tight (mirrors the addTripletIfNeeded behaviour we used
    //   with Eigen).
    struct MatRowAccumulator
    {
        Mat A_;
        PetscInt row_;
        std::vector<PetscInt> cols_;
        std::vector<PetscScalar> vals_;

        MatRowAccumulator(Mat A, const label row)
        :
            A_(A),
            row_(static_cast<PetscInt>(row))
        {
            cols_.reserve(64);
            vals_.reserve(64);
        }

        void add(const label col, const scalar value)
        {
            if (mag(value) > SMALL)
            {
                cols_.push_back(static_cast<PetscInt>(col));
                vals_.push_back(static_cast<PetscScalar>(value));
            }
        }

        void flush()
        {
            if (cols_.empty()) return;
            PetscErrorCode ierr = MatSetValues
            (
                A_,
                1, &row_,
                static_cast<PetscInt>(cols_.size()), cols_.data(),
                vals_.data(),
                ADD_VALUES
            );
            PCHECK(ierr);
            cols_.clear();
            vals_.clear();
        }
    };

    // Diagonal (lumped) mass matrix with constant coefficient on every cell.
    // Pre-allocated with exactly 1 non-zero per row.
    void assembleDiagonalMassMatrix
    (
        const fvMesh& mesh,
        const scalar coefficient,
        Mat& M
    )
    {
        const label n = mesh.nCells();
        if (M) MatDestroy(&M);
        M = createPetscMat(n, 1);

        PetscErrorCode ierr;
        forAll(mesh.C(), cellI)
        {
            const PetscInt idx = static_cast<PetscInt>(cellI);
            const PetscScalar v = coefficient;
            ierr = MatSetValues(M, 1, &idx, 1, &idx, &v, INSERT_VALUES);
            PCHECK(ierr);
        }
        ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); PCHECK(ierr);
        ierr = MatAssemblyEnd  (M, MAT_FINAL_ASSEMBLY); PCHECK(ierr);
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
        Mat& M
    )
    {
        const bool twoD = mesh.nGeometricD() == 2;
        const vectorField& C = mesh.C();

        // Stencils and quadrature data from the LRE object.
        const labelListList& stencils = LREInterp.globalCellStencils();
        const CompactListList<point>& cellQP = LREInterp.cellQuadPoints();
        const CompactListList<scalar>& cellQW = LREInterp.cellQuadWeight();
        const CompactListList<vector>& gradCoeffs = LREInterp.QRGradCoeffs();
        const CompactListList<symmTensor>& hessCoeffs =
            LREInterp.cellHessianCoeffs();
        const CompactListList<LRE::symmTensor3Order>& thirdCoeffs =
            LREInterp.cellThirdDerivCoeffs();

        const label n = mesh.nCells();

        // Pre-allocate non-zeros per row using the actual stencil size.
        // Each row has the stencil size + 1 self entry.
        std::vector<PetscInt> nnzPerRow(n, 0);
        forAll(stencils, cellI)
        {
            nnzPerRow[cellI] = static_cast<PetscInt>(stencils[cellI].size() + 1);
        }

        if (M) MatDestroy(&M);
        PetscErrorCode ierr;
        ierr = MatCreate(PETSC_COMM_SELF, &M); PCHECK(ierr);
        ierr = MatSetSizes(M, n, n, n, n); PCHECK(ierr);
        ierr = MatSetType(M, MATSEQAIJ); PCHECK(ierr);
        ierr = MatSeqAIJSetPreallocation(M, 0, nnzPerRow.data()); PCHECK(ierr);
        ierr = MatSetOption
        (
            M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE
        ); PCHECK(ierr);

        forAll(stencils, cellI)
        {
            MatRowAccumulator row(M, cellI);

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

                    row.add(curStencil[cI], coefficient*w*coeff);
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

                row.add(cellI, coefficient*w*selfCoeff);
            }

            row.flush();
        }

        ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); PCHECK(ierr);
        ierr = MatAssemblyEnd  (M, MAT_FINAL_ASSEMBLY); PCHECK(ierr);
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

    // Standard 5/7-point orthogonal FV Laplacian. Each cell talks only to
    // its face neighbours, so we pre-allocate ~7 entries per row (works for
    // any 1D/2D/3D structured or unstructured mesh up to face count ~ 6).
    void assembleStandardOrthogonalStiffnessMatrix
    (
        const fvMesh& mesh,
        const volTensorField& conductivity,
        Mat& K
    )
    {
        const vectorField& C = mesh.C();
        const scalarField& V = mesh.V();
        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();
        const label n = mesh.nCells();

        if (K) MatDestroy(&K);
        K = createPetscMat(n, 7);

        // Direct MatSetValues per face: cheaper than a per-row accumulator
        // because each face touches only 4 entries (2 owner + 2 neighbour).
        PetscErrorCode ierr;
        auto add = [&](const label r, const label c, const scalar v) -> void
        {
            if (mag(v) <= SMALL) return;
            const PetscInt rr = static_cast<PetscInt>(r);
            const PetscInt cc = static_cast<PetscInt>(c);
            const PetscScalar vv = v;
            ierr = MatSetValues(K, 1, &rr, 1, &cc, &vv, ADD_VALUES);
            PCHECK(ierr);
        };

        forAll(neighbour, faceI)
        {
            const label own = owner[faceI];
            const label nei = neighbour[faceI];

            const tensor Df = 0.5*(conductivity[own] + conductivity[nei]);
            const scalar a =
                orthogonalDiffusionCoeff(mesh.Sf()[faceI], C[nei] - C[own], Df);

            add(own, own, -a/max(V[own], SMALL));
            add(own, nei,  a/max(V[own], SMALL));
            add(nei, own,  a/max(V[nei], SMALL));
            add(nei, nei, -a/max(V[nei], SMALL));
        }

        const surfaceVectorField& Cf = mesh.Cf();

        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            const word bcType =
                patch.lookupPatchField<volScalarField, scalar>("Vm").type();

            // Insulating (zeroGradient) and empty patches contribute zero
            // flux to the matrix; their boundary value enters only via the
            // RHS in the calling code, if at all.
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

                    add(own, own, -a/max(V[own], SMALL));
                }
            }
        }

        ierr = MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY); PCHECK(ierr);
        ierr = MatAssemblyEnd  (K, MAT_FINAL_ASSEMBLY); PCHECK(ierr);
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
        Mat& K
    )
    {
        const CompactListList<point>& faceQP = LREInterp.faceQuadPoints();
        const CompactListList<scalar>& faceQW = LREInterp.faceQuadWeight();
        const labelListList& faceStencils = LREInterp.globalFaceStencils();
        const List<CompactListList<vector>>& faceGradCoeffs =
            LREInterp.QRGradFaceGPCoeffs();

        const scalarField& V = mesh.V();
        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();
        const label n = mesh.nCells();

        // Heuristic preallocation for HO stencils: typical LRE stencils span
        // ~15-50 cells; we reserve a comfortable upper bound (60) per row.
        // PETSc grows dynamically beyond this, but with a one-time penalty.
        if (K) MatDestroy(&K);
        K = createPetscMat(n, 60);

        PetscErrorCode ierr;
        auto add = [&](const label r, const label c, const scalar v) -> void
        {
            if (mag(v) <= SMALL) return;
            const PetscInt rr = static_cast<PetscInt>(r);
            const PetscInt cc = static_cast<PetscInt>(c);
            const PetscScalar vv = v;
            ierr = MatSetValues(K, 1, &rr, 1, &cc, &vv, ADD_VALUES);
            PCHECK(ierr);
        };

        // -------- Internal faces -------------------------------------------
        forAll(neighbour, faceI)
        {
            const label own = owner[faceI];
            const label nei = neighbour[faceI];
            const vector Sf = mesh.Sf()[faceI];
            const scalar area = mag(Sf) + VSMALL;
            const vector nrm = Sf/area;          // unit normal owner -> nei
            const labelList& curStencil = faceStencils[faceI];

            forAll(faceQP[faceI], qpI)
            {
                const scalar w = faceQW[faceI][qpI];

                // For each stencil cell j we accumulate its contribution
                // to the flux through this face / cell volume.
                forAll(curStencil, cI)
                {
                    const label col = curStencil[cI];
                    const vector gCoeff = faceGradCoeffs[faceI][qpI][cI];
                    // Owner-cell conductivity is used at internal faces;
                    // exact for homogeneous tissue. For heterogeneous media
                    // a face-averaged (harmonic mean) conductivity would
                    // restore consistency at material interfaces.
                    const scalar fluxCoeff =
                        area*w*(nrm & (conductivity[own] & gCoeff));

                    // Divergence theorem: +flux to owner, -flux to neighbour.
                    add(own, col,  fluxCoeff/max(V[own], SMALL));
                    add(nei, col, -fluxCoeff/max(V[nei], SMALL));
                }
            }
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
                const vector nrm = Sf/area;
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
                            area*w*(nrm & (conductivity[own] & gCoeff));

                        add(own, col, fluxCoeff/max(V[own], SMALL));
                    }
                }
            }
        }

        ierr = MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY); PCHECK(ierr);
        ierr = MatAssemblyEnd  (K, MAT_FINAL_ASSEMBLY); PCHECK(ierr);
    }

    // (solveSparseSystem removed - replaced by solveLinearKSP further up.)

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
        const scalar samplePower
    )
    {
        const fileName sampleDir
        (
            runTime.path()/"postProcessing"/"highOrderElectroActivationFoamImplicitJfNK_PETSc"
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

        Info<< "P8 activation time = " << p8Time << " s ("
            << 1000.0*p8Time << " ms)" << nl
            << "Activation samples written to " << sampleDir << nl;
    }
}

int main(int argc, char* argv[])
{
    // ----------------------------------------------------------------- //
    // MPI / PETSc lifetime ordering (critical):
    //
    //   1. OpenFOAM's setRootCaseLists.H constructs argList, which calls
    //      MPI_Init (when OpenFOAM is built against MPI).
    //   2. PetscInitialize() then detects MPI is already initialised and
    //      *does not* take ownership of it.
    //   3. All PETSc-using objects live inside an inner { ... } scope so
    //      their destructors (Mat, Vec, RAII wrappers) run BEFORE
    //      PetscFinalize().
    //   4. PetscFinalize() tears down PETSc but leaves MPI alone, because
    //      it didn't initialise it.
    //   5. main returns; argList destructor finally calls MPI_Finalize.
    //
    // Reversing this order (PetscInitialize first) makes PETSc own MPI;
    // then PetscFinalize() calls MPI_Finalize() while OpenFOAM destructors
    // still need it -> "MPI_Comm_get_attr() after MPI_FINALIZE" abort at
    // exit, even though the solver itself ran to completion.
    // ----------------------------------------------------------------- //
    #include "setRootCaseLists.H"

    PetscInitialize(&argc, &argv, nullptr, nullptr);

    // Inner scope: every PETSc object (and every OpenFOAM object that
    // owns/uses PETSc handles via RAII) lives here. Closing brace runs all
    // their destructors before PetscFinalize().
    {
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    const label nStates = ionicModel->nEqns();
    TNNP* tnnpModel = dynamic_cast<TNNP*>(&ionicModel());

    Field<Field<scalar>> states
    (
        totalIionIntegrationPoints,
        Field<scalar>(nStates, 0.0)
    );

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

    Info<< "Running high-order electro activation implicit solver" << nl
        << "Dimension = " << dim << nl
        << "dx = " << dx << nl
        << "requested dt = " << dt << nl
        << "explicit stable dt reference = " << dtExplicitReference << nl
        << "implicitScheme = " << implicitScheme << nl
        << "massMatrix = " << massMatrixMode << nl
        << "useHighOrder_Vm = " << useHighOrder_Vm << nl
        << "useHighOrder_Iion = " << useHighOrder_Iion << nl
        << "Iion integration points = " << totalIionIntegrationPoints << nl
        << "nonlinearMethod = " << nonlinearMethod << nl
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

    // PETSc Mat handles for the global operators. RAII wrappers destroy them
    // automatically at scope exit, even on FatalError.
    PetscMatRAII Mowner;
    PetscMatRAII Kowner;
    PetscMatRAII Lowner;

    if (!useHighOrder_Vm || dim <= 1 || massMatrixMode == "lumped")
    {
        if (massMatrixMode == "consistent" && (!useHighOrder_Vm || dim <= 1))
        {
            Info<< "Consistent mass requested without high-order Vm; "
                << "using the diagonal/lumped mass matrix." << nl;
        }
        assembleDiagonalMassMatrix(mesh, 1.0, Mowner.slot());
    }
    else
    {
        assembleConsistentMassMatrixHO
        (
            mesh, LREInterp_VmPtr(), 1.0, Mowner.slot()
        );
    }

    if (useHighOrder_Vm && dim > 1)
    {
        assembleHighOrderStiffnessMatrix
        (
            mesh, conductivity, LREInterp_VmPtr(), Kowner.slot()
        );
    }
    else
    {
        assembleStandardOrthogonalStiffnessMatrix
        (
            mesh, conductivity, Kowner.slot()
        );
    }

    Mat M = Mowner.raw();
    Mat K = Kowner.raw();

    // L = K / (chi * Cm). MatDuplicate copies the sparsity AND values; then
    // MatScale applies the scalar factor in place on the new matrix.
    {
        Mat tmpL = nullptr;
        PetscErrorCode ierr =
            MatDuplicate(K, MAT_COPY_VALUES, &tmpL); PCHECK(ierr);
        ierr = MatScale(tmpL, 1.0/(chi.value()*Cm.value())); PCHECK(ierr);
        Lowner.adopt(tmpL);
    }
    Mat L = Lowner.raw();

    autoPtr<OFstream> nonlinearResidualFilePtr;
    if (writeNonlinearResiduals)
    {
        const fileName residualDir
        (
            runTime.path()/"postProcessing"/"highOrderElectroActivationFoamImplicitJfNK_PETSc"
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
        // Build the two theta-scheme operator matrices for this time step.
        // We duplicate M then apply AXPY with L to produce A and B.
        // PETSc requires SAME_NONZERO_PATTERN when M and L share sparsity
        // (the consistent HO case); SUBSET_NONZERO_PATTERN otherwise, since
        // diagonal-lumped M is a strict subset of the full-stencil L.
        PetscMatRAII AImplicitOwner;
        PetscMatRAII BImplicitOwner;
        {
            Mat tmpA = nullptr;
            Mat tmpB = nullptr;
            PetscErrorCode ierr;

            // tmpA = M/dt
            ierr = MatDuplicate(M, MAT_COPY_VALUES, &tmpA); PCHECK(ierr);
            ierr = MatScale(tmpA, 1.0/dt); PCHECK(ierr);
            // tmpA += -theta * L   →   tmpA = M/dt - theta*L
            ierr = MatAXPY
            (
                tmpA, -theta, L, DIFFERENT_NONZERO_PATTERN
            ); PCHECK(ierr);
            AImplicitOwner.adopt(tmpA);

            // tmpB = M/dt
            ierr = MatDuplicate(M, MAT_COPY_VALUES, &tmpB); PCHECK(ierr);
            ierr = MatScale(tmpB, 1.0/dt); PCHECK(ierr);
            if (theta < 1.0 - SMALL)
            {
                // CN (or any theta<1) needs the explicit half of the
                // diffusion operator on the RHS. For pure BE (theta=1) this
                // term vanishes and we keep B = M/dt only.
                ierr = MatAXPY
                (
                    tmpB, 1.0 - theta, L, DIFFERENT_NONZERO_PATTERN
                ); PCHECK(ierr);
            }
            BImplicitOwner.adopt(tmpB);
        }
        Mat AImplicit = AImplicitOwner.raw();
        Mat BImplicit = BImplicitOwner.raw();

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

            ionicModel->calculateCurrent
            (
                t0,
                dt,
                VmIntegrationPoints,
                IionIntegrationPoints,
                statesOld
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

        scalarField sourceOld(mesh.nCells(), 0.0);
        scalarField sourceGuess(mesh.nCells(), 0.0);
        forAll(sourceOld, cellI)
        {
            sourceOld[cellI] =
               -IionOld[cellI]
              + externalStimulusCurrent[cellI]/(chi.value()*Cm.value());
        }

        // PETSc Vec mirror of VmOld. Stays constant for the rest of the
        // time-step (the RHS at t^n is rebuilt against it every iteration).
        PetscVecRAII VnOwner; VnOwner.adopt(createPetscVec(mesh.nCells()));
        Vec Vn = VnOwner.raw();
        fieldToPetscVec(VmOld, Vn);

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
        //
        // During the upstroke this puts Newton ~50x closer to Vm^{n+1}
        // than the lagged Vm_n. The clamping prevents wild overshoots
        // when the previous step happened to include the I_Na peak.
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
                // matrix-free JFNK Jacobian. VmIntegrationPoints is the
                // working buffer for the ODE call so we clamp it in place;
                // this does not affect VmCandidate or any persistent field.
                if (jfnkClampODEInput)
                {
                    clampPhysical
                    (
                        VmIntegrationPoints,
                        jfnkInitGuessVmMin,
                        jfnkInitGuessVmMax
                    );
                }

                ionicModel->solveODE
                (
                    t0,
                    dt,
                    VmIntegrationPoints,
                    IionIntegrationPoints,
                    statesCandidate
                );

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
        // Re-usable PETSc Vec for the right-hand side; lives for the entire
        // nonlinear iteration and is overwritten in place each call.
        PetscVecRAII rhsBufferOwner;
        rhsBufferOwner.adopt(createPetscVec(mesh.nCells()));
        Vec rhsBuffer = rhsBufferOwner.raw();

        // Build the right-hand side of the implicit system for a given Iion
        // candidate into the provided PETSc Vec.
        //
        //   rhs = BImplicit * Vm^n  +  theta*src^{n+1}  +  (1-theta)*src^n
        //
        // No explicit dt factor is needed on the source terms: after writing
        // (Vm^{n+1} - Vm^n)/dt = ... and dividing both sides by dt, the
        // sources naturally appear with weights theta and (1-theta) only.
        // The reaction source is  src = -Iion/(chi*Cm) + Iext/(chi*Cm).
        auto rhsFromIion = [&]
        (
            const volScalarField& IionCandidate,
            Vec out
        )
        {
            PetscErrorCode ierr;
            // out <- BImplicit * Vn
            ierr = MatMult(BImplicit, Vn, out); PCHECK(ierr);

            // Add the cell-wise blended source in a single GetArray pass.
            PetscScalar* outArr = nullptr;
            ierr = VecGetArray(out, &outArr); PCHECK(ierr);
            const scalar invChiCm = 1.0/(chi.value()*Cm.value());
            forAll(sourceGuess, cellI)
            {
                const scalar sourceNp1 =
                   -IionCandidate[cellI]
                  + externalStimulusCurrent[cellI]*invChiCm;
                outArr[cellI] +=
                    theta*sourceNp1
                  + (1.0 - theta)*sourceOld[cellI];
            }
            ierr = VecRestoreArray(out, &outArr); PCHECK(ierr);
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
            const label lineSearchIters = -1  // -1 = not applicable (Picard / early exit)
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
                // Only print the line search count when we actually ran it
                // (i.e. in the JFNK Newton branch). For Picard / early-exit
                // paths the field is suppressed to keep the log tidy.
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

        // ----------------------------------------------------------------- //
        // Nonlinear solve via manual Newton + PETSc KSP (matrix-free).
        //
        // We retain the explicit Newton loop here (rather than handing the
        // whole problem to PETSc SNES) for two practical reasons:
        //   - Convergence diagnostics are written per Newton iteration with
        //     Vm / Iion / per-state residuals; SNES surfaces only one
        //     residual norm by default.
        //   - The ionic ODE inside formResidual() has side-effects on
        //     statesGuess and IionGuess which are easier to control from
        //     the host loop than from a SNES callback context struct.
        //
        // The expensive piece is the inner linear solve   J * delta = -R
        // with finite-difference Jacobian-vector products. We delegate that
        // entirely to PETSc:
        //   - A MatShell wraps a user matvec lambda implementing
        //     J*v ≈ (F(x + eps*v) - F(x)) / eps
        //   - KSP (GMRES with modified Gram-Schmidt, restart, ...) drives
        //     the Krylov iteration with restart, preconditioning, etc.
        //
        // For non-JFNK methods (Picard / diagonalIion) we assemble the
        // sparse Jacobian directly into a PETSc Mat and call solveLinearKSP.
        // ----------------------------------------------------------------- //

        // Persistent PETSc scratch vectors used inside the Newton loop.
        // Re-used across iterations to avoid create/destroy churn.
        PetscVecRAII xVecOwner;    xVecOwner.adopt(createPetscVec(mesh.nCells()));
        PetscVecRAII Rcurr;        Rcurr.adopt(createPetscVec(mesh.nCells()));
        PetscVecRAII RatX;         RatX.adopt(createPetscVec(mesh.nCells()));
        PetscVecRAII deltaVec;     deltaVec.adopt(createPetscVec(mesh.nCells()));
        PetscVecRAII negRcurr;     negRcurr.adopt(createPetscVec(mesh.nCells()));
        PetscVecRAII rhsVec;       rhsVec.adopt(createPetscVec(mesh.nCells()));
        PetscVecRAII matVecOutVec; matVecOutVec.adopt(createPetscVec(mesh.nCells()));
        // Backup of x used by the Armijo line search to roll back failed
        // trial steps without recomputing the Krylov solve.
        PetscVecRAII xBackup;      xBackup.adopt(createPetscVec(mesh.nCells()));
        Vec x        = xVecOwner.raw();
        Vec Rcurrent = Rcurr.raw();

        fieldToPetscVec(VmGuess, x);

        // Build F(x) into out: copies x → VmCandidate, integrates the ionic
        // ODE, computes  F = AImplicit * x - rhs(Iion(VmCandidate)).
        // When updateGuess is true the persistent VmGuess/statesGuess/
        // IionGuess fields are updated (used at the start of each Newton
        // iteration). Otherwise local copies are used so the matvec finite-
        // difference probes do not pollute the persistent state.
        auto computeResidual =
        [&]
        (
            Vec xCandidate,
            const bool updateGuess,
            Vec out
        )
        {
            if (updateGuess)
            {
                petscVecToField(xCandidate, VmGuess);
                VmGuess.correctBoundaryConditions();
                evaluateNonlinearFields(VmGuess, statesGuess, IionGuess, true);
                rhsFromIion(IionGuess, rhsVec.raw());
                PetscErrorCode ierr = MatMult(AImplicit, xCandidate, out); PCHECK(ierr);
                ierr = VecAXPY(out, -1.0, rhsVec.raw()); PCHECK(ierr);
            }
            else
            {
                // Probe call: do not modify persistent fields. We need a
                // temporary volScalarField for VmCandidate so we can call
                // evaluateNonlinearFields; statesTmp starts from statesOld.
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

                petscVecToField(xCandidate, VmTmp);
                VmTmp.correctBoundaryConditions();
                evaluateNonlinearFields(VmTmp, statesTmp, IionTmp, true);
                rhsFromIion(IionTmp, rhsVec.raw());
                PetscErrorCode ierr = MatMult(AImplicit, xCandidate, out); PCHECK(ierr);
                ierr = VecAXPY(out, -1.0, rhsVec.raw()); PCHECK(ierr);
            }
        };

        if (useJFNK)
        {
            // ----- Matrix-free shell for the Jacobian-vector product ------
            //
            // PETSc MatShell wraps an arbitrary matvec into a Mat-like
            // object. KSP can then use it as the system operator without
            // ever assembling J explicitly. We register a C-style trampoline
            // that forwards into a std::function held in the context.
            struct JfNKShellCtx
            {
                Vec x;                    // current Newton iterate
                Vec R;                    // F(x)
                Vec probeOut;             // scratch for F(x + eps*v)
                scalar epsScale;
                std::function<void(Vec, bool, Vec)> computeResidual;
            };
            JfNKShellCtx shellCtx;
            shellCtx.x = x;
            shellCtx.R = Rcurrent;
            shellCtx.probeOut = matVecOutVec.raw();
            shellCtx.computeResidual = computeResidual;

            auto trampoline = [](Mat shellMat, Vec v, Vec out) -> PetscErrorCode
            {
                JfNKShellCtx* c = nullptr;
                PetscErrorCode ierr = MatShellGetContext(shellMat, (void**)&c);
                if (ierr) return ierr;

                PetscReal vNorm = 0.0;
                ierr = VecNorm(v, NORM_2, &vNorm); CHKERRQ(ierr);
                const scalar eps =
                    c->epsScale / Foam::max(scalar(vNorm), SMALL);

                // x_perturbed = x + eps * v   (stored back in c->probeOut)
                ierr = VecWAXPY(c->probeOut, eps, v, c->x); CHKERRQ(ierr);
                // out = F(x + eps*v)         (probe call: no global side-effect)
                c->computeResidual(c->probeOut, false, out);
                // out = (F(x + eps*v) - F(x)) / eps
                ierr = VecAXPY(out, -1.0, c->R); CHKERRQ(ierr);
                ierr = VecScale(out, 1.0/eps); CHKERRQ(ierr);
                return 0;
            };

            // MatShell uses C linkage; the trampoline is a stateless lambda
            // so we can cast it to a plain function pointer.
            using TrampolineFn = PetscErrorCode (*)(Mat, Vec, Vec);
            TrampolineFn matMulFn = trampoline;

            Mat Jshell;
            PetscErrorCode ierr;
            ierr = MatCreateShell
            (
                PETSC_COMM_SELF,
                mesh.nCells(), mesh.nCells(),
                mesh.nCells(), mesh.nCells(),
                &shellCtx,
                &Jshell
            ); PCHECK(ierr);
            ierr = MatShellSetOperation
            (
                Jshell, MATOP_MULT, reinterpret_cast<void(*)()>(matMulFn)
            ); PCHECK(ierr);

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

                // R = F(x)  (also updates VmGuess, statesGuess, IionGuess)
                computeResidual(x, true, Rcurrent);

                const scalar coupledResidual = relativeL2Norm(Rcurrent, x);
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

                // Early-exit: if the current iterate already satisfies the
                // tolerances we skip the linear solve. Often fires on the
                // first Newton step when the previous time-step solution is
                // already close to the implicit solution.
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
                        corr + 1, 0, 0.0,
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

                // Update the matvec context so the trampoline sees the
                // current Newton iterate and residual.
                PetscReal xNorm = 0.0;
                ierr = VecNorm(x, NORM_2, &xNorm); PCHECK(ierr);
                shellCtx.epsScale = jfnkEpsilon*(1.0 + scalar(xNorm));

                // negRcurr = -Rcurrent   (right-hand side of J*delta = -R)
                ierr = VecCopy(Rcurrent, negRcurr.raw()); PCHECK(ierr);
                ierr = VecScale(negRcurr.raw(), -1.0); PCHECK(ierr);

                // KSP solve  J * delta = -R  with the MatShell as operator.
                // PCNONE because we have no explicit matrix to factor; one
                // could plug a Jacobi/SOR/AMG preconditioner here if a true
                // operator were assembled. KSPGMRES with MGS + restart is
                // the safe default for stiff ionic Jacobians.
                KSP ksp;
                ierr = KSPCreate(PETSC_COMM_SELF, &ksp); PCHECK(ierr);
                ierr = KSPSetOperators(ksp, Jshell, Jshell); PCHECK(ierr);
                ierr = KSPSetType(ksp, KSPGMRES); PCHECK(ierr);
                ierr = KSPGMRESSetRestart
                (
                    ksp, static_cast<PetscInt>(jfnkMaxKrylovIterations)
                ); PCHECK(ierr);
                ierr = KSPGMRESSetOrthogonalization
                (
                    ksp,
                    KSPGMRESModifiedGramSchmidtOrthogonalization
                ); PCHECK(ierr);

                PC pc;
                ierr = KSPGetPC(ksp, &pc); PCHECK(ierr);
                ierr = PCSetType(pc, PCNONE); PCHECK(ierr);

                const PetscInt totalIters = static_cast<PetscInt>
                (
                    jfnkMaxKrylovIterations*(jfnkMaxRestarts + 1)
                );
                ierr = KSPSetTolerances
                (
                    ksp,
                    jfnkLinearTolerance,    // relative
                    PETSC_DEFAULT,
                    PETSC_DEFAULT,
                    totalIters
                ); PCHECK(ierr);
                ierr = KSPSetFromOptions(ksp); PCHECK(ierr);

                ierr = KSPSolve(ksp, negRcurr.raw(), deltaVec.raw());
                PCHECK(ierr);

                PetscInt kspIts = 0;
                PetscReal kspRnorm = 0.0;
                ierr = KSPGetIterationNumber(ksp, &kspIts); PCHECK(ierr);
                ierr = KSPGetResidualNorm(ksp, &kspRnorm); PCHECK(ierr);
                ierr = KSPDestroy(&ksp); PCHECK(ierr);

                const label gmresIterations = static_cast<label>(kspIts);
                const scalar gmresError = scalar(kspRnorm);

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
                scalar coupledResidualNew = GREAT;

                if (jfnkLineSearch)
                {
                    // Snapshot of the current iterate; restored on each
                    // backtracking trial.
                    ierr = VecCopy(x, xBackup.raw()); PCHECK(ierr);

                    // Squared norm of F at the *current* iterate, against
                    // which we test the Armijo condition. Rcurrent was
                    // computed at the top of this corr iteration.
                    PetscReal RoldNorm = 0.0;
                    ierr = VecNorm(Rcurrent, NORM_2, &RoldNorm); PCHECK(ierr);
                    const PetscReal Rold2 = RoldNorm*RoldNorm;

                    PetscReal alpha = nonlinearRelaxation;

                    while (true)
                    {
                        // x = xBackup + alpha * delta
                        ierr = VecCopy(xBackup.raw(), x); PCHECK(ierr);
                        ierr = VecAXPY(x, alpha, deltaVec.raw());
                        PCHECK(ierr);

                        // F(x_new); commits VmGuess/statesGuess/IionGuess
                        computeResidual(x, true, RatX.raw());

                        PetscReal RnewNorm = 0.0;
                        ierr = VecNorm
                        (
                            RatX.raw(), NORM_2, &RnewNorm
                        ); PCHECK(ierr);
                        const PetscReal Rnew2 = RnewNorm*RnewNorm;

                        // Armijo sufficient-decrease test
                        if
                        (
                            Rnew2
                         <= (1.0 - 2.0*jfnkArmijoC*alpha)*Rold2
                        )
                        {
                            break;  // accept trial
                        }

                        ++lineSearchIters;

                        if (lineSearchIters >= jfnkLineSearchMaxIter)
                        {
                            break;  // budget exhausted; outer loop will
                                     // notice non-convergence and the
                                     // step rollback will catch it
                        }

                        alpha *= 0.5;
                        if (alpha < jfnkLineSearchAlphaMin)
                        {
                            // alpha floored; do one final trial with the
                            // floor and accept whatever comes out.
                            alpha = jfnkLineSearchAlphaMin;
                            ierr = VecCopy(xBackup.raw(), x); PCHECK(ierr);
                            ierr = VecAXPY(x, alpha, deltaVec.raw());
                            PCHECK(ierr);
                            computeResidual(x, true, RatX.raw());
                            break;
                        }
                    }
                    coupledResidualNew = relativeL2Norm(RatX.raw(), x);
                }
                else
                {
                    // Legacy path: unconditional damped step (preserves
                    // reproducibility with earlier runs when explicitly
                    // requested via jfnkLineSearch false).
                    ierr = VecAXPY(x, nonlinearRelaxation, deltaVec.raw());
                    PCHECK(ierr);
                    computeResidual(x, true, RatX.raw());
                    coupledResidualNew = relativeL2Norm(RatX.raw(), x);
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
                    gmresIterations, gmresError,
                    coupledResidualNew,
                    VmResidual,
                    IionResidual,
                    maxStateResidual,
                    stateResiduals,
                    converged,
                    lineSearchIters    // exposes Armijo backtracks
                );

                nonlinearIters = corr + 1;
                if (converged)
                {
                    nonlinearConverged = true;
                    break;
                }
            }

            ierr = MatDestroy(&Jshell); PCHECK(ierr);
        }
        else
        {
            // -------------- Picard / diagonalIion branch ------------------
            // Solve the implicit system at each iterate by direct KSP on
            // AImplicit (optionally augmented with a diagonal Jacobian of
            // Iion to accelerate convergence on stiff sources).
            PetscVecRAII VsolVec;  VsolVec.adopt(createPetscVec(mesh.nCells()));
            PetscVecRAII VprevVec; VprevVec.adopt(createPetscVec(mesh.nCells()));
            PetscVecRAII nonlinResVec;
            nonlinResVec.adopt(createPetscVec(mesh.nCells()));
            PetscVecRAII rhsTmpVec;
            rhsTmpVec.adopt(createPetscVec(mesh.nCells()));

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
                rhsFromIion(IionGuess, rhsTmpVec.raw());

                // Duplicate AImplicit so we can locally add the diagonal
                // Iion linearisation without polluting the operator used by
                // the next iteration.
                Mat ACurrent = nullptr;
                PetscErrorCode ierr;
                ierr = MatDuplicate
                (
                    AImplicit, MAT_COPY_VALUES, &ACurrent
                ); PCHECK(ierr);

                if (useDiagonalIion)
                {
                    // Augment ACurrent with -theta * d(source)/dVm on the
                    // diagonal, and subtract the corresponding RHS shift.
                    computeSourceDerivative(VmGuess, IionGuess);
                    addDiagonalToMatrix
                    (
                        sourceDerivative.primitiveField(),
                       -theta,
                        ACurrent
                    );

                    PetscScalar* rhsArr = nullptr;
                    ierr = VecGetArray(rhsTmpVec.raw(), &rhsArr); PCHECK(ierr);
                    const scalarField& dF = sourceDerivative.primitiveField();
                    forAll(dF, cellI)
                    {
                        rhsArr[cellI] -=
                            theta*dF[cellI]*VmGuess.primitiveField()[cellI];
                    }
                    ierr = VecRestoreArray(rhsTmpVec.raw(), &rhsArr);
                    PCHECK(ierr);
                }

                label linearIterations = 0;
                scalar linearError = GREAT;
                solveLinearKSP
                (
                    ACurrent,
                    rhsTmpVec.raw(),
                    VsolVec.raw(),
                    implicitLinearSolver,
                    implicitPreconditioner,
                    implicitMaxIterations,  // krylov dim ~ maxIter for Picard
                    0,                       // no restart for direct LU
                    implicitTolerance,
                    implicitMaxIterations,
                    linearIterations,
                    linearError
                );
                ierr = MatDestroy(&ACurrent); PCHECK(ierr);

                // Damped update:  x_new = (1-omega) * x_old + omega * Vsol
                fieldToPetscVec(VmGuess, VprevVec.raw());
                ierr = VecCopy(VprevVec.raw(), x); PCHECK(ierr);
                ierr = VecAXPBY
                (
                    x,
                    nonlinearRelaxation,           // alpha (on Vsol)
                    1.0 - nonlinearRelaxation,     // beta  (on x = Vprev)
                    VsolVec.raw()
                ); PCHECK(ierr);

                petscVecToField(x, VmGuess);
                VmGuess.correctBoundaryConditions();
                evaluateNonlinearFields(VmGuess, statesGuess, IionGuess, true);

                // Coupled nonlinear residual = AImplicit*VmGuess - rhs
                rhsFromIion(IionGuess, rhsTmpVec.raw());
                ierr = MatMult(AImplicit, x, nonlinResVec.raw()); PCHECK(ierr);
                ierr = VecAXPY(nonlinResVec.raw(), -1.0, rhsTmpVec.raw());
                PCHECK(ierr);
                const scalar coupledResidual =
                    relativeL2Norm(nonlinResVec.raw(), rhsTmpVec.raw());

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
                    linearIterations, linearError,
                    coupledResidual,
                    VmResidual, IionResidual,
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
            if (useHighOrder_Iion)
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

    runTime.printExecutionTime(Info);

    Info<< "End" << nl << endl;

    }  // <-- inner scope end: Time/mesh/fields/RAII destructors run here.
       //                      MPI is still alive (OpenFOAM owns it).

    // Tear down PETSc. Because OpenFOAM initialised MPI (step 1 above),
    // PETSc tracks that it did not own MPI and PetscFinalize() will NOT
    // call MPI_Finalize. argList's destructor (after main returns) is
    // therefore free to perform its own MPI cleanup.
    PetscFinalize();

    return 0;
}
