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

#include "fvCFD.H"
#include "ionicModel.H"
#include "TNNP.H"
#include "LRE.H"
#include "Field.H"
#include "volFields.H"
#include <cmath>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <vector>

namespace
{
    using SpMat = Eigen::SparseMatrix<scalar, Eigen::RowMajor>;
    using Triplet = Eigen::Triplet<scalar>;
    using EigVec = Eigen::Matrix<scalar, Eigen::Dynamic, 1>;

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
                // Happy breakdown (H(j+1,j) ≈ 0) means the Krylov space is
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
        triplets.reserve(mesh.nCells()*40);

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

    void assembleStandardOrthogonalStiffnessMatrix
    (
        const fvMesh& mesh,
        const volTensorField& conductivity,
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
            const scalar a =
                orthogonalDiffusionCoeff(mesh.Sf()[faceI], C[nei] - C[own], Df);

            addTripletIfNeeded(triplets, own, own, -a/max(V[own], SMALL));
            addTripletIfNeeded(triplets, own, nei,  a/max(V[own], SMALL));
            addTripletIfNeeded(triplets, nei, own,  a/max(V[nei], SMALL));
            addTripletIfNeeded(triplets, nei, nei, -a/max(V[nei], SMALL));
        }

        const surfaceVectorField& Cf = mesh.Cf();

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
                }
            }
        }

        K.resize(mesh.nCells(), mesh.nCells());
        K.setFromTriplets(triplets.begin(), triplets.end());
        K.makeCompressed();
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
        SpMat& K
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

        std::vector<Triplet> triplets;
        triplets.reserve(mesh.nFaces()*40);

        // -------- Internal faces -------------------------------------------
        forAll(neighbour, faceI)
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

                // For each stencil cell j we accumulate its contribution
                // to the flux through this face / cell volume.
                forAll(curStencil, cI)
                {
                    const label col = curStencil[cI];
                    const vector gCoeff = faceGradCoeffs[faceI][qpI][cI];
                    // For homogeneous tissue conductivity[own] == conductivity[nei];
                    // for heterogeneous media a face-averaged conductivity
                    // would be more accurate.
                    const scalar fluxCoeff =
                        area*w*(n & (conductivity[own] & gCoeff));

                    // Divergence theorem: +flux to owner, -flux to neighbour.
                    addTripletIfNeeded(triplets, own, col, fluxCoeff/max(V[own], SMALL));
                    addTripletIfNeeded(triplets, nei, col, -fluxCoeff/max(V[nei], SMALL));
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
            }
        }

        K.resize(mesh.nCells(), mesh.nCells());
        K.setFromTriplets(triplets.begin(), triplets.end());
        K.makeCompressed();
    }

    EigVec solveSparseSystem
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
            runTime.path()/"postProcessing"/"highOrderElectroActivationFoamImplicitJfNK"
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
    #include "setRootCaseLists.H"
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

    SpMat M;
    SpMat K;
    SpMat L;

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
        assembleConsistentMassMatrixHO(mesh, LREInterp_VmPtr(), 1.0, M);
    }

    if (useHighOrder_Vm && dim > 1)
    {
        assembleHighOrderStiffnessMatrix(mesh, conductivity, LREInterp_VmPtr(), K);
    }
    else
    {
        assembleStandardOrthogonalStiffnessMatrix(mesh, conductivity, K);
    }

    L = K;
    L *= 1.0/(chi.value()*Cm.value());

    autoPtr<OFstream> nonlinearResidualFilePtr;
    if (writeNonlinearResiduals)
    {
        const fileName residualDir
        (
            runTime.path()/"postProcessing"/"highOrderElectroActivationFoamImplicitJfNK"
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
                // Standard JfNK residual: F(Vm) = A·Vm - rhs(Iion(Vm))
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
                // product:    J*v ≈ (F(x + eps*v) - F(x)) / eps
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
                const EigVec delta = solveGMRES
                (
                    matVec,
                    -R,
                    jfnkMaxKrylovIterations,
                    jfnkMaxRestarts,
                    jfnkLinearTolerance,
                    gmresIterations,
                    gmresError
                );

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
                    implicitLinearSolver,
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

    return 0;
}
