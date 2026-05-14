#include "fvCFD.H"
#include "LRE.H"
#include <cmath>
#include <chrono>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <vector>


namespace
{
    using SpMat = Eigen::SparseMatrix<scalar, Eigen::RowMajor>;
    using Triplet = Eigen::Triplet<scalar>;
    using EigVec = Eigen::Matrix<scalar, Eigen::Dynamic, 1>;

    struct MMSConfig
    {
        word type;      // example_0 | example_1 | example_2 | example_3 | example_4 | example_5
        scalar A;
        scalar beta;    // used by example_0 | example_1 | example_2 | example_3
        scalar gamma;   // used by example_4
        scalar omega;   // used by example_4
        label dim;
    };

    label activeDim(const MMSConfig& mms)
    {
        return max(mms.dim, label(1));
    }

    scalar zFactor
    (
        const point& p,
        const scalar k,
        const bool useCos,
        const MMSConfig& mms
    )
    {
        if (activeDim(mms) < 3)
        {
            return 1.0;
        }

        return useCos ? std::cos(k*p.z()) : std::sin(k*p.z());
    }

    scalar electroF(const point& p, const MMSConfig& mms)
    {
        const scalar pi = constant::mathematical::pi;

        scalar f = std::cos(pi*p.x());

        if (activeDim(mms) >= 2)
        {
            f *= std::cos(2.0*pi*p.y());
        }
        if (activeDim(mms) >= 3)
        {
            f *= std::cos(3.0*pi*p.z());
        }

        return f;
    }

    scalar electroG(const point& p, const MMSConfig& mms)
    {
        if (activeDim(mms) == 1)
        {
            return 1.0 + p.x();
        }
        else if (activeDim(mms) == 2)
        {
            return 1.0 + p.x()*sqr(p.y());
        }

        return 1.0 + p.x()*sqr(p.y())*pow3(p.z());
    }

    scalar electroBeta
    (
        const scalar alphaX,
        const scalar alphaY,
        const scalar alphaZ,
        const MMSConfig& mms
    )
    {
        const scalar pi2 = sqr(constant::mathematical::pi);

        if (activeDim(mms) == 1)
        {
            return -pi2*alphaX;
        }
        else if (activeDim(mms) == 2)
        {
            return -pi2*(alphaX + 4.0*alphaY);
        }

        return -pi2*(alphaX + 4.0*alphaY + 9.0*alphaZ);
    }

    scalar exactElectroIion
    (
        const point& p,
        const scalar t,
        const scalar alphaX,
        const scalar alphaY,
        const scalar alphaZ,
        const MMSConfig& mms
    )
    {
        const scalar T = std::sqrt(1.0 + t)*electroF(p, mms);
        const scalar G = electroG(p, mms);
        const scalar u1 = (1.0 + t)*G + T;
        const scalar u2 = 1.0/((1.0 + t)*std::sqrt(G));
        const scalar u3 = 0.0;

        return
            -0.5*(u1 + u3 - T)*sqr(u2)*(T - u3)
          + electroBeta(alphaX, alphaY, alphaZ, mms)*(T - u3);
    }

    scalar exactTemperature
    (
        const point& p,
        const scalar t,
        const MMSConfig& mms
    )
    {
        const scalar x = p.x();
        const scalar y = p.y();
        const scalar pi = constant::mathematical::pi;

        if (mms.type == "example_0")
        {
            return
                mms.A
                *std::sin(mms.beta*pi*x)
                *std::sin(mms.beta*pi*y)
                *zFactor(p, mms.beta*pi, false, mms)
                *std::sin(mms.beta*pi*t);
        }
        else if (mms.type == "example_1")
        {
            return
                mms.A
                *std::cos(mms.beta*pi*x)
                *std::sin(mms.beta*pi*y)
                *zFactor(p, mms.beta*pi, false, mms)
                *std::sin(mms.beta*pi*t);
        }
        else if (mms.type == "example_2")
        {
            return
                mms.A
                *std::cos(mms.beta*pi*x)
                *std::cos(mms.beta*pi*y)
                *zFactor(p, mms.beta*pi, false, mms)
                *std::sin(mms.beta*pi*t);
        }
        else if (mms.type == "example_3")
        {
            return
                mms.A
                *std::cos(mms.beta*pi*x)
                *std::cos(mms.beta*pi*y)
                *zFactor(p, mms.beta*pi, true, mms)
                *std::sin(mms.beta*pi*t);
        }
        else if (mms.type == "example_4")
        {
            const scalar k = 0.5*mms.gamma*pi;

            return
                mms.A
                *std::sin(k*x)
                *std::sin(k*y)
                *zFactor(p, k, false, mms)
                *std::sin(mms.omega*pi*t);
        }
        else if (mms.type == "example_5")
        {
            return mms.A*std::sqrt(1.0 + t)*electroF(p, mms);
        }

        FatalErrorInFunction
            << "Unknown mmsType: " << mms.type << nl
            << "Valid options: example_0 | example_1 | example_2 | example_3 | example_4 | example_5"
            << exit(FatalError);

        return 0.0;
    }

    scalar exactSource
    (
        const point& p,
        const scalar t,
        const scalar alphaX,
        const scalar alphaY,
        const scalar alphaZ,
        const scalar rhoCp,
        const MMSConfig& mms
    )
    {
        const scalar x = p.x();
        const scalar y = p.y();
        const scalar pi = constant::mathematical::pi;

        if (mms.type == "example_0")
        {
            const scalar a = rhoCp*mms.beta*pi*std::cos(mms.beta*pi*t);
            const scalar b =
                (alphaX + alphaY + (activeDim(mms) >= 3 ? alphaZ : 0.0))
               *pow(mms.beta*pi,2)
               *std::sin(mms.beta*pi*t);

            return
                mms.A
                *std::sin(mms.beta*pi*x)
                *std::sin(mms.beta*pi*y)
                *zFactor(p, mms.beta*pi, false, mms)
                *(a + b);
        }
        else if (mms.type == "example_1")
        {
            const scalar a = rhoCp*mms.beta*pi*std::cos(mms.beta*pi*t);
            const scalar b =
                (alphaX + alphaY + (activeDim(mms) >= 3 ? alphaZ : 0.0))
               *pow(mms.beta*pi,2)
               *std::sin(mms.beta*pi*t);

            return
                mms.A
                *std::cos(mms.beta*pi*x)
                *std::sin(mms.beta*pi*y)
                *zFactor(p, mms.beta*pi, false, mms)
                *(a + b);
        }
        else if (mms.type == "example_2")
        {
            const scalar a = rhoCp*mms.beta*pi*std::cos(mms.beta*pi*t);
            const scalar b =
                (alphaX + alphaY + (activeDim(mms) >= 3 ? alphaZ : 0.0))
               *pow(mms.beta*pi,2)
               *std::sin(mms.beta*pi*t);

            return
                mms.A
                *std::cos(mms.beta*pi*x)
                *std::cos(mms.beta*pi*y)
                *zFactor(p, mms.beta*pi, false, mms)
                *(a + b);
        }
        else if (mms.type == "example_3")
        {
            const scalar a = rhoCp*mms.beta*pi*std::cos(mms.beta*pi*t);
            const scalar b =
                (alphaX + alphaY + (activeDim(mms) >= 3 ? alphaZ : 0.0))
               *pow(mms.beta*pi,2)
               *std::sin(mms.beta*pi*t);

            return
                mms.A
                *std::cos(mms.beta*pi*x)
                *std::cos(mms.beta*pi*y)
                *zFactor(p, mms.beta*pi, true, mms)
                *(a + b);
        }
        else if (mms.type == "example_4")
        {
            const scalar k = 0.5*mms.gamma*pi;
            const scalar spatial =
                std::sin(k*x)*std::sin(k*y)*zFactor(p, k, false, mms);

            const scalar T = mms.A*spatial*std::sin(mms.omega*pi*t);
            const scalar dTdt = mms.A*spatial*mms.omega*pi*std::cos(mms.omega*pi*t);

            const scalar divKGradT =
                - (alphaX + alphaY + (activeDim(mms) >= 3 ? alphaZ : 0.0))
                *pow(k,2)*T;

            return rhoCp*dTdt - divKGradT;
        }
        else if (mms.type == "example_5")
        {
            return
                mms.A
               *exactElectroIion
                (
                    p,
                    t,
                    alphaX,
                    alphaY,
                    alphaZ,
                    mms
                );
        }

        FatalErrorInFunction
            << "Unknown mmsType: " << mms.type << nl
            << "Valid options: example_0 | example_1 | example_2 | example_3 | example_4 | example_5"
            << exit(FatalError);

        return 0.0;
    }

    label estimatedN(const fvMesh& mesh)
    {
        const label nCellsGlobal =
            returnReduce(mesh.nCells(), sumOp<label>());

        const scalar dim = max(scalar(mesh.nGeometricD()), 1.0);
        return label(std::round(std::pow(scalar(nCellsGlobal), 1.0/dim)));
    }

    scalar characteristicDx(const fvMesh& mesh)
    {
        const boundBox bb(mesh.points());
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
            scalar(returnReduce(mesh.nCells(), sumOp<label>()));

        return std::pow
        (
            activeMeasure/max(nCellsGlobal, SMALL),
            1.0/max(scalar(mesh.nGeometricD()), 1.0)
        );
    }

    scalar volumeWeightedL1(const volScalarField& fld)
    {
        const scalarField& f = fld.primitiveField();
        const scalarField& V = fld.mesh().V();

        return gSum(mag(f)*V)/gSum(V);
    }

    scalar volumeWeightedL2(const volScalarField& fld)
    {
        const scalarField& f = fld.primitiveField();
        const scalarField& V = fld.mesh().V();

        return std::sqrt(gSum(sqr(f)*V)/gSum(V));
    }

    scalar linfNorm(const volScalarField& fld)
    {
        const scalarField magFld(mag(fld.primitiveField()));
        return gMax(magFld);
    }

    scalar errorTCellCentred(const volScalarField& fld)
    {
        const scalarField& f = fld.primitiveField();
        return std::sqrt(gSum(sqr(f)));
    }

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
            << "Valid options are lumped (or diagonal) and consistent (or consistentHO)"
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

    EigVec sourceToEigVec(const volScalarField& sourceQ)
    {
        return fieldToEigVec(sourceQ);
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
        const scalar rhoCp,
        SpMat& M
    )
    {
        std::vector<Triplet> triplets;
        triplets.reserve(mesh.nCells());

        forAll(mesh.C(), cellI)
        {
            triplets.emplace_back(cellI, cellI, rhoCp);
        }

        M.resize(mesh.nCells(), mesh.nCells());
        M.setFromTriplets(triplets.begin(), triplets.end());
        M.makeCompressed();
    }

    void assembleConsistentMassMatrixHO
    (
        const fvMesh& mesh,
        const LRE& LREInterp_T,
        const scalar rhoCp,
        SpMat& M
    )
    {
        const bool twoD = mesh.nGeometricD() == 2;
        const vectorField& C = mesh.C();

        const labelListList& stencils = LREInterp_T.globalCellStencils();
        const CompactListList<point>& cellQP = LREInterp_T.cellQuadPoints();
        const CompactListList<scalar>& cellQW = LREInterp_T.cellQuadWeight();
        const CompactListList<vector>& gradCoeffs = LREInterp_T.QRGradCoeffs();
        const CompactListList<symmTensor>& hessCoeffs =
            LREInterp_T.cellHessianCoeffs();
        const CompactListList<LRE::symmTensor3Order>& thirdCoeffs =
            LREInterp_T.cellThirdDerivCoeffs();

        std::vector<Triplet> triplets;
        triplets.reserve(mesh.nCells()*40);

        forAll(stencils, cellI)
        {
            const labelList& curStencil = stencils[cellI];
            const label selfCoeffI = curStencil.size();

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

                    if (LREInterp_T.order() >= 2)
                    {
                        coeff += 0.5*quadraticForm(hessCoeffs[cellI][cI], d);
                    }

                    if (LREInterp_T.order() >= 3)
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
                        rhoCp*w*coeff
                    );
                }

                scalar selfCoeff = 1.0 + (gradCoeffs[cellI][selfCoeffI] & d);

                if (LREInterp_T.order() >= 2)
                {
                    selfCoeff +=
                        0.5*quadraticForm(hessCoeffs[cellI][selfCoeffI], d);
                }

                if (LREInterp_T.order() >= 3)
                {
                    selfCoeff +=
                        (1.0/6.0)
                       *LRE::cubicForm(thirdCoeffs[cellI][selfCoeffI], d, twoD);
                }

                addTripletIfNeeded
                (
                    triplets,
                    cellI,
                    cellI,
                    rhoCp*w*selfCoeff
                );
            }
        }

        M.resize(mesh.nCells(), mesh.nCells());
        M.setFromTriplets(triplets.begin(), triplets.end());
        M.makeCompressed();
    }

    void assembleHighOrderStiffnessMatrix
    (
        const fvMesh& mesh,
        const volTensorField& conductivity,
        const LRE& LREInterp_T,
        SpMat& K
    )
    {
        const CompactListList<point>& faceQP = LREInterp_T.faceQuadPoints();
        const CompactListList<scalar>& faceQW = LREInterp_T.faceQuadWeight();
        const labelListList& faceStencils = LREInterp_T.globalFaceStencils();
        const List<CompactListList<vector>>& faceGradCoeffs =
            LREInterp_T.QRGradFaceGPCoeffs();

        const scalarField& V = mesh.V();
        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();

        std::vector<Triplet> triplets;
        triplets.reserve(mesh.nFaces()*40);

        forAll(owner, faceI)
        {
            const label own = owner[faceI];
            const label nei = neighbour[faceI];

            const vector Sf = mesh.Sf()[faceI];
            const scalar area = mag(Sf) + VSMALL;
            const vector n = Sf/area;

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

        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            const word bcType = patch.lookupPatchField<volScalarField, scalar>("T").type();

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

                        addTripletIfNeeded
                        (
                            triplets,
                            own,
                            col,
                            fluxCoeff/max(V[own], SMALL)
                        );
                    }
                }
            }
        }

        K.resize(mesh.nCells(), mesh.nCells());
        K.setFromTriplets(triplets.begin(), triplets.end());
        K.makeCompressed();
    }

    EigVec assembleHighOrderBoundaryVector
    (
        const fvMesh& mesh,
        const volTensorField& conductivity,
        const LRE& LREInterp_T,
        const scalar t,
        const MMSConfig& mms
    )
    {
        EigVec b = EigVec::Zero(mesh.nCells());

        const CompactListList<point>& faceQP = LREInterp_T.faceQuadPoints();
        const CompactListList<scalar>& faceQW = LREInterp_T.faceQuadWeight();
        const labelListList& faceStencils = LREInterp_T.globalFaceStencils();
        const List<CompactListList<vector>>& faceGradCoeffs =
            LREInterp_T.QRGradFaceGPCoeffs();

        const scalarField& V = mesh.V();
        const labelUList& owner = mesh.owner();

        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            const word bcType = patch.lookupPatchField<volScalarField, scalar>("T").type();

            if (bcType != "fixedValue" && bcType != "fixedVoltage")
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
                    const scalar Tbc =
                        exactTemperature(faceQP[globalFaceI][qpI], t, mms);

                    const vector gGhost =
                        faceGradCoeffs[globalFaceI][qpI][ghostID];

                    const scalar fluxCoeff =
                        area*w*(n & (conductivity[own] & gGhost));

                    b[own] += fluxCoeff/max(V[own], SMALL)*Tbc;
                }
            }
        }

        return b;
    }

    EigVec solveSparseSystem
    (
        const SpMat& A,
        const EigVec& b,
        const word& linearSolver,
        const scalar tol,
        const label maxIter
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

            Info<< "Implicit HO solve: iterations = "
                << solver.iterations()
                << ", estimated error = "
                << solver.error() << nl << endl;

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

    void applyExactBoundaryValues
    (
        volScalarField& T,
        const scalar t,
        const MMSConfig& mms
    )
    {
        const fvMesh& mesh = T.mesh();
        const surfaceVectorField& Cf = mesh.Cf();

        forAll(T.boundaryField(), patchI)
        {
            const word bcType = T.boundaryField()[patchI].type();

            if (bcType == "fixedValue" || bcType == "fixedVoltage")
            {
                fvPatchScalarField& Tp = T.boundaryFieldRef()[patchI];
                const vectorField& CfPatch = Cf.boundaryField()[patchI];

                forAll(Tp, faceI)
                {
                    Tp[faceI] = exactTemperature(CfPatch[faceI], t, mms);
                }
            }
        }

        T.correctBoundaryConditions();
    }

    void eigVecToField(const EigVec& v, volScalarField& fld)
    {
        scalarField& f = fld.primitiveFieldRef();

        forAll(f, cellI)
        {
            f[cellI] = v[cellI];
        }
    }

    void fillExactFields
    (
        volScalarField& TExact,
        volScalarField& sourceQ,
        const scalar t,
        const scalar alphaX,
        const scalar alphaY,
        const scalar alphaZ,
        const scalar rhoCp,
        const MMSConfig& mms,
        const LRE& LREInterp_T
    )
    {
        const fvMesh& mesh = TExact.mesh();
        const vectorField& C = mesh.C();
        const CompactListList<point>& cellQP = LREInterp_T.cellQuadPoints();
        const CompactListList<scalar>& cellQW = LREInterp_T.cellQuadWeight();

        forAll(C, cellI)
        {
            TExact[cellI] = exactTemperature(C[cellI], t, mms);

            scalar qBar = 0.0;
            scalar wSum = 0.0;

            forAll(cellQP[cellI], qpI)
            {
                qBar +=
                    cellQW[cellI][qpI]
                   *exactSource
                    (
                        cellQP[cellI][qpI],
                        t,
                        alphaX,
                        alphaY,
                        alphaZ,
                        rhoCp,
                        mms
                    );

                wSum += cellQW[cellI][qpI];
            }

            sourceQ[cellI] = qBar/max(wSum, SMALL);
        }

        const surfaceVectorField& Cf = mesh.Cf();

        forAll(TExact.boundaryField(), patchI)
        {
            fvPatchScalarField& Tp = TExact.boundaryFieldRef()[patchI];
            const vectorField& CfPatch = Cf.boundaryField()[patchI];
            const word bcType = Tp.type();

            if (bcType == "fixedValue" || bcType == "fixedVoltage")
            {
                forAll(Tp, faceI)
                {
                    Tp[faceI] = exactTemperature(CfPatch[faceI], t, mms);
                }
            }
        }

        forAll(sourceQ.boundaryField(), patchI)
        {
            fvPatchScalarField& qp = sourceQ.boundaryFieldRef()[patchI];

            forAll(qp, faceI)
            {
                qp[faceI] = 0.0;
            }
        }
    }

    void computeHighOrderLaplacian
    (
        const volScalarField& T,
        const volTensorField& conductivity,
        const LRE& LREInterp_T,
        surfaceScalarField& fluxT_HO,
        volScalarField& lapTHO
    )
    {
        const fvMesh& mesh = T.mesh();

        autoPtr<List<List<vector>>> gradTQuadPtr =
            LREInterp_T.gradScalarFaceQuad(T);

        List<List<vector>>& gradTQuad = gradTQuadPtr.ref();
        const CompactListList<scalar>& faceQuadW =
            LREInterp_T.faceQuadWeight();

        const surfaceVectorField nHat(mesh.Sf()/mesh.magSf());
        const scalarField& magSfInternal = mesh.magSf().internalField();

        for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
        {
            const label owner = mesh.owner()[faceI];
            const vector& faceNormal = nHat[faceI];
            const scalar faceArea = magSfInternal[faceI];

            fluxT_HO[faceI] = 0.0;

            forAll(gradTQuad[faceI], pI)
            {
                const vector Kg = conductivity[owner] & gradTQuad[faceI][pI];

                fluxT_HO[faceI] +=
                    faceArea*(faceNormal & Kg)*faceQuadW[faceI][pI];
            }
        }

        forAll(fluxT_HO.boundaryField(), patchI)
        {
            scalarField& patchFlux = fluxT_HO.boundaryFieldRef()[patchI];

            if (patchFlux.size() == 0)
            {
                continue;
            }

            const word bcType = T.boundaryField()[patchI].type();
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

                forAll(gradTQuad[globalFaceI], pI)
                {
                    const vector Kg =
                        conductivity[owner] & gradTQuad[globalFaceI][pI];

                    patchFlux[faceI] +=
                        faceArea*(faceNormal & Kg)*faceQuadW[globalFaceI][pI];
                }
            }
        }

        lapTHO = fvc::div(fluxT_HO);
    }

    struct GaussErrorSummary
    {
        scalar L1;
        scalar L2;
        scalar Linf;
        scalar normL2;
        scalar relL2; // %
    };

    GaussErrorSummary computeGaussReconstructedError
    (
        const volScalarField& T,
        const scalar t,
        const MMSConfig& mms,
        const LRE& LREInterp_T
    )
    {
        const fvMesh& mesh = T.mesh();
        const bool twoD = mesh.nGeometricD() == 2;

        const CompactListList<point>& cellQP = LREInterp_T.cellQuadPoints();
        const CompactListList<scalar>& cellQW = LREInterp_T.cellQuadWeight();

        tmp<volVectorField> tGradT = LREInterp_T.grad(T);
        const volVectorField& gradT = tGradT();

        tmp<volSymmTensorField> tHessT;
        const volSymmTensorField* hessTPtr = nullptr;
        if (LREInterp_T.order() >= 2)
        {
            tHessT = LREInterp_T.hessian(T);
            hessTPtr = &tHessT();
        }

        autoPtr<List<LRE::symmTensor3Order>> thirdDerivPtr;
        const List<LRE::symmTensor3Order>* thirdTPtr = nullptr;
        if (LREInterp_T.order() >= 3)
        {
            thirdDerivPtr = LREInterp_T.thirdDeriv(T);
            thirdTPtr = &thirdDerivPtr();
        }

        const vectorField& C = mesh.C();
        const scalarField& V = mesh.V();

        scalar e1 = 0.0, e2 = 0.0, n2 = 0.0, volTot = 0.0;
        scalar eInf = 0.0;

        forAll(T.internalField(), cellI)
        {
            const scalar cellV = V[cellI];

            forAll(cellQP[cellI], qpI)
            {
                const point& gp = cellQP[cellI][qpI];
                const scalar w = cellQW[cellI][qpI];
                const vector d = gp - C[cellI];

                scalar TNum = T[cellI] + (gradT[cellI] & d);

                if (hessTPtr) TNum += 0.5*quadraticForm((*hessTPtr)[cellI], d);
                if (thirdTPtr)
                {
                    TNum += (1.0/6.0)*LRE::cubicForm((*thirdTPtr)[cellI], d, twoD);
                }

                const scalar TAnal = exactTemperature(gp, t, mms);
                const scalar eAbs = mag(TAnal - TNum);

                e1 += cellV*w*eAbs;
                e2 += cellV*w*sqr(eAbs);
                n2 += cellV*w*sqr(TAnal);
                volTot += cellV*w;
                eInf = max(eInf, eAbs);
            }
        }

        GaussErrorSummary out;
        out.L1 = e1/max(volTot, SMALL);
        out.L2 = std::sqrt(e2/max(volTot, SMALL));
        out.Linf = eInf;
        out.normL2 = std::sqrt(n2/max(volTot, SMALL));
        out.relL2 = 100.0*out.L2/max(out.normL2, SMALL);
        return out;
    }

    void writeSummary
    (
        const volScalarField& TError,
        const volScalarField& rhsField,
        const Time& runTime,
        const label nTimeSteps,
        const scalar dt,
        const scalar errorT,
        const scalar gaussL1,
        const scalar gaussL2,
        const scalar gaussLinf,
        const scalar gaussNormL2,
        const scalar gaussRelL2,
        const word& implicitScheme,
        const word& massMatrixType,
        const bool useHighOrderT,
        const label lreN,
        const label lreNn,
        const scalar lreK,
        const label lreMaxStencilSize,
        const scalar setupWallTime,
        const scalar timeLoopWallTime,
        const scalar postProcessWallTime,
        const scalar totalWallTime
    )
    {
        const fvMesh& mesh = TError.mesh();

        const scalar L1_T   = volumeWeightedL1(TError);
        const scalar L2_T   = volumeWeightedL2(TError);
        const scalar Linf_T = linfNorm(TError);

        const scalar Linf_RHS = linfNorm(rhsField);

        const label N = estimatedN(mesh);
        const scalar dx = characteristicDx(mesh);
        const word dimName = name(mesh.nGeometricD()) + "D";

        const fileName outFile =
            runTime.path()
        / (dimName + "_" + name(N) + "_cells_transient.dat");

        OFstream os(outFile);

        os  << "Manufactured-solution error summary (cell-centred):" << nl
            << "Field     L1-error       L2-error       Linf-error" << nl
            << "T      " << L1_T << "   " << L2_T << "   " << Linf_T << nl
            << "-------------------------------------------------" << nl << nl

            << "MATLAB-like cell-centred error:" << nl
            << "Field     error_T" << nl
            << "TC      " << errorT << nl
            << "-------------------------------------------------" << nl << nl

            << "Gauss-reconstructed error summary:" << nl
            << "Field     L1-error       L2-error       Linf-error       normL2          relL2(%)" << nl
            << "TG      " << gaussL1 << "   " << gaussL2 << "   " << gaussLinf
            << "   " << gaussNormL2 << "   " << gaussRelL2 << nl
            << "-------------------------------------------------" << nl << nl

            << "RHS summary at final time:" << nl
            << "Field     Linf-RHS" << nl
            << "R      " << Linf_RHS << nl
            << "-------------------------------------------------" << nl << nl

            << "Simulation summary:" << nl
            << "-------------------" << nl
            << "Final time            = " << runTime.value() << nl
            << "Number of cells (N)   = " << N << nl
            << "Dimension             = " << dimName << nl
            << "Grid spacing (dx)     = " << dx << nl
            << "Time step (dt)        = " << dt << nl
            << "Number of steps       = " << nTimeSteps << nl
            << "Implicit scheme       = " << implicitScheme << nl
            << "Mass matrix           = " << massMatrixType << nl
            << "useHighOrder_T        = " << (useHighOrderT ? "true" : "false") << nl
            << "LRE N                 = " << lreN << nl
            << "LRE Nn                = " << lreNn << nl
            << "LRE k                 = " << lreK << nl
            << "LRE maxStencilSize    = " << lreMaxStencilSize << nl
            << "-------------------" << nl << nl

            << "Wall-clock times [s]:" << nl
            << "setup                 = " << setupWallTime << nl
            << "timeLoop              = " << timeLoopWallTime << nl
            << "postProcess           = " << postProcessWallTime << nl
            << "total                 = " << totalWallTime << nl
            << "-------------------" << nl;

        Info<< "Wrote summary to " << outFile << nl
            << "Cell-centred error_T = " << errorT << nl
            << "T error: L1 = " << L1_T
            << ", L2 = " << L2_T
            << ", Linf = " << Linf_T << nl
            << "Gauss error: L1 = " << gaussL1
            << ", L2 = " << gaussL2
            << ", Linf = " << gaussLinf
            << ", relL2 = " << gaussRelL2 << "%" << nl
            << "Final-time RHS Linf = " << Linf_RHS << nl
            << "Timing [s]: setup = " << setupWallTime
            << ", loop = " << timeLoopWallTime
            << ", post = " << postProcessWallTime
            << ", total = " << totalWallTime << endl;
    }
}

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    const auto tStartTotal = std::chrono::steady_clock::now();
    const auto tStartSetup = tStartTotal;

    MMSConfig mms;
    mms.type  = exactSolutionProperties.lookupOrDefault<word>("mmsType", "example_1");
    mms.A     = exactSolutionProperties.lookupOrDefault<scalar>("amplitude", 1.0);
    mms.beta  = exactSolutionProperties.lookupOrDefault<scalar>("beta", 1.0);
    mms.gamma = exactSolutionProperties.lookupOrDefault<scalar>("gamma", 3.0);
    mms.omega = exactSolutionProperties.lookupOrDefault<scalar>("omega", 4.0);
    mms.dim   = mesh.nGeometricD();

    const scalar alphaXValue = alphaX.value();
    const scalar alphaYValue = alphaY.value();
    const scalar alphaZValue = alphaZ.value();

    const scalar rhoCp = rho.value()*Cp.value();
    const scalar dt = runTime.deltaTValue();

    const scalar theta = thetaFromScheme(implicitScheme);
    const word massMatrixMode = normalizedMassMatrixType(massMatrixType);
    const scalar sourceSign = (mms.type == "example_5") ? -1.0 : 1.0;

    SpMat M;
    SpMat K;
    SpMat AImplicit;
    SpMat BImplicit;

    if (useHighOrder_T)
    {
        if (massMatrixMode == "lumped")
        {
            assembleDiagonalMassMatrix(mesh, rhoCp, M);
        }
        else
        {
            assembleConsistentMassMatrixHO(mesh, LREInterp_T, rhoCp, M);
        }

        assembleHighOrderStiffnessMatrix(mesh, conductivity, LREInterp_T, K);

        AImplicit = M;
        AImplicit *= (1.0/dt);
        AImplicit -= theta*K;

        BImplicit = M;
        BImplicit *= (1.0/dt);

        if (theta < 1.0 - SMALL)
        {
            BImplicit += (1.0 - theta)*K;
        }
    }

    fillExactFields
    (
        TExact,
        sourceQ,
        runTime.value(),
        alphaXValue,
        alphaYValue,
        alphaZValue,
        rhoCp,
        mms,
        LREInterp_T
    );

    T = TExact;
    applyExactBoundaryValues(T, runTime.value(), mms);

    label nTimeSteps = 0;

    const auto tEndSetup = std::chrono::steady_clock::now();
    const auto tStartLoop = tEndSetup;

    while (runTime.value() < runTime.endTime().value() - SMALL)
    {
        const scalar t = runTime.value();

        fillExactFields
        (
            TExact,
            sourceQ,
            t,
            alphaXValue,
            alphaYValue,
            alphaZValue,
            rhoCp,
            mms,
            LREInterp_T
        );

        volScalarField sourceQn("sourceQn", sourceQ);
        sourceQn *= sourceSign;

        fillExactFields
        (
            TExact,
            sourceQ,
            t + dt,
            alphaXValue,
            alphaYValue,
            alphaZValue,
            rhoCp,
            mms,
            LREInterp_T
        );

        volScalarField sourceQnp1("sourceQnp1", sourceQ);
        sourceQnp1 *= sourceSign;

        if (useHighOrder_T)
        {
            const EigVec Tn = fieldToEigVec(T);
            const EigVec qn = sourceToEigVec(sourceQn);

            const EigVec bcn =
                assembleHighOrderBoundaryVector
                (
                    mesh,
                    conductivity,
                    LREInterp_T,
                    t,
                    mms
                );

            const EigVec qnp1 = sourceToEigVec(sourceQnp1);

            const EigVec bcnp1 =
                assembleHighOrderBoundaryVector
                (
                    mesh,
                    conductivity,
                    LREInterp_T,
                    t + dt,
                    mms
                );

            const EigVec rhs =
                BImplicit*Tn
              + theta*(qnp1 + bcnp1)
              + (1.0 - theta)*(qn + bcn);

            const EigVec Tnp1 =
                solveSparseSystem
                (
                    AImplicit,
                    rhs,
                    implicitLinearSolver,
                    implicitTolerance,
                    implicitMaxIterations
                );

            eigVecToField(Tnp1, T);

            applyExactBoundaryValues(T, t + dt, mms);
        }
        else
        {
            applyExactBoundaryValues(T, t, mms);

            const volScalarField Tn("Tn", T);
            const volScalarField lapTn(fvc::laplacian(conductivity, T));

            applyExactBoundaryValues(T, t + dt, mms);

            const dimensionedScalar rhoCpOverDt
            (
                "rhoCpOverDt",
                dimless/sqr(dimLength),
                rhoCp/dt
            );

            fvScalarMatrix TEqn
            (
                fvm::Sp(rhoCpOverDt, T)
              - theta*fvm::laplacian(conductivity, T)
             ==
                rhoCpOverDt*Tn
              + (1.0 - theta)*lapTn
              + theta*sourceQnp1
              + (1.0 - theta)*sourceQn
            );

            TEqn.solve();

            applyExactBoundaryValues(T, t + dt, mms);
        }

        if (useHighOrder_T)
        {
            computeHighOrderLaplacian
            (
                T,
                conductivity,
                LREInterp_T,
                fluxT_HO,
                lapTHO
            );
        }
        else
        {
            lapTHO = fvc::laplacian(conductivity, T);
        }

        if (mms.type == "example_5")
        {
            residualField = lapTHO - sourceQ;
        }
        else
        {
            residualField = lapTHO + sourceQ;
        }

        ++nTimeSteps;
        ++runTime;

        if (nTimeSteps % 50 == 0 || nTimeSteps <= 5)
        {
            Info<< "Time step " << nTimeSteps
                << " : t = " << runTime.value()
                << " , RHS Linf = " << linfNorm(residualField)
                << endl;
        }

        if (runTime.outputTime())
        {
            fillExactFields
            (
                TExact,
                sourceQ,
                runTime.value(),
                alphaXValue,
                alphaYValue,
                alphaZValue,
                rhoCp,
                mms,
                LREInterp_T
            );

            TError = T - TExact;
            operatorError = residualField;
            runTime.write();
        }
    }

    const auto tEndLoop = std::chrono::steady_clock::now();
    const auto tStartPost = tEndLoop;

    fillExactFields
    (
        TExact,
        sourceQ,
        runTime.value(),
        alphaXValue,
        alphaYValue,
        alphaZValue,
        rhoCp,
        mms,
        LREInterp_T
    );

    if (useHighOrder_T)
    {
        computeHighOrderLaplacian
        (
            T,
            conductivity,
            LREInterp_T,
            fluxT_HO,
            lapTHO
        );
    }
    else
    {
        lapTHO = fvc::laplacian(conductivity, T);
    }

    if (mms.type == "example_5")
    {
        residualField = lapTHO - sourceQ;
    }
    else
    {
        residualField = lapTHO + sourceQ;
    }
    TError = T - TExact;
    operatorError = residualField;

    const scalar errorT = errorTCellCentred(TError);

    GaussErrorSummary gaussErr;
    if (useHighOrder_T)
    {
        gaussErr =
            computeGaussReconstructedError(T, runTime.value(), mms, LREInterp_T);
    }
    else
    {
        gaussErr.L1 = volumeWeightedL1(TError);
        gaussErr.L2 = volumeWeightedL2(TError);
        gaussErr.Linf = linfNorm(TError);
        gaussErr.normL2 = volumeWeightedL2(TExact);
        gaussErr.relL2 = 100.0*gaussErr.L2/max(gaussErr.normL2, SMALL);
    }

    const auto tEndPost = std::chrono::steady_clock::now();

    const scalar setupWallTime =
        std::chrono::duration_cast<std::chrono::duration<scalar>>(tEndSetup - tStartSetup).count();

    const scalar timeLoopWallTime =
        std::chrono::duration_cast<std::chrono::duration<scalar>>(tEndLoop - tStartLoop).count();

    const scalar postProcessWallTime =
        std::chrono::duration_cast<std::chrono::duration<scalar>>(tEndPost - tStartPost).count();

    const scalar totalWallTime =
        std::chrono::duration_cast<std::chrono::duration<scalar>>(tEndPost - tStartTotal).count();

    runTime.write();

    writeSummary
    (
        TError,
        residualField,
        runTime,
        nTimeSteps,
        dt,
        errorT,
        gaussErr.L1,
        gaussErr.L2,
        gaussErr.Linf,
        gaussErr.normL2,
        gaussErr.relL2,
        implicitScheme,
        massMatrixMode,
        useHighOrder_T,
        lreN,
        lreNn,
        lreK,
        lreMaxStencil,
        setupWallTime,
        timeLoopWallTime,
        postProcessWallTime,
        totalWallTime
    );

    Info<< "End" << nl << endl;
    return 0;
}
