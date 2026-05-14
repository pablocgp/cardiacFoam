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

    enum ManufacturedFieldID
    {
        mfVm = 0,
        mfU1 = 1,
        mfU2 = 2
    };

    struct FieldErrorSummary
    {
        scalar L1;
        scalar L2;
        scalar Linf;
        scalar errorCell;
    };

    struct ReconstructionErrorSummary
    {
        scalar L1;
        scalar L2;
        scalar Linf;
        scalar normL2;
        scalar relL2;
    };

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

    scalar exactVm(const point& p, const scalar t, const label dim)
    {
        return std::sqrt(1.0 + t)*computeF(p, dim);
    }

    scalar exactU1(const point& p, const scalar t, const label dim)
    {
        return (1.0 + t)*computeG(p, dim) + exactVm(p, t, dim);
    }

    scalar exactU2(const point& p, const scalar t, const label dim)
    {
        return 1.0/((1.0 + t)*std::sqrt(computeG(p, dim)));
    }

    scalar exactU3(const point&, const scalar, const label)
    {
        return 0.0;
    }

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

    scalar vmSplitSourcePDE
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
            0.5*(u1 + u3 - V)*sqr(u2)*(V - u3)
          + (beta/(chiVal*CmVal))*u3;
    }

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
            returnReduce(mesh.nCells(), sumOp<label>());

        return std::pow
        (
            activeMeasure/max(nCellsGlobal, scalar(1.0)),
            1.0/max(scalar(mesh.nGeometricD()), scalar(1.0))
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

    scalar cellErrorNorm(const volScalarField& fld)
    {
        return std::sqrt(gSum(sqr(fld.primitiveField())));
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

    FieldErrorSummary computeFieldErrorSummary(const volScalarField& fld)
    {
        FieldErrorSummary s;
        s.L1 = volumeWeightedL1(fld);
        s.L2 = volumeWeightedL2(fld);
        s.Linf = linfNorm(fld);
        s.errorCell = cellErrorNorm(fld);
        return s;
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

    EigVec sourceToEigVec(const volScalarField& fld)
    {
        return fieldToEigVec(fld);
    }

    void eigVecToField(const EigVec& v, volScalarField& fld)
    {
        scalarField& f = fld.primitiveFieldRef();

        forAll(f, cellI)
        {
            f[cellI] = v[cellI];
        }
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

        const labelListList& stencils = LREInterp.globalCellStencils();
        const CompactListList<point>& cellQP = LREInterp.cellQuadPoints();
        const CompactListList<scalar>& cellQW = LREInterp.cellQuadWeight();
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

                addTripletIfNeeded
                (
                    triplets,
                    cellI,
                    cellI,
                    coefficient*w*selfCoeff
                );
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
            const scalar a = orthogonalDiffusionCoeff(mesh.Sf()[faceI], C[nei] - C[own], Df);

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

            if (bcType == "empty" || bcType == zeroGradientFvPatchScalarField::typeName)
            {
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

                    addTripletIfNeeded(triplets, own, own, -a/max(V[own], SMALL));
                }
            }
        }

        K.resize(mesh.nCells(), mesh.nCells());
        K.setFromTriplets(triplets.begin(), triplets.end());
        K.makeCompressed();
    }

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

        forAll(neighbour, faceI)
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
            const word bcType =
                patch.lookupPatchField<volScalarField, scalar>("Vm").type();

            if (bcType == "empty")
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
                        vector gCoeff =
                            faceGradCoeffs[globalFaceI][qpI][cI];

                        if (bcType == "zeroGradient")
                        {
                            gCoeff -= (gCoeff & n)*n;
                        }

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

    EigVec assembleStandardOrthogonalBoundaryVector
    (
        const fvMesh& mesh,
        const volTensorField& conductivity,
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
            }
        }

        return b;
    }

    EigVec assembleHighOrderBoundaryVector
    (
        const fvMesh& mesh,
        const volTensorField& conductivity,
        const LRE& LREInterp,
        const scalar t,
        const label dim
    )
    {
        EigVec b = EigVec::Zero(mesh.nCells());

        const CompactListList<point>& faceQP = LREInterp.faceQuadPoints();
        const CompactListList<scalar>& faceQW = LREInterp.faceQuadWeight();
        const labelListList& faceStencils = LREInterp.globalFaceStencils();
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

                    const scalar fluxCoeff =
                        area*w*(n & (conductivity[own] & gGhost));

                    b[own] += fluxCoeff/max(V[own], SMALL)*Vbc;
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

            Info<< "Implicit FDA HO solve: iterations = "
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

    void cellQuadraturePointsAndWeights
    (
        const fvMesh& mesh,
        const label cellI,
        List<point>& qPoints,
        scalarField& qWeights
    )
    {
        static const scalar xi[4] =
        {
            -0.8611363115940526,
            -0.3399810435848563,
             0.3399810435848563,
             0.8611363115940526
        };

        static const scalar wi[4] =
        {
            0.1739274225687269,
            0.3260725774312731,
            0.3260725774312731,
            0.1739274225687269
        };

        const label dim = mesh.nGeometricD();
        const label n1D = 4;

        label nQ = 1;
        for (label d = 0; d < dim; ++d)
        {
            nQ *= n1D;
        }

        qPoints.setSize(nQ);
        qWeights.setSize(nQ);

        point pMin(GREAT, GREAT, GREAT);
        point pMax(-GREAT, -GREAT, -GREAT);

        const labelList& cellPts = mesh.cellPoints()[cellI];
        const pointField& pts = mesh.points();

        forAll(cellPts, i)
        {
            const point& p = pts[cellPts[i]];

            pMin.x() = min(pMin.x(), p.x());
            pMin.y() = min(pMin.y(), p.y());
            pMin.z() = min(pMin.z(), p.z());

            pMax.x() = max(pMax.x(), p.x());
            pMax.y() = max(pMax.y(), p.y());
            pMax.z() = max(pMax.z(), p.z());
        }

        const point& c = mesh.C()[cellI];

        for (label qI = 0; qI < nQ; ++qI)
        {
            label idx = qI;
            point p = c;
            scalar w = 1.0;

            for (label d = 0; d < dim; ++d)
            {
                const label i = idx % n1D;
                idx /= n1D;

                const scalar lo = pMin[d];
                const scalar hi = pMax[d];
                const scalar mid = 0.5*(lo + hi);
                const scalar half = 0.5*(hi - lo);

                p[d] = mid + half*xi[i];
                w *= wi[i];
            }

            qPoints[qI] = p;
            qWeights[qI] = w;
        }
    }

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

        autoPtr<List<List<vector>>> gradVmQuadPtr =
            LREInterp_Vm.gradScalarFaceQuad(Vm);

        List<List<vector>>& gradVmQuad = gradVmQuadPtr.ref();
        const CompactListList<scalar>& faceQuadW =
            LREInterp_Vm.faceQuadWeight();

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

                fluxVm_HO[faceI] +=
                    faceArea*(faceNormal & Dg)*faceQuadW[faceI][pI];
            }
        }

        forAll(fluxVm_HO.boundaryField(), patchI)
        {
            scalarField& patchFlux = fluxVm_HO.boundaryFieldRef()[patchI];

            if (patchFlux.size() == 0)
            {
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

                    patchFlux[faceI] +=
                        faceArea*(faceNormal & Dg)*faceQuadW[globalFaceI][pI];
                }
            }
        }

        lapVm = fvc::div(fluxVm_HO);
    }

    void computeVolumeAveragedIion
    (
        const volScalarField& Vm,
        const volScalarField& u1,
        const volScalarField& u2,
        const volScalarField& u3,
        const scalar beta,
        const scalar chiVal,
        const scalar CmVal,
        const bool useHighOrderVm,
        const bool useHighOrderStates,
        const LRE& LREInterp_Vm,
        const LRE& LREInterp_states,
        volScalarField& Iion
    )
    {
        const fvMesh& mesh = Vm.mesh();

        if (!useHighOrderStates || mesh.nGeometricD() == 1)
        {
            forAll(Vm.internalField(), cellI)
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
            return;
        }

        const bool twoD = mesh.nGeometricD() == 2;

        tmp<volVectorField> tGradVm =
                useHighOrderVm
            ? LREInterp_Vm.grad(Vm)
            : fvc::grad(Vm);

        const volVectorField& gradVm = tGradVm();


        tmp<volVectorField> tGradU1 = LREInterp_states.grad(u1);
        const volVectorField& gradU1 = tGradU1();

        tmp<volVectorField> tGradU2 = LREInterp_states.grad(u2);
        const volVectorField& gradU2 = tGradU2();

        tmp<volSymmTensorField> tHessVm;
        const volSymmTensorField* hessVmPtr = nullptr;
        if (useHighOrderVm && LREInterp_Vm.order() >= 2)
        {
            tHessVm = LREInterp_Vm.hessian(Vm);
            hessVmPtr = &tHessVm();
        }

        tmp<volSymmTensorField> tHessU1;
        const volSymmTensorField* hessU1Ptr = nullptr;
        if (LREInterp_states.order() >= 2)
        {
            tHessU1 = LREInterp_states.hessian(u1);
            hessU1Ptr = &tHessU1();
        }

        tmp<volSymmTensorField> tHessU2;
        const volSymmTensorField* hessU2Ptr = nullptr;
        if (LREInterp_states.order() >= 2)
        {
            tHessU2 = LREInterp_states.hessian(u2);
            hessU2Ptr = &tHessU2();
        }

        autoPtr<List<LRE::symmTensor3Order>> thirdVmPtr;
        const List<LRE::symmTensor3Order>* thirdVmList = nullptr;
        if (useHighOrderVm && LREInterp_Vm.order() >= 3)
        {
            thirdVmPtr = LREInterp_Vm.thirdDeriv(Vm);
            thirdVmList = &thirdVmPtr();
        }

        autoPtr<List<LRE::symmTensor3Order>> thirdU1Ptr;
        const List<LRE::symmTensor3Order>* thirdU1List = nullptr;
        if (LREInterp_states.order() >= 3)
        {
            thirdU1Ptr = LREInterp_states.thirdDeriv(u1);
            thirdU1List = &thirdU1Ptr();
        }

        autoPtr<List<LRE::symmTensor3Order>> thirdU2Ptr;
        const List<LRE::symmTensor3Order>* thirdU2List = nullptr;
        if (LREInterp_states.order() >= 3)
        {
            thirdU2Ptr = LREInterp_states.thirdDeriv(u2);
            thirdU2List = &thirdU2Ptr();
        }

        const vectorField& C = mesh.C();

        List<point> qPoints;
        scalarField qWeights;

        forAll(Vm.internalField(), cellI)
        {
            cellQuadraturePointsAndWeights(mesh, cellI, qPoints, qWeights);

            scalar iBar = 0.0;
            scalar wSum = 0.0;

            forAll(qPoints, qpI)
            {
                const vector d = qPoints[qpI] - C[cellI];
                const scalar w = qWeights[qpI];

                const symmTensor* HVm =
                    hessVmPtr ? &((*hessVmPtr)[cellI]) : nullptr;
                const symmTensor* HU1 =
                    hessU1Ptr ? &((*hessU1Ptr)[cellI]) : nullptr;
                const symmTensor* HU2 =
                    hessU2Ptr ? &((*hessU2Ptr)[cellI]) : nullptr;

                const LRE::symmTensor3Order* TVm =
                    thirdVmList ? &((*thirdVmList)[cellI]) : nullptr;
                const LRE::symmTensor3Order* TU1 =
                    thirdU1List ? &((*thirdU1List)[cellI]) : nullptr;
                const LRE::symmTensor3Order* TU2 =
                    thirdU2List ? &((*thirdU2List)[cellI]) : nullptr;

                const scalar Vg =
                    reconstructFromTaylor
                    (
                        Vm[cellI],
                        gradVm[cellI],
                        HVm,
                        TVm,
                        d,
                        twoD
                    );

                const scalar u1g =
                    reconstructFromTaylor
                    (
                        u1[cellI],
                        gradU1[cellI],
                        HU1,
                        TU1,
                        d,
                        twoD
                    );

                const scalar u2g =
                    reconstructFromTaylor
                    (
                        u2[cellI],
                        gradU2[cellI],
                        HU2,
                        TU2,
                        d,
                        twoD
                    );

                iBar +=
                    w
                   *ionicCurrentPDE
                    (
                        Vg,
                        u1g,
                        u2g,
                        0.0,
                        beta,
                        chiVal,
                        CmVal
                    );

                wSum += w;
            }

            Iion[cellI] = iBar/max(wSum, SMALL);
        }

        Iion.correctBoundaryConditions();
    }

    void computeVolumeAveragedVmSource
    (
        const volScalarField& Vm,
        const volScalarField& u1,
        const volScalarField& u2,
        const volScalarField& u3,
        const scalar beta,
        const scalar chiVal,
        const scalar CmVal,
        const bool useHighOrderVm,
        const bool useHighOrderStates,
        const LRE& LREInterp_Vm,
        const LRE& LREInterp_states,
        volScalarField& sourceVm
    )
    {
        const fvMesh& mesh = Vm.mesh();

        if (!useHighOrderStates || mesh.nGeometricD() == 1)
        {
            forAll(Vm.internalField(), cellI)
            {
                sourceVm[cellI] =
                    vmSplitSourcePDE
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
            return;
        }

        const bool twoD = mesh.nGeometricD() == 2;

        tmp<volVectorField> tGradVm =
                useHighOrderVm
            ? LREInterp_Vm.grad(Vm)
            : fvc::grad(Vm);

        const volVectorField& gradVm = tGradVm();

        tmp<volVectorField> tGradU1 = LREInterp_states.grad(u1);
        const volVectorField& gradU1 = tGradU1();

        tmp<volVectorField> tGradU2 = LREInterp_states.grad(u2);
        const volVectorField& gradU2 = tGradU2();

        tmp<volSymmTensorField> tHessVm;
        const volSymmTensorField* hessVmPtr = nullptr;
        if (useHighOrderVm && LREInterp_Vm.order() >= 2)
        {
            tHessVm = LREInterp_Vm.hessian(Vm);
            hessVmPtr = &tHessVm();
        }

        tmp<volSymmTensorField> tHessU1;
        const volSymmTensorField* hessU1Ptr = nullptr;
        if (LREInterp_states.order() >= 2)
        {
            tHessU1 = LREInterp_states.hessian(u1);
            hessU1Ptr = &tHessU1();
        }

        tmp<volSymmTensorField> tHessU2;
        const volSymmTensorField* hessU2Ptr = nullptr;
        if (LREInterp_states.order() >= 2)
        {
            tHessU2 = LREInterp_states.hessian(u2);
            hessU2Ptr = &tHessU2();
        }

        autoPtr<List<LRE::symmTensor3Order>> thirdVmPtr;
        const List<LRE::symmTensor3Order>* thirdVmList = nullptr;
        if (useHighOrderVm && LREInterp_Vm.order() >= 3)
        {
            thirdVmPtr = LREInterp_Vm.thirdDeriv(Vm);
            thirdVmList = &thirdVmPtr();
        }

        autoPtr<List<LRE::symmTensor3Order>> thirdU1Ptr;
        const List<LRE::symmTensor3Order>* thirdU1List = nullptr;
        if (LREInterp_states.order() >= 3)
        {
            thirdU1Ptr = LREInterp_states.thirdDeriv(u1);
            thirdU1List = &thirdU1Ptr();
        }

        autoPtr<List<LRE::symmTensor3Order>> thirdU2Ptr;
        const List<LRE::symmTensor3Order>* thirdU2List = nullptr;
        if (LREInterp_states.order() >= 3)
        {
            thirdU2Ptr = LREInterp_states.thirdDeriv(u2);
            thirdU2List = &thirdU2Ptr();
        }

        const vectorField& C = mesh.C();

        List<point> qPoints;
        scalarField qWeights;

        forAll(Vm.internalField(), cellI)
        {
            cellQuadraturePointsAndWeights(mesh, cellI, qPoints, qWeights);

            scalar sourceBar = 0.0;
            scalar wSum = 0.0;

            forAll(qPoints, qpI)
            {
                const vector d = qPoints[qpI] - C[cellI];
                const scalar w = qWeights[qpI];

                const symmTensor* HVm =
                    hessVmPtr ? &((*hessVmPtr)[cellI]) : nullptr;
                const symmTensor* HU1 =
                    hessU1Ptr ? &((*hessU1Ptr)[cellI]) : nullptr;
                const symmTensor* HU2 =
                    hessU2Ptr ? &((*hessU2Ptr)[cellI]) : nullptr;

                const LRE::symmTensor3Order* TVm =
                    thirdVmList ? &((*thirdVmList)[cellI]) : nullptr;
                const LRE::symmTensor3Order* TU1 =
                    thirdU1List ? &((*thirdU1List)[cellI]) : nullptr;
                const LRE::symmTensor3Order* TU2 =
                    thirdU2List ? &((*thirdU2List)[cellI]) : nullptr;

                const scalar Vg =
                    reconstructFromTaylor
                    (
                        Vm[cellI],
                        gradVm[cellI],
                        HVm,
                        TVm,
                        d,
                        twoD
                    );

                const scalar u1g =
                    reconstructFromTaylor
                    (
                        u1[cellI],
                        gradU1[cellI],
                        HU1,
                        TU1,
                        d,
                        twoD
                    );

                const scalar u2g =
                    reconstructFromTaylor
                    (
                        u2[cellI],
                        gradU2[cellI],
                        HU2,
                        TU2,
                        d,
                        twoD
                    );

                sourceBar +=
                    w
                   *vmSplitSourcePDE
                    (
                        Vg,
                        u1g,
                        u2g,
                        u3[cellI],
                        beta,
                        chiVal,
                        CmVal
                    );

                wSum += w;
            }

            sourceVm[cellI] = sourceBar/max(wSum, SMALL);
        }

        sourceVm.correctBoundaryConditions();
    }

    void updateStateFieldsTheta
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
        const scalar theta,
        const label nIterations
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

        u1NI = u1O;
        u2NI = u2O;
        u3NI = u3O;

        for (label iter = 0; iter < max(nIterations, label(1)); ++iter)
        {
            forAll(u1NI, cellI)
            {
                scalar du1Old = 0.0;
                scalar du2Old = 0.0;
                scalar du3Old = 0.0;

                reactionRates
                (
                    VO[cellI],
                    u1O[cellI],
                    u2O[cellI],
                    u3O[cellI],
                    du1Old,
                    du2Old,
                    du3Old
                );

                scalar du1New = 0.0;
                scalar du2New = 0.0;
                scalar du3New = 0.0;

                reactionRates
                (
                    VN[cellI],
                    u1NI[cellI],
                    u2NI[cellI],
                    u3NI[cellI],
                    du1New,
                    du2New,
                    du3New
                );

                u1NI[cellI] =
                    u1O[cellI]
                  + dt*((1.0 - theta)*du1Old + theta*du1New);

                u2NI[cellI] =
                    u2O[cellI]
                  + dt*((1.0 - theta)*du2Old + theta*du2New);

                u3NI[cellI] =
                    u3O[cellI]
                  + dt*((1.0 - theta)*du3Old + theta*du3New);
            }
        }
    }

    ReconstructionErrorSummary computeGaussReconstructedError
    (
        const volScalarField& fld,
        const scalar t,
        const ManufacturedFieldID fieldID,
        const label dim,
        const LRE& LREInterp
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

        autoPtr<List<LRE::symmTensor3Order>> thirdPtr;
        const List<LRE::symmTensor3Order>* thirdList = nullptr;
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
                    const LRE::symmTensor3Order* T3 =
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

    void writeSummary
    (
        const Time& runTime,
        const scalar beta,
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
        const label lreN,
        const label lreNn,
        const scalar lreK,
        const label lreMaxStencilSize,
        const label nonlinearIterations,
        const label stateIterations,
        const scalar setupWallTime,
        const scalar timeLoopWallTime,
        const scalar postProcessWallTime,
        const scalar totalWallTime
    )
    {
        const label N = estimatedN(mesh);
        const scalar dx = characteristicDx(mesh);
        const word dimName = name(mesh.nGeometricD()) + "D";

        const fileName outFile =
            runTime.path()
          / (dimName + "_" + name(N) + "_cells_transient.dat");

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
            << "LRE N                 = " << lreN << nl
            << "LRE Nn                = " << lreNn << nl
            << "LRE k                 = " << lreK << nl
            << "LRE maxStencilSize    = " << lreMaxStencilSize << nl
            << "Nonlinear iterations  = " << nonlinearIterations << nl
            << "State iterations      = " << stateIterations << nl
            << "beta                  = " << beta << nl
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
    #include "setRootCaseLists.H"
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
    const scalar betaScale = beta*lapScale;
    const scalar theta = thetaFromScheme(implicitScheme);
    const word massMatrixMode = normalizedMassMatrixType(massMatrixType);

    Info<< "Running high-order manufactured FDA implicit solver" << nl
        << "Dimension = " << dim << nl
        << "beta = " << beta << nl
        << "dt = " << dt << nl
        << "implicitScheme = " << implicitScheme << nl
        << "massMatrix = " << massMatrixMode << endl;

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

    SpMat M;
    SpMat K;
    SpMat L;
    SpMat AImplicit;
    SpMat BImplicit;

    if (!useHighOrder_Vm || massMatrixMode == "lumped")
    {
        assembleDiagonalMassMatrix(mesh, 1.0, M);
    }
    else
    {
        assembleConsistentMassMatrixHO(mesh, LREInterp_Vm, 1.0, M);
    }

    if (useHighOrder_Vm)
    {
        assembleHighOrderStiffnessMatrix(mesh, conductivity, LREInterp_Vm, K);
    }
    else
    {
        assembleStandardOrthogonalStiffnessMatrix(mesh, conductivity, K);
    }

    L = K;
    L *= lapScale;
    SpMat betaReaction = M;
    betaReaction *= betaScale;
    L -= betaReaction;

    AImplicit = M;
    AImplicit *= (1.0/dt);
    AImplicit -= theta*L;

    BImplicit = M;
    BImplicit *= (1.0/dt);

    if (theta < 1.0 - SMALL)
    {
        BImplicit += (1.0 - theta)*L;
    }

    const auto tEndSetup = std::chrono::steady_clock::now();
    const auto tStartLoop = tEndSetup;

    label nSteps = 0;

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

        computeVolumeAveragedVmSource
        (
            VmOld,
            u1Old,
            u2Old,
            u3Old,
            beta,
            chiVal,
            CmVal,
            useHighOrder_Vm,
            useHighOrder_states,
            LREInterp_Vm,
            LREInterp_states,
            sourceVm
        );

        const EigVec Vn = fieldToEigVec(VmOld);
        const EigVec sourceN = sourceToEigVec(sourceVm);

        EigVec bcN =
                useHighOrder_Vm
            ? assembleHighOrderBoundaryVector(mesh, conductivity, LREInterp_Vm, t, dim)
            : assembleStandardOrthogonalBoundaryVector(mesh, conductivity, t, dim);

            EigVec bcNp1 =
                useHighOrder_Vm
            ? assembleHighOrderBoundaryVector(mesh, conductivity, LREInterp_Vm, t + dt, dim)
            : assembleStandardOrthogonalBoundaryVector(mesh, conductivity, t + dt, dim);


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

        for
        (
            label corr = 0;
            corr < max(implicitNonlinearIterations, label(1));
            ++corr
        )
        {
            computeVolumeAveragedVmSource
            (
                VmGuess,
                u1Guess,
                u2Guess,
                u3Guess,
                beta,
                chiVal,
                CmVal,
                useHighOrder_Vm,
                useHighOrder_states,
                LREInterp_Vm,
                LREInterp_states,
                sourceVm
            );

            const EigVec sourceNp1 = sourceToEigVec(sourceVm);

            const EigVec rhs =
                BImplicit*Vn
              + theta*(sourceNp1 + bcNp1)
              + (1.0 - theta)*(sourceN + bcN);

            const EigVec Vnp1 =
                solveSparseSystem
                (
                    AImplicit,
                    rhs,
                    implicitLinearSolver,
                    implicitTolerance,
                    implicitMaxIterations
                );

            eigVecToField(Vnp1, VmGuess);
            applyExactVmBoundaryValues(VmGuess, t + dt, dim);

            updateStateFieldsTheta
            (
                VmOld,
                VmGuess,
                u1Old,
                u2Old,
                u3Old,
                u1Guess,
                u2Guess,
                u3Guess,
                dt,
                theta,
                stateImplicitIterations
            );

            updateStateBoundaryValues(u1Guess, u2Guess, u3Guess, t + dt, dim);
            u1Guess.correctBoundaryConditions();
            u2Guess.correctBoundaryConditions();
            u3Guess.correctBoundaryConditions();
        }

        Vm.primitiveFieldRef() = VmGuess.primitiveField();
        u1.primitiveFieldRef() = u1Guess.primitiveField();
        u2.primitiveFieldRef() = u2Guess.primitiveField();
        u3.primitiveFieldRef() = u3Guess.primitiveField();

        applyExactVmBoundaryValues(Vm, t + dt, dim);
        updateStateBoundaryValues(u1, u2, u3, t + dt, dim);
        u1.correctBoundaryConditions();
        u2.correctBoundaryConditions();
        u3.correctBoundaryConditions();

        computeVolumeAveragedIion
        (
            Vm,
            u1,
            u2,
            u3,
            beta,
            chiVal,
            CmVal,
            useHighOrder_Vm,
            useHighOrder_states,
            LREInterp_Vm,
            LREInterp_states,
            Iion
        );

        if (useHighOrder_Vm)
        {
            computeHighOrderLaplacian(Vm, conductivity, LREInterp_Vm, fluxVm_HO, lapVm);
        }
        else
        {
            lapVm = fvc::laplacian(conductivity, Vm);
        }

        rhsVm = lapVm/(chi*Cm) - Iion;

        ++nSteps;
        ++runTime;

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

    computeVolumeAveragedIion
    (
        Vm,
        u1,
        u2,
        u3,
        beta,
        chiVal,
        CmVal,
        useHighOrder_Vm,
        useHighOrder_states,
        LREInterp_Vm,
        LREInterp_states,
        Iion
    );

    if (useHighOrder_Vm)
    {
        computeHighOrderLaplacian(Vm, conductivity, LREInterp_Vm, fluxVm_HO, lapVm);
    }
    else
    {
        lapVm = fvc::laplacian(conductivity, Vm);
    }

    rhsVm = lapVm/(chi*Cm) - Iion;

    VmError = Vm - VmExact;
    u1Error = u1 - u1Exact;
    u2Error = u2 - u2Exact;

    const FieldErrorSummary VmCell = computeFieldErrorSummary(VmError);
    const FieldErrorSummary u1Cell = computeFieldErrorSummary(u1Error);
    const FieldErrorSummary u2Cell = computeFieldErrorSummary(u2Error);

    const ReconstructionErrorSummary VmHO =
            useHighOrder_Vm
        ? computeGaussReconstructedError(Vm, runTime.value(), mfVm, dim, LREInterp_Vm)
        : cellAsReconstructionError(VmCell, VmExact);

    const ReconstructionErrorSummary u1HO =
            useHighOrder_states
        ? computeGaussReconstructedError(u1, runTime.value(), mfU1, dim, LREInterp_states)
        : cellAsReconstructionError(u1Cell, u1Exact);

    const ReconstructionErrorSummary u2HO =
            useHighOrder_states
        ? computeGaussReconstructedError(u2, runTime.value(), mfU2, dim, LREInterp_states)
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

    runTime.write();

    writeSummary
    (
        runTime,
        beta,
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
        lreN,
        lreNn,
        lreK,
        lreMaxStencil,
        implicitNonlinearIterations,
        stateImplicitIterations,
        setupWallTime,
        timeLoopWallTime,
        postProcessWallTime,
        totalWallTime
    );

    Info<< "End" << nl << endl;
    return 0;
}
