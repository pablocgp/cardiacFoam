#include "fvCFD.H"
#include "LRE.H"
#include <cmath>

namespace
{
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

        tmp<volVectorField> tGradVm = LREInterp_Vm.grad(Vm);
        const volVectorField& gradVm = tGradVm();

        tmp<volVectorField> tGradU1 = LREInterp_states.grad(u1);
        const volVectorField& gradU1 = tGradU1();

        tmp<volVectorField> tGradU2 = LREInterp_states.grad(u2);
        const volVectorField& gradU2 = tGradU2();

        tmp<volSymmTensorField> tHessVm;
        const volSymmTensorField* hessVmPtr = nullptr;
        if (LREInterp_Vm.order() >= 2)
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
        if (LREInterp_Vm.order() >= 3)
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
        const fvMesh& mesh
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
            << "beta                  = " << beta << nl
            << "-------------------" << nl;

        Info<< "Wrote summary to " << outFile << nl
            << "Vm error: L1 = " << VmCell.L1
            << ", L2 = " << VmCell.L2
            << ", Linf = " << VmCell.Linf << nl
            << "Vm Gauss error: L1 = " << VmHO.L1
            << ", L2 = " << VmHO.L2
            << ", Linf = " << VmHO.Linf
            << ", Relative = " << VmHO.relL2 << "%" << nl
            << "Final-time RHS Linf = " << linfNorm(rhsVm) << endl;
    }
}

int main(int argc, char* argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    const label dim = mesh.nGeometricD();
    const scalar dt = runTime.deltaTValue();
    const scalar beta = computeBeta(conductivity, dim);
    const scalar chiVal = chi.value();
    const scalar CmVal = Cm.value();

    Info<< "Running high-order manufactured FDA explicit solver" << nl
        << "Dimension = " << dim << nl
        << "beta = " << beta << nl
        << "dt = " << dt << endl;

    if (dim == 1 && (useHighOrder_Vm || useHighOrder_states))
    {
        WarningInFunction
            << "1D LRE face/cell quadrature is not implemented in the current "
            << "solids4foam LRE class. The 1D run uses the standard explicit "
            << "FV operator and cell-centred ionic current." << endl;
    }

    fillExactFields(VmExact, u1Exact, u2Exact, runTime.value(), dim);

    Vm.primitiveFieldRef() = VmExact.primitiveField();
    u1.primitiveFieldRef() = u1Exact.primitiveField();
    u2.primitiveFieldRef() = u2Exact.primitiveField();
    u3 = dimensionedScalar("zero", dimless, 0.0);

    Vm.correctBoundaryConditions();
    updateStateBoundaryValues(u1, u2, u3, runTime.value(), dim);
    u1.correctBoundaryConditions();
    u2.correctBoundaryConditions();
    u3.correctBoundaryConditions();

    label nSteps = 0;

    while (runTime.value() < runTime.endTime().value() - SMALL)
    {
        const scalar t = runTime.value();

        computeVolumeAveragedIion
        (
            Vm,
            u1,
            u2,
            u3,
            beta,
            chiVal,
            CmVal,
            useHighOrder_states,
            LREInterp_Vm,
            LREInterp_states,
            Iion
        );

        if (useHighOrder_Vm && dim > 1)
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
        else
        {
            lapVm = fvc::laplacian(conductivity, Vm);
        }

        rhsVm = lapVm/(chi*Cm) - Iion;

        scalarField& VmI = Vm.primitiveFieldRef();
        scalarField& u1I = u1.primitiveFieldRef();
        scalarField& u2I = u2.primitiveFieldRef();
        scalarField& u3I = u3.primitiveFieldRef();

        forAll(VmI, cellI)
        {
            scalar du1dt = 0.0;
            scalar du2dt = 0.0;
            scalar du3dt = 0.0;

            reactionRates
            (
                VmI[cellI],
                u1I[cellI],
                u2I[cellI],
                u3I[cellI],
                du1dt,
                du2dt,
                du3dt
            );

            VmI[cellI] += dt*rhsVm[cellI];
            u1I[cellI] += dt*du1dt;
            u2I[cellI] += dt*du2dt;
            u3I[cellI] += dt*du3dt;
        }

        Vm.correctBoundaryConditions();
        updateStateBoundaryValues(u1, u2, u3, t + dt, dim);
        u1.correctBoundaryConditions();
        u2.correctBoundaryConditions();
        u3.correctBoundaryConditions();

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
        useHighOrder_states,
        LREInterp_Vm,
        LREInterp_states,
        Iion
    );

    if (useHighOrder_Vm && dim > 1)
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

    const ReconsstructionErrorSummary u2HO =
            useHighOrder_states
        ? computeGaussReconstructedError(u2, runTime.value(), mfU2, dim, LREInterp_states)
        : cellAsReconstructionError(u2Cell, u2Exact);

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
        mesh
    );

    Info<< "End" << nl << endl;
    return 0;
}
