#include "fvCFD.H"
#include "LRE.H"
#include <cmath>
#include <chrono>

namespace
{
    enum ManufacturedFieldID
    {
        mfVm = 0,
        mfU1 = 1,
        mfU2 = 2
    };

    struct CellErrorSummary
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

    struct TimingSummary
    {
        scalar setupWallTime;
        scalar timeLoopWallTime;
        scalar postProcessWallTime;
        scalar solverWallTime;
    };

    struct QuadratureState
    {
        List<scalarField> u1;
        List<scalarField> u2;
        List<scalarField> u3;
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
        else
        {
            return
                std::cos(pi*p.x())
               *std::cos(2.0*pi*p.y())
               *std::cos(3.0*pi*p.z());
        }
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
        else
        {
            return 1.0 + p.x()*sqr(p.y())*pow3(p.z());
        }
    }

    scalar exactVm(const point& p, const scalar t, const label dim)
    {
        return std::sqrt(1.0 + t)*computeF(p, dim);
    }

    scalar exactU1(const point& p, const scalar t, const label dim)
    {
        return (1.0 + t)*computeG(p, dim) + std::sqrt(1.0 + t)*computeF(p, dim);
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
        else
        {
            return -pi2*(D.xx() + 4.0*D.yy() + 9.0*D.zz());
        }
    }

    label estimatedN(const fvMesh& mesh)
    {
        const label nCellsGlobal =
            returnReduce(mesh.nCells(), sumOp<label>());

        return label(std::round(std::cbrt(static_cast<double>(nCellsGlobal))));
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
        const scalarField& f = fld.primitiveField();
        return std::sqrt(gSum(sqr(f)));
    }

    CellErrorSummary computeCellErrorSummary(const volScalarField& fld)
    {
        CellErrorSummary s;
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

    List<scalarField> makeQuadratureStorage
    (
        const CompactListList<point>& cellQP,
        const scalar initValue
    )
    {
        List<scalarField> storage(cellQP.size());

        forAll(cellQP, cellI)
        {
            storage[cellI].setSize(cellQP[cellI].size());

            forAll(storage[cellI], qpI)
            {
                storage[cellI][qpI] = initValue;
            }
        }

        return storage;
    }

    QuadratureState makeQuadratureStateStorage
    (
        const CompactListList<point>& cellQP
    )
    {
        QuadratureState qState;
        qState.u1 = makeQuadratureStorage(cellQP, 0.0);
        qState.u2 = makeQuadratureStorage(cellQP, 0.0);
        qState.u3 = makeQuadratureStorage(cellQP, 0.0);
        return qState;
    }

    void initializeQuadratureStatesExact
    (
        QuadratureState& qState,
        const CompactListList<point>& cellQP,
        const scalar t,
        const label dim
    )
    {
        forAll(cellQP, cellI)
        {
            forAll(cellQP[cellI], qpI)
            {
                const point& gp = cellQP[cellI][qpI];

                qState.u1[cellI][qpI] = exactU1(gp, t, dim);
                qState.u2[cellI][qpI] = exactU2(gp, t, dim);
                qState.u3[cellI][qpI] = exactU3(gp, t, dim);
            }
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
            const vectorField& CfPatch = Cf.boundaryField()[patchI];

            forAll(u1p, faceI)
            {
                u1p[faceI] = exactU1(CfPatch[faceI], t, dim);
                u2p[faceI] = exactU2(CfPatch[faceI], t, dim);
                u3p[faceI] = exactU3(CfPatch[faceI], t, dim);
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

        const fvMesh& mesh = u1Exact.mesh();
        const surfaceVectorField& Cf = mesh.Cf();

        forAll(u1Exact.boundaryField(), patchI)
        {
            fvPatchScalarField& u1p = u1Exact.boundaryFieldRef()[patchI];
            fvPatchScalarField& u2p = u2Exact.boundaryFieldRef()[patchI];
            const vectorField& CfPatch = Cf.boundaryField()[patchI];

            forAll(u1p, faceI)
            {
                u1p[faceI] = exactU1(CfPatch[faceI], t, dim);
                u2p[faceI] = exactU2(CfPatch[faceI], t, dim);
            }
        }

        VmExact.correctBoundaryConditions();
        u1Exact.correctBoundaryConditions();
        u2Exact.correctBoundaryConditions();
    }

    void projectQuadratureStatesToCells
    (
        const QuadratureState& qState,
        const CompactListList<scalar>& cellQW,
        volScalarField& u1,
        volScalarField& u2,
        volScalarField& u3
    )
    {
        scalarField& u1I = u1.primitiveFieldRef();
        scalarField& u2I = u2.primitiveFieldRef();
        scalarField& u3I = u3.primitiveFieldRef();

        forAll(qState.u1, cellI)
        {
            scalar u1Bar = 0.0;
            scalar u2Bar = 0.0;
            scalar u3Bar = 0.0;
            scalar wSum = 0.0;

            forAll(qState.u1[cellI], qpI)
            {
                const scalar w = cellQW[cellI][qpI];

                u1Bar += w*qState.u1[cellI][qpI];
                u2Bar += w*qState.u2[cellI][qpI];
                u3Bar += w*qState.u3[cellI][qpI];
                wSum += w;
            }

            u1I[cellI] = u1Bar/max(wSum, SMALL);
            u2I[cellI] = u2Bar/max(wSum, SMALL);
            u3I[cellI] = u3Bar/max(wSum, SMALL);
        }
    }

    void reconstructScalarAtCellQuadrature
    (
        const volScalarField& fld,
        const LRE& LREInterp,
        const CompactListList<point>& cellQP,
        List<scalarField>& fldQP
    )
    {
        const fvMesh& mesh = fld.mesh();
        const bool twoD = mesh.nGeometricD() == 2;
        const vectorField& C = mesh.C();

        if (fldQP.size() != cellQP.size())
        {
            fldQP.setSize(cellQP.size());
        }

        tmp<volVectorField> tGrad = LREInterp.grad(fld);
        const volVectorField& grad = tGrad();

        tmp<volSymmTensorField> tHess;
        const volSymmTensorField* hessPtr = nullptr;
        if (LREInterp.order() >= 2)
        {
            tHess = LREInterp.hessian(fld);
            hessPtr = &tHess();
        }

        autoPtr<List<LRE::symmTensor3Order>> thirdPtr;
        const List<LRE::symmTensor3Order>* thirdList = nullptr;
        if (LREInterp.order() >= 3)
        {
            thirdPtr = LREInterp.thirdDeriv(fld);
            thirdList = &thirdPtr();
        }

        forAll(cellQP, cellI)
        {
            if (fldQP[cellI].size() != cellQP[cellI].size())
            {
                fldQP[cellI].setSize(cellQP[cellI].size());
            }

            const symmTensor* H =
                hessPtr ? &((*hessPtr)[cellI]) : nullptr;
            const LRE::symmTensor3Order* T3 =
                thirdList ? &((*thirdList)[cellI]) : nullptr;

            forAll(cellQP[cellI], qpI)
            {
                const vector d = cellQP[cellI][qpI] - C[cellI];

                fldQP[cellI][qpI] =
                    reconstructFromTaylor
                    (
                        fld[cellI],
                        grad[cellI],
                        H,
                        T3,
                        d,
                        twoD
                    );
            }
        }
    }

    void computeHighOrderLaplacian
    (
        const volScalarField& Vm,
        const volTensorField& conductivity,
        const LRE& LREInterp_Vm,
        volScalarField& lapVm
    )
    {
        const fvMesh& mesh = Vm.mesh();

        surfaceScalarField fluxVm_HO
        (
            IOobject
            (
                "fluxVm_HO",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", lapVm.dimensions()*dimVolume, 0.0)
        );

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

            forAll(patchFlux, faceI)
            {
                patchFlux[faceI] = 0.0;
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
        if (!useHighOrderStates)
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

        const fvMesh& mesh = Vm.mesh();
        const bool twoD = mesh.nGeometricD() == 2;

        const CompactListList<point>& cellQP = LREInterp_states.cellQuadPoints();
        const CompactListList<scalar>& cellQW = LREInterp_states.cellQuadWeight();

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

        forAll(Vm.internalField(), cellI)
        {
            scalar iBar = 0.0;
            scalar wSum = 0.0;

            forAll(cellQP[cellI], qpI)
            {
                const point& gp = cellQP[cellI][qpI];
                const vector d = gp - C[cellI];
                const scalar w = cellQW[cellI][qpI];

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

    void computeVolumeAveragedIionFromQuadratureStates
    (
        const List<scalarField>& VmQP,
        const QuadratureState& qState,
        const CompactListList<scalar>& cellQW,
        const scalar beta,
        const scalar chiVal,
        const scalar CmVal,
        volScalarField& Iion
    )
    {
        forAll(VmQP, cellI)
        {
            scalar iBar = 0.0;
            scalar wSum = 0.0;

            forAll(VmQP[cellI], qpI)
            {
                const scalar w = cellQW[cellI][qpI];

                iBar +=
                    w
                   *ionicCurrentPDE
                    (
                        VmQP[cellI][qpI],
                        qState.u1[cellI][qpI],
                        qState.u2[cellI][qpI],
                        qState.u3[cellI][qpI],
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

    void advanceQuadratureStatesExplicit
    (
        const List<scalarField>& VmQP,
        QuadratureState& qState,
        const scalar dt
    )
    {
        forAll(VmQP, cellI)
        {
            forAll(VmQP[cellI], qpI)
            {
                scalar du1dt = 0.0;
                scalar du2dt = 0.0;
                scalar du3dt = 0.0;

                reactionRates
                (
                    VmQP[cellI][qpI],
                    qState.u1[cellI][qpI],
                    qState.u2[cellI][qpI],
                    qState.u3[cellI][qpI],
                    du1dt,
                    du2dt,
                    du3dt
                );

                qState.u1[cellI][qpI] += dt*du1dt;
                qState.u2[cellI][qpI] += dt*du2dt;
                qState.u3[cellI][qpI] += dt*du3dt;
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
        const bool twoD = mesh.nGeometricD() == 2;

        const CompactListList<point>& cellQP = LREInterp.cellQuadPoints();
        const CompactListList<scalar>& cellQW = LREInterp.cellQuadWeight();

        tmp<volVectorField> tGrad = LREInterp.grad(fld);
        const volVectorField& grad = tGrad();

        tmp<volSymmTensorField> tHess;
        const volSymmTensorField* hessPtr = nullptr;
        if (LREInterp.order() >= 2)
        {
            tHess = LREInterp.hessian(fld);
            hessPtr = &tHess();
        }

        autoPtr<List<LRE::symmTensor3Order>> thirdPtr;
        const List<LRE::symmTensor3Order>* thirdList = nullptr;
        if (LREInterp.order() >= 3)
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

        forAll(fld.internalField(), cellI)
        {
            const scalar cellV = V[cellI];

            forAll(cellQP[cellI], qpI)
            {
                const point& gp = cellQP[cellI][qpI];
                const scalar w = cellQW[cellI][qpI];
                const scalar meas = cellV*w;
                const vector d = gp - C[cellI];

                const symmTensor* H =
                    hessPtr ? &((*hessPtr)[cellI]) : nullptr;
                const LRE::symmTensor3Order* T3 =
                    thirdList ? &((*thirdList)[cellI]) : nullptr;

                const scalar num =
                    reconstructFromTaylor
                    (
                        fld[cellI],
                        grad[cellI],
                        H,
                        T3,
                        d,
                        twoD
                    );

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

    ReconstructionErrorSummary computeQuadratureFieldError
    (
        const List<scalarField>& fldQP,
        const fvMesh& mesh,
        const CompactListList<point>& cellQP,
        const CompactListList<scalar>& cellQW,
        const ManufacturedFieldID fieldID,
        const scalar t,
        const label dim
    )
    {
        const scalarField& V = mesh.V();

        scalar e1 = 0.0;
        scalar e2 = 0.0;
        scalar eInf = 0.0;
        scalar n2 = 0.0;
        scalar volTot = 0.0;

        forAll(fldQP, cellI)
        {
            const scalar cellV = V[cellI];

            forAll(fldQP[cellI], qpI)
            {
                const scalar meas = cellV*cellQW[cellI][qpI];
                const scalar ex = exactFieldValue(fieldID, cellQP[cellI][qpI], t, dim);
                const scalar err = fldQP[cellI][qpI] - ex;

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

    void writeSummary
    (
        const Time& runTime,
        const scalar beta,
        const label nSteps,
        const scalar dt,
        const bool useHighOrderODE,
        const CellErrorSummary& VmCell,
        const CellErrorSummary& u1Cell,
        const CellErrorSummary& u2Cell,
        const ReconstructionErrorSummary& VmHO,
        const ReconstructionErrorSummary& u1HO,
        const ReconstructionErrorSummary& u2HO,
        const TimingSummary& timing,
        const fvMesh& mesh
    )
    {
        const label N = estimatedN(mesh);
        const scalar dx = std::cbrt(mesh.V().average().value());
        const word dimName = name(mesh.nGeometricD()) + "D";

        const fileName outFile =
            runTime.path()
          / (dimName + "_" + name(N) + "_cells_explicit.dat");

        OFstream os(outFile);

        os  << "Manufactured-solution error summary (t = "
            << runTime.value() << "):" << nl
            << "Cell-centred error summary:" << nl
            << "Field     L1-error       L2-error       Linf-error      error_cell" << nl
            << "CVm      " << VmCell.L1 << "   " << VmCell.L2 << "   "
            << VmCell.Linf << "   " << VmCell.errorCell << nl
            << "Cu1      " << u1Cell.L1 << "   " << u1Cell.L2 << "   "
            << u1Cell.Linf << "   " << u1Cell.errorCell << nl
            << "Cu2      " << u2Cell.L1 << "   " << u2Cell.L2 << "   "
            << u2Cell.Linf << "   " << u2Cell.errorCell << nl
            << "-------------------------------------------------" << nl << nl
            << "High-order state/error summary:" << nl
            << "Field     L1_HO          L2_HO          Linf_HO" << nl;

        if (useHighOrderODE)
        {
            os  << "Note: Hu1/Hu2 and Gu1/Gu2 are evaluated directly from persistent "
                << "cell quadrature-point ODE states." << nl;
        }
        else
        {
            os  << "Note: Hu1/Hu2 and Gu1/Gu2 are evaluated from high-order "
                << "reconstruction of cell-centred states." << nl;
        }

        os  << "HVm      " << VmHO.L1 << "   " << VmHO.L2 << "   "
            << VmHO.Linf << nl
            << "Hu1      " << u1HO.L1 << "   " << u1HO.L2 << "   "
            << u1HO.Linf << nl
            << "Hu2      " << u2HO.L1 << "   " << u2HO.L2 << "   "
            << u2HO.Linf << nl
            << "-------------------------------------------------" << nl << nl
            << "Gauss-reconstructed L2 summary:" << nl
            << "Field     ea_mesh_gauss   norm_mesh_gauss   er_mesh_gauss(%)" << nl
            << "GVm      " << VmHO.L2 << "   " << VmHO.normL2
            << "   " << VmHO.relL2 << nl
            << "Gu1      " << u1HO.L2 << "   " << u1HO.normL2
            << "   " << u1HO.relL2 << nl
            << "Gu2      " << u2HO.L2 << "   " << u2HO.normL2
            << "   " << u2HO.relL2 << nl
            << "-------------------------------------------------" << nl << nl
            << "Simulation summary:" << nl
            << "-------------------" << nl
            << "State integration mode     = "
            << (useHighOrderODE ? "quadrature-point ODE states (Approach B)"
                                : "cell-centred ODE states (Approach A)") << nl
            << "Wall-clock timing measured inside this solver:" << nl
            << "setup wall time (s)       = " << timing.setupWallTime << nl
            << "    Time spent before the explicit time loop starts." << nl
            << "    This includes field initialisation, exact manufactured fields," << nl
            << "    quadrature-state initialisation, and boundary-condition preparation." << nl
            << "time loop wall time (s)   = " << timing.timeLoopWallTime << nl
            << "    Time spent inside the transient explicit marching loop." << nl
            << "    This is the main algorithmic cost of the solver." << nl
            << "postprocess wall time (s) = " << timing.postProcessWallTime << nl
            << "    Time spent after the time loop finishes." << nl
            << "    This includes final error evaluation, summaries, and writes." << nl
            << "solver wall time (s)      = " << timing.solverWallTime << nl
            << "    Total wall-clock time measured inside the solver." << nl
            << "    It is approximately setup + time loop + postprocess." << nl
            << "Number of cells (N)      = " << N << nl
            << "Grid spacing (dx)        = " << dx << nl
            << "Time step (dt)           = " << dt << nl
            << "Number of steps          = " << nSteps << nl
            << "Final simulation time    = " << runTime.value() << nl
            << "beta                     = " << beta << nl
            << "-------------------" << nl;

        Info<< "Wrote summary to " << outFile << nl
            << "Vm:  Linf(cell) = " << VmCell.Linf
            << ", L2_HO = " << VmHO.L2
            << ", Linf_HO = " << VmHO.Linf << nl
            << "u1:  Linf(cell) = " << u1Cell.Linf
            << ", L2_HO = " << u1HO.L2
            << ", Linf_HO = " << u1HO.Linf << nl
            << "u2:  Linf(cell) = " << u2Cell.Linf
            << ", L2_HO = " << u2HO.L2
            << ", Linf_HO = " << u2HO.Linf << nl
            << "Timing: setup = " << timing.setupWallTime
            << " s, loop = " << timing.timeLoopWallTime
            << " s, post = " << timing.postProcessWallTime
            << " s, total = " << timing.solverWallTime
            << " s" << endl;
    }
}

int main(int argc, char* argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    using Clock = std::chrono::steady_clock;

    const auto solverWallStart = Clock::now();
    const auto setupWallStart = Clock::now();

    const label dim = mesh.nGeometricD();
    const scalar dt = runTime.deltaTValue();
    const scalar beta = computeBeta(conductivity, dim);
    const scalar chiVal = chi.value();
    const scalar CmVal = Cm.value();
    const scalar chiCmVal = chiVal*CmVal;

    const CompactListList<point>& stateQP = LREInterp_states.cellQuadPoints();
    const CompactListList<scalar>& stateQW = LREInterp_states.cellQuadWeight();

    QuadratureState qState;
    List<scalarField> VmQP;

    if (useHighOrder_ODE)
    {
        qState = makeQuadratureStateStorage(stateQP);
        VmQP = makeQuadratureStorage(stateQP, 0.0);
    }

    Info<< "Running manufactured explicit solver" << nl
        << "Dimension = " << dim << nl
        << "beta = " << beta << nl
        << "dt = " << dt << nl
        << "useHighOrder_ODE = " << useHighOrder_ODE << endl;

    fillExactFields(VmExact, u1Exact, u2Exact, runTime.value(), dim);

    Vm.primitiveFieldRef() = VmExact.primitiveField();

    if (useHighOrder_ODE)
    {
        initializeQuadratureStatesExact(qState, stateQP, runTime.value(), dim);
        projectQuadratureStatesToCells(qState, stateQW, u1, u2, u3);
    }
    else
    {
        u1.primitiveFieldRef() = u1Exact.primitiveField();
        u2.primitiveFieldRef() = u2Exact.primitiveField();
        u3 = dimensionedScalar("zero", dimless, 0.0);
    }

    Vm.correctBoundaryConditions();
    updateStateBoundaryValues(u1, u2, u3, runTime.value(), dim);
    u1.correctBoundaryConditions();
    u2.correctBoundaryConditions();
    u3.correctBoundaryConditions();

    const scalar setupWallTime =
        std::chrono::duration<scalar>(Clock::now() - setupWallStart).count();

    const auto timeLoopWallStart = Clock::now();

    label nSteps = 0;

    while (runTime.value() < runTime.endTime().value() - SMALL)
    {
        const scalar t = runTime.value();
        scalarField& VmI = Vm.primitiveFieldRef();

        if (useHighOrder_ODE)
        {
            reconstructScalarAtCellQuadrature
            (
                Vm,
                LREInterp_Vm,
                stateQP,
                VmQP
            );

            computeVolumeAveragedIionFromQuadratureStates
            (
                VmQP,
                qState,
                stateQW,
                beta,
                chiVal,
                CmVal,
                Iion
            );
        }
        else
        {
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
        }

        if (useHighOrder_Vm)
        {
            computeHighOrderLaplacian
            (
                Vm,
                conductivity,
                LREInterp_Vm,
                lapVm
            );
        }
        else
        {
            lapVm = fvc::laplacian(conductivity, Vm);
        }

        if (useHighOrder_ODE)
        {
            forAll(VmI, cellI)
            {
                const scalar dVdt = lapVm[cellI]/chiCmVal - Iion[cellI];
                VmI[cellI] += dt*dVdt;
            }

            advanceQuadratureStatesExplicit(VmQP, qState, dt);
            projectQuadratureStatesToCells(qState, stateQW, u1, u2, u3);
        }
        else
        {
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

                const scalar dVdt = lapVm[cellI]/chiCmVal - Iion[cellI];

                VmI[cellI] += dt*dVdt;
                u1I[cellI] += dt*du1dt;
                u2I[cellI] += dt*du2dt;
                u3I[cellI] += dt*du3dt;
            }
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
                << " , Linf(Vm) = " << linfNorm(Vm)
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

    const scalar timeLoopWallTime =
        std::chrono::duration<scalar>(Clock::now() - timeLoopWallStart).count();

    const auto postProcessWallStart = Clock::now();

    fillExactFields(VmExact, u1Exact, u2Exact, runTime.value(), dim);

    if (useHighOrder_ODE)
    {
        reconstructScalarAtCellQuadrature
        (
            Vm,
            LREInterp_Vm,
            stateQP,
            VmQP
        );

        computeVolumeAveragedIionFromQuadratureStates
        (
            VmQP,
            qState,
            stateQW,
            beta,
            chiVal,
            CmVal,
            Iion
        );
    }
    else
    {
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
    }

    if (useHighOrder_Vm)
    {
        computeHighOrderLaplacian
        (
            Vm,
            conductivity,
            LREInterp_Vm,
            lapVm
        );
    }
    else
    {
        lapVm = fvc::laplacian(conductivity, Vm);
    }

    VmError = Vm - VmExact;
    u1Error = u1 - u1Exact;
    u2Error = u2 - u2Exact;

    const CellErrorSummary VmCell = computeCellErrorSummary(VmError);
    const CellErrorSummary u1Cell = computeCellErrorSummary(u1Error);
    const CellErrorSummary u2Cell = computeCellErrorSummary(u2Error);

    const ReconstructionErrorSummary VmHO =
        computeGaussReconstructedError(Vm, runTime.value(), mfVm, dim, LREInterp_Vm);

    ReconstructionErrorSummary u1HO{};
    ReconstructionErrorSummary u2HO{};

    if (useHighOrder_ODE)
    {
        u1HO =
            computeQuadratureFieldError
            (
                qState.u1,
                mesh,
                stateQP,
                stateQW,
                mfU1,
                runTime.value(),
                dim
            );

        u2HO =
            computeQuadratureFieldError
            (
                qState.u2,
                mesh,
                stateQP,
                stateQW,
                mfU2,
                runTime.value(),
                dim
            );
    }
    else
    {
        u1HO =
            computeGaussReconstructedError(u1, runTime.value(), mfU1, dim, LREInterp_states);

        u2HO =
            computeGaussReconstructedError(u2, runTime.value(), mfU2, dim, LREInterp_states);
    }

    runTime.write();

    const scalar postProcessWallTime =
        std::chrono::duration<scalar>(Clock::now() - postProcessWallStart).count();

    const scalar solverWallTime =
        std::chrono::duration<scalar>(Clock::now() - solverWallStart).count();

    const TimingSummary timing
    {
        setupWallTime,
        timeLoopWallTime,
        postProcessWallTime,
        solverWallTime
    };

    writeSummary
    (
        runTime,
        beta,
        nSteps,
        dt,
        useHighOrder_ODE,
        VmCell,
        u1Cell,
        u2Cell,
        VmHO,
        u1HO,
        u2HO,
        timing,
        mesh
    );

    Info<< "End" << nl << endl;
    return 0;
}
