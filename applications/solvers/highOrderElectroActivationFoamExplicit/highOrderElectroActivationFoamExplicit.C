/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

Solver
    highOrderElectroActivationFoamExplicit

Description
    Explicit monodomain solver for cardiac electrophysiology. The ionic model is
    selected at run time and advanced with an explicit operator-splitting
    strategy. The diffusion operator can use either the standard finite-volume
    laplacian or the LRE high-order face-quadrature reconstruction. The ionic
    source can also be volume-integrated at LRE cell quadrature points, with
    the ionic states evolved directly at those quadrature points.

    The solver records activationTime as the first crossing of the prescribed
    activation threshold, default 0 V. It also writes point-sampled activation
    data at P8 and along the P1-P8 diagonal for the Niederer benchmark.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ionicModel.H"
#include "LRE.H"
#include "Field.H"
#include "volFields.H"
#include <cmath>

namespace
{
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

            IionCells[cellI] = iBar/max(wSum, SMALL);
        }

        Iion.correctBoundaryConditions();
    }

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
            runTime.path()/"postProcessing"/"highOrderElectroActivationFoamExplicit"
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
    const scalar dtStable = computeStableDeltaT
    (
        conductivity,
        chi.value(),
        Cm.value(),
        CFL,
        dx,
        dim
    );

    scalar dt = min(runTime.deltaTValue(), dtStable);
    runTime.setDeltaT(dt);

    Info<< "Running high-order electro activation explicit solver" << nl
        << "Dimension = " << dim << nl
        << "dx = " << dx << nl
        << "CFL = " << CFL << nl
        << "stable dt = " << dtStable << nl
        << "using dt = " << dt << nl
        << "useHighOrder_Vm = " << useHighOrder_Vm << nl
        << "useHighOrder_Iion = " << useHighOrder_Iion << nl
        << "Iion integration points = " << totalIionIntegrationPoints << nl
        << "activationThreshold = " << activationThreshold << " V" << nl
        << "stopAfterPointActivation = " << stopAfterPointActivation << nl
        << "stopActivationPoint = " << stopActivationPoint << nl
        << "stopDelayAfterActivation = " << stopDelayAfterActivation << " s"
        << nl << endl;

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

    label nSteps = 0;

    while (runTime.value() < effectiveEndTime - SMALL)
    {
        const scalar t0 = runTime.value();
        const scalar remaining = effectiveEndTime - t0;

        if (remaining <= SMALL)
        {
            break;
        }

        if (remaining < dt)
        {
            dt = remaining;
            runTime.setDeltaT(dt);
        }

        const scalarField VmOld(Vm.primitiveField());

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
                Vm,
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
                states
            );

            averageIionIntegrationPointsToCells
            (
                IionIntegrationPoints,
                LREInterp_IionPtr(),
                Iion
            );
        }
        else
        {
            ionicModel->calculateCurrent
            (
                t0,
                dt,
                Vm.internalField(),
                Iion,
                states
            );
            Iion.correctBoundaryConditions();
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

        scalarField& VmI = Vm.primitiveFieldRef();
        const scalarField& rhsVmI = rhsVm.primitiveField();
        forAll(VmI, cellI)
        {
            VmI[cellI] += dt*rhsVmI[cellI];
        }
        Vm.correctBoundaryConditions();

        if (useHighOrder_Iion)
        {
            reconstructVmAtIionIntegrationPoints
            (
                Vm,
                useHighOrder_Vm,
                useHighOrder_Vm ? &LREInterp_VmPtr() : nullptr,
                LREInterp_IionPtr(),
                VmIntegrationPoints
            );

            ionicModel->solveODE
            (
                t0,
                dt,
                VmIntegrationPoints,
                IionIntegrationPoints,
                states
            );

            averageIionIntegrationPointsToCells
            (
                IionIntegrationPoints,
                LREInterp_IionPtr(),
                Iion
            );
        }
        else
        {
            ionicModel->solveODE
            (
                t0,
                dt,
                Vm.internalField(),
                Iion,
                states
            );
            Iion.correctBoundaryConditions();
        }

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

        updateActivationTimes
        (
            VmOld,
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
