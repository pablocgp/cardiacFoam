#include "fvCFD.H"
#include "LRE.H"
#include <cmath>

namespace
{
    scalar exactTemperature
    (
        const point& p,
        const scalar t,
        const scalar A,
        const scalar beta
    )
    {
        const scalar x = p.x();
        return A*(x + std::exp(-beta*t)*pow(x, 5)*(1.0 - x));
    }

    scalar exactSource
    (
        const point& p,
        const scalar t,
        const scalar A,
        const scalar alpha,
        const scalar rhoCp,
        const scalar beta
    )
    {
        const scalar x = p.x();
        const scalar g = std::exp(-beta*t);
        const scalar f = pow(x, 5)*(1.0 - x);
        const scalar d2f = 20.0*pow(x, 3) - 30.0*pow(x, 4);

        return -A*g*(rhoCp*beta*f + alpha*d2f);
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

    void fillExactFields
    (
        volScalarField& TExact,
        volScalarField& sourceQ,
        const scalar t,
        const scalar A,
        const scalar alpha,
        const scalar rhoCp,
        const scalar beta,
        const LRE& LREInterp_T
    )
    {
        const fvMesh& mesh = TExact.mesh();
        const vectorField& C = mesh.C();
        const CompactListList<point>& cellQP = LREInterp_T.cellQuadPoints();
        const CompactListList<scalar>& cellQW = LREInterp_T.cellQuadWeight();

        forAll(C, cellI)
        {
            TExact[cellI] = exactTemperature(C[cellI], t, A, beta);

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
                        A,
                        alpha,
                        rhoCp,
                        beta
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
                    Tp[faceI] = exactTemperature(CfPatch[faceI], t, A, beta);
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

    void computeGaussReconstructedError
    (
        const volScalarField& T,
        const scalar t,
        const scalar A,
        const scalar beta,
        const LRE& LREInterp_T,
        scalar& eaTGauss,
        scalar& normTGauss,
        scalar& erTGauss
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

        scalar e2 = 0.0;
        scalar n2 = 0.0;
        scalar volTot = 0.0;

        forAll(T.internalField(), cellI)
        {
            const scalar cellV = V[cellI];

            forAll(cellQP[cellI], qpI)
            {
                const point& gp = cellQP[cellI][qpI];
                const scalar w = cellQW[cellI][qpI];
                const vector d = gp - C[cellI];

                scalar TNum = T[cellI] + (gradT[cellI] & d);

                if (hessTPtr)
                {
                    TNum += 0.5*quadraticForm((*hessTPtr)[cellI], d);
                }

                if (thirdTPtr)
                {
                    TNum +=
                        (1.0/6.0)
                       *LRE::cubicForm((*thirdTPtr)[cellI], d, twoD);
                }

                const scalar TAnal = exactTemperature(gp, t, A, beta);

                e2 += cellV*w*sqr(TAnal - TNum);
                n2 += cellV*w*sqr(TAnal);
                volTot += cellV*w;
            }
        }

        eaTGauss = std::sqrt(e2/max(volTot, SMALL));
        normTGauss = std::sqrt(n2/max(volTot, SMALL));
        erTGauss = 100.0*eaTGauss/max(normTGauss, SMALL);
    }

    void writeSummary
    (
        const volScalarField& TError,
        const volScalarField& rhsField,
        const Time& runTime,
        const label nTimeSteps,
        const scalar dt,
        const scalar errorT,
        const scalar eaTGauss,
        const scalar normTGauss,
        const scalar erTGauss
    )
    {
        const fvMesh& mesh = TError.mesh();

        const scalar L1_T   = volumeWeightedL1(TError);
        const scalar L2_T   = volumeWeightedL2(TError);
        const scalar Linf_T = linfNorm(TError);

        const scalar Linf_RHS = linfNorm(rhsField);

        const label N = estimatedN(mesh);
        const scalar dx = std::cbrt(mesh.V().average().value());
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
            << "Field     ea_T_mesh_gauss   norm_T_mesh_gauss   er_T_mesh_gauss(%)" << nl
            << "TG      " << eaTGauss << "   " << normTGauss << "   " << erTGauss << nl
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
            << "-------------------" << nl;

        Info<< "Wrote summary to " << outFile << nl
            << "Cell-centred error_T = " << errorT << nl
            << "T error: L1 = " << L1_T
            << ", L2 = " << L2_T
            << ", Linf = " << Linf_T << nl
            << "Gauss error: ea_T_mesh_gauss = " << eaTGauss
            << ", er_T_mesh_gauss = " << erTGauss << "%" << nl
            << "Final-time RHS Linf = " << Linf_RHS << endl;
    }
}

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    const scalar amplitude =
        exactSolutionProperties.lookupOrDefault<scalar>("amplitude", 1.0);

    const scalar beta =
        exactSolutionProperties.lookupOrDefault<scalar>("beta", 1.0);

    const scalar rhoCp = rho.value()*Cp.value();
    const scalar dt = runTime.deltaTValue();

    fillExactFields
    (
        TExact,
        sourceQ,
        runTime.value(),
        amplitude,
        alpha.value(),
        rhoCp,
        beta,
        LREInterp_T
    );

    T = TExact;
    T.correctBoundaryConditions();

    label nTimeSteps = 0;

    while (runTime.value() < runTime.endTime().value() - SMALL)
    {
        const scalar t = runTime.value();

        fillExactFields
        (
            TExact,
            sourceQ,
            t,
            amplitude,
            alpha.value(),
            rhoCp,
            beta,
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

        residualField = lapTHO + sourceQ;

        forAll(T.internalField(), cellI)
        {
            T[cellI] += (dt/rhoCp)*residualField[cellI];
        }

        T.correctBoundaryConditions();

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
                amplitude,
                alpha.value(),
                rhoCp,
                beta,
                LREInterp_T
            );

            TError = T - TExact;
            operatorError = residualField;
            runTime.write();
        }
    }

    fillExactFields
    (
        TExact,
        sourceQ,
        runTime.value(),
        amplitude,
        alpha.value(),
        rhoCp,
        beta,
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

    residualField = lapTHO + sourceQ;
    TError = T - TExact;
    operatorError = residualField;

    scalar eaTGauss = 0.0;
    scalar normTGauss = 0.0;
    scalar erTGauss = 0.0;

    computeGaussReconstructedError
    (
        T,
        runTime.value(),
        amplitude,
        beta,
        LREInterp_T,
        eaTGauss,
        normTGauss,
        erTGauss
    );

    const scalar errorT = errorTCellCentred(TError);

    runTime.write();

    writeSummary
    (
        TError,
        residualField,
        runTime,
        nTimeSteps,
        dt,
        errorT,
        eaTGauss,
        normTGauss,
        erTGauss
    );

    Info<< "End" << nl << endl;
    return 0;
}
