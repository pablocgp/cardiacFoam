IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary spatialIntegrationProperties
(
    IOobject
    (
        "spatialIntegrationProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary exactSolutionProperties
(
    IOobject
    (
        "exactSolutionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar alphaX
(
    "alphaX",
    dimless,
    transportProperties
);

const dimensionedScalar alphaY
(
    "alphaY",
    dimless,
    transportProperties
);

const dimensionedScalar alphaZ
(
    "alphaZ",
    dimless,
    transportProperties
);

const dimensionedScalar rho
(
    "rho",
    dimless,
    transportProperties
);

const dimensionedScalar Cp
(
    "Cp",
    dimless,
    transportProperties
);


volTensorField conductivity
(
    IOobject
    (
        "conductivity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "conductivity",
        dimless,
        tensor
        (
            alphaX.value(), 0.0,            0.0,
            0.0,            alphaY.value(), 0.0,
            0.0,            0.0,            alphaZ.value()
        )
    )
);

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("T", dimless, 0.0)
);

volScalarField TExact
(
    IOobject
    (
        "TExact",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    T
);


volScalarField sourceQ
(
    IOobject
    (
        "sourceQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sourceQ", dimless/sqr(dimLength), 0.0),
    T.boundaryField().types()
);



volScalarField lapTHO
(
    IOobject
    (
        "lapTHO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("lapTHO", dimless/sqr(dimLength), 0.0)
);

volScalarField residualField
(
    IOobject
    (
        "residualField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("residualField", dimless/sqr(dimLength), 0.0)
);

volScalarField TError
(
    IOobject
    (
        "TError",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("TError", dimless, 0.0)
);

surfaceScalarField fluxT_HO
(
    IOobject
    (
        "fluxT_HO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimLength, 0.0)
);

const Switch useHighOrder_T
(
    spatialIntegrationProperties.lookup("useHighOrder_T")
);

const word implicitScheme
(
    spatialIntegrationProperties.lookupOrDefault<word>
    (
        "implicitScheme",
        "backwardEuler"
    )
);

const word massMatrixType
(
    spatialIntegrationProperties.lookupOrDefault<word>
    (
        "massMatrix",
        "lumped"
    )
);

const dictionary implicitHOCoeffs
(
    spatialIntegrationProperties.subDict("implicitHOCoeffs")
);

const word implicitLinearSolver
(
    implicitHOCoeffs.lookupOrDefault<word>
    (
        "linearSolver",
        "BiCGSTAB"
    )
);

const scalar implicitTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>
    (
        "tolerance",
        1.0e-12
    )
);

const label implicitMaxIterations
(
    implicitHOCoeffs.lookupOrDefault<label>
    (
        "maxIterations",
        2000
    )
);

boolList includePatchInStencils(mesh.boundaryMesh().size(), false);
forAll(includePatchInStencils, patchI)
{
    const word bcType = T.boundaryField()[patchI].type();

    if (bcType == "fixedValue" || bcType == "fixedVoltage" || bcType == "zeroGradient")
    {
        includePatchInStencils[patchI] = true;
    }
}

const dictionary lreCoeffsT
(
    spatialIntegrationProperties.subDict("highOrderCoeffs").subDict("LRECoeffs_T")
);

const label lreN = lreCoeffsT.lookupOrDefault<label>("N", -1);
const label lreNn = lreCoeffsT.lookupOrDefault<label>("Nn", -1);
const label lreMaxStencil = lreCoeffsT.lookupOrDefault<label>("maxStencilSize", -1);
const scalar lreK = lreCoeffsT.lookupOrDefault<scalar>("k", -1.0);

LRE LREInterp_T(mesh, includePatchInStencils, lreCoeffsT);

volScalarField operatorError
(
    IOobject
    (
        "operatorError",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("operatorError", dimless/sqr(dimLength), 0.0)
);
