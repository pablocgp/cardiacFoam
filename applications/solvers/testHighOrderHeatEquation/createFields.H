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

const dimensionedScalar alpha
(
    "alpha",
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
    dimensionedTensor("conductivity", dimless, alpha.value()*tensor::I)
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

boolList includePatchInStencils(mesh.boundaryMesh().size(), false);
forAll(includePatchInStencils, patchI)
{
    const word bcType = T.boundaryField()[patchI].type();

    if (bcType == "fixedValue" || bcType == "fixedVoltage")
    {
        includePatchInStencils[patchI] = true;
    }
}

LRE LREInterp_T
(
    mesh,
    includePatchInStencils,
    spatialIntegrationProperties.subDict("highOrderCoeffs").subDict("LRECoeffs_T")
);

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
