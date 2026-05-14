const dimensionSet dimVoltage(dimMass*dimArea/(pow3(dimTime)*dimCurrent));

Info<< "Reading cardiac properties\n" << endl;
IOdictionary cardiacProperties
(
    IOobject
    (
        "cardiacProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading spatial integration properties\n" << endl;
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

Info<< "Reading exact solution properties\n" << endl;
IOdictionary exactSolutionProperties
(
    IOobject
    (
        "exactSolutionProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar chi
(
    "chi",
    dimArea/dimVolume,
    cardiacProperties
);

const dimensionedScalar Cm
(
    "Cm",
    dimCurrent*dimTime/(dimVoltage*dimArea),
    cardiacProperties
);

volTensorField conductivity
(
    IOobject
    (
        "conductivity",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "zero",
        pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
        tensor::zero
    )
);

if (!conductivity.headerOk())
{
    Info<< "conductivity not found on disk, using constant value from cardiacProperties\n" << endl;

    conductivity =
        dimensionedTensor
        (
            dimensionedSymmTensor
            (
                "conductivity",
                pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
                cardiacProperties
            ) & tensor(I)
        );
}

volScalarField Vm
(
    IOobject
    (
        "Vm",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Vm", dimVoltage, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField u1
(
    IOobject
    (
        "u1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("u1", dimless, 0.0),
    fixedValueFvPatchScalarField::typeName
);

volScalarField u2
(
    IOobject
    (
        "u2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("u2", dimless, 0.0),
    fixedValueFvPatchScalarField::typeName
);

volScalarField u3
(
    IOobject
    (
        "u3",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("u3", dimless, 0.0),
    fixedValueFvPatchScalarField::typeName
);

volScalarField Iion
(
    IOobject
    (
        "ionicCurrent",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("ionicCurrent", dimVoltage/dimTime, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField lapVm
(
    IOobject
    (
        "lapVm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "lapVm",
        pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume)
       *(dimVoltage/dimLength)*dimArea/dimVolume,
        0.0
    )
);

volScalarField rhsVm
(
    IOobject
    (
        "rhsVm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhsVm", dimVoltage/dimTime, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField VmExact
(
    IOobject
    (
        "VmExact",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Vm
);

volScalarField u1Exact
(
    IOobject
    (
        "u1Exact",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    u1
);

volScalarField u2Exact
(
    IOobject
    (
        "u2Exact",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    u2
);

volScalarField VmError
(
    IOobject
    (
        "VmError",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("VmError", dimVoltage, 0.0),
    Vm.boundaryField().types()
);

volScalarField u1Error
(
    IOobject
    (
        "u1Error",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("u1Error", dimless, 0.0),
    u1.boundaryField().types()
);

volScalarField u2Error
(
    IOobject
    (
        "u2Error",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("u2Error", dimless, 0.0),
    u2.boundaryField().types()
);

surfaceScalarField fluxVm_HO
(
    IOobject
    (
        "fluxVm_HO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", lapVm.dimensions()*dimVolume, 0.0)
);

const Switch useHighOrder_Vm
(
    spatialIntegrationProperties.lookupOrDefault<Switch>("useHighOrder_Vm", true)
);

const Switch useHighOrder_Iion
(
    spatialIntegrationProperties.found("useHighOrder_Iion")
  ? spatialIntegrationProperties.lookupOrDefault<Switch>("useHighOrder_Iion", true)
  : spatialIntegrationProperties.lookupOrDefault<Switch>("useHighOrder_states", true)
);

boolList includePatchInVmStencils(mesh.boundaryMesh().size(), false);
forAll(includePatchInVmStencils, patchI)
{
    const word bcType = Vm.boundaryField()[patchI].type();

    if
    (
        bcType == fixedValueFvPatchScalarField::typeName
     || bcType == zeroGradientFvPatchScalarField::typeName
     || bcType == "fixedVoltage"
    )
    {
        includePatchInVmStencils[patchI] = true;
    }
}

boolList includePatchInIionStencils(mesh.boundaryMesh().size(), false);
forAll(includePatchInIionStencils, patchI)
{
    const word bcType = u1.boundaryField()[patchI].type();

    if
    (
        bcType == fixedValueFvPatchScalarField::typeName
     || bcType == zeroGradientFvPatchScalarField::typeName
    )
    {
        includePatchInIionStencils[patchI] = true;
    }
}

const dictionary& highOrderCoeffs =
    spatialIntegrationProperties.subDict("highOrderCoeffs");

word IionCoeffsName("LRECoeffs_Iion");
if (!highOrderCoeffs.found(IionCoeffsName))
{
    if (highOrderCoeffs.found("LRECoeffs_states"))
    {
        IionCoeffsName = "LRECoeffs_states";
        WarningInFunction
            << "Using legacy highOrderCoeffs/LRECoeffs_states for "
            << "Iion quadrature. Prefer LRECoeffs_Iion." << nl;
    }
    else
    {
        FatalErrorInFunction
            << "highOrderCoeffs contains neither LRECoeffs_Iion nor "
            << "legacy LRECoeffs_states." << nl
            << abort(FatalError);
    }
}

LRE LREInterp_Vm
(
    mesh,
    includePatchInVmStencils,
    highOrderCoeffs.subDict("LRECoeffs_Vm")
);

LRE LREInterp_Iion
(
    mesh,
    includePatchInIionStencils,
    highOrderCoeffs.subDict(IionCoeffsName)
);

label totalIionIntegrationPoints = mesh.nCells();
if (useHighOrder_Iion)
{
    totalIionIntegrationPoints = 0;
    const CompactListList<point>& cellIionQuadP =
        LREInterp_Iion.cellQuadPoints();

    forAll(mesh.cells(), cellI)
    {
        totalIionIntegrationPoints += cellIionQuadP[cellI].size();
    }
}

Info<< "Iion integration points: " << totalIionIntegrationPoints << nl << endl;
