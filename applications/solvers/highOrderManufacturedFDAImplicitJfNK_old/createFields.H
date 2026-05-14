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

Info<< "Reading time integration properties (ODE solver controls)\n" << endl;
IOdictionary timeIntegrationProperties
(
    IOobject
    (
        "timeIntegrationProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

Info<< "Reading exact solution properties (optional)\n" << endl;
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

volScalarField manufacturedSourceVm
(
    IOobject
    (
        "manufacturedSourceVm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("manufacturedSourceVm", dimVoltage/dimTime, 0.0),
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

// ----------------------------------------------------------------------- //
// Solver controls (mirror highOrderElectroActivationFoamImplicitJfNK).
// ----------------------------------------------------------------------- //

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
    implicitHOCoeffs.lookupOrDefault<word>("linearSolver", "BiCGSTAB")
);

const scalar implicitTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>("tolerance", 1.0e-12)
);

const label implicitMaxIterations
(
    implicitHOCoeffs.lookupOrDefault<label>("maxIterations", 2000)
);

const label implicitNonlinearIterations
(
    implicitHOCoeffs.lookupOrDefault<label>("nonlinearIterations", 3)
);

const label implicitMinNonlinearIterations
(
    implicitHOCoeffs.lookupOrDefault<label>("minNonlinearIterations", 1)
);

const word nonlinearMethod
(
    implicitHOCoeffs.lookupOrDefault<word>("nonlinearMethod", "picard")
);

const scalar nonlinearTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>("nonlinearTolerance", 1.0e-8)
);

const scalar nonlinearVmTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>("nonlinearVmTolerance", nonlinearTolerance)
);

const scalar nonlinearStatesTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>("nonlinearStatesTolerance", nonlinearTolerance)
);

const scalar nonlinearIionTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>("nonlinearIionTolerance", nonlinearTolerance)
);

const Switch nonlinearRequireStatesConvergence
(
    implicitHOCoeffs.lookupOrDefault<Switch>("nonlinearRequireStatesConvergence", false)
);

// If the nonlinear solver exhausts its iteration budget without converging,
// roll back Vm/Iion/states to the values at the beginning of the time step
// (default = false). Setting this to true accepts the unconverged iterate;
// useful only if the user knows the residual is "good enough".
const Switch nonlinearAcceptUnconverged
(
    implicitHOCoeffs.lookupOrDefault<Switch>("nonlinearAcceptUnconverged", false)
);

const scalar nonlinearRelaxation =
    max
    (
        min
        (
            implicitHOCoeffs.lookupOrDefault<scalar>("nonlinearRelaxation", 1.0),
            scalar(1.0)
        ),
        scalar(SMALL)
    );

// ----------------------------------------------------------------------- //
// JFNK controls (GMRES(m) + restart + Armijo line search + init guess).
// ----------------------------------------------------------------------- //

const label jfnkMaxKrylovIterations
(
    implicitHOCoeffs.lookupOrDefault<label>("jfnkMaxKrylovIterations", 30)
);

const label jfnkMaxRestarts
(
    implicitHOCoeffs.lookupOrDefault<label>("jfnkMaxRestarts", 3)
);

const scalar jfnkLinearTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>("jfnkLinearTolerance", 1.0e-4)
);

const scalar jfnkEpsilon
(
    implicitHOCoeffs.lookupOrDefault<scalar>("jfnkEpsilon", 1.0e-7)
);

// Initial-guess strategy at the start of each time step:
//   0: x_0 = Vm_n  (constant extrapolation, legacy)
//   1: x_0 = Vm_n + (Vm_n - Vm_{n-1}) * (dt / dt_prev)  (linear extrapolation)
// On the first time step there is no Vm_{n-1}, so we fall back to order 0.
const label jfnkInitGuessOrder
(
    implicitHOCoeffs.lookupOrDefault<label>("jfnkInitGuessOrder", 1)
);

// Safety bracket on the linear-extrapolation initial guess and (optionally)
// on Vm fed into the ionic ODE solver. Defaults are wide because MMS Vm is
// not in physiological volts (it is O(1) by construction).
const scalar jfnkInitGuessVmMin
(
    implicitHOCoeffs.lookupOrDefault<scalar>("jfnkInitGuessVmMin", -1.0e30)
);

const scalar jfnkInitGuessVmMax
(
    implicitHOCoeffs.lookupOrDefault<scalar>("jfnkInitGuessVmMax",  1.0e30)
);

// When true, Vm is clamped to [jfnkInitGuessVmMin, jfnkInitGuessVmMax] before
// being passed to ionicModel->solveODE. Off by default for MMS (the exact
// Vm is bounded by construction and the clamp is unnecessary).
const Switch jfnkClampODEInput
(
    implicitHOCoeffs.lookupOrDefault<Switch>("jfnkClampODEInput", false)
);

// Armijo backtracking line search around the damped Newton step.
const Switch jfnkLineSearch
(
    implicitHOCoeffs.lookupOrDefault<Switch>("jfnkLineSearch", true)
);

const label jfnkLineSearchMaxIter
(
    implicitHOCoeffs.lookupOrDefault<label>("jfnkLineSearchMaxIter", 10)
);

const scalar jfnkLineSearchAlphaMin
(
    implicitHOCoeffs.lookupOrDefault<scalar>("jfnkLineSearchAlphaMin", 1.0e-4)
);

const scalar jfnkArmijoC
(
    implicitHOCoeffs.lookupOrDefault<scalar>("jfnkArmijoC", 1.0e-4)
);

// FD step used by the diagonalIion linearisation (finite-difference Jacobian
// of Iion w.r.t. Vm, matched to solver2).
const scalar diagonalIionEpsilon
(
    implicitHOCoeffs.lookupOrDefault<scalar>("diagonalIionEpsilon", 1.0e-6)
);

const Switch writeNonlinearResiduals
(
    implicitHOCoeffs.lookupOrDefault<Switch>("writeNonlinearResiduals", true)
);

// ----------------------------------------------------------------------- //
// LRE high-order interpolators for Vm and Iion.
// ----------------------------------------------------------------------- //

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
    const word bcType = Vm.boundaryField()[patchI].type();

    if
    (
        bcType == fixedValueFvPatchScalarField::typeName
     || bcType == zeroGradientFvPatchScalarField::typeName
     || bcType == "fixedVoltage"
    )
    {
        includePatchInIionStencils[patchI] = true;
    }
}

const dictionary& highOrderCoeffs =
    spatialIntegrationProperties.subDict("highOrderCoeffs");

const dictionary lreCoeffsVm(highOrderCoeffs.subDict("LRECoeffs_Vm"));

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

const label lreN = lreCoeffsVm.lookupOrDefault<label>("N", -1);
const label lreNn = lreCoeffsVm.lookupOrDefault<label>("Nn", -1);
const label lreMaxStencil = lreCoeffsVm.lookupOrDefault<label>("maxStencilSize", -1);
const scalar lreK = lreCoeffsVm.lookupOrDefault<scalar>("k", -1.0);

LRE LREInterp_Vm
(
    mesh,
    includePatchInVmStencils,
    lreCoeffsVm
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

// ----------------------------------------------------------------------- //
// Polymorphic ionicModel (manufacturedFDA used here). Same factory call as
// in highOrderManufacturedFDAImplicit / highOrderElectroActivationFoamImplicit.
// ----------------------------------------------------------------------- //

dictionary solutionVariablesMemory(cardiacProperties);
solutionVariablesMemory.merge(timeIntegrationProperties);

autoPtr<ionicModel> ionicModel =
    ionicModel::New
    (
        solutionVariablesMemory,
        totalIionIntegrationPoints,
        runTime.deltaTValue()
    );

// State output fields (u1, u2, u3 from manufacturedFDA's exportedFieldNames).
// They are diagnostic projections of `states` to cell centres; the live state
// data lives in the `states` Field<Field<scalar>> built in main().
const wordList exportedNames = ionicModel->exportedFieldNames();

PtrList<volScalarField> outFields(exportedNames.size());
PtrList<volScalarField> outFieldsExact(exportedNames.size());
PtrList<volScalarField> outFieldsError(exportedNames.size());

forAll(exportedNames, i)
{
    outFields.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                exportedNames[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(exportedNames[i], dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    outFieldsExact.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                exportedNames[i] + "Exact",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(exportedNames[i] + "Exact", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    outFieldsError.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                exportedNames[i] + "Error",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(exportedNames[i] + "Error", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}
