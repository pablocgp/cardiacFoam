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

volScalarField sourceVm
(
    IOobject
    (
        "sourceVm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sourceVm", dimVoltage/dimTime, 0.0),
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
    implicitHOCoeffs.lookupOrDefault<word>("linearSolver", "GMRES")
);

const word linearSolverBackend
(
    implicitHOCoeffs.lookupOrDefault<word>("linearSolverBackend", "PETSc")
);

const scalar implicitTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>("tolerance", 1.0e-12)
);

const label implicitMaxIterations
(
    implicitHOCoeffs.lookupOrDefault<label>("maxIterations", 2000)
);

const word petscLinearKspType
(
    implicitHOCoeffs.lookupOrDefault<word>
    (
        "petscLinearKspType",
        implicitLinearSolver
    )
);

const word petscLinearPcType
(
    implicitHOCoeffs.lookupOrDefault<word>("petscLinearPcType", "ilu")
);

const label petscLinearRestart
(
    implicitHOCoeffs.lookupOrDefault<label>("petscLinearRestart", 50)
);

const word petscLinearOptionsPrefix
(
    implicitHOCoeffs.lookupOrDefault<word>
    (
        "petscLinearOptionsPrefix",
        "fdalinear_"
    )
);

const Switch petscUseOptions
(
    implicitHOCoeffs.lookupOrDefault<Switch>("petscUseOptions", true)
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
    implicitHOCoeffs.lookupOrDefault<word>("nonlinearMethod", "Picard")
);

// ----------------------------------------------------------------------- //
// State (u1,u2,u3) integration strategy for the high-order Iion path.
//
//   cellCentredReconstruct (default): the reaction ODE is integrated once per
//       cell at the cell centre (updateStateFieldsODE); the states are then
//       reconstructed at the Iion Gauss points with a dedicated high-order LRE
//       (LREInterp_states, order p1/p2/p3 set by highOrderCoeffs/LRECoeffs_states).
//       Iion is evaluated at the Gauss points (using the reconstructed Vm and
//       states) and integrated to the cell by quadrature. Cost: 1 ODE/cell.
//
//   gaussPointODE: legacy path. The ODE is integrated independently at every
//       Iion Gauss point (updateIntegrationPointStatesODE). Cost: Nq ODE/cell.
//
// Only relevant when useHighOrder_Iion is true and dim > 1; otherwise the
// cell-centred path is used regardless.
const word stateIntegrationMode
(
    implicitHOCoeffs.lookupOrDefault<word>
    (
        "stateIntegrationMode",
        "cellCentredReconstruct"
    )
);

if
(
    stateIntegrationMode != "cellCentredReconstruct"
 && stateIntegrationMode != "gaussPointODE"
)
{
    FatalErrorInFunction
        << "Unknown stateIntegrationMode: " << stateIntegrationMode << nl
        << "Valid options are cellCentredReconstruct (default) and gaussPointODE"
        << exit(FatalError);
}

const bool reconstructStatesFromCellCentres =
    (stateIntegrationMode == "cellCentredReconstruct");

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

const word jfnkLinearSolverBackend
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkLinearSolverBackend", "PETSc")
);

const word jfnkPetscKspType
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkPetscKspType", "gmres")
);

const word jfnkPetscPcType
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkPetscPcType", "none")
);

const label jfnkPetscRestart
(
    implicitHOCoeffs.lookupOrDefault<label>
    (
        "jfnkPetscRestart",
        jfnkMaxKrylovIterations
    )
);

const word jfnkPetscOptionsPrefix
(
    implicitHOCoeffs.lookupOrDefault<word>
    (
        "jfnkPetscOptionsPrefix",
        "jfnk_"
    )
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
// being passed to the manufactured state ODE integrator. Off by default for
// MMS because the exact Vm is bounded by construction.
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

// JFNK left preconditioner. The default uses the same local source
// linearisation as diagonalIion, assembled as
// P ~= AImplicit - theta*d(source)/d(Vm).
const word jfnkPreconditioner
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkPreconditioner", "diagonalIion")
);

const scalar jfnkPreconditionerDropTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>
    (
        "jfnkPreconditionerDropTolerance",
        1.0e-4
    )
);

const label jfnkPreconditionerFillFactor
(
    implicitHOCoeffs.lookupOrDefault<label>
    (
        "jfnkPreconditionerFillFactor",
        10
    )
);

const label jfnkPreconditionerUpdateFrequency
(
    max
    (
        implicitHOCoeffs.lookupOrDefault<label>
        (
            "jfnkPreconditionerUpdateFrequency",
            1
        ),
        label(1)
    )
);

const word jfnkPreconditionerKspType
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkPreconditionerKspType", "preonly")
);

const word jfnkPreconditionerPcType
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkPreconditionerPcType", "ilu")
);

const scalar jfnkPreconditionerTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>
    (
        "jfnkPreconditionerTolerance",
        1.0e-10
    )
);

const label jfnkPreconditionerMaxIterations
(
    implicitHOCoeffs.lookupOrDefault<label>
    (
        "jfnkPreconditionerMaxIterations",
        1
    )
);

const word jfnkPreconditionerOptionsPrefix
(
    implicitHOCoeffs.lookupOrDefault<word>
    (
        "jfnkPreconditionerOptionsPrefix",
        "jfnkpc_"
    )
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

const Switch profileTimings
(
    implicitHOCoeffs.lookupOrDefault<Switch>("profileTimings", true)
);

const word memoryOptimization
(
    implicitHOCoeffs.lookupOrDefault<word>("memoryOptimization", "auto")
);

const label memoryOptimizationCellThreshold
(
    max
    (
        implicitHOCoeffs.lookupOrDefault<label>
        (
            "memoryOptimizationCellThreshold",
            8000
        ),
        label(0)
    )
);

const Switch memoryOptimizationAdaptiveTripletReserve
(
    implicitHOCoeffs.lookupOrDefault<Switch>
    (
        "memoryOptimizationAdaptiveTripletReserve",
        true
    )
);

const label memoryOptimizationFluxTripletsPerFace
(
    max
    (
        implicitHOCoeffs.lookupOrDefault<label>
        (
            "memoryOptimizationFluxTripletsPerFace",
            200
        ),
        label(16)
    )
);

const label memoryOptimizationStabilisedTripletsPerFace
(
    max
    (
        implicitHOCoeffs.lookupOrDefault<label>
        (
            "memoryOptimizationStabilisedTripletsPerFace",
            800
        ),
        label(16)
    )
);

const Switch memoryOptimizationTrimHeap
(
    implicitHOCoeffs.lookupOrDefault<Switch>("memoryOptimizationTrimHeap", true)
);

const Switch memoryOptimizationCompactMassAssembly
(
    implicitHOCoeffs.lookupOrDefault<Switch>
    (
        "memoryOptimizationCompactMassAssembly",
        true
    )
);

const word stateODESolver
(
    implicitHOCoeffs.lookupOrDefault<word>
    (
        "stateODESolver",
        timeIntegrationProperties.lookupOrDefault<word>("solver", "RKF45")
    )
);

const scalar stateODEInitialStep
(
    implicitHOCoeffs.lookupOrDefault<scalar>
    (
        "stateODEInitialStep",
        timeIntegrationProperties.lookupOrDefault<scalar>("initialODEStep", 1.0e-5)
    )
);

const scalar stateODEAbsTol
(
    implicitHOCoeffs.lookupOrDefault<scalar>
    (
        "stateODEAbsTol",
        timeIntegrationProperties.lookupOrDefault<scalar>("absTol", SMALL)
    )
);

const scalar stateODERelTol
(
    implicitHOCoeffs.lookupOrDefault<scalar>
    (
        "stateODERelTol",
        timeIntegrationProperties.lookupOrDefault<scalar>
        (
            "relTol",
            implicitHOCoeffs.lookupOrDefault<scalar>("stateODETolerance", 1.0e-4)
        )
    )
);

const label stateODEMaxSteps
(
    implicitHOCoeffs.lookupOrDefault<label>
    (
        "stateODEMaxSteps",
        timeIntegrationProperties.lookupOrDefault<label>
        (
            "maxSteps",
            implicitHOCoeffs.lookupOrDefault<label>("stateODEMaxSubSteps", 1000000000)
        )
    )
);

const Switch stateODEUseOpenMP
(
    implicitHOCoeffs.lookupOrDefault<Switch>
    (
        "stateODEUseOpenMP",
        true
    )
);

const label stateODEOpenMPThreshold
(
    implicitHOCoeffs.lookupOrDefault<label>
    (
        "stateODEOpenMPThreshold",
        256
    )
);

const label stateODENumThreads
(
    implicitHOCoeffs.lookupOrDefault<label>
    (
        "stateODENumThreads",
        0
    )
);

// ----------------------------------------------------------------------- //
// LRE high-order interpolators for Vm and Iion.
// ----------------------------------------------------------------------- //

// Classify each boundary patch into one of three roles for the LRE:
//   - includePatchIn{Vm,Iion}Stencils  (Dirichlet → ghost point with value)
//   - mirrorPatchNames{Vm,Iion}        (zeroGradient → mirror-ghost stencil)
//   - skipped                          (everything else, e.g. empty)
//
// The mirror-ghost path was added in LRE to recover full high-order
// convergence at zeroGradient boundaries (see LRE constructor in
// $SOLIDS4FOAM_DIR/src/solids4FoamModels/numerics/LRE/LRE.C). A patch must
// not appear simultaneously in includePatchIn*Stencils and in
// mirrorPatchNames* — the LRE constructor enforces that with a FatalError.

boolList includePatchInVmStencils(mesh.boundaryMesh().size(), false);
boolList includePatchInIionStencils(mesh.boundaryMesh().size(), false);
wordList mirrorPatchNamesVm;
wordList mirrorPatchNamesIion;

forAll(mesh.boundaryMesh(), patchI)
{
    const word bcType = Vm.boundaryField()[patchI].type();
    const word& patchName = mesh.boundaryMesh()[patchI].name();

    if
    (
        bcType == fixedValueFvPatchScalarField::typeName
     || bcType == "fixedVoltage"
    )
    {
        includePatchInVmStencils[patchI] = true;
        includePatchInIionStencils[patchI] = true;
    }
    else if (bcType == zeroGradientFvPatchScalarField::typeName)
    {
        mirrorPatchNamesVm.append(patchName);
        mirrorPatchNamesIion.append(patchName);
    }
    // empty / processor / anything else: leave both lists at false; the
    // LRE will treat the patch with its default boundary handling (skip
    // for empty, NotImplemented for unsupported types).
}

const dictionary& highOrderCoeffs =
    spatialIntegrationProperties.subDict("highOrderCoeffs");

const scalar stabilisationAlpha
(
    max
    (
        highOrderCoeffs.lookupOrDefault<scalar>("stabilisationAlpha", 0.0),
        scalar(0.0)
    )
);

dictionary lreCoeffsVm(highOrderCoeffs.subDict("LRECoeffs_Vm"));
// Inject the mirror-patch list before constructing the LRE so that
// calcGlobalQRFaceGPCoeffs builds the doubled (mirror) stencil for those
// patches. The entry is optional and harmless when empty.
if (!mirrorPatchNamesVm.empty())
{
    lreCoeffsVm.set("mirrorPatches", mirrorPatchNamesVm);
}

// Iion quadrature uses LRECoeffs_Iion exclusively. The legacy fallback that
// accepted LRECoeffs_states for Iion has been removed: LRECoeffs_states now has
// a single, distinct meaning (the state-reconstruction interpolator used by
// stateIntegrationMode = cellCentredReconstruct).
const word IionCoeffsName("LRECoeffs_Iion");
if (!highOrderCoeffs.found(IionCoeffsName))
{
    FatalErrorInFunction
        << "highOrderCoeffs does not contain LRECoeffs_Iion (required for "
        << "the Iion quadrature)." << nl
        << abort(FatalError);
}

const label lreN = lreCoeffsVm.lookupOrDefault<label>("N", -1);
const label lreNn = lreCoeffsVm.lookupOrDefault<label>("Nn", -1);
const label lreMaxStencil = lreCoeffsVm.lookupOrDefault<label>("maxStencilSize", -1);
const scalar lreK = lreCoeffsVm.lookupOrDefault<scalar>("k", -1.0);

dictionary lreCoeffsIion(highOrderCoeffs.subDict(IionCoeffsName));
if (!mirrorPatchNamesIion.empty())
{
    lreCoeffsIion.set("mirrorPatches", mirrorPatchNamesIion);
}
const label lreIionN = lreCoeffsIion.lookupOrDefault<label>("N", -1);
const label lreIionNn = lreCoeffsIion.lookupOrDefault<label>("Nn", -1);
const label lreIionMaxStencil =
    lreCoeffsIion.lookupOrDefault<label>("maxStencilSize", -1);
const scalar lreIionK = lreCoeffsIion.lookupOrDefault<scalar>("k", -1.0);

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
    lreCoeffsIion
);

// ----------------------------------------------------------------------- //
// LREInterp_states: high-order interpolator used by
// stateIntegrationMode = cellCentredReconstruct to reconstruct the
// cell-centred states (u1,u2,u3) at the Iion Gauss points. Its order
// (p1/p2/p3) is set independently via highOrderCoeffs/LRECoeffs_states.
//
// Only built when actually needed (mode active, high-order Iion, dim > 1);
// otherwise left null so the legacy gaussPointODE path pays no extra cost.
// ----------------------------------------------------------------------- //
const bool buildStateInterp =
    reconstructStatesFromCellCentres
 && useHighOrder_Iion
 && mesh.nGeometricD() > 1;

// Build the include-list for the state stencils from u1's boundary types
// (the three states share fixedValue Dirichlet BCs in the MMS setup). Mirror
// the Vm classification: fixedValue -> Dirichlet ghost; zeroGradient -> mirror.
boolList includePatchInStateStencils(mesh.boundaryMesh().size(), false);
wordList mirrorPatchNamesStates;
forAll(mesh.boundaryMesh(), patchI)
{
    const word bcType = u1.boundaryField()[patchI].type();
    const word& patchName = mesh.boundaryMesh()[patchI].name();

    if
    (
        bcType == fixedValueFvPatchScalarField::typeName
     || bcType == "fixedVoltage"
    )
    {
        includePatchInStateStencils[patchI] = true;
    }
    else if (bcType == zeroGradientFvPatchScalarField::typeName)
    {
        mirrorPatchNamesStates.append(patchName);
    }
}

autoPtr<LRE> LREInterp_statesPtr;
if (buildStateInterp)
{
    if (!highOrderCoeffs.found("LRECoeffs_states"))
    {
        FatalErrorInFunction
            << "stateIntegrationMode = cellCentredReconstruct requires "
            << "highOrderCoeffs/LRECoeffs_states (state-reconstruction "
            << "interpolator)." << nl
            << abort(FatalError);
    }

    dictionary lreCoeffsStates(highOrderCoeffs.subDict("LRECoeffs_states"));
    if (!mirrorPatchNamesStates.empty())
    {
        lreCoeffsStates.set("mirrorPatches", mirrorPatchNamesStates);
    }

    LREInterp_statesPtr.reset
    (
        new LRE(mesh, includePatchInStateStencils, lreCoeffsStates)
    );

    Info<< "LREInterp_states built (order " << LREInterp_statesPtr().order()
        << ") for stateIntegrationMode = cellCentredReconstruct" << nl << endl;
}

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
