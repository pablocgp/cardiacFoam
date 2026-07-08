// Define voltage dimensions for convenience
const dimensionSet dimVoltage(dimMass*dimArea/(pow3(dimTime)*dimCurrent));

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
    dimensionedScalar("Vm", dimVoltage, -0.084),
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
    Vm.boundaryField().types()
);

volScalarField externalStimulusCurrent
(
    IOobject
    (
        "externalStimulusCurrent",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("externalStimulusCurrent", dimCurrent/dimVolume, 0.0),
    Vm.boundaryField().types()
);

Info<< "Reading stimulus protocol" << endl;
IOdictionary stimulusProtocol
(
    IOobject
    (
        "stimulusProtocol",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

List<boundBox> stimulusBoxes;
List<scalar> stimulusStartTimes;

if
(
    stimulusProtocol.found("stimulusLocationMinList")
 && stimulusProtocol.found("stimulusLocationMaxList")
)
{
    List<point> mins;
    List<point> maxs;
    stimulusProtocol.lookup("stimulusLocationMinList") >> mins;
    stimulusProtocol.lookup("stimulusLocationMaxList") >> maxs;

    if (mins.size() != maxs.size())
    {
        FatalErrorInFunction
            << "stimulusLocationMinList and stimulusLocationMaxList must have "
            << "the same size." << abort(FatalError);
    }

    stimulusBoxes.setSize(mins.size());
    forAll(mins, i)
    {
        stimulusBoxes[i] = boundBox(mins[i], maxs[i]);
    }
}
else
{
    stimulusBoxes.setSize(1);
    stimulusBoxes[0] = boundBox
    (
        point(stimulusProtocol.lookup("stimulusLocationMin")),
        point(stimulusProtocol.lookup("stimulusLocationMax"))
    );
}

if (stimulusProtocol.found("stimulusStartTimeList"))
{
    stimulusProtocol.lookup("stimulusStartTimeList") >> stimulusStartTimes;
    if (stimulusStartTimes.size() != stimulusBoxes.size())
    {
        FatalErrorInFunction
            << "stimulusStartTimeList must have the same size as the stimulus "
            << "box list." << abort(FatalError);
    }
}
else
{
    stimulusStartTimes.setSize(stimulusBoxes.size());
    const scalar startTime =
        stimulusProtocol.lookupOrDefault<scalar>("stimulusStartTime", 0.0);

    forAll(stimulusStartTimes, i)
    {
        stimulusStartTimes[i] = startTime;
    }
}

List<labelHashSet> stimCellSets(stimulusBoxes.size());
labelHashSet allStimulusCells;

forAll(mesh.C(), cellI)
{
    const point& c = mesh.C()[cellI];
    forAll(stimulusBoxes, bI)
    {
        if (stimulusBoxes[bI].contains(c))
        {
            stimCellSets[bI].insert(cellI);
            allStimulusCells.insert(cellI);
            break;
        }
    }
}

List<labelList> stimulusCellIDsList(stimulusBoxes.size());
forAll(stimulusBoxes, bI)
{
    stimulusCellIDsList[bI] = stimCellSets[bI].toc();
}

const dimensionedScalar stimulusDuration
(
    "stimulusDuration",
    dimTime,
    stimulusProtocol
);

const dimensionedScalar stimulusIntensity
(
    "stimulusIntensity",
    dimCurrent/dimVolume,
    stimulusProtocol
);

Info<< "Stimulus boxes: " << stimulusBoxes.size() << nl
    << "Cells in stimulus region: " << allStimulusCells.size() << nl
    << "Stimulus start time(s): " << stimulusStartTimes << nl
    << "Stimulus duration: " << stimulusDuration.value() << nl
    << "Stimulus intensity: " << stimulusIntensity.value() << nl << endl;

Info<< "Reading numerical solver properties" << endl;
IOdictionary timeIntegrationProperties
(
    IOobject
    (
        "timeIntegrationProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading cardiac properties" << endl;
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

Info<< "chi = " << chi.value() << " " << chi.dimensions() << nl
    << "Cm = " << Cm.value() << " " << Cm.dimensions() << nl << endl;

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
    Info<< "conductivity not found on disk, using constant value from "
        << cardiacProperties.name() << endl;

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
else
{
    Info<< "conductivity field read from " << runTime.timeName() << endl;
}

Info<< "Reading optional spatial integration properties" << endl;
IOdictionary spatialIntegrationProperties
(
    IOobject
    (
        "spatialIntegrationProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

const Switch useHighOrder_Vm
(
    spatialIntegrationProperties.lookupOrDefault<Switch>("useHighOrder_Vm", false)
);

const Switch useHighOrder_Iion
(
    spatialIntegrationProperties.found("useHighOrder_Iion")
  ? spatialIntegrationProperties.lookupOrDefault<Switch>("useHighOrder_Iion", false)
  : spatialIntegrationProperties.lookupOrDefault<Switch>("useHighOrder_states", false)
);

// State (u/gating/concentration) handling for the high-order Iion path:
//   cellCentredReconstruct (default): the reaction ODE is integrated once per
//       cell at the cell centre; the states are then reconstructed at the Iion
//       Gauss points with a dedicated high-order LRE (LREInterp_states, order
//       set by highOrderCoeffs/LRECoeffs_states) and Iion is evaluated there.
//   gaussPointODE: legacy path. The ODE is integrated independently at every
//       Iion Gauss point. Cost: Nq ODE/cell.
// Only relevant when useHighOrder_Iion is true and dim > 1; otherwise the
// cell-centred path is used regardless.
const word stateIntegrationMode
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<word>
    (
        "stateIntegrationMode",
        "cellCentredReconstruct"
    )
  : word("cellCentredReconstruct")
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

// Front-aware hybrid: when enabled (orthogonal to stateIntegrationMode), the
// reaction ODE is integrated at the Iion Gauss points only in cells where the
// membrane potential varies strongly within the cell (the activation front and
// the stimulus), and at the cell centre + LRE reconstruction elsewhere. This
// recovers the sub-cell resolution of the ionic kinetics at the sharp upstroke
// (where the cell-average path fails to ignite/conduct) while keeping most of
// the cell-centred speed-up. When false, the solver behaves exactly as before
// (per stateIntegrationMode). Only active when useHighOrder_Iion and dim > 1.
const Switch frontAwareHybrid
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<Switch>
    (
        "frontAwareHybrid",
        false
    )
  : Switch(false)
);

// A cell is treated as a "front" cell (Gauss-point ODEs) when the sub-cell range
// of Vm exceeds this threshold [V] (default 5 mV), or when it is stimulated.
const scalar hybridFrontVmRange
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "hybridFrontVmRange",
        5.0e-3
    )
  : scalar(5.0e-3)
);

// Number of cell layers by which the (per-time-step frozen) front band is
// dilated, so the moving front stays inside the Gauss-resolved zone and the
// front/smooth transition is smooth. Default 2.
const label hybridFrontDilation
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
    (
        "hybridFrontDilation",
        2
    )
  : label(2)
);

// Barth--Jespersen-style limiter on the state reconstruction: bound each
// reconstructed Gauss value to the [min,max] of the cell + its face neighbours.
// Prevents high-order (p3) overshoot of the stiff TNNP states at under-resolved
// fronts from driving the ionic algebra to overflow (SIGFPE). Default on;
// applies to cellCentredReconstruct and the front-aware hybrid (smooth cells).
const Switch stateReconstructLimiter
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<Switch>
    (
        "stateReconstructLimiter",
        true
    )
  : Switch(true)
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

const word implicitLinearSolver
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<word>
    (
        "linearSolver",
        "GMRES"
    )
  : word("GMRES")
);

const word linearSolverBackend
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<word>
    (
        "linearSolverBackend",
        "PETSc"
    )
  : word("PETSc")
);

const scalar implicitTolerance =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "tolerance",
        1e-10
    )
  : scalar(1e-10);

const label implicitMaxIterations =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
    (
        "maxIterations",
        2000
    )
  : label(2000);

const word petscLinearKspType
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<word>
    (
        "petscLinearKspType",
        implicitLinearSolver
    )
  : implicitLinearSolver
);

const word petscLinearPcType
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<word>
    (
        "petscLinearPcType",
        "ilu"
    )
  : word("ilu")
);

const label petscLinearRestart
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
    (
        "petscLinearRestart",
        50
    )
  : label(50)
);

const word petscLinearOptionsPrefix
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<word>
    (
        "petscLinearOptionsPrefix",
        "activation_linear_"
    )
  : word("activation_linear_")
);

const Switch petscUseOptions
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<Switch>
    (
        "petscUseOptions",
        true
    )
  : Switch(true)
);

const label implicitNonlinearIterations =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
    (
        "nonlinearIterations",
        1
    )
  : label(1);

const label implicitMinNonlinearIterations =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
    (
        "minNonlinearIterations",
        1
    )
  : label(1);

const word nonlinearMethod
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<word>
    (
        "nonlinearMethod",
        "picard"
    )
  : word("picard")
);

const scalar nonlinearTolerance =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "nonlinearTolerance",
        1e-6
    )
  : scalar(1e-6);

const scalar nonlinearVmTolerance =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "nonlinearVmTolerance",
        nonlinearTolerance
    )
  : nonlinearTolerance;

const scalar nonlinearStatesTolerance =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "nonlinearStatesTolerance",
        nonlinearTolerance
    )
  : nonlinearTolerance;

const scalar nonlinearIionTolerance =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "nonlinearIionTolerance",
        nonlinearTolerance
    )
  : nonlinearTolerance;

const Switch nonlinearRequireStatesConvergence
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<Switch>
    (
        "nonlinearRequireStatesConvergence",
        false
    )
  : Switch(false)
);

// If the nonlinear solver exhausts its iteration budget without converging,
// roll back Vm/Iion/states to the values at the beginning of the time step
// (default = false). Setting this to true accepts the unconverged iterate;
// useful only if the user knows the residual is "good enough".
const Switch nonlinearAcceptUnconverged
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<Switch>
    (
        "nonlinearAcceptUnconverged",
        false
    )
  : Switch(false)
);

const scalar nonlinearRelaxation =
    max
    (
        min
        (
            spatialIntegrationProperties.found("implicitHOCoeffs")
          ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
            (
                "nonlinearRelaxation",
                1.0
            )
          : scalar(1.0),
            scalar(1.0)
        ),
        scalar(SMALL)
    );

// Maximum dimension of each inner Krylov subspace in GMRES(m).
// Larger values give better convergence but use more memory (O(m*nCells)).
const label jfnkMaxKrylovIterations =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
    (
        "jfnkMaxKrylovIterations",
        30
    )
  : label(30);

// Maximum number of outer GMRES restarts. After krylovDim inner Arnoldi steps
// the subspace is discarded and rebuilt from the current residual, allowing
// up to (krylovDim * (maxRestarts + 1)) total matrix-vector products.
const label jfnkMaxRestarts =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
    (
        "jfnkMaxRestarts",
        3
    )
  : label(3);

const scalar jfnkLinearTolerance =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "jfnkLinearTolerance",
        1e-4
    )
  : scalar(1e-4);

const scalar jfnkEpsilon =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "jfnkEpsilon",
        1e-7
    )
  : scalar(1e-7);

const word jfnkLinearSolverBackend
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<word>
    (
        "jfnkLinearSolverBackend",
        "PETSc"
    )
  : word("PETSc")
);

const word jfnkPetscKspType
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<word>
    (
        "jfnkPetscKspType",
        "gmres"
    )
  : word("gmres")
);

const word jfnkPetscPcType
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<word>
    (
        "jfnkPetscPcType",
        "none"
    )
  : word("none")
);

const label jfnkPetscRestart
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
    (
        "jfnkPetscRestart",
        jfnkMaxKrylovIterations
    )
  : jfnkMaxKrylovIterations
);

const word jfnkPetscOptionsPrefix
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<word>
    (
        "jfnkPetscOptionsPrefix",
        "activation_jfnk_"
    )
  : word("activation_jfnk_")
);

// ----------------------------------------------------------------------- //
// Initial-guess strategy for Newton at the start of each time step.
//
//   0: x_0 = Vm_n  (constant extrapolation; legacy behaviour)
//   1: x_0 = Vm_n + (Vm_n - Vm_{n-1}) * (dt / dt_prev)  (linear extrapolation)
//
// During the cardiac AP upstroke dVm/dt is large but smooth on the time
// scale of dt, so a linear extrapolation puts Newton MUCH closer to
// Vm^{n+1} than the lagged Vm_n. This drastically reduces the size of
// the Newton step (and therefore the FD-Jacobian noise sensitivity).
// On the first time step there is no Vm_{n-1}, so we fall back to order 0.
// ----------------------------------------------------------------------- //
const label jfnkInitGuessOrder =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
    (
        "jfnkInitGuessOrder",
        1
    )
  : label(1);

// Safety bracket on the extrapolated initial guess. If the linear
// extrapolation overshoots beyond [vmExtrapMin, vmExtrapMax] (in volts)
// we clamp; this protects Newton from initialising in unphysical territory.
const scalar jfnkInitGuessVmMin =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "jfnkInitGuessVmMin",
        -0.100
    )
  : scalar(-0.100);

const scalar jfnkInitGuessVmMax =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "jfnkInitGuessVmMax",
         0.060
    )
  : scalar(0.060);

// ----------------------------------------------------------------------- //
// Protective clamp on Vm passed to the ionic ODE solver.
//
// JFNK's matrix-free Jacobian probes  F(x + eps * v)  along arbitrary
// Krylov directions. When dt is large (or the Newton iterate is far from
// the solution) those probes can push Vm outside the physiologically
// valid range, where the ionic model is no longer calibrated and RKF45
// stalls trying to integrate ultra-stiff dynamics.
//
// This clamps Vm to [jfnkInitGuessVmMin, jfnkInitGuessVmMax] only on the
// path going into ionicModel->solveODE(...). Outside the bracket F(x)
// becomes slightly non-smooth, but the actual fixed-point solution lives
// inside the bracket so the clamp does not change the converged result.
// ----------------------------------------------------------------------- //
const Switch jfnkClampODEInput
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<Switch>
    (
        "jfnkClampODEInput",
        true
    )
  : Switch(true)
);

// ----------------------------------------------------------------------- //
// Armijo backtracking line search (globalisation of Newton).
//
// At each Newton iteration the update  x <- x + alpha * delta  is rejected
// if  ||F(x_new)||^2 > (1 - 2*armijoC*alpha) * ||F(x_old)||^2  (insufficient
// decrease). alpha is then halved and retried, until the condition holds,
// the iteration budget is exhausted, or alpha falls below the floor.
// ----------------------------------------------------------------------- //
const Switch jfnkLineSearch
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<Switch>
    (
        "jfnkLineSearch",
        true
    )
  : Switch(true)
);

const label jfnkLineSearchMaxIter =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
    (
        "jfnkLineSearchMaxIter",
        10
    )
  : label(10);

const scalar jfnkLineSearchAlphaMin =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "jfnkLineSearchAlphaMin",
        1e-4
    )
  : scalar(1e-4);

const scalar jfnkArmijoC =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "jfnkArmijoC",
        1e-4
    )
  : scalar(1e-4);

const scalar diagonalIionEpsilon =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "diagonalIionEpsilon",
        1e-6
    )
  : scalar(1e-6);

const Switch writeNonlinearResiduals
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<Switch>
    (
        "writeNonlinearResiduals",
        true
    )
  : Switch(true)
);

const Switch profileTimings
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<Switch>
    (
        "profileTimings",
        true
    )
  : Switch(true)
);

const word memoryOptimization
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<word>
    (
        "memoryOptimization",
        "auto"
    )
  : word("auto")
);

const label memoryOptimizationCellThreshold
(
    max
    (
        spatialIntegrationProperties.found("implicitHOCoeffs")
      ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
        (
            "memoryOptimizationCellThreshold",
            8000
        )
      : label(8000),
        label(0)
    )
);

const Switch memoryOptimizationAdaptiveTripletReserve
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<Switch>
    (
        "memoryOptimizationAdaptiveTripletReserve",
        true
    )
  : Switch(true)
);

const label memoryOptimizationFluxTripletsPerFace
(
    max
    (
        spatialIntegrationProperties.found("implicitHOCoeffs")
      ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
        (
            "memoryOptimizationFluxTripletsPerFace",
            200
        )
      : label(200),
        label(16)
    )
);

const label memoryOptimizationStabilisedTripletsPerFace
(
    max
    (
        spatialIntegrationProperties.found("implicitHOCoeffs")
      ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
        (
            "memoryOptimizationStabilisedTripletsPerFace",
            800
        )
      : label(800),
        label(16)
    )
);

const Switch memoryOptimizationTrimHeap
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<Switch>
    (
        "memoryOptimizationTrimHeap",
        true
    )
  : Switch(true)
);

const Switch memoryOptimizationCompactMassAssembly
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<Switch>
    (
        "memoryOptimizationCompactMassAssembly",
        true
    )
  : Switch(true)
);

const scalar activationThreshold =
    spatialIntegrationProperties.lookupOrDefault<scalar>("activationThreshold", 0.0);

const point diagonalStart
(
    spatialIntegrationProperties.lookupOrDefault<point>
    (
        "diagonalStart",
        point(0.0, 0.0, 7.0e-3)
    )
);

const point diagonalEnd
(
    spatialIntegrationProperties.lookupOrDefault<point>
    (
        "diagonalEnd",
        point(1.9999e-2, 3.0e-3, 0.0)
    )
);

const point p8Point
(
    spatialIntegrationProperties.lookupOrDefault<point>
    (
        "P8Point",
        diagonalEnd
    )
);

const Switch stopAfterPointActivation
(
    spatialIntegrationProperties.lookupOrDefault<Switch>
    (
        "stopAfterPointActivation",
        false
    )
);

const point stopActivationPoint
(
    spatialIntegrationProperties.lookupOrDefault<point>
    (
        "stopActivationPoint",
        p8Point
    )
);

const scalar stopDelayAfterActivation =
    max
    (
        spatialIntegrationProperties.lookupOrDefault<scalar>
        (
            "stopDelayAfterActivation",
            0.0
        ),
        0.0
    );

// Early-abort the run (cleanly, exit 0) when the activation fails, so a
// non-igniting / conduction-blocked configuration is not run to endTime.
// The reason is recorded in postProcessing/.../activationStatus.dat.
//   - no ignition: nothing has activated by activationIgnitionTimeout.
//   - conduction block: no new activation for activationStallTime while the
//     stop point has not been reached (needs stopAfterPointActivation).
const Switch abortOnActivationFailure
(
    spatialIntegrationProperties.lookupOrDefault<Switch>
    (
        "abortOnActivationFailure",
        false
    )
);

const scalar activationIgnitionTimeout
(
    spatialIntegrationProperties.lookupOrDefault<scalar>
    (
        "activationIgnitionTimeout",
        0.02
    )
);

const scalar activationStallTime
(
    spatialIntegrationProperties.lookupOrDefault<scalar>
    (
        "activationStallTime",
        0.02
    )
);

const label nDiagonalSamples =
    spatialIntegrationProperties.lookupOrDefault<label>("nDiagonalSamples", 20);

const label sampleNeighbours =
    spatialIntegrationProperties.lookupOrDefault<label>("sampleNeighbours", 8);

const scalar samplePower =
    spatialIntegrationProperties.lookupOrDefault<scalar>("samplePower", 2.0);

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

const scalar stabilisationAlpha
(
    spatialIntegrationProperties.found("highOrderCoeffs")
  ? max
    (
        spatialIntegrationProperties.subDict("highOrderCoeffs").lookupOrDefault<scalar>
        (
            "stabilisationAlpha",
            0.0
        ),
        scalar(0.0)
    )
  : scalar(0.0)
);

autoPtr<LRE> LREInterp_VmPtr;

const bool needLREForVm = useHighOrder_Vm || stabilisationAlpha > SMALL;

if (needLREForVm)
{
    if
    (
        !spatialIntegrationProperties.found("highOrderCoeffs")
     || !spatialIntegrationProperties.subDict("highOrderCoeffs").found("LRECoeffs_Vm")
    )
    {
        FatalErrorInFunction
            << "Vm LRE coefficients are required when useHighOrder_Vm is "
            << "true or stabilisationAlpha > 0, but "
            << "constant/spatialIntegrationProperties does not contain "
            << "highOrderCoeffs/LRECoeffs_Vm." << nl
            << abort(FatalError);
    }

    const dictionary& highOrderCoeffs =
        spatialIntegrationProperties.subDict("highOrderCoeffs");

    boolList includePatchInVmStencils(mesh.boundaryMesh().size(), false);
    wordList mirrorPatchNamesVm;

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
        }
        else if (bcType == zeroGradientFvPatchScalarField::typeName)
        {
            mirrorPatchNamesVm.append(patchName);
        }
    }

    dictionary lreCoeffsVm(highOrderCoeffs.subDict("LRECoeffs_Vm"));
    if (!mirrorPatchNamesVm.empty())
    {
        lreCoeffsVm.set("mirrorPatches", mirrorPatchNamesVm);
    }

    LREInterp_VmPtr.reset
    (
        new LRE
        (
            mesh,
            includePatchInVmStencils,
            lreCoeffsVm
        )
    );
}

autoPtr<LRE> LREInterp_IionPtr;

// LREInterp_states: used by stateIntegrationMode = cellCentredReconstruct to
// reconstruct the cell-centred ionic states at the Iion Gauss points. Order
// (p1/p2/p3) set independently via highOrderCoeffs/LRECoeffs_states.
autoPtr<LRE> LREInterp_statesPtr;

if (useHighOrder_Iion)
{
    if (!spatialIntegrationProperties.found("highOrderCoeffs"))
    {
        FatalErrorInFunction
            << "useHighOrder_Iion is true, but constant/spatialIntegrationProperties "
            << "does not contain highOrderCoeffs." << nl
            << abort(FatalError);
    }

    const dictionary& highOrderCoeffs =
        spatialIntegrationProperties.subDict("highOrderCoeffs");

    // Iion quadrature uses LRECoeffs_Iion exclusively. The legacy fallback that
    // accepted LRECoeffs_states for the Iion quadrature has been removed:
    // LRECoeffs_states now means the state-reconstruction interpolator.
    const word IionCoeffsName("LRECoeffs_Iion");
    if (!highOrderCoeffs.found(IionCoeffsName))
    {
        FatalErrorInFunction
            << "useHighOrder_Iion is true, but highOrderCoeffs does not "
            << "contain LRECoeffs_Iion." << nl
            << abort(FatalError);
    }

    boolList includePatchInIionStencils(mesh.boundaryMesh().size(), false);
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
            includePatchInIionStencils[patchI] = true;
        }
        else if (bcType == zeroGradientFvPatchScalarField::typeName)
        {
            mirrorPatchNamesIion.append(patchName);
        }
    }

    dictionary lreCoeffsIion(highOrderCoeffs.subDict(IionCoeffsName));
    if (!mirrorPatchNamesIion.empty())
    {
        lreCoeffsIion.set("mirrorPatches", mirrorPatchNamesIion);
    }

    LREInterp_IionPtr.reset
    (
        new LRE
        (
            mesh,
            includePatchInIionStencils,
            lreCoeffsIion
        )
    );

    // Build the state-reconstruction interpolator only when the new mode is
    // active and Gauss points exist. The states are not volScalarFields (they
    // live in a Field<Field<scalar>>) but behave like Vm at the boundary, so
    // the Iion include-list / mirror-patch classification is reused verbatim.
    // Both cellCentredReconstruct and the front-aware hybrid reconstruct the
    // cell-centred states at the Gauss points (the hybrid only in the smooth
    // cells), so the interpolator is built for either.
    if
    (
        (reconstructStatesFromCellCentres || frontAwareHybrid)
     && mesh.nGeometricD() > 1
    )
    {
        if (!highOrderCoeffs.found("LRECoeffs_states"))
        {
            FatalErrorInFunction
                << "stateIntegrationMode = cellCentredReconstruct or "
                << "frontAwareHybrid requires highOrderCoeffs/LRECoeffs_states "
                << "(state-reconstruction interpolator)." << nl
                << abort(FatalError);
        }

        dictionary lreCoeffsStates(highOrderCoeffs.subDict("LRECoeffs_states"));
        if (!mirrorPatchNamesIion.empty())
        {
            lreCoeffsStates.set("mirrorPatches", mirrorPatchNamesIion);
        }

        LREInterp_statesPtr.reset
        (
            new LRE
            (
                mesh,
                includePatchInIionStencils,
                lreCoeffsStates
            )
        );

        Info<< "LREInterp_states built (order "
            << LREInterp_statesPtr().order()
            << ") for "
            << (frontAwareHybrid ? "frontAwareHybrid" : "cellCentredReconstruct")
            << nl << endl;
    }
}

label totalIionIntegrationPoints = mesh.nCells();
if (useHighOrder_Iion)
{
    totalIionIntegrationPoints = 0;
    const CompactListList<point>& cellIionQuadP =
        LREInterp_IionPtr().cellQuadPoints();

    forAll(mesh.cells(), cellI)
    {
        totalIionIntegrationPoints += cellIionQuadP[cellI].size();
    }
}

Info<< "Iion integration points: " << totalIionIntegrationPoints << nl << endl;

// Merge the dictionaries and create a memory object to read all model properties.
dictionary solutionVariablesMemory(cardiacProperties);
solutionVariablesMemory.merge(timeIntegrationProperties);
solutionVariablesMemory.merge(stimulusProtocol);

autoPtr<EmbeddedTNNPModel> ionicModel
(
    new EmbeddedTNNPModel
    (
        solutionVariablesMemory,
        totalIionIntegrationPoints,
        runTime.deltaTValue()
    )
);

volScalarField activationTime
(
    IOobject
    (
        "activationTime",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("activationTime", dimTime, 0.0),
    Vm.boundaryField().types()
);

boolList calculateActivationTime(mesh.nCells(), true);

forAll(calculateActivationTime, cellI)
{
    if (activationTime[cellI] > SMALL)
    {
        calculateActivationTime[cellI] = false;
    }
}

volVectorField activationVelocity
(
    IOobject
    (
        "activationVelocity",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(1.0/(activationTime + dimensionedScalar("SMALL", dimTime, SMALL)))
);

PtrList<volScalarField> outFields;
const wordList exportedNames = ionicModel->exportedFieldNames();
outFields.setSize(exportedNames.size());

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
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            Vm.boundaryField().types()
        )
    );
}
