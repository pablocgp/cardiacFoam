// Define voltage dimensions for convenience
const dimensionSet dimVoltage(dimMass*dimArea/(pow3(dimTime)*dimCurrent));

// Transmembrane potential, the single PDE unknown. Resting value -84 mV;
// zeroGradient everywhere unless the case overrides it, which is the
// no-flux (electrically insulated) tissue boundary. READ_IF_PRESENT so a
// run can be restarted from a written field.
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

// Ionic current density produced by the cell model, the nonlinear source of
// the monodomain equation. Written for post-processing only; it inherits
// Vm's boundary types so the two fields stay compatible in expressions.
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

// Externally applied stimulus current, non-zero only inside the stimulus
// regions and only during their time windows (see applyStimulus).
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

// Stimulus protocol: geometric regions, their start times, duration and
// intensity. Read from constant/stimulusProtocol.
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

// One bounding box and one start time per stimulus site. The plural
// (...List) keys allow several sites, e.g. an S1-S2 protocol; the singular
// keys are kept as the one-site fallback.
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

// Cells inside each stimulus box, resolved once at start-up because the
// boxes do not move. A cell belongs to the first box that contains it
// (break), so overlapping boxes do not stimulate the same cell twice.
// In parallel each rank collects only its own cells, which is exactly what
// applyStimulus needs.
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

// The same sets as plain lists, which is the form the time loop iterates.
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

// Time integration and ODE controls: constant/timeIntegrationProperties.
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

// Tissue and cell-model properties: surface-to-volume ratio chi, membrane
// capacitance Cm, conductivity tensor and the ionic model name.
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

// Conductivity tensor D of the monodomain equation. Read from disk when a
// field is present - which is how anisotropic, fibre-oriented conductivity
// enters - and otherwise built as an isotropic tensor from the scalar in
// cardiacProperties.
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

// Spatial discretisation, nonlinear solver and linear solver controls:
// constant/spatialIntegrationProperties. Optional - every key below has a
// default that reproduces the low-order, single-backend behaviour.
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

// --------------------------------------------------------------------
// High-order switches
// --------------------------------------------------------------------
// Diffusion operator: false = two-point orthogonal finite volume,
// true = LRE reconstruction over the stencil in LRECoeffs_Vm.
const Switch useHighOrder_Vm
(
    spatialIntegrationProperties.lookupOrDefault<Switch>("useHighOrder_Vm", false)
);

// Ionic source: false = one evaluation at the cell centre, true =
// quadrature over the cell (which is what makes the source term match the
// order of the diffusion term). useHighOrder_states is the old key name.
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

// Cached form of the mode above, used at every point where the two paths
// diverge.
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


// --------------------------------------------------------------------
// Time discretisation and mass operator
// --------------------------------------------------------------------
// Theta scheme: backwardEuler (theta = 1, first order, unconditionally
// stable) or crankNicolson (theta = 1/2, second order).
const word implicitScheme
(
    spatialIntegrationProperties.lookupOrDefault<word>
    (
        "implicitScheme",
        "backwardEuler"
    )
);

// lumped puts the cell volume on the diagonal; consistent assembles the
// high-order mass matrix over the LRE stencil, which is what makes the
// temporal accuracy match the spatial one. Forced to lumped when
// useHighOrder_Vm is false.
const word massMatrixType
(
    spatialIntegrationProperties.lookupOrDefault<word>
    (
        "massMatrix",
        "lumped"
    )
);

// --------------------------------------------------------------------
// Linear solver for the assembled branches (Picard / diagonalIion)
// --------------------------------------------------------------------
// Solver name in backend-neutral spelling (GMRES, BiCGSTAB, SparseLU);
// each backend maps it onto its own.
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

// Eigen or PETSc. Only PETSc works at more than one rank: Eigen has no
// distributed matrix, so its path is serial by construction.
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

// Relative residual at which the linear solve stops.
const scalar implicitTolerance =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "tolerance",
        1e-10
    )
  : scalar(1e-10);

// Iteration cap for the linear solve.
const label implicitMaxIterations =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
    (
        "maxIterations",
        2000
    )
  : label(2000);

// PETSc-specific overrides of the above. Defaulting kspType to the
// backend-neutral name keeps a case that only sets linearSolver working.
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

// Preconditioner. NOTE for parallel comparisons: PETSc wraps ilu and lu in
// block-Jacobi at more than one rank, which makes them partition-dependent
// by construction. Use jacobi when a run has to be compared
// iteration-by-iteration against a serial one.
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

// Restart length m of GMRES(m); ignored by the other KSP types.
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

// Options prefix, so this solve can be tuned from the command line or
// PETSC_OPTIONS without affecting the JFNK ones.
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

// Let PETSc's options database override the settings above
// (KSPSetFromOptions).
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

// Phase B: when the assembled matrix only changes on its diagonal each
// nonlinear iteration (diagonalIion), optionally reuse the existing
// preconditioner instead of recomputing it (approximation; saves the ILU/LU
// refactorisation). Default false = exact behaviour (recompute every update).
const Switch petscReusePreconditioner
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<Switch>
    (
        "petscReusePreconditioner",
        false
    )
  : Switch(false)
);

// --------------------------------------------------------------------
// Nonlinear PDE/source coupling
// --------------------------------------------------------------------
// Iteration budget for the nonlinear solve. The default of 1 makes the
// scheme semi-implicit: one source evaluation and one linear solve per
// time step, with no convergence test.
const label implicitNonlinearIterations =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
    (
        "nonlinearIterations",
        1
    )
  : label(1);

// Iterations to perform before the stopping test is consulted at all,
// which prevents an accidental exit on a residual that has not moved yet.
const label implicitMinNonlinearIterations =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
    (
        "minNonlinearIterations",
        1
    )
  : label(1);

// picard (fixed point on the source), diagonalIion (adds a local
// dIion/dVm to the matrix diagonal) or JFNK (Newton with a matrix-free
// Jacobian). All three converge to the same solution; they differ in how
// many iterations and how much work per iteration.
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

// Default relative tolerance for every nonlinear residual; the three keys
// below override it per quantity.
const scalar nonlinearTolerance =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "nonlinearTolerance",
        1e-6
    )
  : scalar(1e-6);

// Stopping tolerance on the Vm increment.
const scalar nonlinearVmTolerance =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "nonlinearVmTolerance",
        nonlinearTolerance
    )
  : nonlinearTolerance;

// Stopping tolerance on the state (gating/concentration) increments.
const scalar nonlinearStatesTolerance =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "nonlinearStatesTolerance",
        nonlinearTolerance
    )
  : nonlinearTolerance;

// Stopping tolerance on the ionic-current increment.
const scalar nonlinearIionTolerance =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "nonlinearIionTolerance",
        nonlinearTolerance
    )
  : nonlinearTolerance;

// Whether the state residuals must also be converged, or only Vm and Iion.
// The states are far stiffer, so requiring them is much more expensive.
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

// Under-relaxation of the nonlinear update, clipped to (0, 1]. 1 is no
// relaxation and is the default.
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

// Relative tolerance of the JFNK inner Krylov solve. Loose on purpose:
// the inexact-Newton idea is to solve the Jacobian system only as well as
// the current Newton residual warrants.
const scalar jfnkLinearTolerance =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "jfnkLinearTolerance",
        1e-4
    )
  : scalar(1e-4);

// Step of the finite difference that provides the Jacobian-vector product,
// J v ~ (F(u + eps v) - F(u))/eps. Too large and the product is
// inaccurate, too small and it is cancellation noise.
const scalar jfnkEpsilon =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "jfnkEpsilon",
        1e-7
    )
  : scalar(1e-7);

// Backend for the JFNK inner Krylov solve, independent of the one used by
// the assembled branches. PETSc drives a KSP over a matrix-free MatShell;
// Eigen uses the in-house restarted GMRES.
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

// KSP type of the inner solve. GMRES is the usual choice: the matrix-free
// Jacobian is not symmetric.
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

// PETSc PC for the inner solve. Usually none, because the real
// preconditioning is done by jfnkPreconditioner below (a PCShell).
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

// Restart length of the inner GMRES; defaults to the Krylov dimension so
// the two cannot silently disagree.
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

// Options prefix of the inner solve, kept separate from the PDE solve's.
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

// Phase A3: preconditioner for the matrix-free JFNK Jacobian.
//   "none"      : no preconditioner (legacy behaviour, unchanged results path).
//   "diffusion" : right-preconditioning with M^{-1}, where M = AImplicit is the
//                 (constant, SPD) implicit diffusion operator = the linear part
//                 of the Jacobian. M is LU-factorised once (reused across steps),
//                 so each Krylov step costs one cheap triangular solve; this
//                 sharply cuts the number of matrix-free matvecs (each of which
//                 re-integrates the ionic ODE, the dominant cost).
const word jfnkPreconditioner
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<word>
    (
        "jfnkPreconditioner",
        "none"
    )
  : word("none")
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

// Upper bracket of the extrapolated initial guess. Together with
// ...VmMin it keeps the guess inside a physical range, so that a bad
// extrapolation cannot hand the ionic model a voltage it cannot evaluate.
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

// Maximum backtracking steps in the Armijo search.
const label jfnkLineSearchMaxIter =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
    (
        "jfnkLineSearchMaxIter",
        10
    )
  : label(10);

// Smallest step the search will accept before giving up.
const scalar jfnkLineSearchAlphaMin =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "jfnkLineSearchAlphaMin",
        1e-4
    )
  : scalar(1e-4);

// Constant c of the Armijo sufficient-decrease condition
// ||F(u + a d)|| <= (1 - c a) ||F(u)||.
const scalar jfnkArmijoC =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "jfnkArmijoC",
        1e-4
    )
  : scalar(1e-4);

// Finite-difference step used by the diagonalIion method to build the
// local derivative dIion/dVm that it adds to the matrix diagonal.
const scalar diagonalIionEpsilon =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "diagonalIionEpsilon",
        1e-6
    )
  : scalar(1e-6);

// Write the per-iteration nonlinear residual history, for diagnosing
// convergence.
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

// Print the setup/loop/post timing breakdown and the memory checkpoints.
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

// --------------------------------------------------------------------
// Memory controls for large high-order 3-D cases
// --------------------------------------------------------------------
// auto keeps the fast path for small meshes and enables the lower-RSS
// choices above the cell threshold below; on/off force either.
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

// Cell count above which auto switches the low-memory choices on.
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

// Size the assembly's triplet buffer from the actual operator instead of
// the worst case, which on a stabilised tet mesh is several times larger.
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

// Added for cardiacFoam: how often the per-step "Step N ... min(Vm) max(Vm)"
// line is printed. Was a hard-coded `nSteps % 100`.
//
// It is exposed because that line is the only time-resolved output the solver
// produces, and at dt = 1e-4 s a cadence of 100 samples every 10 ms. The TNNP
// upstroke is ~1 ms wide, so a 0-D single-cell test reading dV/dt from it
// measures ~12 V/s against a true ~300 V/s - a sampling artefact that looks
// like a broken ionic model. Setting this to 1 gives a sample per time step.
//
// Default 100: unchanged behaviour for every existing case.
const label stepLogInterval
(
    max
    (
        spatialIntegrationProperties.found("implicitHOCoeffs")
      ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
        (
            "stepLogInterval",
            100
        )
      : label(100),
        label(1)
    )
);

// Added for cardiacFoam: the chunk size of the chunked high-order stiffness
// assembly, which used to be a hard-coded 50000 local.
//
// It is a capacity control - K is assembled as sum_i K_chunk_i and
// addTripletIfNeeded filters individual contributions, not accumulated ones, so
// the chunking is transparent to the result. It is exposed because that claim
// was untestable while the value was baked into the binary: verifying it needs
// the SAME mesh assembled with more than one chunk and with exactly one, and
// with 50000 fixed any given mesh gives only one of the two. With the key, the
// regression test is a 640-cell case at chunk = 100 against the default.
//
// Note the chunk size does change the floating-point summation ORDER of K, so
// compare with a tolerance, not bitwise.
const label assemblyFaceChunk
(
    max
    (
        spatialIntegrationProperties.found("implicitHOCoeffs")
      ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
        (
            "assemblyFaceChunk",
            50000
        )
      : label(50000),
        label(1)
    )
);

// Triplets reserved per face for a flux-only (unstabilised) assembly.
// Capacity only: it never changes the assembled result.
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

// As above for a stabilised assembly, where each face also carries the
// reconstruction rows of both cells.
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

// Return the heap fragmentation left by the per-quadrature-point QR
// factorisations to the OS (malloc_trim) once setup is done.
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

// Assemble the consistent mass matrix in chunks rather than accumulating
// every triplet first.
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

// --------------------------------------------------------------------
// Activation tracking and benchmark sampling
// --------------------------------------------------------------------
// Vm level whose first crossing defines the activation time of a cell.
const scalar activationThreshold =
    spatialIntegrationProperties.lookupOrDefault<scalar>("activationThreshold", 0.0);

// Ends of the diagonal along which activation times are sampled, and the
// P8 corner point, following the Niederer benchmark. Defaults are that
// benchmark's geometry.
const point diagonalStart
(
    spatialIntegrationProperties.lookupOrDefault<point>
    (
        "diagonalStart",
        point(0.0, 0.0, 7.0e-3)
    )
);

// See diagonalStart.
const point diagonalEnd
(
    spatialIntegrationProperties.lookupOrDefault<point>
    (
        "diagonalEnd",
        point(1.9999e-2, 3.0e-3, 0.0)
    )
);

// Benchmark probe point; defaults to the far end of the diagonal.
const point p8Point
(
    spatialIntegrationProperties.lookupOrDefault<point>
    (
        "P8Point",
        diagonalEnd
    )
);

// Stop the run once the probe point has activated (plus the delay below),
// instead of integrating to endTime. What the benchmark measures is the
// arrival time, so there is nothing to gain from continuing.
const Switch stopAfterPointActivation
(
    spatialIntegrationProperties.lookupOrDefault<Switch>
    (
        "stopAfterPointActivation",
        false
    )
);

// Point whose activation triggers that stop.
const point stopActivationPoint
(
    spatialIntegrationProperties.lookupOrDefault<point>
    (
        "stopActivationPoint",
        p8Point
    )
);

// Extra time integrated after that activation, e.g. to capture the
// repolarisation that follows it.
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

// Abort if nothing has activated by this time: the configuration failed to
// ignite.
const scalar activationIgnitionTimeout
(
    spatialIntegrationProperties.lookupOrDefault<scalar>
    (
        "activationIgnitionTimeout",
        0.02
    )
);

// Abort if no new cell activates for this long while the stop point has
// not been reached: conduction block.
const scalar activationStallTime
(
    spatialIntegrationProperties.lookupOrDefault<scalar>
    (
        "activationStallTime",
        0.02
    )
);

// Number of sample points along the diagonal.
const label nDiagonalSamples =
    spatialIntegrationProperties.lookupOrDefault<label>("nDiagonalSamples", 20);

// Cells used by the inverse-distance interpolation at a sample point.
const label sampleNeighbours =
    spatialIntegrationProperties.lookupOrDefault<label>("sampleNeighbours", 8);

// Exponent of that inverse-distance weighting.
const scalar samplePower =
    spatialIntegrationProperties.lookupOrDefault<scalar>("samplePower", 2.0);

// --------------------------------------------------------------------
// Diagnostic fields, written for post-processing only
// --------------------------------------------------------------------
// Diffusion term div(D grad Vm) as evaluated by computeHighOrderLaplacian.
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

// Right-hand side of the semi-discrete equation, i.e. what dVm/dt equals.
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

// Face fluxes n . (D grad Vm) from the high-order reconstruction, whose
// divergence is lapVm.
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

// --------------------------------------------------------------------
// High-order reconstruction (LRE / moving least squares)
// --------------------------------------------------------------------
// Weight of the face-jump stabilisation. Zero (the default, and what the
// benchmarks use) switches it off entirely. It penalises the difference
// between the two cell values across a face and the value their
// reconstructions predict, which damps the odd-even decoupling the wide
// stencil can otherwise support; it vanishes for a field the
// reconstruction represents exactly, so it does not change the formal
// order.
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

// Reconstruction used by the Vm diffusion operator and by the consistent
// mass matrix. Held by pointer because it is built only when needed - its
// construction is the memory and time peak of the setup.
autoPtr<highOrderInterp> LREInterp_VmPtr;

// The stabilisation needs the reconstruction even when the operator itself
// is low order, hence the second condition.
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
        new highOrderInterp
        (
            mesh,
            includePatchInVmStencils,
            lreCoeffsVm
        )
    );
}

// Reconstruction used to evaluate Vm at the ionic-source quadrature points
// of each cell.
autoPtr<highOrderInterp> LREInterp_IionPtr;

// LREInterp_states: used by stateIntegrationMode = cellCentredReconstruct to
// reconstruct the cell-centred ionic states at the Iion Gauss points. Order
// (p1/p2/p3) set independently via highOrderCoeffs/LRECoeffs_states.
autoPtr<highOrderInterp> LREInterp_statesPtr;

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
        new highOrderInterp
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
            new highOrderInterp
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

// Number of points at which the ionic model has to be evaluated: one per
// cell in the low-order path, and the sum of the per-cell quadrature
// points in the high-order one. This sizes every per-point state array of
// the ionic model, so it is computed before the model is constructed.
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

// Modified for cardiacFoam: mesh.nCells() is passed as the size of the
// cell-centre scratch bank. The front-aware hybrid integrates Gauss points and
// cell centres in two separate solveODE() calls per iteration, and each set
// needs scratch slots of its own.
autoPtr<EmbeddedTNNPModel> ionicModel
(
    new EmbeddedTNNPModel
    (
        solutionVariablesMemory,
        totalIionIntegrationPoints,
        mesh.nCells(),
        runTime.deltaTValue()
    )
);

// Time at which each cell first crossed activationThreshold; zero means
// "not yet activated". READ_IF_PRESENT so a restart keeps the activation
// map already computed.
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

// Cells still waiting to activate. Cells that already carry an activation
// time (from a restart) are switched off so their value is not overwritten
// by a later crossing.
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
// State and algebraic quantities the ionic model exposes for writing; one
// volScalarField is created per name so they can be inspected in ParaView.
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
