// Voltage dimensions, not among OpenFOAM's predefined sets.
const dimensionSet dimVoltage(dimMass*dimArea/(pow3(dimTime)*dimCurrent));

// --------------------------------------------------------------------
// Dictionaries
// --------------------------------------------------------------------
// chi, Cm and the conductivity tensor of the manufactured problem.
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

// High-order switches, nonlinear and linear solver controls, LRE stencil
// definitions and the memory controls. The bulk of what follows.
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

// State-ODE integrator and its tolerances.
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

// Optional overrides of the manufactured solution's parameters.
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

// Conductivity tensor D. Read from disk if a field is present - which is
// how a non-uniform or anisotropic D enters - and otherwise built as an
// isotropic tensor from the scalar in cardiacProperties. Note computeBeta()
// reads cell 0, so the manufactured source assumes it is uniform.
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

// --------------------------------------------------------------------
// Solution fields
// --------------------------------------------------------------------
// Transmembrane potential: the PDE unknown. Its boundary condition is the
// manufactured one, imposed as fixedValue and refreshed every time step.
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

// State variables of the reaction system. They are LOCAL ODEs with no
// spatial coupling, so their boundary values matter only where a
// reconstruction stencil reaches a boundary.
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

// See u1.
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

// See u1. Its exact solution is identically zero; it is carried so the code
// path matches the physiological solver's.
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

// Ionic current, the nonlinear source. Written for post-processing.
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

// --------------------------------------------------------------------
// Diagnostic fields, written for post-processing
// --------------------------------------------------------------------
// Diffusion term div(D grad Vm) as evaluated explicitly.
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
    zeroGradientFvPatchScalarField::typeName
);

// Source term of the Vm equation, -Iion plus the manufactured forcing.
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

// Exact manufactured solution sampled at the cell centres, and the
// corresponding error fields. These are what the error norms are computed
// from and what makes a wrong result visible in ParaView rather than only
// in a table.
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

// See VmExact.
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

// See VmExact.
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

// Vm minus its exact value.
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

// See VmError.
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

// See VmError.
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

// ----------------------------------------------------------------------- //
// Solver controls (mirror highOrderElectroActivationFoamImplicitJfNK).
// ----------------------------------------------------------------------- //

const Switch useHighOrder_Vm
(
    spatialIntegrationProperties.lookupOrDefault<Switch>("useHighOrder_Vm", true)
);

// Ionic source: false = one evaluation at the cell centre, true = quadrature
// over the cell, which is what makes the source term match the order of the
// diffusion term. Evaluating a nonlinear function at the centre is only
// second-order accurate however good the operator is.
const Switch useHighOrder_Iion
(
    spatialIntegrationProperties.found("useHighOrder_Iion")
  ? spatialIntegrationProperties.lookupOrDefault<Switch>("useHighOrder_Iion", true)
  : spatialIntegrationProperties.lookupOrDefault<Switch>("useHighOrder_states", true)
);

// --------------------------------------------------------------------
// Time discretisation and mass operator
// --------------------------------------------------------------------
// Theta scheme: backwardEuler (theta = 1, O(dt)) or crankNicolson
// (theta = 1/2, O(dt^2)).
const word implicitScheme
(
    spatialIntegrationProperties.lookupOrDefault<word>
    (
        "implicitScheme",
        "backwardEuler"
    )
);

// lumped puts the cell volume on the diagonal; consistent assembles the
// high-order mass matrix over the LRE stencil. Only consistent lets the
// temporal accuracy match the spatial one.
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

// --------------------------------------------------------------------
// Linear solver for the assembled branches (Picard / diagonalIion)
// --------------------------------------------------------------------
// Solver name in backend-neutral spelling (GMRES, BiCGSTAB, SparseLU).
const word implicitLinearSolver
(
    implicitHOCoeffs.lookupOrDefault<word>("linearSolver", "GMRES")
);

// Eigen or PETSc. Only PETSc works at more than one rank: Eigen has no
// distributed matrix, so its path is serial by construction.
const word linearSolverBackend
(
    implicitHOCoeffs.lookupOrDefault<word>("linearSolverBackend", "PETSc")
);

// Relative residual at which the linear solve stops. Tight here because a
// verification run must not have the linear solve above the discretisation
// error.
const scalar implicitTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>("tolerance", 1.0e-12)
);

// Iteration cap for the linear solve.
const label implicitMaxIterations
(
    implicitHOCoeffs.lookupOrDefault<label>("maxIterations", 2000)
);

// PETSc-specific overrides of the above.
const word petscLinearKspType
(
    implicitHOCoeffs.lookupOrDefault<word>
    (
        "petscLinearKspType",
        implicitLinearSolver
    )
);

// Preconditioner. NOTE for parallel comparisons: PETSc wraps ilu and lu in
// block-Jacobi at more than one rank, making them partition-dependent by
// construction. Use jacobi when comparing a parallel run against a serial
// one iteration by iteration.
const word petscLinearPcType
(
    implicitHOCoeffs.lookupOrDefault<word>("petscLinearPcType", "ilu")
);

// Restart length m of GMRES(m); ignored by the other KSP types.
const label petscLinearRestart
(
    implicitHOCoeffs.lookupOrDefault<label>("petscLinearRestart", 50)
);

// Options prefix, so this solve can be tuned from PETSC_OPTIONS without
// affecting the JFNK ones.
const word petscLinearOptionsPrefix
(
    implicitHOCoeffs.lookupOrDefault<word>
    (
        "petscLinearOptionsPrefix",
        "fdalinear_"
    )
);

// Let PETSc's options database override the settings above.
const Switch petscUseOptions
(
    implicitHOCoeffs.lookupOrDefault<Switch>("petscUseOptions", true)
);

// --------------------------------------------------------------------
// Nonlinear PDE/source coupling
// --------------------------------------------------------------------
// Iteration budget per time step.
const label implicitNonlinearIterations
(
    implicitHOCoeffs.lookupOrDefault<label>("nonlinearIterations", 3)
);

// Iterations performed before the stopping test is consulted, which
// prevents an exit on a residual that has not moved yet.
const label implicitMinNonlinearIterations
(
    implicitHOCoeffs.lookupOrDefault<label>("minNonlinearIterations", 1)
);

// Picard (fixed point on the source), diagonalIion (adds a local
// dIion/dVm to the matrix diagonal) or JFNK (Newton with a matrix-free
// Jacobian). All three converge to the same solution and differ in
// iterations versus work per iteration - comparing them on this
// manufactured problem is one of the purposes of the solver.
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

// True when the states are integrated once per cell and reconstructed at
// the quadrature points, false when an independent ODE is integrated at
// every quadrature point.
const bool reconstructStatesFromCellCentres =
    (stateIntegrationMode == "cellCentredReconstruct");

// Default relative tolerance for the nonlinear residuals; the three keys
// below override it per quantity.
const scalar nonlinearTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>("nonlinearTolerance", 1.0e-8)
);

// Stopping tolerance on the Vm increment.
const scalar nonlinearVmTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>("nonlinearVmTolerance", nonlinearTolerance)
);

// Stopping tolerance on the state increments.
const scalar nonlinearStatesTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>("nonlinearStatesTolerance", nonlinearTolerance)
);

// Stopping tolerance on the ionic-current increment.
const scalar nonlinearIionTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>("nonlinearIionTolerance", nonlinearTolerance)
);

// Whether the state residuals must also be converged, or only Vm and Iion.
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

// Under-relaxation of the nonlinear update, clipped to (0, 1].
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

// --------------------------------------------------------------------
// JFNK
// --------------------------------------------------------------------
// Outer restarts of the inner GMRES.
const label jfnkMaxRestarts
(
    implicitHOCoeffs.lookupOrDefault<label>("jfnkMaxRestarts", 3)
);

// Relative tolerance of the inner Krylov solve. Loose on purpose: the
// inexact-Newton idea is to solve the Jacobian system only as well as the
// current Newton residual warrants.
const scalar jfnkLinearTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>("jfnkLinearTolerance", 1.0e-4)
);

// Step of the finite difference providing the Jacobian-vector product,
// J v ~ (F(u + eps v) - F(u))/eps. Too large and the product is inaccurate,
// too small and it is cancellation noise; the optimum is set by the noise
// floor of the residual, i.e. by the state-ODE tolerance.
const scalar jfnkEpsilon
(
    implicitHOCoeffs.lookupOrDefault<scalar>("jfnkEpsilon", 1.0e-7)
);

// Backend for the inner Krylov solve, independent of the assembled one.
const word jfnkLinearSolverBackend
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkLinearSolverBackend", "PETSc")
);

// KSP type of the inner solve; GMRES, since the Jacobian is not symmetric.
const word jfnkPetscKspType
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkPetscKspType", "gmres")
);

// PETSc PC for the inner solve; usually none, because the preconditioning
// is done by the shell below.
const word jfnkPetscPcType
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkPetscPcType", "none")
);

// Restart length of the inner GMRES.
const label jfnkPetscRestart
(
    implicitHOCoeffs.lookupOrDefault<label>
    (
        "jfnkPetscRestart",
        jfnkMaxKrylovIterations
    )
);

// Options prefix of the inner solve.
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

// Upper bracket of the extrapolated initial guess; with ...VmMin it keeps
// the guess physical, so a bad extrapolation cannot hand the reaction terms
// a value they cannot evaluate.
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

// Maximum backtracking steps in the Armijo search.
const label jfnkLineSearchMaxIter
(
    implicitHOCoeffs.lookupOrDefault<label>("jfnkLineSearchMaxIter", 10)
);

// Smallest step the search will accept before giving up.
const scalar jfnkLineSearchAlphaMin
(
    implicitHOCoeffs.lookupOrDefault<scalar>("jfnkLineSearchAlphaMin", 1.0e-4)
);

// Constant c of the Armijo condition
// ||F(u + a d)|| <= (1 - c a) ||F(u)||.
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

// ILUT drop tolerance of the JFNK preconditioner: entries below it are
// discarded, trading accuracy of M^{-1} for fill-in.
const scalar jfnkPreconditionerDropTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>
    (
        "jfnkPreconditionerDropTolerance",
        1.0e-4
    )
);

// ILUT fill factor: bound on the entries kept per row.
const label jfnkPreconditionerFillFactor
(
    implicitHOCoeffs.lookupOrDefault<label>
    (
        "jfnkPreconditionerFillFactor",
        10
    )
);

// How often the preconditioner is rebuilt, in nonlinear iterations. The
// linear part of the Jacobian is constant, so a stale preconditioner is
// usually still a good one.
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

// KSP/PC used to APPLY the preconditioner; preonly + ilu means "apply the
// factorisation once", which is the intended use.
const word jfnkPreconditionerKspType
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkPreconditionerKspType", "preonly")
);

// See jfnkPreconditionerKspType.
const word jfnkPreconditionerPcType
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkPreconditionerPcType", "ilu")
);

// Tolerance and iteration cap of that inner application; both are nominal
// when the KSP is preonly.
const scalar jfnkPreconditionerTolerance
(
    implicitHOCoeffs.lookupOrDefault<scalar>
    (
        "jfnkPreconditionerTolerance",
        1.0e-10
    )
);

// See jfnkPreconditionerTolerance.
const label jfnkPreconditionerMaxIterations
(
    implicitHOCoeffs.lookupOrDefault<label>
    (
        "jfnkPreconditionerMaxIterations",
        1
    )
);

// Options prefix of the preconditioner solve.
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

// Write the per-iteration nonlinear history (NonlinearConvergenceRecord).
const Switch writeNonlinearResiduals
(
    implicitHOCoeffs.lookupOrDefault<Switch>("writeNonlinearResiduals", true)
);

// Print the setup/loop/post timing breakdown and the memory checkpoints.
const Switch profileTimings
(
    implicitHOCoeffs.lookupOrDefault<Switch>("profileTimings", true)
);

// --------------------------------------------------------------------
// Memory controls for large high-order 3-D cases
// --------------------------------------------------------------------
// auto keeps the fast path for small meshes and enables the lower-RSS
// choices above the cell threshold below; on/off force either.
const word memoryOptimization
(
    implicitHOCoeffs.lookupOrDefault<word>("memoryOptimization", "auto")
);

// Cell count above which auto switches the low-memory choices on.
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

// Size the assembly's triplet buffer from the actual operator instead of
// the worst case, which on a stabilised tet mesh is several times larger.
const Switch memoryOptimizationAdaptiveTripletReserve
(
    implicitHOCoeffs.lookupOrDefault<Switch>
    (
        "memoryOptimizationAdaptiveTripletReserve",
        true
    )
);

// Triplets reserved per face for a flux-only assembly. Capacity only: it
// never changes the assembled result.
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

// As above for a stabilised assembly, where each face also carries the
// reconstruction rows of both of its cells.
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

// Return the heap fragmentation left by the per-quadrature-point QR
// factorisations to the OS (malloc_trim) once setup is done.
const Switch memoryOptimizationTrimHeap
(
    implicitHOCoeffs.lookupOrDefault<Switch>("memoryOptimizationTrimHeap", true)
);

// Assemble the consistent mass matrix in chunks rather than accumulating
// every triplet first.
const Switch memoryOptimizationCompactMassAssembly
(
    implicitHOCoeffs.lookupOrDefault<Switch>
    (
        "memoryOptimizationCompactMassAssembly",
        true
    )
);

// --------------------------------------------------------------------
// State-ODE integration
// --------------------------------------------------------------------
// Integrator for the local reaction ODEs: RKF45 (adaptive, embedded 4(5)
// error estimate) or RK4 (fixed step).
const word stateODESolver
(
    implicitHOCoeffs.lookupOrDefault<word>
    (
        "stateODESolver",
        timeIntegrationProperties.lookupOrDefault<word>("solver", "RKF45")
    )
);

// First sub-step attempted; the adaptive controller takes over from there.
const scalar stateODEInitialStep
(
    implicitHOCoeffs.lookupOrDefault<scalar>
    (
        "stateODEInitialStep",
        timeIntegrationProperties.lookupOrDefault<scalar>("initialODEStep", 1.0e-5)
    )
);

// Absolute and relative tolerances of the adaptive sub-stepping. These set
// the NOISE FLOOR of the nonlinear residual, which is what limits how small
// jfnkEpsilon can usefully be and how tight nonlinearTolerance can be
// asked to be.
const scalar stateODEAbsTol
(
    implicitHOCoeffs.lookupOrDefault<scalar>
    (
        "stateODEAbsTol",
        timeIntegrationProperties.lookupOrDefault<scalar>("absTol", SMALL)
    )
);

// See stateODEAbsTol.
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

// Cap on sub-steps per PDE step, so a stiff point cannot hang the run.
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

// Integrate the per-point ODEs with OpenMP. Numerically identical either
// way (each point is an independent ODE); it pays off only where the ODE
// dominates the cost.
const Switch stateODEUseOpenMP
(
    implicitHOCoeffs.lookupOrDefault<Switch>
    (
        "stateODEUseOpenMP",
        true
    )
);

// Minimum number of points on a rank before the OpenMP region is entered
// at all: below it the fork/join costs more than the work.
const label stateODEOpenMPThreshold
(
    implicitHOCoeffs.lookupOrDefault<label>
    (
        "stateODEOpenMPThreshold",
        256
    )
);

// Thread count for that region; 0 or less leaves it to OpenMP.
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
// --------------------------------------------------------------------
// LRE reconstruction stencils
// --------------------------------------------------------------------
// How each boundary patch enters the reconstruction stencils, classified
// from Vm's boundary condition:
//   fixedValue  -> a Dirichlet ghost value is added to the stencil, which
//                  is what makes the reconstruction exact up to the
//                  boundary for the manufactured solution;
//   zeroGradient-> the stencil is mirrored across the patch;
//   empty/processor/other -> neither; the library applies its own default.
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

// Sub-dictionary holding the stabilisation weight and the three LRE stencil
// definitions (LRECoeffs_Vm, _Iion and _states).
const dictionary& highOrderCoeffs =
    spatialIntegrationProperties.subDict("highOrderCoeffs");

// Weight of the face-jump stabilisation. Zero disables it. It penalises
// the difference between the two cell values across a face and what their
// reconstructions predict, which damps the odd-even decoupling the wide
// stencil can otherwise support; it vanishes for a field the
// reconstruction represents exactly, so it does not change the formal
// order.
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

// Stencil parameters echoed to the log so a run records what it was built
// with. Only N (polynomial order), Nn (extra stencil points) and k (weight
// exponent) reach the library; maxStencilSize is inert, the adapter drops
// it.
const label lreN = lreCoeffsVm.lookupOrDefault<label>("N", -1);
const label lreNn = lreCoeffsVm.lookupOrDefault<label>("Nn", -1);
const label lreMaxStencil = lreCoeffsVm.lookupOrDefault<label>("maxStencilSize", -1);
const scalar lreK = lreCoeffsVm.lookupOrDefault<scalar>("k", -1.0);

dictionary lreCoeffsIion(highOrderCoeffs.subDict(IionCoeffsName));
if (!mirrorPatchNamesIion.empty())
{
    lreCoeffsIion.set("mirrorPatches", mirrorPatchNamesIion);
}
// Same echo for the Iion interpolator.
const label lreIionN = lreCoeffsIion.lookupOrDefault<label>("N", -1);
const label lreIionNn = lreCoeffsIion.lookupOrDefault<label>("Nn", -1);
const label lreIionMaxStencil =
    lreCoeffsIion.lookupOrDefault<label>("maxStencilSize", -1);
const scalar lreIionK = lreCoeffsIion.lookupOrDefault<scalar>("k", -1.0);

// Reconstruction used by the Vm diffusion operator and the consistent mass
// matrix. Constructing it is the memory and time peak of the setup: it
// builds every stencil and factorises one least-squares problem per cell
// and per face quadrature point.
highOrderInterp LREInterp_Vm
(
    mesh,
    includePatchInVmStencils,
    lreCoeffsVm
);

// Reconstruction used to evaluate Vm at the ionic quadrature points. Kept
// separate from LREInterp_Vm so the source term can be integrated at a
// different order from the diffusion operator, which is what the (Vm,
// states, Iion) triads of the convergence driver sweep.
highOrderInterp LREInterp_Iion
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

// Reconstruction of the cell-centred states at the ionic quadrature points.
// Built only for stateIntegrationMode = cellCentredReconstruct, where the
// ODE is solved once per cell instead of at every quadrature point.
autoPtr<highOrderInterp> LREInterp_statesPtr;
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
        new highOrderInterp(mesh, includePatchInStateStencils, lreCoeffsStates)
    );

    Info<< "LREInterp_states built (order " << LREInterp_statesPtr().order()
        << ") for stateIntegrationMode = cellCentredReconstruct" << nl << endl;
}

// Number of points at which the reaction terms are evaluated: one per cell
// in the low-order path, the sum of the per-cell quadrature points in the
// high-order one. Sizes every per-point state array.
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
