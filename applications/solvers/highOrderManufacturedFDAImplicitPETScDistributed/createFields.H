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

// All it supplied here were second NAMES for the state-ODE controls
// (solver, initialODEStep, absTol, relTol, maxSteps), which are read as
// stateODE* from implicitHOCoeffs and always won when both were present - two
// names for one setting. The launchers still write the file, for the electro
// solver's benefit; this solver ignores it.

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
// Solver controls.
// ----------------------------------------------------------------------- //

const Switch useHighOrder_Vm
(
    spatialIntegrationProperties.lookupOrDefault<Switch>("useHighOrder_Vm", true)
);

// Ionic source: false = one evaluation at the cell centre, true = quadrature
// over the cell, which is what makes the source term match the order of the
// diffusion term. Evaluating a nonlinear function at the centre is only
// second-order accurate however good the operator is.
// useHighOrder_states was the previous name of this key and is no longer 
// accepted. A dictionary that still carries it is refused rather than ignored: 
// ignoring it would quietly fall back to the default and change the order the 
// case was written to measure. The old name was wrong anyway - this key is the 
// quadrature of Iion, while how the states are treated is stateIntegrationMode 
// plus LRECoeffs_states.
if (spatialIntegrationProperties.found("useHighOrder_states"))
{
    FatalErrorInFunction
        << "useHighOrder_states was renamed to useHighOrder_Iion." << nl
        << "It selects the quadrature of the ionic source, not the treatment"
        << " of the states: that is stateIntegrationMode."
        << exit(FatalError);
}

const Switch useHighOrder_Iion
(
    spatialIntegrationProperties.lookupOrDefault<Switch>("useHighOrder_Iion", true)
);

// --------------------------------------------------------------------
// Time discretisation and mass operator
// --------------------------------------------------------------------
// Two families:
//
//   theta scheme  backwardEuler (theta = 1, O(dt)) or crankNicolson
//                 (theta = 1/2, O(dt^2));
//   BDF-k         BDF1, BDF2, BDF3, BDF4 - fully implicit, O(dt^k).
//
// BDF1 and backwardEuler are the same scheme and are asserted to produce
// byte-identical operators.
//
// BDF5 and BDF6 are deliberately absent. Dahlquist's second barrier says no
// A-stable linear multistep method has order above 2, so everything past BDF2
// is only A(alpha)-stable: alpha = 86.0 deg for BDF3 and 73.4 for BDF4, which a
// diffusion-dominated spectrum stays well inside, but 51.8 for BDF5 and 17.8
// for BDF6, which the reaction Jacobian of a depolarisation front does not
// respect. BDF7 and beyond are zero-unstable and do not exist.
//
// Worth knowing when choosing between the families: Crank-Nicolson is A-stable
// but NOT L-stable - its amplification factor tends to -1 rather than 0 as
// lambda*dt -> -infinity, so the stiffest modes are not damped, they ring.
// BDF2 is second order AND L-stable, so on a stiff front it is the better of
// the two second-order options.
const word implicitScheme
(
    spatialIntegrationProperties.lookupOrDefault<word>
    (
        "implicitScheme",
        "backwardEuler"
    )
);

// How the k-1 extra levels a BDF-k step needs are supplied before the first
// step, when no history exists yet.
//
//   ESDIRK3   take the first k-1 steps with ESDIRK3 and rotate the results
//             into the history. DEFAULT, and the only mode that preserves the
//             order of BDF3 and BDF4.
//   exact     evaluate the manufactured solution at t0 - dt, t0 - 2 dt, ...
//             Looks like the natural choice for a verification solver and is
//             not: it is O(tau_h dt) wrong with respect to the semi-discrete
//             system and caps the observed order at 1 (see below). Needs an
//             exact solution, which a physical solver never has.
//   constant  repeat V^0. O(dt) wrong, so it also caps the order at 1, with a
//             constant that does not shrink with the mesh.
//
// exact and constant are kept ONLY so the cost of a bad startup stays
// measurable (suite test M20), not as usable options.
//
// ESDIRK3 is the only mode here that supplies starting values consistent with
// the SEMI-DISCRETE system rather than with the PDE, which is exactly what the
// other two get wrong. It is also self-starting, so it needs no exact solution
// and ports to a physical solver unchanged - which matters, because a solver
// with no manufactured solution has no other option, and because a multistep
// method must be RESTARTED at any discontinuity in the source (a stimulus
// switching on or off invalidates the history it interpolates through).
//
// ESDIRK3 is sufficient for every BDF offered: a one-step method of order p
// leaves starting values accurate to O(dt^(p+1)), so it covers BDF up to k = 4.
// ESDIRK3 on the current binary BDF4 fits 3.89 by Cauchy self-convergence
// (route A of suite test M02) and 3.94 against the analytic solution (route B).
//
// The solver REFUSES a startup of insufficient order rather than running it:
// that failure produces no error and no warning, only a quietly capped order,
// which is the trap the exact/constant modes fell into and cost a full
// diagnosis to find.
//
// Why that distinction decides the observed order, measured rather than
// argued: a k-step method needs starting values accurate to O(dt^k). Seeding
// from the exact PDE solution is accurate to O(tau_h dt), because the exact
// solution of the PDE is NOT the solution of the semi-discrete system and the
// mismatch scales with the span of the history. That pins the global order at
// 1 with a coefficient proportional to the spatial truncation error, and the
// cap lifts only when tau_h <~ dt^(k-1).
//
// Measured on a 2-D hexa MMS case by suite test M20, on BDF3, isolating the
// startup as e(dt) = ||V_startup(dt) - V_ESDIRK3(dt)|| at the same dt and the
// same mesh, so that BDF3's own truncation error cancels and what is left is
// what the startup introduced:
//
//   exact      e ~ dt^1.01 at N=40 and dt^1.13 at N=80, with e falling 15.5x
//              (= 2^3.95) per mesh refinement - the coefficient is tau_h.
//   constant   e ~ dt^1.00 on both meshes, with e UNCHANGED by refinement
//              (ratio 1.00) - the coefficient is dVm*/dt, which does not
//              depend on h.
//
// That h-dependence is the whole difference between the two modes; at N=40
// constant's e is 1.28e4 times exact's, which is the ratio of those two
// coefficients. Judging this on raw Cauchy differences instead does NOT work:
// they mix the startup term with BDF3's own O(dt^3) error, which is identical
// on both meshes, and the mesh ratio comes out diluted.
//
// A ramp over BDF orders (BDF1, then BDF2, ...) is deliberately NOT offered:
// it changes a_0 at each of the first k-1 steps, so it needs AImplicit rebuilt
// and every cached factorisation invalidated, and its own first step is only
// O(dt^2) accurate anyway. ESDIRK3 costs one extra constant matrix and is
// third order, so it dominates that option.
const word bdfStartup
(
    spatialIntegrationProperties.lookupOrDefault<word>
    (
        "bdfStartup",
        "ESDIRK3"
    )
);

// Which manufactured solution the run is verified against.
//
//   sines    (default) V = sqrt(1+t) cos(pi x) cos(2 pi y) ..., with the
//            states built on the polynomial G. Smooth and non-polynomial, so
//            no reconstruction here is exact on it and the measured SPATIAL
//            order is the scheme's.
//   uniform  every spatial factor is 1, so V = sqrt(1+t) and the states are
//            functions of t alone. Every reconstruction, quadrature and flux
//            in this solver is exact on a constant, so the spatial error is
//            zero to roundoff and the error against the exact solution is
//            purely TEMPORAL. That is what makes a temporal order measurable
//            against the analytic solution rather than against a second
//            numerical one.
//
// The catch, and the reason sines stays the default: a constant has zero
// Laplacian, so under uniform the diffusion operator contributes nothing and
// only the time integration of the reaction system is being verified.
// Diffusion is covered on the sines profile by Cauchy self-convergence
// (route A of suite test M02); the uniform profile is its route B.
//
// Why not a profile that is exact AND diffusive: it would have to be a
// discrete eigenfunction of the LRE operator, which no closed form gives, and
// the zeroGradient boundaries rule out the linear profile that would
// otherwise work.
const word manufacturedSolution
(
    spatialIntegrationProperties.lookupOrDefault<word>
    (
        "manufacturedSolution",
        "sines"
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
        "consistent"
    )
);

const dictionary implicitHOCoeffs
(
    spatialIntegrationProperties.subDict("implicitHOCoeffs")
);

// --------------------------------------------------------------------
// Linear solver for the assembled branches (Picard / diagonalIion)
// --------------------------------------------------------------------
// Solver name in backend-neutral spelling, exactly one of GMRES, BiCGSTAB and
// SparseLU. Every option key below has a single accepted spelling, compared
// case-sensitively, and is validated here at start-up - see requireOneOf.
const word implicitLinearSolver
(
    implicitHOCoeffs.lookupOrDefault<word>("linearSolver", "GMRES")
);

requireOneOf
(
    "linearSolver", implicitLinearSolver, {"GMRES", "BiCGSTAB", "SparseLU"}
);

// Eigen or PETSc. Only PETSc works at more than one rank: Eigen has no
// distributed matrix, so its path is serial by construction.
const word linearSolverBackend
(
    implicitHOCoeffs.lookupOrDefault<word>("linearSolverBackend", "PETSc")
);

requireOneOf("linearSolverBackend", linearSolverBackend, {"PETSc", "Eigen"});

// GMRES exists only on the PETSc side, and it is the default of linearSolver. 
// Without this a case that switched the backend to Eigen and nothing else ran 
// through the whole set-up and died in the first linear solve.
if (linearSolverBackend == "Eigen" && implicitLinearSolver == "GMRES")
{
    FatalErrorInFunction
        << "linearSolver GMRES is not available on the Eigen backend." << nl
        << "With linearSolverBackend Eigen use BiCGSTAB or SparseLU."
        << exit(FatalError);
}

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

// PETSc KSP of the assembled solve, in PETSc's own spelling. Its default is the
// TRANSLATION of linearSolver (GMRES -> gmres, BiCGSTAB -> bcgs, SparseLU ->
// preonly), not the raw word: the default used to be linearSolver itself, a
// name from the neutral vocabulary sitting inside a key of PETSc's. Set it only
// to choose a KSP linearSolver has no name for, which in practice means cg.
const word petscLinearKspType
(
    implicitHOCoeffs.lookupOrDefault<word>
    (
        "petscLinearKspType",
        petscKspTypeForLinearSolver(implicitLinearSolver)
    )
);

requireOneOf
(
    "petscLinearKspType", petscLinearKspType, {"gmres", "bcgs", "cg", "preonly"}
);

// Preconditioner. NOTE for parallel comparisons: PETSc wraps ilu and lu in
// block-Jacobi at more than one rank, making them partition-dependent by
// construction. Use jacobi when comparing a parallel run against a serial
// one iteration by iteration.
const word petscLinearPcType
(
    implicitHOCoeffs.lookupOrDefault<word>("petscLinearPcType", "ilu")
);

// A closed list: the preconditioners the suite uses and tests (M07, M18). Any
// other PETSc PC can still be tried through the options database, e.g.
// -fdalinear_pc_type asm, which does not go through this dictionary.
requireOneOf
(
    "petscLinearPcType",
    petscLinearPcType,
    {"none", "jacobi", "ilu", "lu", "gamg", "hypre"}
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

// Whether PETSc's options database may override everything configured above,
// through each solve's prefix. OFF by default: the dictionary is then the only
// source, and a run cannot be redirected by something in the environment
// without leaving a trace in the case.
//
// Turning it on is also what makes PETSc's own diagnostics reachable, since
// they are options too: PETSC_OPTIONS="-fdalinear_ksp_view" prints the KSP and
// PC actually in use, which is not what the banner prints - the banner shows
// the names as requested, before the database is consulted.
const Switch petscUseOptions
(
    implicitHOCoeffs.lookupOrDefault<Switch>("petscUseOptions", false)
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
//   gaussPointODE: the ODE is integrated independently at every Iion Gauss 
//       point (updateIntegrationPointStatesODE). Cost: Nq ODE/cell.
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
// floor of the residual, i.e. by the state-ODE tolerance (stateODERelTol).
const scalar jfnkEpsilon
(
    implicitHOCoeffs.lookupOrDefault<scalar>("jfnkEpsilon", 1.0e-7)
);

// Backend for the inner Krylov solve, independent of the assembled one.
const word jfnkLinearSolverBackend
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkLinearSolverBackend", "PETSc")
);

requireOneOf
(
    "jfnkLinearSolverBackend", jfnkLinearSolverBackend, {"PETSc", "Eigen"}
);

// KSP type of the inner solve; GMRES, since the Jacobian is not symmetric.
const word jfnkPetscKspType
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkPetscKspType", "gmres")
);

requireOneOf
(
    "jfnkPetscKspType", jfnkPetscKspType, {"gmres", "bcgs", "cg", "preonly"}
);

// PETSc PC for the inner solve; usually none, because the preconditioning
// is done by the shell below.
const word jfnkPetscPcType
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkPetscPcType", "none")
);

requireOneOf
(
    "jfnkPetscPcType",
    jfnkPetscPcType,
    {"none", "jacobi", "ilu", "lu", "gamg", "hypre"}
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
//   0: x_0 = Vm_n  (constant extrapolation)
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

requireOneOf
(
    "jfnkPreconditionerKspType",
    jfnkPreconditionerKspType,
    {"gmres", "bcgs", "cg", "preonly"}
);

// See jfnkPreconditionerKspType.
const word jfnkPreconditionerPcType
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkPreconditionerPcType", "ilu")
);

requireOneOf
(
    "jfnkPreconditionerPcType",
    jfnkPreconditionerPcType,
    {"none", "jacobi", "ilu", "lu", "gamg", "hypre"}
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
// of Iion w.r.t. Vm).
const scalar diagonalIionEpsilon
(
    implicitHOCoeffs.lookupOrDefault<scalar>("diagonalIionEpsilon", 1.0e-6)
);

// Write the per-iteration nonlinear history (NonlinearConvergenceRecord).
const Switch writeNonlinearResiduals
(
    implicitHOCoeffs.lookupOrDefault<Switch>("writeNonlinearResiduals", true)
);

// Write Vm, the states and the exact/error fields at the LAST instant, to open
// them in ParaView or compare them against another run. AUTO_WRITE on the
// fields is not enough: runTime.write() only writes at a write time, and with
// the writeInterval the launchers set that never happens. See the writeNow() at
// the end of the solver.
const Switch writeFinalFields
(
    implicitHOCoeffs.lookupOrDefault<Switch>("writeFinalFields", false)
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

// how often the per-step "Time step N : t = ... , RHS Linf = ..." line is 
// printed. Was a hard-coded `nSteps % 50`.
//
// Note this also gates updateNumericalLaplacian(), because the printed
// RHS Linf cannot be formed without it - the cadence and the diagnostic are
// coupled by construction, so setting this to 1 makes every step pay for the
// Laplacian. That is the intended trade when a test needs a time-resolved
// trace; it is the wrong default, hence 50.
const label stepLogInterval
(
    max
    (
        implicitHOCoeffs.lookupOrDefault<label>("stepLogInterval", 50),
        label(1)
    )
);

// the chunk sizes of the two chunked assemblies, which used to be hard-coded 
// 50000 locals.
//
// They are capacity controls - K is assembled as sum_i K_chunk_i and
// addTripletIfNeeded filters individual contributions, not accumulated ones, so
// the chunking is transparent to the result. They are exposed because that
// claim was untestable while the value was baked into the binary: verifying it
// needs the SAME mesh assembled with more than one chunk and with exactly one,
// and with 50000 fixed any given mesh gives only one of the two. With the key,
// the regression test is a 900-cell case at chunk = 100 against the default.
//
// Note the chunk size does change the floating-point summation ORDER of K, so
// compare with a tolerance, not bitwise.
const label assemblyFaceChunk
(
    max
    (
        implicitHOCoeffs.lookupOrDefault<label>("assemblyFaceChunk", 50000),
        label(1)
    )
);

const label assemblyCellChunk
(
    max
    (
        implicitHOCoeffs.lookupOrDefault<label>("assemblyCellChunk", 50000),
        label(1)
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
// error estimate), RK4 (fixed step) or Euler. Validated here: it used to be
// checked only inside the ODE advance, i.e. on the first step of the time loop,
// after the whole set-up had run.
const word stateODESolver
(
    implicitHOCoeffs.lookupOrDefault<word>("stateODESolver", "RKF45")
);

requireOneOf("stateODESolver", stateODESolver, {"RKF45", "RK4", "Euler"});

// Which polynomial in time drives the state ODEs across a PDE step.
//
//   history  (default) the interpolant through the Vm levels the scheme has
//            already stored: the k+1 levels of BDF-k, or the stage values of
//            the ESDIRK startup. Costs nothing but the interpolation itself.
//   linear   the two-node ramp between V^n and the candidate, which is what.
//            Kept ONLY so its cost stays measurable:
//            it caps every BDF at order 2, and nothing reports it. The ramp
//            is O(dt^2) wrong pointwise, O(dt^3) once integrated over the
//            step, so O(dt^2) globally.
//
//   Suite test M21: BDF4 fits 2.04 with linear and 3.89 with history.
//
// backwardEuler, crankNicolson and BDF1 store no level older than V^n, so they
// take the two-node ramp under either value and are unaffected. So is the
// gaussPointODE state mode, which always uses the ramp: raising it would mean
// storing every level at every quadrature point.
const word stateODEDriver
(
    implicitHOCoeffs.lookupOrDefault<word>("stateODEDriver", "history")
);

// First sub-step attempted; the adaptive controller takes over from there.
const scalar stateODEInitialStep
(
    implicitHOCoeffs.lookupOrDefault<scalar>("stateODEInitialStep", 1.0e-5)
);

// Absolute and relative tolerances of the adaptive sub-stepping. These set
// the NOISE FLOOR of the nonlinear residual, which is what limits how small
// jfnkEpsilon can usefully be and how tight nonlinearTolerance can be
// asked to be.
const scalar stateODEAbsTol
(
    implicitHOCoeffs.lookupOrDefault<scalar>("stateODEAbsTol", SMALL)
);

// See stateODEAbsTol.
const scalar stateODERelTol
(
    implicitHOCoeffs.lookupOrDefault<scalar>("stateODERelTol", 1.0e-4)
);

// Cap on sub-steps per PDE step, so a stiff point cannot hang the run.
const label stateODEMaxSteps
(
    implicitHOCoeffs.lookupOrDefault<label>("stateODEMaxSteps", 1000000000)
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
// High-order reconstruction stencils for Vm and Iion.
// ----------------------------------------------------------------------- //
// How each boundary patch enters the reconstruction stencils, classified
// from Vm's boundary condition:
//   fixedValue / fixedVoltage -> the boundary value is added to the stencil
//                  as a Dirichlet ghost point (includePatchIn{Vm,Iion}Stencils),
//                  which keeps the reconstruction exact up to the boundary for
//                  the manufactured solution;
//   anything else -> not added. On a zeroGradient patch the boundary faces use
//                  a one-sided stencil, and the Neumann condition enters as a
//                  zero boundary flux: the finite-volume treatment, and exact
//                  here, since the manufactured solution has zero normal
//                  derivative on every boundary of the unit square.
//
// No mirrored stencil is built for zeroGradient patches. 
boolList includePatchInVmStencils(mesh.boundaryMesh().size(), false);
boolList includePatchInIionStencils(mesh.boundaryMesh().size(), false);

forAll(mesh.boundaryMesh(), patchI)
{
    const word bcType = Vm.boundaryField()[patchI].type();

    if
    (
        bcType == fixedValueFvPatchScalarField::typeName
     || bcType == "fixedVoltage"
    )
    {
        includePatchInVmStencils[patchI] = true;
        includePatchInIionStencils[patchI] = true;
    }
}

// Sub-dictionary holding the stabilisation weight and the three LRE stencil
// definitions (LRECoeffs_Vm, _Iion and _states).
const dictionary& highOrderCoeffs =
    spatialIntegrationProperties.subDict("highOrderCoeffs");

dictionary lreCoeffsVm(highOrderCoeffs.subDict("LRECoeffs_Vm"));

// Stencil parameters echoed to the log so a run records what it was built
// with. Only N, Nn and k reach the library; maxStencilSize is inert, the adapter
// drops it.
//
//   N   the polynomial DEGREE of the local reconstruction.
//   Nn  the stencil cells ADDED on top of the minimum that degree needs, NOT
//       the total. The minimum is the number of terms in the Taylor polynomial,
//       (N+1)(N+2)/2 in 2-D and (N+1)(N+2)(N+3)/6 in 3-D, which is the number of
//       unknowns the local least-squares fit determines. So
//
//           stencil size = minimum + Nn
//
//       and Nn is the redundancy of that fit. With N 3 in 2-D the minimum is 10,
//       so Nn 10 gives a 20-cell stencil, two points per unknown.
//   k   exponent of the radial weight function.
const label lreN = lreCoeffsVm.lookupOrDefault<label>("N", -1);
const label lreNn = lreCoeffsVm.lookupOrDefault<label>("Nn", -1);
const label lreMaxStencil = lreCoeffsVm.lookupOrDefault<label>("maxStencilSize", -1);
const scalar lreK = lreCoeffsVm.lookupOrDefault<scalar>("k", -1.0);

// Iion quadrature.
const word IionCoeffsName("LRECoeffs_Iion");
if (!highOrderCoeffs.found(IionCoeffsName))
{
    FatalErrorInFunction
        << "highOrderCoeffs does not contain LRECoeffs_Iion (required for "
        << "the Iion quadrature)." << nl
        << exit(FatalError);
}

dictionary lreCoeffsIion(highOrderCoeffs.subDict(IionCoeffsName));
// Same echo for the Iion interpolator.
const label lreIionN = lreCoeffsIion.lookupOrDefault<label>("N", -1);
const label lreIionNn = lreCoeffsIion.lookupOrDefault<label>("Nn", -1);
const label lreIionMaxStencil =
    lreCoeffsIion.lookupOrDefault<label>("maxStencilSize", -1);
const scalar lreIionK = lreCoeffsIion.lookupOrDefault<scalar>("k", -1.0);

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

// There is deliberately NO key for the non-orthogonal correction of the
// low-order operator. It is always on.
//
// A bare two-point flux resolves only the part of grad(Vm) along the line joining
// the two cell centres, so where that line is not parallel to the face normal it
// drops the tangential remainder and the operator is INCONSISTENT: order 0.81 on
// unstructured triangles instead of 2, with an error 230x larger than the same run
// on hexahedra. No study wants that, so it is not selectable - the alternative
// would be a switch whose only effect is to make the solver wrong. On an
// orthogonal mesh with isotropic conductivity the correction is identically zero
// and costs nothing, so there is no performance argument for a switch either.
// See nonOrthFluxVector() in the solver for the decomposition.

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
// otherwise left null so the gaussPointODE path pays no extra cost.
// ----------------------------------------------------------------------- //
const bool buildStateInterp =
    reconstructStatesFromCellCentres
 && useHighOrder_Iion
 && mesh.nGeometricD() > 1;

// Build the include-list for the state stencils from u1's boundary types, with
// the same classification as Vm: a fixedValue patch adds its boundary value to
// the stencil, anything else does not. In the MMS set-up the three states share
// fixedValue boundary conditions.
boolList includePatchInStateStencils(mesh.boundaryMesh().size(), false);
forAll(mesh.boundaryMesh(), patchI)
{
    const word bcType = u1.boundaryField()[patchI].type();

    if
    (
        bcType == fixedValueFvPatchScalarField::typeName
     || bcType == "fixedVoltage"
    )
    {
        includePatchInStateStencils[patchI] = true;
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
            << exit(FatalError);
    }

    dictionary lreCoeffsStates(highOrderCoeffs.subDict("LRECoeffs_states"));

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
