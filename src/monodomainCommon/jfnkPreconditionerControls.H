/*---------------------------------------------------------------------------*\
    cardiacFoam - shared code of the monodomain solvers

    Dictionary controls of the JFNK preconditioner. Included from each
    solver's createFields.H, inside main(), where the dictionary
    `implicitHOCoeffs` (spatialIntegrationProperties/implicitHOCoeffs) is in
    scope. Same keys, same defaults and same validation in both solvers.
\*---------------------------------------------------------------------------*/

// JFNK preconditioner (applied on the RIGHT, see PetscShellKspSolver).
// The default uses the same local source
// linearisation as diagonalIion, assembled as
// P ~= AImplicit - theta*d(source)/d(Vm).
const word jfnkPreconditioner
(
    implicitHOCoeffs.lookupOrDefault<word>("jfnkPreconditioner", "diagonalIion")
);

requireOneOf
(
    "jfnkPreconditioner", jfnkPreconditioner, {"none", "AImplicit", "diagonalIion"}
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
