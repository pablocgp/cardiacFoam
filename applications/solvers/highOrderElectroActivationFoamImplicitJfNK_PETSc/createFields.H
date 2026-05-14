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

// KSP type: GMRES, FGMRES, BiCGSTAB (BCGS), PREONLY (use the PC alone,
// e.g. PC=LU for a direct solve). Default GMRES for general non-symmetric
// matrices typical of theta-scheme operators.
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

// PETSc preconditioner: ILU, JACOBI, BJACOBI, LU (sparse direct), NONE, SOR, ASM.
// PCILU is a strong default for serial structured/unstructured meshes.
// For the legacy SparseLU behaviour use linearSolver=PREONLY + preconditioner=LU.
const word implicitPreconditioner
(
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<word>
    (
        "preconditioner",
        "ILU"
    )
  : word("ILU")
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
// we clamp; this protects Newton from initialising in unphysical
// territory (e.g. Vm = +1 V because dVm/dt at the previous step was
// dominated by I_Na rather than the smooth diffusion baseline).
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
// valid range (e.g. -107 mV or +1 V), where the ionic model is no longer
// calibrated and RKF45 stalls trying to integrate ultra-stiff dynamics.
//
// This option clamps Vm to [jfnkInitGuessVmMin, jfnkInitGuessVmMax]
// *only* on the path going into ionicModel->solveODE(...). The clamp
// affects the residual F(x) outside the physical bracket (introducing a
// kink at the boundary), but does not change the actual fixed-point
// solution which is well inside the bracket. The kink makes F slightly
// non-smooth in the unphysical region, which Newton/GMRES handle fine
// via line search.
//
// Effect: failed JFNK probes no longer drag the ionic ODE into a
// pathological regime; the worst case is that Newton sees a flat F there
// and rejects the step via Armijo. Net: large dt is much more robust.
//
// Reuses the same bracket as jfnkInitGuessVmMin/Max to keep tuning
// minimal. Disable by setting jfnkClampODEInput=false.
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
//
// This prevents Newton from overshooting during the cardiac AP upstroke
// (where the linearisation at Vm_old fails to capture the I_Na sigmoid),
// allowing convergence at much larger time steps (dt ~ 1e-3 s).
//
// Set jfnkLineSearch to false to fall back to the legacy fixed-damping
// update for reproducibility with earlier runs.
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

// Maximum number of backtracking trials per Newton iteration before
// accepting the current alpha unconditionally.
const label jfnkLineSearchMaxIter =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<label>
    (
        "jfnkLineSearchMaxIter",
        10
    )
  : label(10);

// Floor on alpha. If backtracking reaches this value without satisfying
// the Armijo condition we accept the current iterate and let the outer
// Newton loop / time-step rollback handle the failure.
const scalar jfnkLineSearchAlphaMin =
    spatialIntegrationProperties.found("implicitHOCoeffs")
  ? spatialIntegrationProperties.subDict("implicitHOCoeffs").lookupOrDefault<scalar>
    (
        "jfnkLineSearchAlphaMin",
        1e-4
    )
  : scalar(1e-4);

// Sufficient-decrease constant c in the Armijo condition. The textbook
// default is 1e-4; smaller c is more permissive (fewer backtracks),
// larger c is stricter.
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

autoPtr<LRE> LREInterp_VmPtr;

if (useHighOrder_Vm)
{
    if
    (
        !spatialIntegrationProperties.found("highOrderCoeffs")
     || !spatialIntegrationProperties.subDict("highOrderCoeffs").found("LRECoeffs_Vm")
    )
    {
        FatalErrorInFunction
            << "useHighOrder_Vm is true, but constant/spatialIntegrationProperties "
            << "does not contain highOrderCoeffs/LRECoeffs_Vm." << nl
            << abort(FatalError);
    }

    boolList includePatchInVmStencils(mesh.boundaryMesh().size(), false);
    forAll(includePatchInVmStencils, patchI)
    {
        const word bcType = Vm.boundaryField()[patchI].type();

        if
        (
            bcType == fixedValueFvPatchScalarField::typeName
         || bcType == "fixedVoltage"
        )
        {
            includePatchInVmStencils[patchI] = true;
        }
    }

    LREInterp_VmPtr.reset
    (
        new LRE
        (
            mesh,
            includePatchInVmStencils,
            spatialIntegrationProperties.subDict("highOrderCoeffs").subDict("LRECoeffs_Vm")
        )
    );
}

autoPtr<LRE> LREInterp_IionPtr;

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
                << "useHighOrder_Iion is true, but highOrderCoeffs contains "
                << "neither LRECoeffs_Iion nor legacy LRECoeffs_states." << nl
                << abort(FatalError);
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

    LREInterp_IionPtr.reset
    (
        new LRE
        (
            mesh,
            includePatchInIionStencils,
            highOrderCoeffs.subDict(IionCoeffsName)
        )
    );
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

autoPtr<ionicModel> ionicModel =
    ionicModel::New
    (
        solutionVariablesMemory,
        totalIionIntegrationPoints,
        runTime.deltaTValue()
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
