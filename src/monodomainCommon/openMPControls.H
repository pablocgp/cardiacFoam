/*---------------------------------------------------------------------------*\
    cardiacFoam - shared code of the monodomain solvers

    OpenMP controls. Included from each solver's createFields.H, inside
    main(), where the dictionary `implicitHOCoeffs` is in scope. The same keys
    drive the same threaded regions in both solvers: the ionic-state ODEs (at
    cells and at quadrature points) and the reconstructions at the quadrature
    points (Vm and the states), plus Eigen's own threaded kernels.
    stateODENumThreads sets the thread count of all of them.
\*---------------------------------------------------------------------------*/

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

#ifdef _OPENMP
if (stateODENumThreads > 0)
{
    omp_set_num_threads(stateODENumThreads);
}

// Eigen parallelises some of its own kernels with OpenMP when the code is
// compiled with it - row-major sparse * dense products above a size threshold,
// i.e. the operator applications on large meshes (K*x, M*x, the JFNK mat-vec)
// - and sizes them from OMP_NUM_THREADS, not from the switch above. Tie Eigen
// to the same switch: false means serial; true follows stateODENumThreads
// (0 = OpenMP's default). (The extra CPU of E33's serial control turned out to
// be OpenBLAS threads under PETSc, not Eigen; that is set from the environment,
// OPENBLAS_NUM_THREADS, which the test harness fixes to 1.)
Eigen::setNbThreads(stateODEUseOpenMP ? max(stateODENumThreads, label(0)) : 1);
#endif
