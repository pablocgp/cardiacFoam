/*---------------------------------------------------------------------------*\
    cardiacFoam - shared code of the monodomain solvers
    (highOrderElectroActivationFoamImplicitPETScDistributed and
    highOrderManufacturedFDAImplicitPETScDistributed)

    Time schemes: theta family and BDF1-4 coefficients, the
    ESDIRK3 tableau, and the Vm driver of the state ODEs.

    INCLUDED INSIDE THE SOLVER'S ANONYMOUS NAMESPACE, after its system,
    OpenFOAM, Eigen and PETSc includes, through monodomainCommon.H. It
    therefore includes nothing itself and has internal linkage in each
    executable. Code here is used verbatim by both solvers: a change made here
    changes both, which is the point.
\*---------------------------------------------------------------------------*/

    // Added for cardiacFoam, ported from the MMS solver: coefficients of the
    // time scheme. With M the mass operator and L the scaled diffusion
    // operator, the two families are
    //
    //   theta family   M (V^{n+1} - V^n)/dt
    //                    = L (theta V^{n+1} + (1-theta) V^n)
    //                    + theta s^{n+1} + (1-theta) s^n
    //
    //   BDF-k          (M/dt) sum_{j=0}^{k} a_j V^{n+1-j}  =  L V^{n+1} + s^{n+1}
    //
    //     A = (a_0/dt) M - theta L        (theta = 1 for BDF)
    //     B = M/dt + (1-theta) L          (= M/dt for BDF)
    //
    // and the step solves A V^{n+1} = B H + source, where the history vector H
    // is V^n for the theta family and -sum_{j>=1} a_j V^{n+1-j} for BDF. So a
    // BDF step is still ONE mass-matrix product: the combination of the stored
    // levels is formed first, as a field. BDF1 and backwardEuler are the same
    // scheme: a_0 = 1, a_1 = -1 gives H = V^n and A = M/dt - L, bit for bit.
    struct TimeSchemeCoeffs
    {
        scalar theta = 1.0;
        //- 0 for the theta family, k for BDF-k
        label bdfOrder = 0;
        //- a_0 .. a_k, trailing entries zero. a_0 is what scales M in A.
        FixedList<scalar, 5> a;

        TimeSchemeCoeffs()
        {
            a[0] = 1.0;
            a[1] = -1.0;
            a[2] = 0.0;
            a[3] = 0.0;
            a[4] = 0.0;
        }
    };

    // Added for cardiacFoam, ported from the MMS solver: Butcher tableau of an
    // ESDIRK - Explicit first stage (raises the stage order to 2, against order
    // reduction on stiff problems), Singly Diagonally Implicit (every implicit
    // stage shares gamma, so ONE matrix M/(gamma dt) - L serves all of them) and
    // stiffly accurate (b equals the last row, so V^{n+1} = Y_s and the method is
    // L-stable). Here it only supplies what BDF-k cannot do by itself: the first
    // k-1 steps, the k-1 steps after a restart, and a step of non-nominal length.
    struct EsdirkTableau
    {
        label nStages = 0;
        label order = 0;
        scalar gamma = 0.0;
        //- Lower-triangular including the diagonal; entries above it are zero
        FixedList<FixedList<scalar, 6>, 6> a;
        //- Abscissae. c[nStages - 1] == 1 for a stiffly accurate method
        FixedList<scalar, 6> c;
    };

    // ESDIRK3(2)4L[2]SA: four stages, classical order 3, L-stable, stiffly
    // accurate. gamma is the root of gamma^3 - 3 gamma^2 + (3/2) gamma - 1/6 = 0
    // in (1/4, 1/2). The coefficients were obtained on the MMS side by solving
    // the order conditions for c = [0, 2 gamma, 3/5, 1], row sums equal to c
    // and b equal to the last row, and checked there: the four third-order
    // conditions hold to ~1e-41 and the fourth-order one misses by 1.4e-3, so
    // the method is exactly third order. Copied verbatim.
    EsdirkTableau esdirk3Tableau()
    {
        EsdirkTableau tb;
        tb.nStages = 4;
        tb.order = 3;
        tb.gamma = 0.435866521508458999416;

        forAll(tb.c, i)
        {
            tb.c[i] = 0.0;
            forAll(tb.a[i], j)
            {
                tb.a[i][j] = 0.0;
            }
        }

        const scalar g = tb.gamma;

        tb.c[0] = 0.0;
        tb.c[1] = 2.0*g;
        tb.c[2] = 0.6;
        tb.c[3] = 1.0;

        tb.a[1][0] = g;
        tb.a[1][1] = g;

        tb.a[2][0] =  0.257648246066427245800;
        tb.a[2][1] = -0.093514767574886245216;
        tb.a[2][2] =  g;

        tb.a[3][0] =  0.187641024346723825161;
        tb.a[3][1] = -0.595297473576954948048;
        tb.a[3][2] =  0.971789927721772123471;
        tb.a[3][3] =  g;

        return tb;
    }

    // Coefficients of the requested scheme. An unknown name is a fatal
    // dictionary error rather than a silent default.
    TimeSchemeCoeffs timeSchemeCoeffs(const word& scheme)
    {
        TimeSchemeCoeffs c;

        if (scheme == "backwardEuler")
        {
            c.theta = 1.0;
            c.bdfOrder = 0;
            return c;
        }

        if (scheme == "crankNicolson")
        {
            c.theta = 0.5;
            c.bdfOrder = 0;
            return c;
        }

        // BDF is fully implicit in K, so theta = 1, and the source is wanted at
        // t^{n+1} only - one source evaluation per step FEWER than
        // Crank-Nicolson needs.
        c.theta = 1.0;

        if (scheme == "BDF1")
        {
            c.bdfOrder = 1;
            c.a[0] =  1.0;
            c.a[1] = -1.0;
        }
        else if (scheme == "BDF2")
        {
            c.bdfOrder = 2;
            c.a[0] =  3.0/2.0;
            c.a[1] = -2.0;
            c.a[2] =  1.0/2.0;
        }
        else if (scheme == "BDF3")
        {
            c.bdfOrder = 3;
            c.a[0] =  11.0/6.0;
            c.a[1] = -3.0;
            c.a[2] =  3.0/2.0;
            c.a[3] = -1.0/3.0;
        }
        else if (scheme == "BDF4")
        {
            c.bdfOrder = 4;
            c.a[0] =  25.0/12.0;
            c.a[1] = -4.0;
            c.a[2] =  3.0;
            c.a[3] = -4.0/3.0;
            c.a[4] =  1.0/4.0;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown implicitScheme: " << scheme << nl
                << "Valid options are backwardEuler, crankNicolson, BDF1,"
                << " BDF2, BDF3 and BDF4." << nl
                << "ESDIRK is not a scheme here: it exists only to supply the"
                << " starting values BDF needs, via the bdfStartup key." << nl
                << "BDF5 and BDF6 are deliberately not offered: past BDF2 no"
                << " linear multistep method is A-stable (Dahlquist), and the"
                << " A(alpha) wedge closes from 86.0 deg at BDF3 and 73.4 at"
                << " BDF4 to 51.8 at BDF5 and 17.8 at BDF6. BDF7 and beyond"
                << " are zero-unstable."
                << exit(FatalError);
        }

        return c;
    }

    // Added for cardiacFoam, ported from the MMS solver: the Vm the state ODEs
    // are driven by, at ONE point, as a set of nodes in normalised time
    // s = tau/interval over [0, 1].
    //
    // The ionic states are integrated over a PDE step with Vm held by the
    // caller, and the ACCURACY OF THAT Vm(t) CAPS THE COUPLING. This solver used
    // to hold it FROZEN at the candidate V^{n+1} for the whole step - a
    // first-order splitting of the Vm-state coupling, whatever the PDE scheme,
    // Crank-Nicolson included (E06 measured a Cauchy ratio of 0.41 ~ 2^-1.3 on
    // CN, where second order gives 0.25). A linear blend between V^n and the
    // candidate is O(dt^2); the MMS measured it capping every scheme above
    // second order (BDF4 fits 2.04 with it and 3.89 with the history
    // interpolant, ESDIRK4 1.88 against 2.99).
    //
    // nNodes == 2 with s = {0, 1} is the linear blend. BDF2-4 pass the k+1
    // levels of their history and the ESDIRK stages their stage values - both
    // already stored, so the higher-order interpolant costs no evaluation.
    // The LAST node is always the candidate, at s = 1.
    struct VmDriverNodes
    {
        label nNodes = 2;
        //- Normalised abscissae; the last one is 1
        FixedList<scalar, 6> s;
        //- Vm [V] at those abscissae
        FixedList<scalar, 6> v;
    };

    // Lagrange interpolation of the driver at normalised time s.
    inline scalar interpolateVmDriver
    (
        const VmDriverNodes& nd,
        const scalar s
    )
    {
        // Written out for two nodes: the path every one-step scheme takes.
        if (nd.nNodes <= 2)
        {
            return (1.0 - s)*nd.v[0] + s*nd.v[1];
        }

        scalar result = 0.0;

        for (label j = 0; j < nd.nNodes; ++j)
        {
            scalar Lj = 1.0;

            for (label m = 0; m < nd.nNodes; ++m)
            {
                if (m != j)
                {
                    Lj *= (s - nd.s[m])/(nd.s[j] - nd.s[m]);
                }
            }

            result += Lj*nd.v[j];
        }

        return result;
    }

    // Added for cardiacFoam, ported from the MMS solver: the nodes of the
    // driver at FIELD level - the same abscissae for every point, one field of
    // values per node, indexed like the Vm array handed to solveODE. The last
    // node is the candidate, which solveODE takes from its own Vm argument, so
    // its field pointer is left null.
    struct VmDriverStages
    {
        label nNodes = 2;
        FixedList<scalar, 6> s;
        FixedList<const scalarField*, 6> field;

        VmDriverStages()
        {
            forAll(s, i)
            {
                s[i] = 0.0;
                field[i] = nullptr;
            }
            s[1] = 1.0;
        }
    };
