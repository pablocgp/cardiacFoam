/*---------------------------------------------------------------------------*\
    cardiacFoam - shared code of the monodomain solvers
    (highOrderElectroActivationFoamImplicitPETScDistributed and
    highOrderManufacturedFDAImplicitPETScDistributed)

    Aggregator: the order below is the dependency order.

    INCLUDED INSIDE THE SOLVER'S ANONYMOUS NAMESPACE, after its system,
    OpenFOAM, Eigen and PETSc includes, through monodomainCommon.H. It
    therefore includes nothing itself and has internal linkage in each
    executable. Code here is used verbatim by both solvers: a change made here
    changes both, which is the point.
\*---------------------------------------------------------------------------*/

#include "monodomainBase.H"
#include "monodomainPetsc.H"
#include "monodomainReconstruction.H"
#include "monodomainAssembly.H"
#include "monodomainKrylov.H"
#include "monodomainTimeSchemes.H"
