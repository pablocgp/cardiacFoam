/*---------------------------------------------------------------------------*\
High-order finite-volume solver for a manufactured reaction-diffusion problem.

    Solves the monodomain equation

        chi Cm dVm/dt = div(D grad Vm) - Iion(Vm, u1, u2, u3),

    coupled to three local reaction ODEs for the states u, against a
    MANUFACTURED solution. Every field has a closed-form exact value at every
    instant, so the error is known rather than estimated and the observed order
    of accuracy can be FITTED. That is what this solver exists for: it measures
    orders.

    High order in space. The diffusion operator, the consistent mass matrix, the
    quadrature of the ionic source and the reconstruction of the states at the
    quadrature points each use a moving-least-squares reconstruction whose degree
    is chosen INDEPENDENTLY. The observed order is the minimum over those three
    ingredients, which is why they are separate keys: the point of the sweeps is
    to find which one is binding.

    High order in time. Two families share one pair of operators: the theta
    methods (backwardEuler, crankNicolson) and BDF1-4. A BDF-k run needs k-1
    starting values accurate with respect to the SEMI-DISCRETE system, which is
    what the ESDIRK3 startup supplies; seeding them from the exact solution caps
    the order at one without any error being reported. The state ODEs are driven
    across the step by the polynomial through the levels the scheme already
    stores, because a two-node ramp would cap the whole solver at order two.

    The nonlinear coupling of Vm and Iion is solved by Picard, by a diagonal
    linearisation of the ionic current, or by JFNK; all three converge to the
    same solution and differ in cost, which is one of the comparisons the solver
    is built to make.

    The linear algebra is distributed: matrices are MPIAIJ on PETSC_COMM_WORLD,
    with each rank owning the rows of its cells and columns numbered globally,
    since a high-order stencil reaches into other ranks. Matrix-vector products
    go through PETSc. A serial Eigen backend is selectable for a single rank.
\*---------------------------------------------------------------------------*/

#include <cmath>
#include <petscksp.h>
#include "fvCFD.H"
#include "highOrderInterp.H"
// Needed by the coupled-patch (processor face)
// stabilisation exchange in the stiffness assembly.
#include "globalIndex.H"
#include "syncTools.H"
#include "PstreamBuffers.H"
#include "processorFvPatch.H"
#include <chrono>
#include <algorithm>
#include <cctype>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>
#include <functional>
#include <sys/resource.h>
#include <string>
#include <vector>
#include <initializer_list>
#ifdef __GLIBC__
#include <malloc.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{

// Shared code of the two monodomain solvers (src/monodomainCommon). Included
// here, inside the anonymous namespace, so it has internal linkage.
#include "monodomainCommon.H"
    // Which manufactured solution is in force.
    //
    //   mmsSines    the default, V = sqrt(1+t) cos(pi x) cos(2 pi y) ...
    //   mmsUniform  V = sqrt(1+t), with every spatial factor equal to one.
    //
    // The uniform profile exists so that temporal order can be measured
    // against the ANALYTIC solution instead of against a second numerical one.
    // Every reconstruction, quadrature and flux in this solver is exact on a
    // uniform field, so the spatial error is zero to roundoff and the whole
    // remaining error is temporal. 
    //
    // Laplacian of a uniform field is zero, so the uniform profile leaves the
    // DIFFUSION operator inert and verifies the time integration of the
    // reaction system alone. Diffusion is covered by the Cauchy self-
    // convergence test on the sines profile. A profile that is both exactly
    // representable AND has non-zero diffusion would have to be a discrete
    // eigenfunction of the LRE operator, which no closed form gives; the
    // boundary conditions are zeroGradient, which also rules out a linear one.
    enum MmsProfile { mmsSines, mmsUniform };

    // Set once in main() from the dictionary, before the first exact value is
    // evaluated, and read-only afterwards - including from the OpenMP loops.
    MmsProfile mmsProfile = mmsSines;


    // Which manufactured field a generic routine is acting on. u3 is absent
    // because its exact solution is identically zero.
    enum ManufacturedFieldID
    {
        mfVm = 0,
        mfU1 = 1,
        mfU2 = 2
    };

    // Error norms of one field against its exact solution, at the CELL CENTRES:
    // volume-weighted L1 and L2, Linf, and the unweighted cell-error norm that
    // the convergence driver fits its rates against.
    struct FieldErrorSummary
    {
        scalar L1;
        scalar L2;
        scalar Linf;
        scalar errorCell;
    };

    // Error norms of one field at the QUADRATURE POINTS, i.e. of the high-order
    // reconstruction rather than of the cell average. normL2 is the norm of the
    // exact field itself and relL2 the resulting percentage, which is what
    // makes the p-refinement comparable across fields of different magnitude.
    struct ReconstructionErrorSummary
    {
        scalar L1;
        scalar L2;
        scalar Linf;
        scalar normL2;
        scalar relL2;
    };

    // One row of the per-time-step nonlinear history: iteration counts, the
    // coupled and per-field residuals, whether the step converged or was rolled
    // back, and the linear-solver and line-search work it took. Written to
    // nonlinearResiduals.dat when writeNonlinearResiduals is on, and the
    // primary evidence when a run has to be compared against a serial one
    // iteration by iteration.
    struct NonlinearConvergenceRecord
    {
        scalar time;
        label step;
        label iterations;
        label maxIterations;
        bool converged;
        bool rolledBack;
        scalar coupledResidual;
        scalar VmResidual;
        scalar u1Residual;
        scalar u2Residual;
        scalar u3Residual;
        scalar maxStateResidual;
        scalar IionResidual;
        label linearIterations;
        scalar linearError;
        label lineSearchIterations;
    };

    // Spatial factor of the manufactured Vm: a product of sines whose
    // wavenumbers increase with the coordinate index, so that no direction is
    // degenerate and the operator is exercised anisotropically. Chosen smooth
    // and non-polynomial, so a p3 reconstruction cannot be exact and the
    // measured order is the scheme's, not the test function's.
    scalar computeF(const point& p, const label dim)
    {
        // Exactly representable by every reconstruction here, so the spatial
        // error vanishes and what is left is the temporal error.
        if (mmsProfile == mmsUniform)
        {
            return 1.0;
        }

        const scalar pi = constant::mathematical::pi;

        if (dim == 1)
        {
            return std::cos(pi*p.x());
        }
        else if (dim == 2)
        {
            return std::cos(pi*p.x())*std::cos(2.0*pi*p.y());
        }

        return
            std::cos(pi*p.x())
           *std::cos(2.0*pi*p.y())
           *std::cos(3.0*pi*p.z());
    }

    // Spatial factor of the manufactured states, polynomial and bounded away
    // from zero (so the square root in u2 is safe).
    scalar computeG(const point& p, const label dim)
    {
        // Uniform too, and it has to be: u2 = 1/((1+t) sqrt(G)) is not a
        // polynomial, so a varying G would put a spatial reconstruction error
        // back into the states and defeat the whole point of the mode.
        if (mmsProfile == mmsUniform)
        {
            return 1.0;
        }

        if (dim == 1)
        {
            return 1.0 + p.x();
        }
        else if (dim == 2)
        {
            return 1.0 + p.x()*sqr(p.y());
        }

        return 1.0 + p.x()*sqr(p.y())*pow3(p.z());
    }

    // Manufactured solution for Vm: sqrt(1+t) * F(x). The square root gives a
    // time dependence that no theta scheme integrates exactly, so the temporal
    // order is measurable too.
    scalar exactVm(const point& p, const scalar t, const label dim)
    {
        return std::sqrt(1.0 + t)*computeF(p, dim);
    }

    // Manufactured u1, deliberately coupled to Vm so the state ODEs cannot be
    // satisfied independently of the PDE.
    scalar exactU1(const point& p, const scalar t, const label dim)
    {
        return (1.0 + t)*computeG(p, dim) + exactVm(p, t, dim);
    }

    // Manufactured u2. Decays like 1/(1+t) and stays positive.
    scalar exactU2(const point& p, const scalar t, const label dim)
    {
        return 1.0/((1.0 + t)*std::sqrt(computeG(p, dim)));
    }

    // Manufactured u3 is identically zero: the third state is carried so that
    // the code path matches the physiological solver's, not to be verified.
    scalar exactU3(const point&, const scalar, const label)
    {
        return 0.0;
    }

    // Exact value of the field named by a ManufacturedFieldID, so the error and
    // boundary routines can be written once for all fields.
    scalar exactFieldValue
    (
        const ManufacturedFieldID fieldID,
        const point& p,
        const scalar t,
        const label dim
    )
    {
        switch (fieldID)
        {
            case mfVm:
                return exactVm(p, t, dim);
            case mfU1:
                return exactU1(p, t, dim);
            case mfU2:
                return exactU2(p, t, dim);
            default:
                return 0.0;
        }
    }

    // Right-hand side of the manufactured state ODEs, du/dt = R(V, u). These
    // are the rates the SOLVER integrates; the manufactured source that makes
    // the exact solution satisfy them is added separately.
    void reactionRates
    (
        const scalar V,
        const scalar u1,
        const scalar u2,
        const scalar u3,
        scalar& du1dt,
        scalar& du2dt,
        scalar& du3dt
    )
    {
        const scalar a = u1 + u3 - V;
        const scalar b = V - u3;

        du1dt = sqr(a)*sqr(u2) + 0.5*a*sqr(u2)*b;
        du2dt = -a*pow(u2, 3);
        du3dt = 0.0;
    }

    // Manufactured ionic current: a nonlinear surrogate of the physiological
    // Iion, cubic in the states and linear in (V - u3). It is nonlinear in Vm,
    // which is what makes the Picard / diagonalIion / JFNK distinction
    // meaningful on this problem.
    scalar ionicCurrentPDE
    (
        const scalar V,
        const scalar u1,
        const scalar u2,
        const scalar u3,
        const scalar beta,
        const scalar chiVal,
        const scalar CmVal
    )
    {
        return
            -0.5*(u1 + u3 - V)*sqr(u2)*(V - u3)
          + (beta/(chiVal*CmVal))*(V - u3);
    }

    // Source term of the Vm equation, i.e. -Iion. Kept as its own function
    // because the sign convention between "ionic current" and "source" is the
    // easiest thing to get wrong when comparing against the physiological
    // solver.
    scalar vmSourcePDE
    (
        const scalar V,
        const scalar u1,
        const scalar u2,
        const scalar u3,
        const scalar beta,
        const scalar chiVal,
        const scalar CmVal
    )
    {
        return -ionicCurrentPDE(V, u1, u2, u3, beta, chiVal, CmVal);
    }

    // Coefficient that makes the manufactured Vm an exact solution of the
    // diffusion operator: for the sine product of computeF, the Laplacian
    // contributes -pi^2 (D_xx + 4 D_yy + 9 D_zz), the wavenumbers squared. Read
    // from cell 0 because the manufactured setup uses a uniform conductivity.
    scalar computeBeta(const volTensorField& conductivity, const label dim)
    {
        // A constant has zero Laplacian, so there is nothing for beta to
        // cancel and the ionic current reduces to the term that supplies
        // dV/dt.
        if (mmsProfile == mmsUniform)
        {
            return 0.0;
        }

        const tensor& D = conductivity[0];
        const scalar pi2 = sqr(constant::mathematical::pi);

        if (dim == 1)
        {
            return -pi2*D.xx();
        }
        else if (dim == 2)
        {
            return -pi2*(D.xx() + 4.0*D.yy());
        }

        return -pi2*(D.xx() + 4.0*D.yy() + 9.0*D.zz());
    }

    // Mesh resolution N recovered from the global cell count, only for
    // labelling output. Reduced, so it does not depend on the decomposition.
    label estimatedN(const fvMesh& mesh)
    {
        const label nCellsGlobal =
            returnReduce(mesh.nCells(), sumOp<label>());

        const scalar dim = max(scalar(mesh.nGeometricD()), 1.0);
        return label(std::round(std::pow(scalar(nCellsGlobal), 1.0/dim)));
    }

    // Representative cell size h = (active domain measure / global cell
    // count)^(1/dim). This is the h the convergence driver fits its orders
    // against, so it must be a property of the mesh and not of the
    // decomposition - hence the reduced bounding box and cell count below.
    scalar characteristicDx(const fvMesh& mesh)
    {
        // Mesh.bounds() is the globally reduced
        // bounding box; boundBox(mesh.points()) is this rank's own points only.
        // The cell count below was already reduced, so mixing the two gave a
        // domain measure that shrank with the rank count and an h that did not
        // match the mesh. run_convergence.py fits the convergence order against
        // this h, so the error would have surfaced as a wrong order rather than
        // as an obviously wrong number.
        const boundBox& bb = mesh.bounds();
        const vector span = bb.max() - bb.min();

        scalar activeMeasure = 1.0;
        if (mesh.nGeometricD() >= 1)
        {
            activeMeasure *= max(span.x(), SMALL);
        }
        if (mesh.nGeometricD() >= 2)
        {
            activeMeasure *= max(span.y(), SMALL);
        }
        if (mesh.nGeometricD() >= 3)
        {
            activeMeasure *= max(span.z(), SMALL);
        }

        const scalar nCellsGlobal =
            returnReduce(mesh.nCells(), sumOp<label>());

        return std::pow
        (
            activeMeasure/max(nCellsGlobal, scalar(1.0)),
            1.0/max(scalar(mesh.nGeometricD()), scalar(1.0))
        );
    }

    // Volume-weighted L1 norm, sum(|f| V)/sum(V), globally reduced. Weighting
    // by volume is what makes the norms comparable between meshes of different
    // cell size, and reduction is what makes them independent of the
    // decomposition.
    scalar volumeWeightedL1(const volScalarField& fld)
    {
        const scalarField& f = fld.primitiveField();
        const scalarField& V = fld.mesh().V();

        return gSum(mag(f)*V)/gSum(V);
    }

    // Volume-weighted L2 norm, sqrt(sum(f^2 V)/sum(V)), globally reduced. This
    // is the Vm_L2 the convergence tables fit their orders against.
    scalar volumeWeightedL2(const volScalarField& fld)
    {
        const scalarField& f = fld.primitiveField();
        const scalarField& V = fld.mesh().V();

        return std::sqrt(gSum(sqr(f)*V)/gSum(V));
    }

    // Global maximum of |f|.
    scalar linfNorm(const volScalarField& fld)
    {
        const scalarField magFld(mag(fld.primitiveField()));
        return gMax(magFld);
    }

    // Unweighted L2 over cells, sqrt(sum(f^2)), globally reduced. Reported
    // alongside the weighted norms because it is the one that does not hide a
    // localised error behind a small cell volume.
    scalar cellErrorNorm(const volScalarField& fld)
    {
        return std::sqrt(gSum(sqr(fld.primitiveField())));
    }

    FieldErrorSummary computeFieldErrorSummary(const volScalarField& fld)
    {
        FieldErrorSummary s;
        s.L1 = volumeWeightedL1(fld);
        s.L2 = volumeWeightedL2(fld);
        s.Linf = linfNorm(fld);
        s.errorCell = cellErrorNorm(fld);
        return s;
    }

    // Same conversion for a source field, kept separate so the two uses read
    // differently at the call site.
    EigVec sourceToEigVec(const volScalarField& fld)
    {
        return fieldToEigVec(fld);
    }



    // Gauss quadrature over one tetrahedron (a,b,c,d), appended to the running
    // lists. Uses a Duffy-style collapsed mapping of the cube onto the
    // tetrahedron, hence the r^2 s Jacobian in the weights.
    void appendTetQuadraturePointsAndWeights
    (
        const point& a,
        const point& b,
        const point& c,
        const point& d,
        DynamicList<point>& qPoints,
        DynamicList<scalar>& qWeights
    )
    {
        static const scalar xi[4] =
        {
            0.06943184420297371,
            0.33000947820757187,
            0.6699905217924281,
            0.9305681557970262
        };

        static const scalar wi[4] =
        {
            0.1739274225687269*0.5,
            0.3260725774312731*0.5,
            0.3260725774312731*0.5,
            0.1739274225687269*0.5
        };

        const scalar tetVol =
            mag((b - a) & ((c - a) ^ (d - a)))/6.0;

        if (tetVol <= VSMALL)
        {
            return;
        }

        for (label i = 0; i < 4; ++i)
        {
            const scalar r = xi[i];
            for (label j = 0; j < 4; ++j)
            {
                const scalar s = xi[j];
                for (label k = 0; k < 4; ++k)
                {
                    const scalar t = xi[k];

                    const scalar l0 = 1.0 - r;
                    const scalar l1 = r*(1.0 - s);
                    const scalar l2 = r*s*(1.0 - t);
                    const scalar l3 = r*s*t;

                    qPoints.append(l0*a + l1*b + l2*c + l3*d);
                    qWeights.append
                    (
                        tetVol*6.0*wi[i]*wi[j]*wi[k]*sqr(r)*s
                    );
                }
            }
        }
    }

    // Quadrature rule for an arbitrary polyhedral cell: the cell is decomposed
    // into tetrahedra on its face triangulation and the cell centre, and the
    // per-tet rules above are concatenated. Weights are PHYSICAL, i.e. they sum
    // to the cell volume.
    void cellQuadraturePointsAndWeights
    (
        const fvMesh& mesh,
        const label cellI,
        List<point>& qPoints,
        scalarField& qWeights
    )
    {
        DynamicList<point> qPointDyn;
        DynamicList<scalar> qWeightDyn;

        const cell& curCell = mesh.cells()[cellI];
        const faceList& faces = mesh.faces();
        const pointField& pts = mesh.points();
        const vectorField& faceCentres = mesh.faceCentres();
        const point& cellCentre = mesh.C()[cellI];

        forAll(curCell, cellFaceI)
        {
            const label faceI = curCell[cellFaceI];
            const face& f = faces[faceI];

            if (f.size() < 3)
            {
                continue;
            }

            const point& faceCentre = faceCentres[faceI];

            forAll(f, fpI)
            {
                const point& p0 = pts[f[fpI]];
                const point& p1 = pts[f[(fpI + 1) % f.size()]];

                appendTetQuadraturePointsAndWeights
                (
                    cellCentre,
                    faceCentre,
                    p0,
                    p1,
                    qPointDyn,
                    qWeightDyn
                );
            }
        }

        if (qPointDyn.size() == 0)
        {
            qPoints.setSize(1);
            qWeights.setSize(1);
            qPoints[0] = cellCentre;
            qWeights[0] = 1.0;
            return;
        }

        scalar wSum = 0.0;
        forAll(qWeightDyn, qI)
        {
            wSum += qWeightDyn[qI];
        }
        wSum = max(wSum, SMALL);

        qPoints.setSize(qPointDyn.size());
        qWeights.setSize(qWeightDyn.size());

        forAll(qPointDyn, qI)
        {
            qPoints[qI] = qPointDyn[qI];
            qWeights[qI] = qWeightDyn[qI]/wSum;
        }
    }

    // Impose the exact state values on Dirichlet patches at time t. The states
    // are local ODEs with no spatial coupling, so this matters only where the
    // reconstruction stencils reach a boundary.
    void updateStateBoundaryValues
    (
        volScalarField& u1,
        volScalarField& u2,
        volScalarField& u3,
        const scalar t,
        const label dim
    )
    {
        const fvMesh& mesh = u1.mesh();
        const surfaceVectorField& Cf = mesh.Cf();

        forAll(u1.boundaryField(), patchI)
        {
            fvPatchScalarField& u1p = u1.boundaryFieldRef()[patchI];
            fvPatchScalarField& u2p = u2.boundaryFieldRef()[patchI];
            fvPatchScalarField& u3p = u3.boundaryFieldRef()[patchI];

            if (u1p.type() == "empty")
            {
                continue;
            }

            const vectorField& CfPatch = Cf.boundaryField()[patchI];

            if (u1p.type() == fixedValueFvPatchScalarField::typeName)
            {
                forAll(u1p, faceI)
                {
                    u1p[faceI] = exactU1(CfPatch[faceI], t, dim);
                    u2p[faceI] = exactU2(CfPatch[faceI], t, dim);
                    u3p[faceI] = exactU3(CfPatch[faceI], t, dim);
                }
            }
        }
    }

    // Evaluate the exact solution at the cell centres into the *Exact fields,
    // for the error measures and for the written output.
    void fillExactFields
    (
        volScalarField& VmExact,
        volScalarField& u1Exact,
        volScalarField& u2Exact,
        const scalar t,
        const label dim
    )
    {
        const vectorField& C = VmExact.mesh().C();

        forAll(C, cellI)
        {
            VmExact[cellI] = exactVm(C[cellI], t, dim);
            u1Exact[cellI] = exactU1(C[cellI], t, dim);
            u2Exact[cellI] = exactU2(C[cellI], t, dim);
        }

        const fvMesh& mesh = VmExact.mesh();
        const surfaceVectorField& Cf = mesh.Cf();

        forAll(VmExact.boundaryField(), patchI)
        {
            fvPatchScalarField& Vp = VmExact.boundaryFieldRef()[patchI];
            fvPatchScalarField& u1p = u1Exact.boundaryFieldRef()[patchI];
            fvPatchScalarField& u2p = u2Exact.boundaryFieldRef()[patchI];

            if (Vp.type() == "empty")
            {
                continue;
            }

            const vectorField& CfPatch = Cf.boundaryField()[patchI];

            if
            (
                Vp.type() == fixedValueFvPatchScalarField::typeName
             || Vp.type() == "fixedVoltage"
            )
            {
                forAll(Vp, faceI)
                {
                    Vp[faceI] = exactVm(CfPatch[faceI], t, dim);
                }
            }

            if (u1p.type() == fixedValueFvPatchScalarField::typeName)
            {
                forAll(u1p, faceI)
                {
                    u1p[faceI] = exactU1(CfPatch[faceI], t, dim);
                    u2p[faceI] = exactU2(CfPatch[faceI], t, dim);
                }
            }
        }

        VmExact.correctBoundaryConditions();
        u1Exact.correctBoundaryConditions();
        u2Exact.correctBoundaryConditions();
    }

    // Ionic current evaluated once per cell at the cell centre: the low-order
    // path, and the reference the quadrature path is compared against.
    void computeCellCentredIion
    (
        const volScalarField& Vm,
        const volScalarField& u1,
        const volScalarField& u2,
        const volScalarField& u3,
        const scalar beta,
        const scalar chiVal,
        const scalar CmVal,
        volScalarField& Iion,
        const bool useOpenMP,
        const label openMPThreshold
    )
    {
        const label nCells = Vm.internalField().size();

        #pragma omp parallel for schedule(static) if(useOpenMP && nCells >= openMPThreshold)
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            Iion[cellI] =
                ionicCurrentPDE
                (
                    Vm[cellI],
                    u1[cellI],
                    u2[cellI],
                    u3[cellI],
                    beta,
                    chiVal,
                    CmVal
                );
        }

        Iion.correctBoundaryConditions();
    }

    // Vm source evaluated once per cell at the cell centre.
    void computeCellCentredVmSource
    (
        const volScalarField& Vm,
        const volScalarField& u1,
        const volScalarField& u2,
        const volScalarField& u3,
        const scalar beta,
        const scalar chiVal,
        const scalar CmVal,
        volScalarField& sourceVm,
        const bool useOpenMP,
        const label openMPThreshold
    )
    {
        const label nCells = Vm.internalField().size();

        #pragma omp parallel for schedule(static) if(useOpenMP && nCells >= openMPThreshold)
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            sourceVm[cellI] =
                vmSourcePDE
                (
                    Vm[cellI],
                    u1[cellI],
                    u2[cellI],
                    u3[cellI],
                    beta,
                    chiVal,
                    CmVal
                );
        }

        sourceVm.correctBoundaryConditions();
    }

    // The three cell-centred states of the manufactured model at the ionic
    // quadrature points, through the shared reconstruction core
    // (monodomainReconstruction.H) - fused, threaded, never limited.
    void reconstructStatesAtIionIntegrationPoints
    (
        const volScalarField& u1,
        const volScalarField& u2,
        const volScalarField& u3,
        const highOrderInterp& LREInterp_states,
        const highOrderInterp& LREInterp_Iion,
        scalarField& u1IntegrationPoints,
        scalarField& u2IntegrationPoints,
        scalarField& u3IntegrationPoints,
        const Switch useOpenMP,
        const label openMPThreshold,
        const Switch profileTimings,
        FineGrainedTimings& fineTimings
    )
    {
        const fvMesh& mesh = u1.mesh();
        const CompactListList<point>& qp = LREInterp_Iion.cellQuadPoints();

        // The three states in one walk of the stencil, as this solver's own
        // loop did before the shared core existed.
        const scalarField* fieldsArr[3] =
            {&u1.primitiveField(), &u2.primitiveField(), &u3.primitiveField()};
        scalarField* outputs[3] =
            {&u1IntegrationPoints, &u2IntegrationPoints, &u3IntegrationPoints};

        auto write = [&outputs](const label ip, const label k, const scalar v)
        {
            (*outputs[k])[ip] = v;
        };

        reconstructCellFieldsAtIntegrationPoints
        (
            mesh,
            UList<const scalarField*>(fieldsArr, 3),
            LREInterp_states,
            qp,
            write,
            nullptr,
            false,
            useOpenMP,
            openMPThreshold,
            profileTimings ? &fineTimings.stateReconDerivWallTime : nullptr,
            profileTimings ? &fineTimings.stateReconExpandWallTime : nullptr
        );

        if (profileTimings)
        {
            ++fineTimings.stateReconCalls;
        }
    }

    // Collapse a quantity defined at the quadrature points back to a cell
    // average, weighting by the physical quadrature weights. The inverse
    // direction of the reconstruction above.
    void averageIntegrationPointFieldToCells
    (
        const scalarField& integrationPointValues,
        const highOrderInterp& LREInterp_Iion,
        volScalarField& field
    )
    {
        const fvMesh& mesh = field.mesh();
        const CompactListList<scalar>& cellIionQuadW =
            LREInterp_Iion.cellQuadWeightPhysical();

        scalarField& cellValues = field.primitiveFieldRef();
        label integrationPointI = 0;

        forAll(mesh.cells(), cellI)
        {
            scalar valueBar = 0.0;
            scalar wSum = 0.0;

            forAll(cellIionQuadW[cellI], qI)
            {
                const scalar w = cellIionQuadW[cellI][qI];
                valueBar += w*integrationPointValues[integrationPointI];
                wSum += w;
                ++integrationPointI;
            }

            cellValues[cellI] = valueBar/max(wSum, SMALL);
        }

        field.correctBoundaryConditions();
    }

    // Seed the per-quadrature-point states from the exact solution at t = 0.
    // Necessary because those states are independent unknowns in the
    // gaussPointODE mode - they are not derived from the cell values.
    void initialiseIionIntegrationPointStates
    (
        const highOrderInterp& LREInterp_Iion,
        const scalar t,
        const label dim,
        scalarField& u1IntegrationPoints,
        scalarField& u2IntegrationPoints,
        scalarField& u3IntegrationPoints
    )
    {
        const CompactListList<point>& cellIionQuadP =
            LREInterp_Iion.cellQuadPoints();

        label integrationPointI = 0;
        forAll(cellIionQuadP, cellI)
        {
            forAll(cellIionQuadP[cellI], qI)
            {
                const point& p = cellIionQuadP[cellI][qI];
                u1IntegrationPoints[integrationPointI] = exactU1(p, t, dim);
                u2IntegrationPoints[integrationPointI] = exactU2(p, t, dim);
                u3IntegrationPoints[integrationPointI] = exactU3(p, t, dim);
                ++integrationPointI;
            }
        }
    }

    // Cell-averaged ionic current obtained by evaluating Iion at every
    // quadrature point and integrating, rather than evaluating it once at the
    // cell centre.
    void computeIionFromIntegrationPoints
    (
        const scalarField& VmIntegrationPoints,
        const scalarField& u1IntegrationPoints,
        const scalarField& u2IntegrationPoints,
        const scalarField& u3IntegrationPoints,
        const scalar beta,
        const scalar chiVal,
        const scalar CmVal,
        scalarField& IionIntegrationPoints,
        const bool useOpenMP,
        const label openMPThreshold
    )
    {
        const label nIntegrationPoints = IionIntegrationPoints.size();

        #pragma omp parallel for schedule(static) if(useOpenMP && nIntegrationPoints >= openMPThreshold)
        for (label integrationPointI = 0; integrationPointI < nIntegrationPoints; ++integrationPointI)
        {
            IionIntegrationPoints[integrationPointI] =
                ionicCurrentPDE
                (
                    VmIntegrationPoints[integrationPointI],
                    u1IntegrationPoints[integrationPointI],
                    u2IntegrationPoints[integrationPointI],
                    u3IntegrationPoints[integrationPointI],
                    beta,
                    chiVal,
                    CmVal
                );
        }
    }

    // Same for the Vm source term.
    void computeVmSourceFromIntegrationPoints
    (
        const scalarField& VmIntegrationPoints,
        const scalarField& u1IntegrationPoints,
        const scalarField& u2IntegrationPoints,
        const scalarField& u3IntegrationPoints,
        const scalar beta,
        const scalar chiVal,
        const scalar CmVal,
        scalarField& sourceIntegrationPoints,
        const bool useOpenMP,
        const label openMPThreshold
    )
    {
        const label nIntegrationPoints = sourceIntegrationPoints.size();

        #pragma omp parallel for schedule(static) if(useOpenMP && nIntegrationPoints >= openMPThreshold)
        for (label integrationPointI = 0; integrationPointI < nIntegrationPoints; ++integrationPointI)
        {
            sourceIntegrationPoints[integrationPointI] =
                vmSourcePDE
                (
                    VmIntegrationPoints[integrationPointI],
                    u1IntegrationPoints[integrationPointI],
                    u2IntegrationPoints[integrationPointI],
                    u3IntegrationPoints[integrationPointI],
                    beta,
                    chiVal,
                    CmVal
                );
        }
    }

    void reactionRatesLinearVm
    (
        const VmDriverNodes& vmNodes,
        const scalar tau,
        const scalar dt,
        const scalar u1,
        const scalar u2,
        const scalar u3,
        scalar& du1dt,
        scalar& du2dt,
        scalar& du3dt
    )
    {
        const scalar alpha =
            dt > SMALL
          ? min(max(tau/dt, scalar(0.0)), scalar(1.0))
          : scalar(1.0);

        const scalar VmTau = interpolateVmDriver(vmNodes, alpha);

        reactionRates(VmTau, u1, u2, u3, du1dt, du2dt, du3dt);
    }

    // One classical fourth-order Runge-Kutta step of the state ODEs at one
    // point. Fixed step: used when the adaptive integrator is not requested.
    void rk4StateStep
    (
        const VmDriverNodes& vmNodes,
        const scalar tau,
        const scalar h,
        const scalar dt,
        const scalar u1,
        const scalar u2,
        const scalar u3,
        scalar& u1New,
        scalar& u2New,
        scalar& u3New
    )
    {
        scalar k11 = 0.0, k12 = 0.0, k13 = 0.0;
        scalar k21 = 0.0, k22 = 0.0, k23 = 0.0;
        scalar k31 = 0.0, k32 = 0.0, k33 = 0.0;
        scalar k41 = 0.0, k42 = 0.0, k43 = 0.0;

        reactionRatesLinearVm
        (
            vmNodes, tau, dt,
            u1, u2, u3,
            k11, k12, k13
        );
        reactionRatesLinearVm
        (
            vmNodes, tau + 0.5*h, dt,
            u1 + 0.5*h*k11,
            u2 + 0.5*h*k12,
            u3 + 0.5*h*k13,
            k21, k22, k23
        );
        reactionRatesLinearVm
        (
            vmNodes, tau + 0.5*h, dt,
            u1 + 0.5*h*k21,
            u2 + 0.5*h*k22,
            u3 + 0.5*h*k23,
            k31, k32, k33
        );
        reactionRatesLinearVm
        (
            vmNodes, tau + h, dt,
            u1 + h*k31,
            u2 + h*k32,
            u3 + h*k33,
            k41, k42, k43
        );

        u1New = u1 + (h/6.0)*(k11 + 2.0*k21 + 2.0*k31 + k41);
        u2New = u2 + (h/6.0)*(k12 + 2.0*k22 + 2.0*k32 + k42);
        u3New = u3 + (h/6.0)*(k13 + 2.0*k23 + 2.0*k33 + k43);
    }

    // One Runge-Kutta-Fehlberg 4(5) step with an embedded error estimate, from
    // which the next step size is chosen. The state ODEs are stiff at the
    // upstroke, so the adaptive step is what keeps the local error below the
    // tolerance without paying for it everywhere.
    void rkf45StateStep
    (
        const VmDriverNodes& vmNodes,
        const scalar tau,
        const scalar h,
        const scalar dt,
        const scalar u1,
        const scalar u2,
        const scalar u3,
        scalar& u1Fifth,
        scalar& u2Fifth,
        scalar& u3Fifth,
        const scalar stateODEAbsTol,
        const scalar stateODERelTol,
        scalar& err
    )
    {
        scalar k11 = 0.0, k12 = 0.0, k13 = 0.0;
        scalar k21 = 0.0, k22 = 0.0, k23 = 0.0;
        scalar k31 = 0.0, k32 = 0.0, k33 = 0.0;
        scalar k41 = 0.0, k42 = 0.0, k43 = 0.0;
        scalar k51 = 0.0, k52 = 0.0, k53 = 0.0;
        scalar k61 = 0.0, k62 = 0.0, k63 = 0.0;

        reactionRatesLinearVm
        (
            vmNodes, tau, dt,
            u1, u2, u3,
            k11, k12, k13
        );
        reactionRatesLinearVm
        (
            vmNodes, tau + h/5.0, dt,
            u1 + h*(1.0/5.0)*k11,
            u2 + h*(1.0/5.0)*k12,
            u3 + h*(1.0/5.0)*k13,
            k21, k22, k23
        );
        reactionRatesLinearVm
        (
            vmNodes, tau + 3.0*h/10.0, dt,
            u1 + h*((3.0/40.0)*k11 + (9.0/40.0)*k21),
            u2 + h*((3.0/40.0)*k12 + (9.0/40.0)*k22),
            u3 + h*((3.0/40.0)*k13 + (9.0/40.0)*k23),
            k31, k32, k33
        );
        reactionRatesLinearVm
        (
            vmNodes, tau + 3.0*h/5.0, dt,
            u1 + h*((3.0/10.0)*k11 - (9.0/10.0)*k21 + (6.0/5.0)*k31),
            u2 + h*((3.0/10.0)*k12 - (9.0/10.0)*k22 + (6.0/5.0)*k32),
            u3 + h*((3.0/10.0)*k13 - (9.0/10.0)*k23 + (6.0/5.0)*k33),
            k41, k42, k43
        );
        reactionRatesLinearVm
        (
            vmNodes, tau + h, dt,
            u1 + h*((-11.0/54.0)*k11 + (5.0/2.0)*k21 - (70.0/27.0)*k31 + (35.0/27.0)*k41),
            u2 + h*((-11.0/54.0)*k12 + (5.0/2.0)*k22 - (70.0/27.0)*k32 + (35.0/27.0)*k42),
            u3 + h*((-11.0/54.0)*k13 + (5.0/2.0)*k23 - (70.0/27.0)*k33 + (35.0/27.0)*k43),
            k51, k52, k53
        );
        reactionRatesLinearVm
        (
            vmNodes, tau + 7.0*h/8.0, dt,
            u1 + h*((1631.0/55296.0)*k11 + (175.0/512.0)*k21 + (575.0/13824.0)*k31 + (44275.0/110592.0)*k41 + (253.0/4096.0)*k51),
            u2 + h*((1631.0/55296.0)*k12 + (175.0/512.0)*k22 + (575.0/13824.0)*k32 + (44275.0/110592.0)*k42 + (253.0/4096.0)*k52),
            u3 + h*((1631.0/55296.0)*k13 + (175.0/512.0)*k23 + (575.0/13824.0)*k33 + (44275.0/110592.0)*k43 + (253.0/4096.0)*k53),
            k61, k62, k63
        );

        u1Fifth = u1 + h*((37.0/378.0)*k11 + (250.0/621.0)*k31 + (125.0/594.0)*k41 + (512.0/1771.0)*k61);
        u2Fifth = u2 + h*((37.0/378.0)*k12 + (250.0/621.0)*k32 + (125.0/594.0)*k42 + (512.0/1771.0)*k62);
        u3Fifth = u3 + h*((37.0/378.0)*k13 + (250.0/621.0)*k33 + (125.0/594.0)*k43 + (512.0/1771.0)*k63);

        const scalar u1Fourth = u1 + h*((2825.0/27648.0)*k11 + (18575.0/48384.0)*k31 + (13525.0/55296.0)*k41 + (277.0/14336.0)*k51 + 0.25*k61);
        const scalar u2Fourth = u2 + h*((2825.0/27648.0)*k12 + (18575.0/48384.0)*k32 + (13525.0/55296.0)*k42 + (277.0/14336.0)*k52 + 0.25*k62);
        const scalar u3Fourth = u3 + h*((2825.0/27648.0)*k13 + (18575.0/48384.0)*k33 + (13525.0/55296.0)*k43 + (277.0/14336.0)*k53 + 0.25*k63);

        err = max
        (
            mag(u1Fifth - u1Fourth)
           /(stateODEAbsTol + stateODERelTol*max(mag(u1Fifth), mag(u1Fourth))),
            max
            (
                mag(u2Fifth - u2Fourth)
               /(stateODEAbsTol + stateODERelTol*max(mag(u2Fifth), mag(u2Fourth))),
                mag(u3Fifth - u3Fourth)
               /(stateODEAbsTol + stateODERelTol*max(mag(u3Fifth), mag(u3Fourth)))
            )
        );
    }

    // Integrate the state ODEs at one point across a full PDE time step,
    // sub-stepping adaptively until the interval is covered. Vm is interpolated
    // linearly across the step (reactionRatesLinearVm), so the sub-steps see a
    // consistent voltage.
    void advanceStateODE
    (
        const VmDriverNodes& vmNodes,
        const scalar u1Old,
        const scalar u2Old,
        const scalar u3Old,
        const scalar dt,
        const word& stateODESolver,
        const scalar stateODEInitialStep,
        const scalar stateODEAbsTol,
        const scalar stateODERelTol,
        const label stateODEMaxSteps,
        scalar& u1New,
        scalar& u2New,
        scalar& u3New
    )
    {
        if (dt <= SMALL)
        {
            u1New = u1Old;
            u2New = u2Old;
            u3New = u3Old;
            return;
        }

        if (stateODESolver == "Euler")
        {
            scalar du1 = 0.0, du2 = 0.0, du3 = 0.0;
            reactionRatesLinearVm
            (
                vmNodes, 0.0, dt,
                u1Old, u2Old, u3Old,
                du1, du2, du3
            );

            u1New = u1Old + dt*du1;
            u2New = u2Old + dt*du2;
            u3New = u3Old + dt*du3;
            return;
        }

        if (stateODESolver == "RK4")
        {
            rk4StateStep
            (
                vmNodes, 0.0, dt, dt,
                u1Old, u2Old, u3Old,
                u1New, u2New, u3New
            );
            return;
        }

        // Validated at start-up in createFields.H, so this is a guard, not a
        // reachable branch.
        if (stateODESolver != "RKF45")
        {
            FatalErrorInFunction
                << "Unknown state ODE solver: " << stateODESolver << nl
                << "Valid options are RKF45, RK4 and Euler"
                << exit(FatalError);
        }

        scalar tau = 0.0;
        scalar h =
            stateODEInitialStep > SMALL
          ? min(dt, stateODEInitialStep)
          : dt;
        scalar u1Cur = u1Old;
        scalar u2Cur = u2Old;
        scalar u3Cur = u3Old;
        const scalar absTol = max(stateODEAbsTol, scalar(SMALL));
        const scalar relTol = max(stateODERelTol, scalar(SMALL));
        const scalar hMin = max(dt*1.0e-12, SMALL);

        for (label subStep = 0; subStep < max(stateODEMaxSteps, label(1)); ++subStep)
        {
            if (tau >= dt - SMALL)
            {
                break;
            }

            h = min(h, dt - tau);

            scalar u1Trial = u1Cur;
            scalar u2Trial = u2Cur;
            scalar u3Trial = u3Cur;
            scalar err = GREAT;

            rkf45StateStep
            (
                vmNodes, tau, h, dt,
                u1Cur, u2Cur, u3Cur,
                u1Trial, u2Trial, u3Trial,
                absTol,
                relTol,
                err
            );

            if (err <= 1.0 || h <= hMin)
            {
                u1Cur = u1Trial;
                u2Cur = u2Trial;
                u3Cur = u3Trial;
                tau += h;
            }

            const scalar factor = min
            (
                scalar(4.0),
                max(scalar(0.1), scalar(0.84)*std::pow(1.0/(err + SMALL), 0.25))
            );
            h = max(hMin, min(dt - tau, h*factor));
        }

        if (tau < dt - 10.0*SMALL)
        {
            WarningInFunction
                << "State RKF45 reached maxSteps before completing dt. "
                << "tau = " << tau << ", dt = " << dt << nl;
        }

        u1New = u1Cur;
        u2New = u2Cur;
        u3New = u3Cur;
    }

    // Advance the states at every ionic quadrature point over one PDE step.
    // Each point is an independent ODE, hence the optional OpenMP loop.
    void updateIntegrationPointStatesODE
    (
        const VmDriverStages& ipDriver,
        const scalarField& VmNewIntegrationPoints,
        const scalarField& u1OldIntegrationPoints,
        const scalarField& u2OldIntegrationPoints,
        const scalarField& u3OldIntegrationPoints,
        scalarField& u1NewIntegrationPoints,
        scalarField& u2NewIntegrationPoints,
        scalarField& u3NewIntegrationPoints,
        const scalar dt,
        const word& stateODESolver,
        const scalar stateODEInitialStep,
        const scalar stateODEAbsTol,
        const scalar stateODERelTol,
        const label stateODEMaxSteps,
        const bool useOpenMP,
        const label openMPThreshold
    )
    {
        const label nIntegrationPoints = u1NewIntegrationPoints.size();

        #pragma omp parallel for schedule(static) if(useOpenMP && nIntegrationPoints >= openMPThreshold)
        for (label integrationPointI = 0; integrationPointI < nIntegrationPoints; ++integrationPointI)
        {
            // The driver's nodes at this point: the two-node linear blend by
            // default, the history polynomial under stateODEGaussPointHistory.
            // A node without a field is the candidate, VmNewIntegrationPoints.
            VmDriverNodes vmNodes;
            vmNodes.nNodes = ipDriver.nNodes;
            for (label j = 0; j < ipDriver.nNodes; ++j)
            {
                vmNodes.s[j] = ipDriver.s[j];
                vmNodes.v[j] =
                    ipDriver.field[j]
                  ? (*ipDriver.field[j])[integrationPointI]
                  : VmNewIntegrationPoints[integrationPointI];
            }

            advanceStateODE
            (
                vmNodes,
                u1OldIntegrationPoints[integrationPointI],
                u2OldIntegrationPoints[integrationPointI],
                u3OldIntegrationPoints[integrationPointI],
                dt,
                stateODESolver,
                stateODEInitialStep,
                stateODEAbsTol,
                stateODERelTol,
                stateODEMaxSteps,
                u1NewIntegrationPoints[integrationPointI],
                u2NewIntegrationPoints[integrationPointI],
                u3NewIntegrationPoints[integrationPointI]
            );
        }
    }

    void updateStateFieldsODE
    (
        const VmDriverStages& vmDriver,
        const volScalarField& u1Old,
        const volScalarField& u2Old,
        const volScalarField& u3Old,
        volScalarField& u1New,
        volScalarField& u2New,
        volScalarField& u3New,
        const scalar dt,
        const word& stateODESolver,
        const scalar stateODEInitialStep,
        const scalar stateODEAbsTol,
        const scalar stateODERelTol,
        const label stateODEMaxSteps,
        const bool useOpenMP,
        const label openMPThreshold
    )
    {
        scalarField& u1NI = u1New.primitiveFieldRef();
        scalarField& u2NI = u2New.primitiveFieldRef();
        scalarField& u3NI = u3New.primitiveFieldRef();

        const scalarField& u1O = u1Old.primitiveField();
        const scalarField& u2O = u2Old.primitiveField();
        const scalarField& u3O = u3Old.primitiveField();

        const label nCells = u1NI.size();
        const label nNodes = vmDriver.nNodes;

        for (label j = 0; j < nNodes; ++j)
        {
            if (!vmDriver.field[j])
            {
                FatalErrorInFunction
                    << "Vm driver node " << j << " of " << nNodes
                    << " was not set. A missing node is a silently wrong ODE"
                    << " driver, not a crash, so it is checked here."
                    << exit(FatalError);
            }
        }

        #pragma omp parallel for schedule(static) if(useOpenMP && nCells >= openMPThreshold)
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            VmDriverNodes vmNodes;
            vmNodes.nNodes = nNodes;

            for (label j = 0; j < nNodes; ++j)
            {
                vmNodes.s[j] = vmDriver.s[j];
                vmNodes.v[j] = (*vmDriver.field[j])[cellI];
            }

            advanceStateODE
            (
                vmNodes,
                u1O[cellI],
                u2O[cellI],
                u3O[cellI],
                dt,
                stateODESolver,
                stateODEInitialStep,
                stateODEAbsTol,
                stateODERelTol,
                stateODEMaxSteps,
                u1NI[cellI],
                u2NI[cellI],
                u3NI[cellI]
            );
        }
    }

    ReconstructionErrorSummary computeGaussReconstructedError
    (
        const volScalarField& fld,
        const scalar t,
        const ManufacturedFieldID fieldID,
        const label dim,
        const highOrderInterp& LREInterp
    )
    {
        const fvMesh& mesh = fld.mesh();
        const bool useTaylor = mesh.nGeometricD() > 1;
        const bool twoD = mesh.nGeometricD() == 2;

        tmp<volVectorField> tGrad;
        const volVectorField* gradPtr = nullptr;
        if (useTaylor)
        {
            tGrad = LREInterp.grad(fld);
            gradPtr = &tGrad();
        }

        tmp<volSymmTensorField> tHess;
        const volSymmTensorField* hessPtr = nullptr;
        if (useTaylor && LREInterp.order() >= 2)
        {
            tHess = LREInterp.hessian(fld);
            hessPtr = &tHess();
        }

        autoPtr<List<highOrderInterp::symmTensor3Order>> thirdPtr;
        const List<highOrderInterp::symmTensor3Order>* thirdList = nullptr;
        if (useTaylor && LREInterp.order() >= 3)
        {
            thirdPtr = LREInterp.thirdDeriv(fld);
            thirdList = &thirdPtr();
        }

        const vectorField& C = mesh.C();
        const scalarField& V = mesh.V();

        scalar e1 = 0.0;
        scalar e2 = 0.0;
        scalar eInf = 0.0;
        scalar n2 = 0.0;
        scalar volTot = 0.0;

        List<point> qPoints;
        scalarField qWeights;

        forAll(fld.internalField(), cellI)
        {
            cellQuadraturePointsAndWeights(mesh, cellI, qPoints, qWeights);

            const scalar cellV = V[cellI];

            forAll(qPoints, qpI)
            {
                const point& gp = qPoints[qpI];
                const scalar w = qWeights[qpI];
                const scalar meas = cellV*w;
                const vector d = gp - C[cellI];

                scalar num = fld[cellI];

                if (gradPtr)
                {
                    const symmTensor* H =
                        hessPtr ? &((*hessPtr)[cellI]) : nullptr;
                    const highOrderInterp::symmTensor3Order* T3 =
                        thirdList ? &((*thirdList)[cellI]) : nullptr;

                    num =
                        reconstructFromTaylor
                        (
                            fld[cellI],
                            (*gradPtr)[cellI],
                            H,
                            T3,
                            d,
                            twoD
                        );
                }

                const scalar ex = exactFieldValue(fieldID, gp, t, dim);
                const scalar err = num - ex;

                e1 += meas*mag(err);
                e2 += meas*sqr(err);
                eInf = max(eInf, mag(err));
                n2 += meas*sqr(ex);
                volTot += meas;
            }
        }

        e1 = returnReduce(e1, sumOp<scalar>());
        e2 = returnReduce(e2, sumOp<scalar>());
        eInf = returnReduce(eInf, maxOp<scalar>());
        n2 = returnReduce(n2, sumOp<scalar>());
        volTot = returnReduce(volTot, sumOp<scalar>());

        ReconstructionErrorSummary s;
        s.L1 = e1/max(volTot, SMALL);
        s.L2 = std::sqrt(e2/max(volTot, SMALL));
        s.Linf = eInf;
        s.normL2 = std::sqrt(n2/max(volTot, SMALL));
        s.relL2 = 100.0*s.L2/max(s.normL2, SMALL);

        return s;
    }

    ReconstructionErrorSummary cellAsReconstructionError
    (
        const FieldErrorSummary& cellErr,
        const volScalarField& exact
    )
    {
        ReconstructionErrorSummary s;
        s.L1 = cellErr.L1;
        s.L2 = cellErr.L2;
        s.Linf = cellErr.Linf;
        s.normL2 = volumeWeightedL2(exact);
        s.relL2 = 100.0*s.L2/max(s.normL2, SMALL);
        return s;
    }

    ReconstructionErrorSummary computeIntegrationPointStateError
    (
        const scalarField& stateIntegrationPoints,
        const scalar t,
        const ManufacturedFieldID fieldID,
        const label dim,
        const highOrderInterp& LREInterp_Iion,
        const fvMesh& mesh
    )
    {
        const CompactListList<point>& cellIionQuadP =
            LREInterp_Iion.cellQuadPoints();
        const CompactListList<scalar>& cellIionQuadW =
            LREInterp_Iion.cellQuadWeightPhysical();
        const scalarField& V = mesh.V();

        scalar e1 = 0.0;
        scalar e2 = 0.0;
        scalar eInf = 0.0;
        scalar n2 = 0.0;
        scalar volTot = 0.0;

        label integrationPointI = 0;
        forAll(mesh.cells(), cellI)
        {
            const scalar cellV = V[cellI];
            scalar wSum = 0.0;

            forAll(cellIionQuadW[cellI], qI)
            {
                wSum += cellIionQuadW[cellI][qI];
            }

            forAll(cellIionQuadP[cellI], qI)
            {
                const point& gp = cellIionQuadP[cellI][qI];
                const scalar meas =
                    cellV*cellIionQuadW[cellI][qI]/max(wSum, SMALL);
                const scalar ex = exactFieldValue(fieldID, gp, t, dim);
                const scalar err = stateIntegrationPoints[integrationPointI] - ex;

                e1 += meas*mag(err);
                e2 += meas*sqr(err);
                eInf = max(eInf, mag(err));
                n2 += meas*sqr(ex);
                volTot += meas;
                ++integrationPointI;
            }
        }

        e1 = returnReduce(e1, sumOp<scalar>());
        e2 = returnReduce(e2, sumOp<scalar>());
        eInf = returnReduce(eInf, maxOp<scalar>());
        n2 = returnReduce(n2, sumOp<scalar>());
        volTot = returnReduce(volTot, sumOp<scalar>());

        ReconstructionErrorSummary s;
        s.L1 = e1/max(volTot, SMALL);
        s.L2 = std::sqrt(e2/max(volTot, SMALL));
        s.Linf = eInf;
        s.normL2 = std::sqrt(n2/max(volTot, SMALL));
        s.relL2 = 100.0*s.L2/max(s.normL2, SMALL);

        return s;
    }

    // Write the run summary consumed by run_convergence.py: the error norms of
    // every field at the cell centres and at the quadrature points, the fitted
    // quantities, and the timing/memory figures. This file is the interface
    // between the solver and the convergence tables, so its layout is part of
    // the driver's contract.
    void writeSummary
    (
        const Time& runTime,
        const label nSteps,
        const scalar dt,
        const FieldErrorSummary& VmCell,
        const FieldErrorSummary& u1Cell,
        const FieldErrorSummary& u2Cell,
        const ReconstructionErrorSummary& VmHO,
        const ReconstructionErrorSummary& u1HO,
        const ReconstructionErrorSummary& u2HO,
        const volScalarField& rhsVm,
        const volScalarField& Iion,
        const fvMesh& mesh,
        const word& implicitScheme,
        const word& bdfStartup,
        const word& stateODEDriver,
        const word& manufacturedSolution,
        const word& massMatrixType,
        const bool useHighOrderVm,
        const bool useHighOrderIion,
        const scalar stabilisationAlpha,
        const word& memoryOptimization,
        const bool memoryOptimizationEffective,
        const label memoryOptimizationCellThreshold,
        const Switch memoryOptimizationAdaptiveTripletReserve,
        const label memoryOptimizationFluxTripletsPerFace,
        const label memoryOptimizationStabilisedTripletsPerFace,
        const Switch memoryOptimizationTrimHeap,
        const Switch memoryOptimizationCompactMassAssembly,
        const bool compactMassAssembly,
        const label stiffnessTripletsPerFaceReserve,
        const label allocatedIionIntegrationPoints,
        const label lreN,
        const label lreNn,
        const scalar lreK,
        const label lreMaxStencilSize,
        const label lreIionN,
        const label lreIionNn,
        const scalar lreIionK,
        const label lreIionMaxStencilSize,
        const word& linearSolverBackend,
        const word& implicitLinearSolver,
        const word& petscLinearKspType,
        const word& petscLinearPcType,
        const scalar implicitTolerance,
        const label implicitMaxIterations,
        const label maxNonlinearIterations,
        const label minNonlinearIterations,
        const scalar nonlinearTolerance,
        const scalar nonlinearVmTolerance,
        const scalar nonlinearStatesTolerance,
        const Switch nonlinearRequireStatesConvergence,
        const Switch nonlinearAcceptUnconverged,
        const scalar nonlinearIionTolerance,
        const scalar nonlinearRelaxation,
        const label jfnkMaxKrylovIterations,
        const label jfnkMaxRestarts,
        const scalar jfnkLinearTolerance,
        const scalar jfnkEpsilon,
        const word& jfnkLinearSolverBackend,
        const word& jfnkPetscKspType,
        const word& jfnkPetscPcType,
        const label jfnkInitGuessOrder,
        const scalar jfnkInitGuessVmMin,
        const scalar jfnkInitGuessVmMax,
        const Switch jfnkClampODEInput,
        const Switch jfnkLineSearch,
        const label jfnkLineSearchMaxIter,
        const scalar jfnkLineSearchAlphaMin,
        const scalar jfnkArmijoC,
        const word& jfnkPreconditioner,
        const scalar jfnkPreconditionerDropTolerance,
        const label jfnkPreconditionerFillFactor,
        const label jfnkPreconditionerUpdateFrequency,
        const scalar diagonalIionEpsilon,
        const word& stateODESolver,
        const scalar stateODEInitialStep,
        const scalar stateODEAbsTol,
        const scalar stateODERelTol,
        const label stateODEMaxSteps,
        const Switch stateODEUseOpenMP,
        const label stateODENumThreads,
        const std::vector<NonlinearConvergenceRecord>& nonlinearHistory,
        const scalar peakMemoryKB,
        const scalar setupWallTime,
        const scalar timeLoopWallTime,
        const scalar postProcessWallTime,
        const scalar totalWallTime,
        const FineGrainedTimings& fineTimings,
        const word& nonlinearMethod
    )
    {
        const label N = estimatedN(mesh);
        const scalar dx = characteristicDx(mesh);
        const word dimName = name(mesh.nGeometricD()) + "D";
        const bool isPicard = nonlinearMethod == "Picard";
        const bool isJFNK = nonlinearMethod == "JFNK";
        const bool isDiagonalIion = nonlinearMethod == "diagonalIion";

        const NonlinearConvergenceRecord* finalNonlinear =
            nonlinearHistory.empty() ? nullptr : &nonlinearHistory.back();

        label maxNonlinearIterationsUsed = 0;
        label nonConvergedNonlinearSteps = 0;
        label rolledBackNonlinearSteps = 0;
        scalar maxCoupledResidual = 0.0;
        scalar maxVmResidual = 0.0;
        scalar maxU1Residual = 0.0;
        scalar maxU2Residual = 0.0;
        scalar maxU3Residual = 0.0;
        scalar maxStateResidual = 0.0;
        scalar maxIionResidual = 0.0;

        for (const NonlinearConvergenceRecord& rec : nonlinearHistory)
        {
            maxNonlinearIterationsUsed =
                max(maxNonlinearIterationsUsed, rec.iterations);
            if (!rec.converged)
            {
                ++nonConvergedNonlinearSteps;
            }
            if (rec.rolledBack)
            {
                ++rolledBackNonlinearSteps;
            }
            maxCoupledResidual = max(maxCoupledResidual, rec.coupledResidual);
            maxVmResidual = max(maxVmResidual, rec.VmResidual);
            maxU1Residual = max(maxU1Residual, rec.u1Residual);
            maxU2Residual = max(maxU2Residual, rec.u2Residual);
            maxU3Residual = max(maxU3Residual, rec.u3Residual);
            maxStateResidual = max(maxStateResidual, rec.maxStateResidual);
            maxIionResidual = max(maxIionResidual, rec.IionResidual);
        }

        // The canonical summary goes to the global
        // case directory, written by the master; the other ranks write a
        // throwaway copy into their own processorN directory.
        //
        // Note the shape of this: every rank must still construct an OFstream.
        // OpenFOAM's file handler can be collective, so a naive
        // "if (!Pstream::master()) return;" here deadlocks the run - the master
        // enters the file operation and the others never do. That cost a hang
        // that looked like a solver problem: the time loop completed, every
        // nonlinear step converged, and then the run simply stopped.
        //
        // runTime.path() points at processorN/ on a decomposed run, which is
        // why the master needs rootPath()/globalCaseName() instead: otherwise
        // the post-processing driver finds no file in the case directory.
        // Caras del borde FISICO. Se calcula aca porque el reduce es colectivo
        // y, como explica la nota de arriba, todos los ranks entran a esta
        // funcion. Los parches empty (zmin/zmax en 2D) aportan 0 solos, porque
        // emptyFvPatch::size() devuelve 0; !coupled() descarta los processor,
        // que son interiores a la malla global. Usar nFaces()-nInternalFaces()
        // contaria ambos y daria un numero distinto en 2D y en cada particion.
        label nBoundaryFacesGlobal = 0;
        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& bPatch = mesh.boundary()[patchI];
            if (!bPatch.coupled())
            {
                nBoundaryFacesGlobal += bPatch.size();
            }
        }
        reduce(nBoundaryFacesGlobal, sumOp<label>());

        const fileName summaryName
        (
            dimName + "_" + name(N) + "_cells_transient.dat"
        );

        const fileName outFile =
            Pstream::master()
          ? runTime.rootPath()/runTime.globalCaseName()/summaryName
          : runTime.path()/summaryName;

        OFstream os(outFile);

        os  << "Manufactured-solution error summary (cell-centred):" << nl
            << "Field     L1-error       L2-error       Linf-error" << nl
            << "Vm      " << VmCell.L1 << "   " << VmCell.L2 << "   "
            << VmCell.Linf << nl
            << "u1      " << u1Cell.L1 << "   " << u1Cell.L2 << "   "
            << u1Cell.Linf << nl
            << "u2      " << u2Cell.L1 << "   " << u2Cell.L2 << "   "
            << u2Cell.Linf << nl
            << "-------------------------------------------------" << nl << nl
            << "MATLAB-like cell-centred error:" << nl
            << "Field     error_cell" << nl
            << "VmC      " << VmCell.errorCell << nl
            << "u1C      " << u1Cell.errorCell << nl
            << "u2C      " << u2Cell.errorCell << nl
            << "-------------------------------------------------" << nl << nl
            << "Gauss-reconstructed error summary:" << nl
            << "Field     L1-error       L2-error       Linf-error       normL2         relL2(%)" << nl
            << "VmG      " << VmHO.L1 << "   " << VmHO.L2 << "   "
            << VmHO.Linf << "   " << VmHO.normL2 << "   "
            << VmHO.relL2 << nl
            << "u1G      " << u1HO.L1 << "   " << u1HO.L2 << "   "
            << u1HO.Linf << "   " << u1HO.normL2 << "   "
            << u1HO.relL2 << nl
            << "u2G      " << u2HO.L1 << "   " << u2HO.L2 << "   "
            << u2HO.Linf << "   " << u2HO.normL2 << "   "
            << u2HO.relL2 << nl
            << "-------------------------------------------------" << nl << nl
            << "RHS summary at final time:" << nl
            << "Field     Linf-RHS" << nl
            << "VmR      " << linfNorm(rhsVm) << nl
            << "IionR    " << linfNorm(Iion) << nl
            << "-------------------------------------------------" << nl << nl
            << "Simulation summary:" << nl
            << "-------------------" << nl
            << "Final time            = " << runTime.value() << nl
            << "Number of cells (N)   = " << N << nl
            // N es el conteo POR DIRECCION. El total no se deriva de N salvo en
            // malla hexaedrica; en triangular/poligonal depende del mallador.
            // returnReduce porque mesh.nCells() es el conteo LOCAL de cada rank.
            << "Number of cells (total)= "
            << returnReduce(mesh.nCells(), sumOp<label>()) << nl
            << "Number of boundary faces = " << nBoundaryFacesGlobal << nl
            << "Dimension             = " << dimName << nl
            << "Grid spacing (dx)     = " << dx << nl
            << "Time step (dt)        = " << dt << nl
            << "Number of steps       = " << nSteps << nl
            << "Implicit scheme       = " << implicitScheme << nl
            // Which startup supplied the BDF history.
            // Recorded even when it is inert (theta schemes and BDF1 ignore
            // it) so that a run is self-describing: the observed temporal
            // order of a BDF run is meaningless without knowing this.
            << "BDF startup           = " << bdfStartup << nl
            // Same reasoning: under linear, every BDF is capped at order 2.
            << "State ODE driver      = " << stateODEDriver << nl
            // Which solution every error in this file is measured against.
            << "Manufactured solution = " << manufacturedSolution << nl
            << "Mass matrix           = " << massMatrixType << nl
            << "useHighOrder_Vm       = " << (useHighOrderVm ? "true" : "false") << nl
            << "Stabilisation alpha   = " << stabilisationAlpha << nl
            << "memoryOptimization    = " << memoryOptimization << nl
            << "memoryOptimization effective = "
            << (memoryOptimizationEffective ? "true" : "false") << nl
            << "memoryOptimization cellThreshold = "
            << memoryOptimizationCellThreshold << nl
            << "adaptive triplet reserve = "
            << (memoryOptimizationAdaptiveTripletReserve ? "true" : "false") << nl
            << "flux triplets/face reserve = "
            << memoryOptimizationFluxTripletsPerFace << nl
            << "stabilised triplets/face reserve = "
            << memoryOptimizationStabilisedTripletsPerFace << nl
            << "memoryOptimization trimHeap = "
            << (memoryOptimizationTrimHeap ? "true" : "false") << nl
            << "memoryOptimization compactMassAssembly = "
            << (memoryOptimizationCompactMassAssembly ? "true" : "false") << nl
            << "compact mass assembly = "
            << (compactMassAssembly ? "true" : "false") << nl
            << "stiffness triplets/face reserve = "
            << stiffnessTripletsPerFaceReserve << nl
            << "Vm LRE N              = " << lreN << nl
            << "Vm LRE Nn             = " << lreNn << nl
            << "Vm LRE k              = " << lreK << nl
            << "Vm LRE maxStencilSize = " << lreMaxStencilSize << nl
            << "useHighOrder_Iion     = " << (useHighOrderIion ? "true" : "false") << nl
            << "allocated Iion IPs    = " << allocatedIionIntegrationPoints << nl
            << "Iion LRE N            = " << lreIionN << nl
            << "Iion LRE Nn           = " << lreIionNn << nl
            << "Iion LRE k            = " << lreIionK << nl
            << "Iion LRE maxStencilSize = " << lreIionMaxStencilSize << nl
            << "Nonlinear method      = " << nonlinearMethod << nl
            << "Max nonlinear iterations = " << maxNonlinearIterations << nl
            << "Min nonlinear iterations = " << minNonlinearIterations << nl
            << "Linear backend        = " << linearSolverBackend << nl
            << "Linear PDE solver     = " << implicitLinearSolver << nl
            << "PETSc linear KSP      = " << petscLinearKspType << nl
            << "PETSc linear PC       = " << petscLinearPcType << nl
            << "Linear PDE tolerance  = " << implicitTolerance << nl
            << "Linear PDE maxIter    = " << implicitMaxIterations << nl
            << "State ODE solver      = " << stateODESolver << nl
            << "State ODE initial step= " << stateODEInitialStep << nl
            << "State ODE absTol      = " << stateODEAbsTol << nl
            << "State ODE relTol      = " << stateODERelTol << nl
            << "State ODE maxSteps    = " << stateODEMaxSteps << nl
            << "State ODE OpenMP      = "
            << (stateODEUseOpenMP ? "true" : "false") << nl
            << "State ODE threads     = " << stateODENumThreads << nl;

        if (isPicard)
        {
            os  << "Vm Picard tolerance   = " << nonlinearVmTolerance << nl
                << "State Picard tolerance= " << nonlinearStatesTolerance << nl
                << "Picard relaxation     = " << nonlinearRelaxation << nl;
        }
        else if (isJFNK)
        {
            os  << "JFNK coupled tolerance= " << nonlinearTolerance << nl
                << "JFNK Vm tolerance     = " << nonlinearVmTolerance << nl
                << "JFNK Iion tolerance   = " << nonlinearIionTolerance << nl
                << "JFNK require states   = "
                << (nonlinearRequireStatesConvergence ? "true" : "false") << nl;
            if (nonlinearRequireStatesConvergence)
            {
                os  << "JFNK state tolerance  = "
                    << nonlinearStatesTolerance << nl;
            }
            os  << "JFNK accept unconverged = "
                << (nonlinearAcceptUnconverged ? "true" : "false") << nl
                << "JFNK Newton relaxation= " << nonlinearRelaxation << nl
                << "JFNK GMRES m          = " << jfnkMaxKrylovIterations << nl
                << "JFNK GMRES restarts   = " << jfnkMaxRestarts << nl
                << "JFNK GMRES tolerance  = " << jfnkLinearTolerance << nl
                << "JFNK epsilon          = " << jfnkEpsilon << nl
                << "JFNK backend          = " << jfnkLinearSolverBackend << nl
                << "JFNK PETSc KSP        = " << jfnkPetscKspType << nl
                << "JFNK PETSc PC         = " << jfnkPetscPcType << nl
                << "JFNK initGuessOrder   = " << jfnkInitGuessOrder << nl
                << "JFNK initGuessVmMin   = " << jfnkInitGuessVmMin << nl
                << "JFNK initGuessVmMax   = " << jfnkInitGuessVmMax << nl
                << "JFNK clampODEInput    = "
                << (jfnkClampODEInput ? "true" : "false") << nl
                << "JFNK lineSearch       = "
                << (jfnkLineSearch ? "true" : "false") << nl;
            if (jfnkLineSearch)
            {
                os  << "JFNK lineSearch maxIt = " << jfnkLineSearchMaxIter << nl
                    << "JFNK lineSearch alphaMin = "
                    << jfnkLineSearchAlphaMin << nl
                    << "JFNK Armijo c         = " << jfnkArmijoC << nl;
            }
            os  << "JFNK preconditioner   = " << jfnkPreconditioner << nl
                << "JFNK PC dropTol       = "
                << jfnkPreconditionerDropTolerance << nl
                << "JFNK PC fillFactor    = "
                << jfnkPreconditionerFillFactor << nl
                << "JFNK PC updateFreq    = "
                << jfnkPreconditionerUpdateFrequency << nl;
        }
        else if (isDiagonalIion)
        {
            os  << "Diagonal coupled tolerance = " << nonlinearTolerance << nl
                << "Diagonal Vm tolerance = " << nonlinearVmTolerance << nl
                << "Diagonal Iion tolerance = " << nonlinearIionTolerance << nl
                << "Diagonal require states = "
                << (nonlinearRequireStatesConvergence ? "true" : "false") << nl;
            if (nonlinearRequireStatesConvergence)
            {
                os  << "Diagonal state tolerance = "
                    << nonlinearStatesTolerance << nl;
            }
            os  << "Diagonal accept unconverged = "
                << (nonlinearAcceptUnconverged ? "true" : "false") << nl
                << "Diagonal relaxation  = " << nonlinearRelaxation << nl
                << "Diagonal Iion epsilon= " << diagonalIionEpsilon << nl;
        }

        os  << "-------------------" << nl << nl
            << nonlinearMethod << " nonlinear convergence summary:" << nl
            << "final it/max          = "
            << (finalNonlinear ? finalNonlinear->iterations : 0)
            << "/" << maxNonlinearIterations << nl
            << "final converged       = "
            << (finalNonlinear && finalNonlinear->converged ? "true" : "false") << nl
            << "final rolledBack      = "
            << (finalNonlinear && finalNonlinear->rolledBack ? "true" : "false") << nl;

        if (!isPicard)
        {
            os  << "final coupled_relL2   = "
                << (finalNonlinear ? finalNonlinear->coupledResidual : GREAT) << nl;
        }

        os  << "final Vm_relL2        = "
            << (finalNonlinear ? finalNonlinear->VmResidual : GREAT) << nl
            << "final u1_relL2        = "
            << (finalNonlinear ? finalNonlinear->u1Residual : GREAT) << nl
            << "final u2_relL2        = "
            << (finalNonlinear ? finalNonlinear->u2Residual : GREAT) << nl
            << "final u3_relL2        = "
            << (finalNonlinear ? finalNonlinear->u3Residual : GREAT) << nl
            << "final maxState_relL2  = "
            << (finalNonlinear ? finalNonlinear->maxStateResidual : GREAT) << nl;

        if (!isPicard)
        {
            os  << "final Iion_relL2      = "
                << (finalNonlinear ? finalNonlinear->IionResidual : GREAT) << nl;
        }

        os  << "max it used           = " << maxNonlinearIterationsUsed << nl
            << "non-converged steps   = " << nonConvergedNonlinearSteps << nl
            << "rolled-back steps     = " << rolledBackNonlinearSteps << nl;

        if (!isPicard)
        {
            os  << "max coupled_relL2     = " << maxCoupledResidual << nl;
        }

        os  << "max Vm_relL2          = " << maxVmResidual << nl
            << "max u1_relL2          = " << maxU1Residual << nl
            << "max u2_relL2          = " << maxU2Residual << nl
            << "max u3_relL2          = " << maxU3Residual << nl
            << "maxState_relL2        = " << maxStateResidual << nl;

        if (!isPicard)
        {
            os  << "max Iion_relL2        = " << maxIionResidual << nl;
        }

        os  << "-------------------------------------------------" << nl << nl
            << "Computational resources:" << nl
            << "peakRSS_kB            = " << peakMemoryKB << nl
            << "peakRSS_MB            = " << peakMemoryKB/1024.0 << nl
            << "-------------------" << nl << nl
            << "Wall-clock times [s]:" << nl
            << "setup                 = " << setupWallTime << nl
            << "timeLoop              = " << timeLoopWallTime << nl
            << "postProcess           = " << postProcessWallTime << nl
            << "total                 = " << totalWallTime << nl
            // Phases inside timeLoop, so these are a SUBSET of it and must not
            // be added to the four above. nonlinearEval against
            // sparseLinearSolve is the measured split between the work OpenMP
            // parallelises and the work MPI does.
            << "nonlinearEval         = "
            << fineTimings.nonlinearEvalWallTime << nl
            << "nonlinearEvalCalls    = "
            << fineTimings.nonlinearEvalCalls << nl
            << "sparseLinearSolve     = "
            << fineTimings.sparseLinearSolveWallTime << nl
            << "sparseLinearSolveCalls = "
            << fineTimings.sparseLinearSolveCalls << nl
            << "jfnkGmres             = " << fineTimings.gmresWallTime << nl
            << "jfnkGmresCalls        = " << fineTimings.gmresCalls << nl
            << "jfnkPcSetup           = "
            << fineTimings.preconditionerSetupWallTime << nl
            << "jfnkPcSetups          = "
            << fineTimings.preconditionerSetups << nl
            << "jfnkPcApply           = "
            << fineTimings.preconditionerApplyWallTime << nl
            << "jfnkPcApplications    = "
            << fineTimings.preconditionerApplications << nl
            // Las dos mitades de la reconstruccion de estados. Subconjunto de
            // nonlinearEval, y solo distintas de cero en el camino
            // cellCentredReconstruct (states_CCpN).
            << "stateReconDeriv       = "
            << fineTimings.stateReconDerivWallTime << nl
            << "stateReconExpand      = "
            << fineTimings.stateReconExpandWallTime << nl
            << "stateReconCalls       = "
            << fineTimings.stateReconCalls << nl
            // Idem para la reconstruccion de Vm, que corre siempre que
            // useHighOrder_Iion este activo y no solo en cellCentredReconstruct.
            << "vmReconDeriv          = "
            << fineTimings.vmReconDerivWallTime << nl
            << "vmReconExpand         = "
            << fineTimings.vmReconExpandWallTime << nl
            << "vmReconCalls          = "
            << fineTimings.vmReconCalls << nl
            << "-------------------" << nl;

        Info<< "Wrote summary to " << outFile << nl
            << "Vm error: L1 = " << VmCell.L1
            << ", L2 = " << VmCell.L2
            << ", Linf = " << VmCell.Linf << nl
            << "Vm Gauss error: L1 = " << VmHO.L1
            << ", L2 = " << VmHO.L2
            << ", Linf = " << VmHO.Linf
            << ", Relative = " << VmHO.relL2 << "%" << nl
            << "Final-time RHS Linf = " << linfNorm(rhsVm) << nl
            << "Timing [s]: setup = " << setupWallTime
            << ", loop = " << timeLoopWallTime
            << ", post = " << postProcessWallTime
            << ", total = " << totalWallTime << endl;
    }
}

int main(int argc, char* argv[])
{
#ifdef __GLIBC__
    // Tighten the glibc allocator before any large allocation. The LRE
    // construction runs ~2.3M QR factorisations on a 3D unstructured tet
    // mesh with N=40 / p3; without these tunables glibc grows the heap by
    // several GB of unreclaimed fragments and the run OOMs on a 16 GB box.
    //   - M_ARENA_MAX=2: prevents glibc from creating ~64 arenas (one set
    //     per ODE OpenMP thread). LRE setup is serial; multiple arenas
    //     just inflate the RSS.
    //   - M_MMAP_THRESHOLD=64 KB: routes mid-sized allocations (Eigen
    //     temporaries during QR) through mmap so they are returned to
    //     the kernel individually when freed.
    //   - M_TRIM_THRESHOLD=64 KB: lets glibc release sbrk-heap padding
    //     back to the OS more aggressively.
    mallopt(M_ARENA_MAX, 2);
    mallopt(M_MMAP_THRESHOLD, 64*1024);
    mallopt(M_TRIM_THRESHOLD, 64*1024);
#endif

    #include "setRootCaseLists.H"
    PetscSession petscSession(argc, argv);
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    const auto tStartTotal = std::chrono::steady_clock::now();
    const auto tStartSetup = tStartTotal;

    const label dim = mesh.nGeometricD();
    const scalar dt = runTime.deltaTValue();

    // Fix the manufactured profile before anything
    // evaluates an exact value - computeBeta on the next line is the first
    // reader. Validated rather than defaulted, because a typo would silently
    // run the other manufactured solution and every error in the run would be
    // measured against the wrong reference.
    if (manufacturedSolution == "uniform")
    {
        mmsProfile = mmsUniform;
    }
    else if (manufacturedSolution != "sines")
    {
        FatalErrorInFunction
            << "Unknown manufacturedSolution: " << manufacturedSolution << nl
            << "Valid options are sines and uniform."
            << exit(FatalError);
    }

    const scalar beta = computeBeta(conductivity, dim);
    const scalar chiVal = chi.value();
    const scalar CmVal = Cm.value();
    const scalar chiCmVal = chiVal*CmVal;
    const scalar lapScale = 1.0/max(chiCmVal, SMALL);
    const TimeSchemeCoeffs timeScheme = timeSchemeCoeffs(implicitScheme);
    const scalar theta = timeScheme.theta;
    const label bdfOrder = timeScheme.bdfOrder;
    const word massMatrixMode = normalizedMassMatrixType(massMatrixType);

    // ESDIRK exists in this solver for ONE purpose: supplying the k-1 starting
    // values a BDF-k step needs. It is deliberately not selectable as a scheme
    // in its own right.
    //
    // The reason is a measurement, not taste. ESDIRK is a one-step method, so
    // it is self-starting and handles variable dt trivially - but it has stage
    // order 2 and therefore suffers order reduction on stiff problems: ESDIRK4
    // measured 3.03 here, capped, while BDF4 reaches 3.89 by Cauchy self-
    // convergence (suite test M02) and 3.94 against the analytic solution
    // (route B of M02). The ESDIRK4 figure predates the time-loop fix and was never
    // re-measured, so read it as the order it showed and not as a comparison
    // digit for digit. BDF has stage order equal to its order and no such cap. And in the regime this solver's
    // physiological sibling lives in, the ODE integrations dominate the runtime
    // (74-99 % measured), where a 4-stage ESDIRK costs about 2.5x the ODE work
    // of a BDF step. So BDF wins on accuracy AND on the cost that matters, and
    // the only thing ESDIRK is uniquely good at here is the startup.
    //
    // ESDIRK3 is what is kept, and it is sufficient: a one-step method of order
    // p leaves starting values accurate to O(dt^(p+1)), so ESDIRK3 covers BDF
    // up to k = 4, which is the highest BDF offered. An ESDIRK4 startup was
    // implemented and measured against it for BDF4 - 3.02 against 3.20, both
    // measured before the time-loop fix and so comparable only to each other,
    // i.e. no gain, because its own order reduction leaves it with the same
    // O(dt^4) local error - so it was removed rather than kept as dead weight.
    const bool esdirkStartup =
        (bdfStartup == "ESDIRK3");

    if
    (
        bdfOrder > 0
     && bdfStartup != "exact"
     && bdfStartup != "constant"
     && !esdirkStartup
    )
    {
        FatalErrorInFunction
            << "Unknown bdfStartup: " << bdfStartup << nl
            << "Valid options are exact, constant and ESDIRK3."
            << exit(FatalError);
    }

    // Validated rather than falling back: a typo would silently select the
    // ramp or the interpolant, and the only symptom is a different order.
    if (stateODEDriver != "history" && stateODEDriver != "linear")
    {
        FatalErrorInFunction
            << "Unknown stateODEDriver: " << stateODEDriver << nl
            << "Valid options are history and linear."
            << exit(FatalError);
    }

    const bool historyStateODEDriver = (stateODEDriver == "history");

    const bool needEsdirkOperator = (bdfOrder > 1 && esdirkStartup);
    const EsdirkTableau esdirkTab = esdirk3Tableau();

    // Invariant rather than a reachable branch today: ESDIRK3 covers every BDF
    // this solver offers. It is kept because adding BDF5 would break it
    // silently - a startup of insufficient order does not fail, it caps the
    // observed order, which is exactly the trap the exact/constant seeding
    // modes fell into and which cost a full diagnosis to find.
    if (bdfOrder > 1 && esdirkStartup && esdirkTab.order + 1 < bdfOrder)
    {
        FatalErrorInFunction
            << "bdfStartup = " << bdfStartup << " supplies starting values of"
            << " order " << (esdirkTab.order + 1) << ", which is below the"
            << " O(dt^" << bdfOrder << ") that " << implicitScheme
            << " requires." << nl
            << "Use a higher-order startup, or the observed order will be"
            << " capped without any error being reported."
            << exit(FatalError);
    }
    // Compared exactly, no longer lower-cased first.
    const word& memoryOptimizationMode = memoryOptimization;

    const bool memoryOptimizationEffective = memoryOptimizationIsEffective
    (
        memoryOptimizationMode,
        useHighOrder_Vm,
        dim,
        mesh,
        memoryOptimizationCellThreshold
    );

    const bool useAdaptiveTripletReserve =
        memoryOptimizationEffective
     && memoryOptimizationAdaptiveTripletReserve;

    const label stiffnessTripletsPerFaceReserve =
            useAdaptiveTripletReserve
        ? (
                stabilisationAlpha > SMALL
              ? memoryOptimizationStabilisedTripletsPerFace
              : memoryOptimizationFluxTripletsPerFace
          )
        : label(800);

    const bool trimHeapAfterLargeSetup =
        memoryOptimizationEffective && memoryOptimizationTrimHeap;

    const bool compactMassAssembly =
        memoryOptimizationEffective && memoryOptimizationCompactMassAssembly;

    const label allocatedIionIntegrationPoints =
            useHighOrder_Iion && dim > 1
        ? totalIionIntegrationPoints
        : label(0);


    scalarField VmOldIntegrationPoints(allocatedIionIntegrationPoints, 0.0);
    scalarField VmGuessIntegrationPoints(allocatedIionIntegrationPoints, 0.0);
    scalarField sourceIntegrationPoints(allocatedIionIntegrationPoints, 0.0);
    scalarField IionIntegrationPoints(allocatedIionIntegrationPoints, 0.0);
    scalarField u1IntegrationPoints(allocatedIionIntegrationPoints, 0.0);
    scalarField u2IntegrationPoints(allocatedIionIntegrationPoints, 0.0);
    scalarField u3IntegrationPoints(allocatedIionIntegrationPoints, 0.0);

    Info<< "Running high-order manufactured FDA implicit PETSc solver" << nl
        << "Dimension = " << dim << nl
        << "beta = " << beta << nl
        << "dt = " << dt << nl
        << "implicitScheme = " << implicitScheme << nl
        << "bdfStartup = " << bdfStartup << nl
        << "stateODEDriver = " << stateODEDriver << nl
        << "stateODEGaussPointHistory = " << stateODEGaussPointHistory << nl
        << "manufacturedSolution = " << manufacturedSolution << nl
        << "massMatrix = " << massMatrixMode << nl
        << "linearSolverBackend = " << linearSolverBackend << nl
        << "PETSc linear KSP = " << petscLinearKspType
        << ", PC = " << petscLinearPcType << nl
        << "nonlinearMethod = " << nonlinearMethod << nl
        << "JFNK linear backend = " << jfnkLinearSolverBackend << nl
        << "JFNK PETSc KSP = " << jfnkPetscKspType
        << ", PC = " << jfnkPetscPcType << nl
        << "stateODESolver = " << stateODESolver << nl
        << "stateODEUseOpenMP = " << stateODEUseOpenMP
        << ", stateODENumThreads = " << stateODENumThreads << nl
        << "nonlinearRelaxation = " << nonlinearRelaxation << nl
        << "useHighOrder_Vm = " << useHighOrder_Vm << nl
        << "useHighOrder_Iion = " << useHighOrder_Iion << nl
        << "stateIntegrationMode = " << stateIntegrationMode << nl
        << "stabilisationAlpha = " << stabilisationAlpha << nl
        << "memoryOptimization = " << memoryOptimization
        << " (effective: "
        << (memoryOptimizationEffective ? "true" : "false") << ")" << nl
        << "stiffnessTripletsPerFaceReserve = "
        << stiffnessTripletsPerFaceReserve << nl
        << "compactMassAssembly = "
        << (compactMassAssembly ? "true" : "false") << nl
        << "jfnkPreconditioner = " << jfnkPreconditioner << nl
        << "jfnkPreconditionerUpdateFrequency = "
        << jfnkPreconditionerUpdateFrequency << nl
        << "Iion integration points = " << totalIionIntegrationPoints
        << " (allocated: " << allocatedIionIntegrationPoints << ")"
        << endl;

    if (reconstructStatesFromCellCentres && useHighOrder_Iion && dim > 1)
    {
        Info<< "LREInterp_states order = " << LREInterp_statesPtr().order()
            << " (states evolved at cell centres, reconstructed at Gauss points)"
            << endl;
    }

    // One spelling per method. picard, jfnk, diagonal
    // and localDiagonal were accepted as aliases and are now refused.
    const bool usePicard = nonlinearMethod == "Picard";
    const bool useDiagonalIion = nonlinearMethod == "diagonalIion";
    const bool useJFNK = nonlinearMethod == "JFNK";

    if (!usePicard && !useDiagonalIion && !useJFNK)
    {
        FatalErrorInFunction
            << "Unknown nonlinearMethod '" << nonlinearMethod << "'. "
            << "Valid options are Picard, JFNK and diagonalIion (exact"
            << " spelling, case-sensitive)."
            << exit(FatalError);
    }

    fillExactFields(VmExact, u1Exact, u2Exact, runTime.value(), dim);

    Vm.primitiveFieldRef() = VmExact.primitiveField();
    u1.primitiveFieldRef() = u1Exact.primitiveField();
    u2.primitiveFieldRef() = u2Exact.primitiveField();
    u3 = dimensionedScalar("zero", dimless, 0.0);

    Vm.correctBoundaryConditions();
    updateStateBoundaryValues(u1, u2, u3, runTime.value(), dim);
    u1.correctBoundaryConditions();
    u2.correctBoundaryConditions();
    u3.correctBoundaryConditions();

    if (useHighOrder_Iion && dim > 1)
    {
        initialiseIionIntegrationPointStates
        (
            LREInterp_Iion,
            runTime.value(),
            dim,
            u1IntegrationPoints,
            u2IntegrationPoints,
            u3IntegrationPoints
        );
    }

    // Publish this rank's global row offset for the
    // linear-solver wrappers. The matrices carry local rows and global columns,
    // so every PETSc insertion has to shift the row index by this amount. Zero
    // in serial, which is why the serial path is unaffected.
    initialiseGlobalCells(mesh);

    Info<< "Distributed linear algebra: rows ["
        << gRowStart << ", " << gRowStart + mesh.nCells() << ") of "
        << LREInterp_Vm.globalCells().totalSize() << " global cells" << endl;

    SpMat M;
    SpMat K;
    SpMat AImplicit;
    SpMat BImplicit;

    // Reject Vm boundary conditions the assembly cannot
    // represent, at startup, instead of producing a silently wrong answer.
    // Mirrors the check in the electro solver.
    //
    // The patch dispatchers key off the Vm patch field's TYPE NAME and only
    // branch on the types below. Anything else used to fall out of the bottom
    // of the `if`, and - measured, not inferred - it did so in two DIFFERENT
    // ways in this solver:
    //
    //   useHighOrder_Vm = true   `fixedGradient` changed the answer, because
    //                            the high-order loop applies the Dirichlet-style
    //                            flux operator with no right-hand side, i.e. a
    //                            disguised homogeneous Dirichlet;
    //   useHighOrder_Vm = false  `fixedGradient` and `calculated` were ignored
    //                            entirely, because the low-order dispatcher has
    //                            no `else` at all, i.e. a disguised zero Neumann.
    //
    // The same case therefore produced two different wrong answers depending on
    // a switch that is supposed to change only the order of accuracy. Both ran
    // to completion, exit 0, no warning.
    //
    // Checked on every rank: `processor` patches only exist in parallel, and a
    // rank-local list would disagree between ranks.
    {
        const wordHashSet supportedVmPatchTypes
        ({
            word("empty"),
            zeroGradientFvPatchScalarField::typeName,
            word("processor")
        });

        wordHashSet offending;
        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            const word bcType =
                patch.lookupPatchField<volScalarField, scalar>("Vm").type();

            if (!supportedVmPatchTypes.found(bcType))
            {
                offending.insert(bcType);
            }
        }

        if (!offending.empty())
        {
            FatalErrorInFunction
                << "Unsupported boundary condition(s) on Vm: "
                << offending.sortedToc() << nl
                << "The assembly only implements " << supportedVmPatchTypes.sortedToc()
                << nl
                << "Any other type is silently ignored by the operator, which"
                << " produces a plausible but wrong answer rather than an error."
                << nl
                << "fixedValue/fixedVoltage are not supported: Dirichlet"
                << " conditions on Vm were removed from this solver and from the"
                << " physiological one alike."
                << exit(FatalError);
        }
    }

    // Report the chunk sizes and how many chunks they
    // actually produce on THIS rank. Without the second number a chunking test
    // cannot tell whether it exercised the multi-chunk path at all - the counts
    // are local, so a mesh that needs two chunks in serial can need one per rank
    // in parallel.
    Info<< "Assembly chunking     : faceChunk = " << assemblyFaceChunk
        << ", cellChunk = " << assemblyCellChunk << endl;
    Pout<< "Assembly chunking     : nInternalFaces = " << mesh.nInternalFaces()
        << " -> " << (mesh.nInternalFaces() + assemblyFaceChunk - 1)
                     /assemblyFaceChunk
        << " face chunk(s); nCells = " << mesh.nCells()
        << " -> " << (mesh.nCells() + assemblyCellChunk - 1)/assemblyCellChunk
        << " cell chunk(s)" << endl;

    if (profileTimings)
    {
        logMemoryCheckpoint("before mass assembly");
    }

    if (!useHighOrder_Vm || massMatrixMode == "lumped")
    {
        assembleDiagonalMassMatrix(mesh, 1.0, M);
    }
    else
    {
        assembleConsistentMassMatrixHO
        (
            mesh,
            LREInterp_Vm,
            1.0,
            compactMassAssembly,
            assemblyCellChunk,
            M
        );
    }

    if (profileTimings)
    {
        logMemoryCheckpoint("after mass assembly");
    }

#ifdef __GLIBC__
    // Reclaim the heap fragmentation left by the ~290k per-cell QR
    // factorisations triggered when LREInterp_Vm populated its cell
    // gradient / Hessian / 3rd-derivative coefficient tables.
    if (trimHeapAfterLargeSetup)
    {
        malloc_trim(0);
    }
#endif

    if (profileTimings)
    {
        logMemoryCheckpoint("before stiffness assembly");
    }

    if (useHighOrder_Vm)
    {
        assembleHighOrderStiffnessMatrix
        (
            mesh,
            conductivity,
            LREInterp_Vm,
            stabilisationAlpha,
            stiffnessTripletsPerFaceReserve,
            assemblyFaceChunk,
            K
        );
    }
    else
    {
        assembleStandardOrthogonalStiffnessMatrix
        (
            mesh,
            conductivity,
            LREInterp_Vm,
            stabilisationAlpha,
            K
        );
    }

    if (profileTimings)
    {
        logMemoryCheckpoint("after stiffness assembly");
    }

    // Consistency check on the assembled diffusion
    // operator. A constant field has zero gradient, so with no-flux boundaries
    // every row of K must sum to zero. This caught a silent failure during the
    // LRE -> movingLeastSquares port: CompactListList::operator[] returns a
    // UList VIEW by value, and binding it to a const labelList& compiled but
    // produced an empty stencil, so K came out with no entries at all. The
    // symptom was a fixed error that did not reduce under refinement, which is
    // easy to mistake for a discretisation problem. Cheap (one pass over the
    // non-zeros) and it doubles as the conservation probe in parallel.
    {
        scalar maxRowSum = 0.0;
        for (label r = 0; r < K.outerSize(); ++r)
        {
            scalar rowSum = 0.0;
            for (SpMat::InnerIterator it(K, r); it; ++it)
            {
                rowSum += it.value();
            }
            maxRowSum = max(maxRowSum, mag(rowSum));
        }
        reduce(maxRowSum, maxOp<scalar>());

        const label nnzGlobal = returnReduce(label(K.nonZeros()), sumOp<label>());

        Info<< "Stiffness matrix: " << nnzGlobal << " non-zeros, "
            << "max |row sum| = " << maxRowSum
            << " (must be ~0: K must annihilate constants)" << endl;

        if (nnzGlobal == 0)
        {
            FatalErrorInFunction
                << "The stiffness matrix is empty. The diffusion operator would"
                << " be absent from the system." << abort(FatalError);
        }
    }

    // With CF_DUMP_STENCILS set, write the MLS stencil
    // MEMBERSHIP - one line per cell, "globalCell : global columns" - so that a
    // serial and a parallel build can be diffed cell by cell. This is one level
    // below CF_DUMP_K: it says whether a differing operator row comes from
    // different stencil members or from different coefficients over the same
    // members.
    if (getenv("CF_DUMP_STENCILS"))
    {
        const globalIndex& gcD = LREInterp_Vm.globalCells();
        const CompactListList<label>& cellSt = LREInterp_Vm.globalCellStencils();

        OFstream os
        (
            "stencil_dump_rank" + Foam::name(Pstream::myProcNo()) + ".dat"
        );
        // Cell centres alongside, at full precision: the selection is a
        // nearest-N by distance, so a centroid that differs in the last bits
        // between a serial and a decomposed mesh is enough to reorder a tie.
        OFstream osc
        (
            "cellCentre_dump_rank" + Foam::name(Pstream::myProcNo()) + ".dat"
        );
        osc.precision(17);
        forAll(mesh.C(), cellI)
        {
            osc << gcD.toGlobal(cellI) << ' ' << mesh.C()[cellI].x() << ' '
                << mesh.C()[cellI].y() << ' ' << mesh.C()[cellI].z() << nl;
        }

        forAll(mesh.C(), cellI)
        {
            os  << gcD.toGlobal(cellI) << " :";
            const UList<label> st = cellSt[cellI];
            forAll(st, i)
            {
                os  << ' ' << st[i];
            }
            os  << nl;
        }

        // Face stencils, keyed by LOCAL face index: the reader maps them with
        // processorN/constant/polyMesh/faceProcAddressing.
        const CompactListList<label>& faceSt = LREInterp_Vm.globalFaceStencils();

        OFstream osf
        (
            "faceStencil_dump_rank" + Foam::name(Pstream::myProcNo()) + ".dat"
        );
        forAll(faceSt, faceI)
        {
            osf << faceI << " :";
            const UList<label> st = faceSt[faceI];
            forAll(st, i)
            {
                osf << ' ' << st[i];
            }
            osf << nl;
        }
    }

    // With CF_DUMP_K set in the environment, write K as
    // (global row, global column, value) triples, one file per rank.
    //
    // This is the check that a row-sum or nnz fingerprint cannot make: mapping
    // the parallel indices back through processorN/constant/polyMesh/
    // cellProcAddressing gives an entry-by-entry comparison against the serial
    // assembly, which is what separates a genuine discretisation difference
    // from partition-dependent preconditioning. Off by default, no cost.
    if (getenv("CF_DUMP_K"))
    {
        OFstream os
        (
            "K_dump_rank" + Foam::name(Pstream::myProcNo()) + ".dat"
        );
        os.precision(17);
        for (label r = 0; r < K.outerSize(); ++r)
        {
            for (SpMat::InnerIterator it(K, r); it; ++it)
            {
                os  << (globalRowStart() + r) << ' ' << it.col() << ' '
                    << it.value() << nl;
            }
        }
    }

#ifdef __GLIBC__
    // Reclaim the heap fragmentation from the ~2.3M per-(face, GP) QR
    // factorisations during LREInterp_Vm.QRGradFaceGPCoeffs() population.
    if (trimHeapAfterLargeSetup)
    {
        malloc_trim(0);
    }
#endif

    // The only thing BDF changes in the operator is
    // the scalar in front of M. a_0 is 1 for the theta family and for BDF1, so
    // backwardEuler and BDF1 build the same matrix down to the bit.
    AImplicit = M;
    AImplicit *= (timeScheme.a[0]/dt);
    AImplicit -= (theta*lapScale)*K;

    // The ESDIRK stage operator. Because the method is
    // SINGLY diagonally implicit, every implicit stage has the same gamma on
    // its diagonal, so one matrix and one factorisation serve all of them -
    // this is the whole reason to insist on "singly". Built here, before M is
    // moved into BImplicit below.
    SpMat AEsdirk;

    if (needEsdirkOperator)
    {
        AEsdirk = M;
        AEsdirk *= (1.0/(esdirkTab.gamma*dt));
        AEsdirk -= lapScale*K;
    }

    // Move M into BImplicit instead of copying: M is not referenced again
    // after this block, so we can hand its storage over and free the
    // ~280 MB it occupied as a separate matrix.
    BImplicit = std::move(M);
    BImplicit *= (1.0/dt);

    if (theta < 1.0 - SMALL)
    {
        BImplicit += ((1.0 - theta)*lapScale)*K;
    }

    // PETSc mirrors of the two implicit operators, used
    // for every matrix-vector product. Eigen cannot do A*x any more, because A
    // has local rows and global columns while x only holds this rank's block;
    // MatMult performs the halo exchange internally. Built once here, since
    // AImplicit and BImplicit are constant while dt is (dt only changes on a
    // final clipped step, which this solver does not take).
    DistributedMatVec applyA;
    DistributedMatVec applyB;
    DistributedMatVec applyK;
    // The ESDIRK stage operator as a product, for the stage residual.
    DistributedMatVec applyEsdirk;
    if (needEsdirkOperator)
    {
        applyEsdirk.reset(AEsdirk);
    }
    applyA.reset(AImplicit);
    applyB.reset(BImplicit);
    applyK.reset(K);

    if (profileTimings)
    {
        logMemoryCheckpoint("after implicit matrix assembly");
    }

#ifdef __GLIBC__
    if (trimHeapAfterLargeSetup)
    {
        malloc_trim(0);
    }
#endif

    // The KSP and PC the two cached PETSc solvers below
    // actually use, resolved once from linearSolver, petscLinearKspType and
    // petscLinearPcType. They must NOT take the raw dictionary values - see
    // resolvePetscKspPc for the silent wrong answer that produced. Logged,
    // because the banner prints the names as requested, before resolution, and
    // those are not what PETSc runs.
    word petscLinearKspTypeResolved;
    word petscLinearPcTypeResolved;
    resolvePetscKspPc
    (
        implicitLinearSolver,
        petscLinearKspType,
        petscLinearPcType,
        petscLinearKspTypeResolved,
        petscLinearPcTypeResolved
    );

    Info<< "PETSc linear solve, resolved: KSP = " << petscLinearKspTypeResolved
        << ", PC = " << petscLinearPcTypeResolved << endl;

    // The persistent solvers of the ESDIRK stage operator, separate from the
    // BDF operator's so that neither factorisation is thrown away. Built on
    // the first stage solve (shared solveSparseSystem); both backends work.
    PetscKspMatrixSolver esdirkPetscSolver;
    PersistentEigenSolvers esdirkEigenSolvers;
    bool esdirkSolverNeedsRebuild = true;

    // The persistent solvers of the assembled solves (Picard and diagonalIion),
    // shared with the physiological solver: the factorisation lives for the
    // whole run and each solve decides whether to rebuild it (first solve),
    // refresh the values (diagonalIion moves the diagonal every iteration) or
    // reuse it untouched (Picard, constant operator). The MMS used to keep a
    // cached solver for Picard on PETSc only and to build a NEW solver - full
    // setup and factorisation - on every diagonalIion iteration.
    PetscKspMatrixSolver persistentPetscSolver;
    PersistentEigenSolvers persistentEigenSolvers;
    bool linearSolverNeedsRebuild = true;

    // The JFNK preconditioner and the controls of the JFNK linear solve,
    // shared with the physiological solver (monodomainKrylov.H).
    JfnkPreconditioner jfnkPc
    (
        jfnkPreconditioner,
        usesPetscBackend(jfnkLinearSolverBackend),
        jfnkPreconditionerKspType,
        jfnkPreconditionerPcType,
        jfnkPreconditionerTolerance,
        jfnkPreconditionerMaxIterations,
        jfnkMaxKrylovIterations,
        jfnkPreconditionerOptionsPrefix,
        petscUseOptions,
        jfnkPreconditionerFillFactor,
        jfnkPreconditionerDropTolerance,
        jfnkPreconditionerUpdateFrequency
    );

    JfnkKrylovControls jfnkKrylov;
    jfnkKrylov.petsc = usesPetscBackend(jfnkLinearSolverBackend);
    jfnkKrylov.kspType = jfnkPetscKspType;
    jfnkKrylov.pcType = jfnkPetscPcType;
    jfnkKrylov.restart = jfnkPetscRestart;
    jfnkKrylov.maxKrylovIterations = jfnkMaxKrylovIterations;
    jfnkKrylov.maxRestarts = jfnkMaxRestarts;
    jfnkKrylov.tolerance = jfnkLinearTolerance;
    jfnkKrylov.optionsPrefix = jfnkPetscOptionsPrefix;
    jfnkKrylov.useOptions = petscUseOptions;

    // Persistent JFNK shell-Krylov solver (PETSc MatShell + KSP). The Krylov
    // subspace and KSP structures are allocated once and reused across every
    // Newton iteration, avoiding the ~25-40 MB alloc/free cycle that was
    // happening once per Newton step in the old solvePetscShellSystem call.
    // The mat-vec / apply-PC callback payloads are refreshed per solve().
    PetscShellKspSolver jfnkPetscShellSolver;

    autoPtr<OFstream> nonlinearResidualFilePtr;
    // Master only, global case directory. See the
    // transient summary above.
    if (writeNonlinearResiduals && Pstream::master())
    {
        const fileName residualDir
        (
            runTime.rootPath()/runTime.globalCaseName()
          / "postProcessing"/"highOrderManufacturedFDAImplicitPETScDistributed"
        );
        mkDir(residualDir);
        nonlinearResidualFilePtr.reset
        (
            new OFstream(residualDir/"nonlinearResiduals.dat")
        );
        nonlinearResidualFilePtr()
            << "# time_s step nonlinearMethod stateODESolver iter linearIterations linearError "
            << "newtonResidual Vm_relL2 u1_relL2 u2_relL2 u3_relL2 Iion_relL2 "
            << "lineSearchIters converged"
            << nl;
    }

    auto updateNumericalLaplacian =
    [&](const scalar evalTime, const bool updateFluxField)
    {
        if (useHighOrder_Vm && updateFluxField)
        {
            computeHighOrderLaplacian
            (
                Vm,
                conductivity,
                LREInterp_Vm,
                fluxVm_HO,
                lapVm
            );
        }

        // Distributed product, see DistributedMatVec.
        EigVec lapNow = applyK(fieldToEigVec(Vm));
        eigVecToField(lapNow, lapVm);
        lapVm.correctBoundaryConditions();
    };

    const auto tEndSetup = std::chrono::steady_clock::now();
    const auto tStartLoop = tEndSetup;

    label nSteps = 0;
    scalarField VmTwoStepsAgo(mesh.nCells(), 0.0);
    scalar dtPrev = 0.0;
    bool hasPrevStep = false;
    std::vector<NonlinearConvergenceRecord> nonlinearHistory;

    // One home for the in-loop timers, so that the summary cannot report a
    // different number from the one the loop accumulated.
    FineGrainedTimings fineTimings;

    // BDF history, i.e. the k-1 levels OLDER than V^n
    // that a BDF-k step needs. VmBdfOlder[j] holds V^{n-1-j}, so index 0 is the
    // most recent of them. V^n itself is NOT duplicated here: it is VmOld,
    // rebuilt at the top of every step. The list is empty for the theta family
    // and for BDF1, neither of which looks back further than one level.
    List<scalarField> VmBdfOlder(max(bdfOrder - 1, label(0)));

    // How many of those levels are valid on the current piece of the solution.
    // A BDF step needs all k-1; until then ESDIRK3 takes the step. Seeded full
    // by the exact/constant startups, empty by ESDIRK3, and reset to zero by a
    // rolled-back step.
    label bdfValidLevels = esdirkStartup ? label(0) : max(bdfOrder - 1, label(0));

    // Sized and zeroed up front so that a level which has not been filled yet
    // is a wrong ANSWER rather than a silently skipped loop over an empty
    // field. Nothing reads them before they are filled: with the ESDIRK
    // startup the first bdfOrder-1 steps are taken by ESDIRK and the rotation
    // at the end of each fills one level.
    forAll(VmBdfOlder, j)
    {
        VmBdfOlder[j].setSize(mesh.nCells(), 0.0);
    }

    if (bdfOrder > 1 && !esdirkStartup)
    {
        const vectorField& Cbdf = mesh.C();
        const scalar t0 = runTime.value();

        forAll(VmBdfOlder, j)
        {
            VmBdfOlder[j].setSize(mesh.nCells());

            if (bdfStartup == "exact")
            {
                // The manufactured solution one further step back per level.
                // sqrt(1+t) is defined for t > -1, so a handful of steps before
                // t0 = 0 is not a problem. This is what isolates the order of
                // the SCHEME from the order of the startup, which is the whole
                // reason a verification solver gets to do this.
                const scalar tj = t0 - scalar(j + 1)*dt;

                forAll(VmBdfOlder[j], cellI)
                {
                    VmBdfOlder[j][cellI] = exactVm(Cbdf[cellI], tj, dim);
                }
            }
            else
            {
                // "constant": repeat the initial field. Needs no exact
                // solution - which is what a physical solver would face - and
                // is wrong by O(dt) per missing level. Offered so that cost is
                // measurable instead of assumed.
                VmBdfOlder[j] = Vm.primitiveField();
            }
        }

        Info<< "BDF history seeded (" << bdfStartup << "): "
            << VmBdfOlder.size() << " level(s) before t = " << t0 << endl;
    }
    else if (bdfOrder > 1)
    {
        Info<< "BDF history will be built by ESDIRK3 over the first "
            << (bdfOrder - 1) << " step(s)" << endl;
    }

    // Stop within half a step of endTime, not within
    // SMALL. The accumulated time is a floating-point sum of dt, and it can land
    // just BELOW endTime by more than SMALL (1e-15): 160 steps of 3.2e-3 sum to
    // 0.51199999999999868, which is 1.3e-15 short of 0.512. The old test then
    // took one extra step and the run finished at 0.5152 with no warning - so a
    // field or an error "at endTime" belonged to a different instant, and a
    // temporal-order study comparing runs at different dt compared fields at
    // different times. Half a step is exact because dt is constant here, and it
    // leaves every run that was already correct with the same step count:
    // checked across 159 151 existing summaries, none of which had an extra step.
    // Nonlinear iterations that stopped at the residual floor instead of the
    // tolerance (nonlinearIterationExit, monodomainBase.H).
    label nResidualFloorExits = 0;

    // Rolled-back steps in a row (nonlinearMaxConsecutiveRollbacks).
    label nConsecutiveRollbacks = 0;

    while (runTime.value() < runTime.endTime().value() - 0.5*dt)
    {
        const scalar t = runTime.value();

        volScalarField VmOld
        (
            IOobject
            (
                "VmOld",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Vm
        );

        volScalarField u1Old
        (
            IOobject
            (
                "u1Old",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            u1
        );

        volScalarField u2Old
        (
            IOobject
            (
                "u2Old",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            u2
        );

        volScalarField u3Old
        (
            IOobject
            (
                "u3Old",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            u3
        );

        scalarField u1OldIntegrationPoints(u1IntegrationPoints);
        scalarField u2OldIntegrationPoints(u2IntegrationPoints);
        scalarField u3OldIntegrationPoints(u3IntegrationPoints);
        scalarField u1GuessIntegrationPoints(u1OldIntegrationPoints);
        scalarField u2GuessIntegrationPoints(u2OldIntegrationPoints);
        scalarField u3GuessIntegrationPoints(u3OldIntegrationPoints);

        if (useHighOrder_Iion && dim > 1)
        {
            reconstructVmAtIionIntegrationPoints
            (
                VmOld,
                useHighOrder_Vm,
                LREInterp_Vm,
                LREInterp_Iion,
                VmOldIntegrationPoints,
                stateODEUseOpenMP,
                stateODEOpenMPThreshold,
                profileTimings,
                fineTimings
            );

            computeVmSourceFromIntegrationPoints
            (
                VmOldIntegrationPoints,
                u1OldIntegrationPoints,
                u2OldIntegrationPoints,
                u3OldIntegrationPoints,
                beta,
                chiVal,
                CmVal,
                sourceIntegrationPoints,
                stateODEUseOpenMP,
                stateODEOpenMPThreshold
            );

            averageIntegrationPointFieldToCells
            (
                sourceIntegrationPoints,
                LREInterp_Iion,
                sourceVm
            );
        }
        else
        {
            computeCellCentredVmSource
            (
                VmOld,
                u1Old,
                u2Old,
                u3Old,
                beta,
                chiVal,
                CmVal,
                sourceVm,
                stateODEUseOpenMP,
                stateODEOpenMPThreshold
            );
        }

        // The history vector handed to B. For the
        // theta family it is V^n; for BDF-k it is -sum_{j>=1} a_j V^{n+1-j},
        // formed HERE as one field so that a BDF step still costs a single
        // mass-matrix product rather than k of them.
        //
        // For BDF1, a_1 = -1 and there are no older levels, so the scaling is
        // by exactly 1.0 and the result is bit-identical to the theta path.
        EigVec Vn = fieldToEigVec(VmOld);

        if (bdfOrder > 0)
        {
            Vn *= -timeScheme.a[1];

            forAll(VmBdfOlder, j)
            {
                const scalar aj = timeScheme.a[j + 2];
                const scalarField& Vj = VmBdfOlder[j];

                forAll(Vj, cellI)
                {
                    Vn[cellI] -= aj*Vj[cellI];
                }
            }
        }

        // Temporary diagnostic, env-gated: is each history level actually the
        // solution at the time level the BDF formula assumes? Prints the max
        // deviation of level j from the exact solution at t - (j+1) dt, which
        // is zero to spatial-error level if the indexing is right and O(dt) if
        // it is off by a step.
        if (bdfOrder > 1 && getenv("CF_DUMP_BDF_HISTORY") && nSteps < 4)
        {
            const vectorField& Cdbg = mesh.C();

            scalar dOld = 0.0;
            forAll(VmOld.primitiveField(), cellI)
            {
                dOld = max
                (
                    dOld,
                    mag(VmOld.primitiveField()[cellI] - exactVm(Cdbg[cellI], t, dim))
                );
            }
            Pout<< "BDFHIST step " << nSteps << " t = " << t
                << " VmOld vs exact(t): " << dOld << endl;

            forAll(VmBdfOlder, j)
            {
                const scalar tj = t - scalar(j + 1)*dt;
                scalar dj = 0.0;
                forAll(VmBdfOlder[j], cellI)
                {
                    dj = max
                    (
                        dj,
                        mag(VmBdfOlder[j][cellI] - exactVm(Cdbg[cellI], tj, dim))
                    );
                }
                Pout<< "BDFHIST step " << nSteps << " level " << j
                    << " expects t = " << tj << " deviation " << dj << endl;
            }
        }
        const EigVec sourceN = sourceToEigVec(sourceVm);

        volScalarField VmGuess
        (
            IOobject
            (
                "VmGuess",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            VmOld
        );

        volScalarField u1Guess
        (
            IOobject
            (
                "u1Guess",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            u1Old
        );

        volScalarField u2Guess
        (
            IOobject
            (
                "u2Guess",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            u2Old
        );

        volScalarField u3Guess
        (
            IOobject
            (
                "u3Guess",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            u3Old
        );

        volScalarField IionGuess
        (
            IOobject
            (
                "IionGuess",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Iion
        );

        if (useJFNK && jfnkInitGuessOrder >= 1 && hasPrevStep && dtPrev > SMALL)
        {
            scalarField& Vg = VmGuess.primitiveFieldRef();
            const scalarField& VnInternal = VmOld.primitiveField();
            const scalar ratio = dt/dtPrev;

            forAll(Vg, cellI)
            {
                scalar v =
                    VnInternal[cellI]
                  + ratio*(VnInternal[cellI] - VmTwoStepsAgo[cellI]);

                v = min(max(v, jfnkInitGuessVmMin), jfnkInitGuessVmMax);
                Vg[cellI] = v;
            }
        }
        VmGuess.correctBoundaryConditions();

        VmTwoStepsAgo = VmOld.primitiveField();
        dtPrev = dt;
        hasPrevStep = true;

        volScalarField sourceDerivative
        (
            IOobject
            (
                "sourceDerivative",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("sourceDerivative", dimless/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );

        // The interval the state ODEs are integrated
        // over, and the time the candidate fields belong to, measured from t.
        //
        // For every one-step scheme this is the whole step and never changes.
        // ESDIRK is the exception: stage i lands at t + c_i dt, so its states
        // must be integrated over c_i dt and its boundary values evaluated
        // there. The ESDIRK block below is the ONLY place that assigns this,
        // and it restores dt before leaving the step - a stale value would
        // silently integrate the states over the wrong interval, which is
        // exactly the kind of error a convergence study reads as a bad order
        // rather than as a bug.
        scalar stageDt = dt;

        // The Vm nodes the state ODEs are driven by. Node 0 is always V^n;
        // the LAST node is always the candidate and is filled in by the
        // evaluator below, because it changes at every nonlinear iteration.
        // The two-node linear blend set here is what backwardEuler,
        // crankNicolson and BDF1 use, bit-identical to what this replaces.
        // BDF2-4 and the ESDIRK startup raise nNodes further down, unless
        // stateODEDriver is linear.
        VmDriverStages vmDriver;
        vmDriver.nNodes = 2;
        vmDriver.s[0] = 0.0;
        vmDriver.s[1] = 1.0;
        vmDriver.field[0] = &VmOld.primitiveField();

        // The same driver at the Gauss points, for the gaussPointODE path: the
        // linear blend from V^n by default, and under stateODEGaussPointHistory
        // the history polynomial, set below with its levels reconstructed at
        // the Gauss points (as in the physiological solver).
        VmDriverStages ipDriver;
        ipDriver.nNodes = 2;
        ipDriver.s[0] = 0.0;
        ipDriver.s[1] = 1.0;
        ipDriver.field[0] = &VmOldIntegrationPoints;
        const bool gaussHistory =
            stateODEGaussPointHistory && historyStateODEDriver
         && useHighOrder_Iion && dim > 1 && !reconstructStatesFromCellCentres;
        autoPtr<volScalarField> VmLevelScratchPtr;
        List<scalarField> VmBdfOlderIntegrationPoints;
        if (gaussHistory)
        {
            VmLevelScratchPtr.reset(new volScalarField("VmLevelScratch", VmOld));
        }

        auto evaluateNonlinearFields =
        [&]
        (
            const volScalarField& VmCandidate,
            volScalarField& u1Candidate,
            volScalarField& u2Candidate,
            volScalarField& u3Candidate,
            scalarField& VmCandidateIntegrationPoints,
            scalarField& u1CandidateIntegrationPoints,
            scalarField& u2CandidateIntegrationPoints,
            scalarField& u3CandidateIntegrationPoints,
            volScalarField& sourceCandidate,
            volScalarField& IionCandidate
        )
        {
            const auto evalStart = std::chrono::steady_clock::now();

            if (useHighOrder_Iion && dim > 1)
            {
                reconstructVmAtIionIntegrationPoints
                (
                    VmCandidate,
                    useHighOrder_Vm,
                    LREInterp_Vm,
                    LREInterp_Iion,
                    VmCandidateIntegrationPoints,
                    stateODEUseOpenMP,
                    stateODEOpenMPThreshold,
                    profileTimings,
                    fineTimings
                );

                if (jfnkClampODEInput)
                {
                    clampPhysical
                    (
                        VmCandidateIntegrationPoints,
                        jfnkInitGuessVmMin,
                        jfnkInitGuessVmMax
                    );
                }

                if (reconstructStatesFromCellCentres)
                {
                    // ODE once per cell (cell-centred), then reconstruct the
                    // states at the Iion Gauss points with LREInterp_states.
                    // The candidate is the LAST driver node and changes every
                    // nonlinear iteration, so it is refreshed here rather than
                    // captured once.
                    vmDriver.field[vmDriver.nNodes - 1] =
                        &VmCandidate.primitiveField();

                    updateStateFieldsODE
                    (
                        vmDriver,
                        u1Old,
                        u2Old,
                        u3Old,
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        stageDt,
                        stateODESolver,
                        stateODEInitialStep,
                        stateODEAbsTol,
                        stateODERelTol,
                        stateODEMaxSteps,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );

                    // Temporary diagnostic, env-gated: overwrite the integrated
                    // states with the exact manufactured ones. This decouples
                    // the PDE from the state ODE completely, so the temporal
                    // order measured with it ON is the order of the PDE step
                    // alone. If BDF3 recovers third order here and not
                    // otherwise, the cap lives in the state path.
                    if (getenv("CF_EXACT_STATES"))
                    {
                        const vectorField& Cxs = mesh.C();
                        scalarField& u1i = u1Candidate.primitiveFieldRef();
                        scalarField& u2i = u2Candidate.primitiveFieldRef();
                        scalarField& u3i = u3Candidate.primitiveFieldRef();

                        forAll(u1i, cellI)
                        {
                            u1i[cellI] = exactU1(Cxs[cellI], t + stageDt, dim);
                            u2i[cellI] = exactU2(Cxs[cellI], t + stageDt, dim);
                            u3i[cellI] = exactU3(Cxs[cellI], t + stageDt, dim);
                        }
                    }

                    updateStateBoundaryValues
                    (
                        u1Candidate, u2Candidate, u3Candidate, t + stageDt, dim
                    );
                    u1Candidate.correctBoundaryConditions();
                    u2Candidate.correctBoundaryConditions();
                    u3Candidate.correctBoundaryConditions();

                    reconstructStatesAtIionIntegrationPoints
                    (
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        LREInterp_statesPtr(),
                        LREInterp_Iion,
                        u1CandidateIntegrationPoints,
                        u2CandidateIntegrationPoints,
                        u3CandidateIntegrationPoints,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold,
                        profileTimings,
                        fineTimings
                    );
                }
                else
                {
                    // Integrate the ODE independently at every Gauss point, 
                    // then average the states back to the cell centres.
                    // Modified for cardiacFoam: stageDt, not dt. Inside an
                    // ESDIRK stage the states are integrated over c_i dt, as
                    // on the cell path; passing dt integrated them over the
                    // whole step at every stage.
                    updateIntegrationPointStatesODE
                    (
                        ipDriver,
                        VmCandidateIntegrationPoints,
                        u1OldIntegrationPoints,
                        u2OldIntegrationPoints,
                        u3OldIntegrationPoints,
                        u1CandidateIntegrationPoints,
                        u2CandidateIntegrationPoints,
                        u3CandidateIntegrationPoints,
                        stageDt,
                        stateODESolver,
                        stateODEInitialStep,
                        stateODEAbsTol,
                        stateODERelTol,
                        stateODEMaxSteps,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );

                    averageIntegrationPointFieldToCells
                    (
                        u1CandidateIntegrationPoints,
                        LREInterp_Iion,
                        u1Candidate
                    );
                    averageIntegrationPointFieldToCells
                    (
                        u2CandidateIntegrationPoints,
                        LREInterp_Iion,
                        u2Candidate
                    );
                    averageIntegrationPointFieldToCells
                    (
                        u3CandidateIntegrationPoints,
                        LREInterp_Iion,
                        u3Candidate
                    );
                }

                computeVmSourceFromIntegrationPoints
                (
                    VmCandidateIntegrationPoints,
                    u1CandidateIntegrationPoints,
                    u2CandidateIntegrationPoints,
                    u3CandidateIntegrationPoints,
                    beta,
                    chiVal,
                    CmVal,
                    sourceIntegrationPoints,
                    stateODEUseOpenMP,
                    stateODEOpenMPThreshold
                );
                averageIntegrationPointFieldToCells
                (
                    sourceIntegrationPoints,
                    LREInterp_Iion,
                    sourceCandidate
                );

                computeIionFromIntegrationPoints
                (
                    VmCandidateIntegrationPoints,
                    u1CandidateIntegrationPoints,
                    u2CandidateIntegrationPoints,
                    u3CandidateIntegrationPoints,
                    beta,
                    chiVal,
                    CmVal,
                    IionIntegrationPoints,
                    stateODEUseOpenMP,
                    stateODEOpenMPThreshold
                );
                averageIntegrationPointFieldToCells
                (
                    IionIntegrationPoints,
                    LREInterp_Iion,
                    IionCandidate
                );
            }
            else
            {
                if (jfnkClampODEInput)
                {
                    scalarField VmClampedInternal(VmCandidate.internalField());
                    clampPhysical
                    (
                        VmClampedInternal,
                        jfnkInitGuessVmMin,
                        jfnkInitGuessVmMax
                    );

                    volScalarField VmClamped
                    (
                        IOobject
                        (
                            "VmClamped",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        VmCandidate
                    );
                    VmClamped.primitiveFieldRef() = VmClampedInternal;
                    VmClamped.correctBoundaryConditions();

                    // The candidate is the LAST driver node and changes every
                    // nonlinear iteration, so it is refreshed here rather than
                    // captured once.
                    vmDriver.field[vmDriver.nNodes - 1] =
                        &VmClamped.primitiveField();

                    updateStateFieldsODE
                    (
                        vmDriver,
                        u1Old,
                        u2Old,
                        u3Old,
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        stageDt,
                        stateODESolver,
                        stateODEInitialStep,
                        stateODEAbsTol,
                        stateODERelTol,
                        stateODEMaxSteps,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );

                    computeCellCentredVmSource
                    (
                        VmClamped,
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        beta,
                        chiVal,
                        CmVal,
                        sourceCandidate,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );

                    computeCellCentredIion
                    (
                        VmClamped,
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        beta,
                        chiVal,
                        CmVal,
                        IionCandidate,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );
                }
                else
                {
                    // The candidate is the LAST driver node and changes every
                    // nonlinear iteration, so it is refreshed here rather than
                    // captured once.
                    vmDriver.field[vmDriver.nNodes - 1] =
                        &VmCandidate.primitiveField();

                    updateStateFieldsODE
                    (
                        vmDriver,
                        u1Old,
                        u2Old,
                        u3Old,
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        stageDt,
                        stateODESolver,
                        stateODEInitialStep,
                        stateODEAbsTol,
                        stateODERelTol,
                        stateODEMaxSteps,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );

                    computeCellCentredVmSource
                    (
                        VmCandidate,
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        beta,
                        chiVal,
                        CmVal,
                        sourceCandidate,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );

                    computeCellCentredIion
                    (
                        VmCandidate,
                        u1Candidate,
                        u2Candidate,
                        u3Candidate,
                        beta,
                        chiVal,
                        CmVal,
                        IionCandidate,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );
                }
            }

            updateStateBoundaryValues(u1Candidate, u2Candidate, u3Candidate, t + stageDt, dim);
            u1Candidate.correctBoundaryConditions();
            u2Candidate.correctBoundaryConditions();
            u3Candidate.correctBoundaryConditions();

            if (profileTimings)
            {
                const auto evalEnd = std::chrono::steady_clock::now();
                fineTimings.nonlinearEvalWallTime +=
                    std::chrono::duration<scalar>(evalEnd - evalStart).count();
                ++fineTimings.nonlinearEvalCalls;
            }
        };

        auto computeSourceDerivative = [&]()
        {
            volScalarField VmPlus
            (
                IOobject
                (
                    "VmPlus",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                VmGuess
            );

            VmPlus.primitiveFieldRef() += diagonalIionEpsilon;
            VmPlus.correctBoundaryConditions();

            volScalarField u1Plus
            (
                IOobject
                (
                    "u1Plus",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                u1Old
            );
            volScalarField u2Plus
            (
                IOobject
                (
                    "u2Plus",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                u2Old
            );
            volScalarField u3Plus
            (
                IOobject
                (
                    "u3Plus",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                u3Old
            );
            volScalarField sourcePlus
            (
                IOobject
                (
                    "sourcePlus",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                sourceVm
            );
            volScalarField IionPlus
            (
                IOobject
                (
                    "IionPlus",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                IionGuess
            );

            scalarField VmPlusIntegrationPoints
            (
                allocatedIionIntegrationPoints,
                0.0
            );
            scalarField u1PlusIntegrationPoints(u1OldIntegrationPoints);
            scalarField u2PlusIntegrationPoints(u2OldIntegrationPoints);
            scalarField u3PlusIntegrationPoints(u3OldIntegrationPoints);

            evaluateNonlinearFields
            (
                VmPlus,
                u1Plus,
                u2Plus,
                u3Plus,
                VmPlusIntegrationPoints,
                u1PlusIntegrationPoints,
                u2PlusIntegrationPoints,
                u3PlusIntegrationPoints,
                sourcePlus,
                IionPlus
            );

            scalarField& deriv = sourceDerivative.primitiveFieldRef();
            forAll(deriv, cellI)
            {
                deriv[cellI] =
                    (sourcePlus[cellI] - sourceVm[cellI])
                   /max(diagonalIionEpsilon, SMALL);
            }
            sourceDerivative.correctBoundaryConditions();
        };

        auto rhsFromSource = [&](const volScalarField& sourceCandidate)
        {
            const EigVec sourceNp1 = sourceToEigVec(sourceCandidate);

            EigVec rhs =
                applyB(Vn)
              + theta*sourceNp1
              + (1.0 - theta)*sourceN;

            return rhs;
        };

        if (useJFNK)
        {
            evaluateNonlinearFields
            (
                VmGuess,
                u1Guess,
                u2Guess,
                u3Guess,
                VmGuessIntegrationPoints,
                u1GuessIntegrationPoints,
                u2GuessIntegrationPoints,
                u3GuessIntegrationPoints,
                sourceVm,
                IionGuess
            );
        }

        bool esdirkHandledStep = false;

        bool nonlinearConverged = false;
        label nonlinearIters = 0;
        scalar finalCoupledResidual = GREAT;
        scalar finalVmResidual = GREAT;
        scalar finalU1Residual = GREAT;
        scalar finalU2Residual = GREAT;
        scalar finalU3Residual = GREAT;
        scalar finalMaxStateResidual = GREAT;
        scalar finalIionResidual = GREAT;
        label finalLinearIterations = 0;
        scalar finalLinearError = GREAT;
        label finalLineSearchIterations = -1;

        auto writeAndPrintResiduals =
        [&]
        (
            const label corr,
            const label linearIterations,
            const scalar linearError,
            const scalar newtonResidual,
            const scalar VmResidual,
            const scalar u1Residual,
            const scalar u2Residual,
            const scalar u3Residual,
            const scalar IionResidual,
            const bool converged,
            const label lineSearchIters = -1
        )
        {
            Info<< "Time = " << t << "," << nl
                << "       Nonlinear method:          " << nonlinearMethod << nl
                << "       Implicit FDA solver:       " << implicitLinearSolver
                << "; iterations = " << linearIterations
                << ", estimated error = " << linearError << nl
                << "       NonLinSolver:              iter = " << corr
                << ", Vm residual = " << VmResidual
                << ", u1 residual = " << u1Residual
                << ", u2 residual = " << u2Residual
                << ", u3 residual = " << u3Residual
                << ", Iion residual = " << IionResidual
                << ", coupled residual = " << newtonResidual;
            if (lineSearchIters >= 0)
            {
                Info<< ", lineSearchIters = " << lineSearchIters;
            }
            Info<< ", converged = " << converged << nl;

            if (writeNonlinearResiduals && nonlinearResidualFilePtr.valid())
            {
                nonlinearResidualFilePtr()
                    << t << ' ' << (nSteps + 1) << ' ' << nonlinearMethod << ' '
                    << stateODESolver << ' '
                    << corr << ' ' << linearIterations << ' ' << linearError << ' '
                    << newtonResidual << ' '
                    << VmResidual << ' ' << u1Residual << ' '
                    << u2Residual << ' ' << u3Residual << ' '
                    << IionResidual << ' '
                    << lineSearchIters << ' ' << converged << nl;
            }

            finalCoupledResidual = newtonResidual;
            finalVmResidual = VmResidual;
            finalU1Residual = u1Residual;
            finalU2Residual = u2Residual;
            finalU3Residual = u3Residual;
            finalMaxStateResidual =
                max(u1Residual, max(u2Residual, u3Residual));
            finalIionResidual = IionResidual;
            finalLinearIterations = linearIterations;
            finalLinearError = linearError;
            finalLineSearchIterations = lineSearchIters;
        };

        // The ESDIRK startup step. This runs ONLY for
        // the first bdfOrder-1 steps of a BDF run, to build the history; it is
        // not reachable as a scheme in its own right (see the note on
        // needEsdirkOperator above for why not).
        //
        // Stage 1 is explicit, Y_1 = V^n. For i >= 2 the stage solves
        //
        //   A_g Y_i = (1/g) [ B V^n + sum_{j<i} a_ij F_j ] + bc(t_i) + s(Y_i)
        //
        // with A_g = M/(g dt) - lapScale K, B = M/dt, t_i = t + c_i dt, and
        //
        //   F_j = lapScale K Y_j + bc(t_j) + s(Y_j).
        //
        // Note bc enters with coefficient 1 while the B and F terms carry 1/g:
        // that falls out of dividing the stage equation by g dt. Getting it
        // wrong would be invisible on this MMS suite, where every case is
        // zeroGradient and bc is identically zero - hence this note.
        //
        // The shape is the same one-step equation the Picard branch below
        // solves, which is why one matrix, one cached factorisation and one
        // fixed-point iteration serve every stage.
        //
        // Stiffly accurate (c_s = 1, b = last row): V^{n+1} = Y_s, the last
        // stage IS the answer and there is no separate combination step.
        //
        // Only Picard is wired up. JFNK and diagonalIion would each need their
        // per-stage Jacobian machinery threaded through, which is separate
        // work; an unknown scheme name is a FatalError upstream, so nothing
        // reaches here by accident.
        // ESDIRK runs only while it is supplying the first bdfOrder-1 starting
        // values. After that the history is full and BDF takes over for the
        // rest of the run - which, for anything but the opening steps, is the
        // whole run.
        // Modified for cardiacFoam: counted in VALID levels rather than in
        // steps, so that a rolled-back step can restart the history (see the
        // rotation at the end of the step), as the physiological solver does.
        const bool esdirkThisStep =
            (bdfOrder > 1 && esdirkStartup && bdfValidLevels < bdfOrder - 1);

        // A BDF step drives its state ODEs with the
        // polynomial through its OWN history levels, for the same reason ESDIRK
        // uses its stage values, and at the same price - the levels are already
        // stored, so the only cost is the interpolation itself.
        //
        // Nodes at s = -(k-1) ... -1, 0, 1 in units of dt measured from t^n,
        // with the candidate V^{n+1} last. The ODEs are integrated over
        // s in [0, 1], which lies inside the node hull, so this interpolates
        // rather than extrapolates.
        //
        // Why it matters, measured (suite test M21): with the two-node blend
        // BDF4 fits order 2.04, with this interpolant 3.89. The blend is
        // O(dt^2) pointwise, O(dt^3) once integrated over the step, so it caps
        // the global order at 2 whatever the scheme. BDF can reach fourth order where
        // ESDIRK4 cannot, because a multistep method has stage order equal to
        // its order and suffers no order reduction.
        //
        // Skipped while ESDIRK is supplying the starting values, since the
        // history is not full yet and that block sets its own driver. Skipped
        // too under stateODEDriver linear, which keeps the two-node default.
        if (bdfOrder > 1 && !esdirkThisStep && historyStateODEDriver)
        {
            vmDriver.nNodes = bdfOrder + 1;

            // Oldest first, so the candidate lands on the last slot, which is
            // the one the evaluator refreshes.
            for (label j = 0; j < bdfOrder - 1; ++j)
            {
                // VmBdfOlder[m] holds V^{n-1-m}, so walk it backwards.
                const label m = bdfOrder - 2 - j;
                vmDriver.s[j] = -scalar(m + 1);
                vmDriver.field[j] = &VmBdfOlder[m];
            }

            vmDriver.s[bdfOrder - 1] = 0.0;
            vmDriver.field[bdfOrder - 1] = &VmOld.primitiveField();

            vmDriver.s[bdfOrder] = 1.0;
            vmDriver.field[bdfOrder] = nullptr;

            if (gaussHistory)
            {
                VmBdfOlderIntegrationPoints.setSize(bdfOrder - 1);
                forAll(VmBdfOlderIntegrationPoints, m)
                {
                    VmBdfOlderIntegrationPoints[m].setSize
                    (
                        VmOldIntegrationPoints.size()
                    );
                    reconstructVmLevelAtIionIntegrationPoints
                    (
                        VmBdfOlder[m],
                        VmLevelScratchPtr(),
                        useHighOrder_Vm,
                        LREInterp_Vm,
                        LREInterp_Iion,
                        VmBdfOlderIntegrationPoints[m],
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold,
                        profileTimings,
                        fineTimings
                    );
                }

                ipDriver.nNodes = bdfOrder + 1;
                for (label j = 0; j < bdfOrder - 1; ++j)
                {
                    const label m = bdfOrder - 2 - j;
                    ipDriver.s[j] = -scalar(m + 1);
                    ipDriver.field[j] = &VmBdfOlderIntegrationPoints[m];
                }
                ipDriver.s[bdfOrder - 1] = 0.0;
                ipDriver.field[bdfOrder - 1] = &VmOldIntegrationPoints;
                ipDriver.s[bdfOrder] = 1.0;
                ipDriver.field[bdfOrder] = nullptr;
            }
        }

        if (esdirkThisStep)
        {
            const label nS = esdirkTab.nStages;
            const scalar gam = esdirkTab.gamma;

            List<EigVec> stageF(nS);
            const EigVec BVn = applyB(fieldToEigVec(VmOld));

            // The stage VALUES, kept so that the state ODEs can be driven by a
            // polynomial through them instead of by a straight line between the
            // two ends of the stage interval. Y_1 = V^n, so index 0 is a copy
            // of VmOld and the rest are filled as the stages are solved.
            //
            // This is what lifts the second-order cap the linear blend imposed:
            // the values are already computed, so the higher-order interpolant
            // costs storage and nothing else.
            List<scalarField> stageVm(nS);
            forAll(stageVm, j)
            {
                // Sized up front: the last node is always overwritten with the
                // candidate before use, but a zero-length field behind a live
                // pointer reads out of bounds instead of failing, so it is not
                // left to that guarantee.
                stageVm[j].setSize(mesh.nCells(), 0.0);
            }
            stageVm[0] = VmOld.primitiveField();

            // The stage values at the Gauss points, for the gaussPointODE path
            // under stateODEGaussPointHistory; each is reconstructed once, when
            // its stage is accepted.
            List<scalarField> stageVmIntegrationPoints(gaussHistory ? nS : 0);
            if (gaussHistory)
            {
                stageVmIntegrationPoints[0] = VmOldIntegrationPoints;
            }

            // Stage 1, explicit: Y_1 = V^n with the states at t^n, which is
            // what stageDt = 0 delivers (advanceStateODE returns u unchanged
            // over a zero interval).
            stageDt = 0.0;
            VmGuess.primitiveFieldRef() = VmOld.primitiveField();
            VmGuess.correctBoundaryConditions();

            evaluateNonlinearFields
            (
                VmGuess,
                u1Guess,
                u2Guess,
                u3Guess,
                VmGuessIntegrationPoints,
                u1GuessIntegrationPoints,
                u2GuessIntegrationPoints,
                u3GuessIntegrationPoints,
                sourceVm,
                IionGuess
            );

            stageF[0] =
                lapScale*applyK(fieldToEigVec(VmGuess))
              + sourceToEigVec(sourceVm);

            nonlinearConverged = true;

            for (label i = 1; i < nS; ++i)
            {
                stageDt = esdirkTab.c[i]*dt;

                // Drive the state ODEs with the polynomial through the stage
                // values already known, in time normalised by THIS stage's
                // interval: s_j = c_j / c_i, so s_0 = 0 and s_i = 1 as the
                // interpolant requires. The last node is the stage being
                // solved, refreshed by the evaluator at every iteration.
                //
                // Note the abscissae need not be increasing - Lagrange does not
                // care, and a six-stage ESDIRK would have c = [0, 1/2, 0.332,
                // 0.62, 0.85, 1], which is not monotone. What matters is that
                // [0, 1] lies inside their hull, which it does because c_1 = 0
                // and c_i are both nodes.
                if (historyStateODEDriver)
                {
                    vmDriver.nNodes = i + 1;
                    for (label j = 0; j <= i; ++j)
                    {
                        vmDriver.s[j] = esdirkTab.c[j]/esdirkTab.c[i];
                        vmDriver.field[j] = &stageVm[j];
                    }

                    if (gaussHistory)
                    {
                        ipDriver.nNodes = i + 1;
                        for (label j = 0; j < i; ++j)
                        {
                            ipDriver.s[j] = esdirkTab.c[j]/esdirkTab.c[i];
                            ipDriver.field[j] = &stageVmIntegrationPoints[j];
                        }
                        ipDriver.s[i] = 1.0;
                        ipDriver.field[i] = nullptr;
                    }
                }
                else
                {
                    // The ramp from stage 0 (= V^n) to the stage being solved.
                    vmDriver.nNodes = 2;
                    vmDriver.s[0] = 0.0;
                    vmDriver.s[1] = 1.0;
                    vmDriver.field[0] = &stageVm[0];
                    vmDriver.field[1] = &stageVm[i];
                }

                EigVec rhsBase = BVn;
                for (label j = 0; j < i; ++j)
                {
                    rhsBase += esdirkTab.a[i][j]*stageF[j];
                }
                rhsBase /= gam;

                // Modified for cardiacFoam: the stage is solved with the same
                // relaxation and the same stopping test as every other
                // nonlinear iteration (nonlinearIterationExit), as in the
                // physiological solver. It used to stop on the Vm increment
                // alone, with no relaxation, so an ESDIRK step could be less
                // converged than the BDF steps it seeds.
                bool stageConverged = false;

                // Fields consistent with the current stage guess and interval.
                evaluateNonlinearFields
                (
                    VmGuess,
                    u1Guess,
                    u2Guess,
                    u3Guess,
                    VmGuessIntegrationPoints,
                    u1GuessIntegrationPoints,
                    u2GuessIntegrationPoints,
                    u3GuessIntegrationPoints,
                    sourceVm,
                    IionGuess
                );

                scalar coupledResidualPrevIter = GREAT;

                for
                (
                    label corr = 0;
                    corr < max(implicitNonlinearIterations, label(1));
                    ++corr
                )
                {
                    const scalarField VmPrevious(VmGuess.primitiveField());
                    const scalarField u1Previous(u1Guess.primitiveField());
                    const scalarField u2Previous(u2Guess.primitiveField());
                    const scalarField u3Previous(u3Guess.primitiveField());
                    const scalarField u1IPPrevious(u1GuessIntegrationPoints);
                    const scalarField u2IPPrevious(u2GuessIntegrationPoints);
                    const scalarField u3IPPrevious(u3GuessIntegrationPoints);
                    const scalarField IionPrevious(IionGuess.primitiveField());

                    const EigVec rhs = rhsBase + sourceToEigVec(sourceVm);

                    label stageLinIters = 0;
                    scalar stageLinError = GREAT;

                    MatrixUpdate stagePolicy = MatrixUpdate::Reuse;
                    if (esdirkSolverNeedsRebuild)
                    {
                        stagePolicy = MatrixUpdate::Rebuild;
                        esdirkSolverNeedsRebuild = false;
                    }

                    const EigVec Ysol = solveSparseSystem
                    (
                        AEsdirk,
                        rhs,
                        linearSolverBackend,
                        implicitLinearSolver,
                        petscLinearKspType,
                        petscLinearPcType,
                        petscLinearRestart,
                        petscLinearOptionsPrefix,
                        petscUseOptions,
                        implicitTolerance,
                        implicitMaxIterations,
                        stageLinIters,
                        stageLinError,
                        esdirkPetscSolver,
                        esdirkEigenSolvers,
                        stagePolicy,
                        petscReusePreconditioner
                    );

                    const EigVec Yprev = fieldToEigVec(VmGuess);
                    eigVecToField
                    (
                        Yprev + nonlinearRelaxation*(Ysol - Yprev), VmGuess
                    );
                    VmGuess.correctBoundaryConditions();

                    evaluateNonlinearFields
                    (
                        VmGuess,
                        u1Guess,
                        u2Guess,
                        u3Guess,
                        VmGuessIntegrationPoints,
                        u1GuessIntegrationPoints,
                        u2GuessIntegrationPoints,
                        u3GuessIntegrationPoints,
                        sourceVm,
                        IionGuess
                    );

                    const EigVec rhsNew = rhsBase + sourceToEigVec(sourceVm);
                    const scalar coupledResidual = relativeL2Norm
                    (
                        applyEsdirk(fieldToEigVec(VmGuess)) - rhsNew, rhsNew
                    );
                    const scalar VmResidual =
                        relativeL2Difference(VmGuess.primitiveField(), VmPrevious);

                    scalar u1Residual = GREAT;
                    scalar u2Residual = GREAT;
                    scalar u3Residual = GREAT;
                    if (useHighOrder_Iion && dim > 1 && !reconstructStatesFromCellCentres)
                    {
                        u1Residual =
                            relativeL2Difference(u1GuessIntegrationPoints, u1IPPrevious);
                        u2Residual =
                            relativeL2Difference(u2GuessIntegrationPoints, u2IPPrevious);
                        u3Residual =
                            relativeL2Difference(u3GuessIntegrationPoints, u3IPPrevious);
                    }
                    else
                    {
                        u1Residual =
                            relativeL2Difference(u1Guess.primitiveField(), u1Previous);
                        u2Residual =
                            relativeL2Difference(u2Guess.primitiveField(), u2Previous);
                        u3Residual =
                            relativeL2Difference(u3Guess.primitiveField(), u3Previous);
                    }
                    const scalar IionResidual =
                        relativeL2Difference(IionGuess.primitiveField(), IionPrevious);

                    const bool converged = nonlinearIterationStops
                    (
                        nonlinearIterationExit
                        (
                            corr,
                            implicitMinNonlinearIterations,
                            coupledResidual,
                            coupledResidualPrevIter,
                            nonlinearTolerance,
                            VmResidual,
                            nonlinearVmTolerance,
                            IionResidual,
                            nonlinearIionTolerance,
                            nonlinearRequireStatesConvergence,
                            max(u1Residual, max(u2Residual, u3Residual)),
                            nonlinearStatesTolerance
                        ),
                        coupledResidual,
                        nonlinearTolerance,
                        nResidualFloorExits
                    );
                    coupledResidualPrevIter = coupledResidual;

                    writeAndPrintResiduals
                    (
                        corr + 1,
                        stageLinIters,
                        stageLinError,
                        coupledResidual,
                        VmResidual,
                        u1Residual,
                        u2Residual,
                        u3Residual,
                        IionResidual,
                        converged
                    );

                    ++nonlinearIters;
                    if (converged)
                    {
                        stageConverged = true;
                        break;
                    }
                }

                if (!stageConverged)
                {
                    nonlinearConverged = false;
                }

                // The fields are already evaluated at the accepted stage value
                // (the last thing the loop did), and F_i feeds every later stage.
                stageF[i] =
                    lapScale*applyK(fieldToEigVec(VmGuess))
                  + sourceToEigVec(sourceVm);

                // Keep the accepted stage value: later stages interpolate their
                // ODE driver through it.
                stageVm[i] = VmGuess.primitiveField();

                if (gaussHistory && i + 1 < nS)
                {
                    stageVmIntegrationPoints[i].setSize
                    (
                        VmOldIntegrationPoints.size()
                    );
                    reconstructVmLevelAtIionIntegrationPoints
                    (
                        stageVm[i],
                        VmLevelScratchPtr(),
                        useHighOrder_Vm,
                        LREInterp_Vm,
                        LREInterp_Iion,
                        stageVmIntegrationPoints[i],
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold,
                        profileTimings,
                        fineTimings
                    );
                }
            }

            // Stiffly accurate: the last stage sits at t + dt and carries both
            // V^{n+1} and the states, so VmGuess/u*Guess already hold what the
            // shared bookkeeping below expects.
            stageDt = dt;

            // Restore the two-node default. stageVm is local to this block, so
            // leaving the driver pointing into it would dangle the moment the
            // block ends - and a dangling read is a wrong answer, not a crash.
            vmDriver.nNodes = 2;
            vmDriver.s[0] = 0.0;
            vmDriver.s[1] = 1.0;
            vmDriver.field[0] = &VmOld.primitiveField();
            vmDriver.field[1] = nullptr;
            ipDriver.nNodes = 2;
            ipDriver.s[0] = 0.0;
            ipDriver.s[1] = 1.0;
            ipDriver.field[0] = &VmOldIntegrationPoints;
            ipDriver.field[1] = nullptr;

            finalCoupledResidual = finalVmResidual;
            esdirkHandledStep = true;
        }

        if (useJFNK && !esdirkHandledStep)
        {
            EigVec x = fieldToEigVec(VmGuess);

            volScalarField VmTmp
            (
                IOobject
                (
                    "VmTmp",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                VmGuess
            );
            volScalarField u1Tmp
            (
                IOobject
                (
                    "u1Tmp",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                u1Guess
            );
            volScalarField u2Tmp
            (
                IOobject
                (
                    "u2Tmp",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                u2Guess
            );
            volScalarField u3Tmp
            (
                IOobject
                (
                    "u3Tmp",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                u3Guess
            );
            volScalarField sourceTmp
            (
                IOobject
                (
                    "sourceTmp",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                sourceVm
            );
            volScalarField IionTmp
            (
                IOobject
                (
                    "IionTmp",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                IionGuess
            );

            scalarField VmTmpIntegrationPoints
            (
                allocatedIionIntegrationPoints,
                0.0
            );
            scalarField u1TmpIntegrationPoints(u1OldIntegrationPoints);
            scalarField u2TmpIntegrationPoints(u2OldIntegrationPoints);
            scalarField u3TmpIntegrationPoints(u3OldIntegrationPoints);

            auto residualFor =
            [&](const EigVec& xCandidate, const bool updateGuess) -> EigVec
            {
                if (updateGuess)
                {
                    eigVecToField(xCandidate, VmGuess);
                    VmGuess.correctBoundaryConditions();
                    evaluateNonlinearFields
                    (
                        VmGuess,
                        u1Guess,
                        u2Guess,
                        u3Guess,
                        VmGuessIntegrationPoints,
                        u1GuessIntegrationPoints,
                        u2GuessIntegrationPoints,
                        u3GuessIntegrationPoints,
                        sourceVm,
                        IionGuess
                    );

                    const EigVec xApplied = fieldToEigVec(VmGuess);
                    return applyA(xApplied) - rhsFromSource(sourceVm);
                }

                VmTmp = VmGuess;
                u1Tmp = u1Guess;
                u2Tmp = u2Guess;
                u3Tmp = u3Guess;
                sourceTmp = sourceVm;
                IionTmp = IionGuess;

                VmTmpIntegrationPoints = 0.0;
                u1TmpIntegrationPoints = u1OldIntegrationPoints;
                u2TmpIntegrationPoints = u2OldIntegrationPoints;
                u3TmpIntegrationPoints = u3OldIntegrationPoints;

                eigVecToField(xCandidate, VmTmp);
                VmTmp.correctBoundaryConditions();

                evaluateNonlinearFields
                (
                    VmTmp,
                    u1Tmp,
                    u2Tmp,
                    u3Tmp,
                    VmTmpIntegrationPoints,
                    u1TmpIntegrationPoints,
                    u2TmpIntegrationPoints,
                    u3TmpIntegrationPoints,
                    sourceTmp,
                    IionTmp
                );

                const EigVec xApplied = fieldToEigVec(VmTmp);
                return applyA(xApplied) - rhsFromSource(sourceTmp);
            };


            jfnkPc.beginStep();

            // The increments measured at the END of the previous corr, which are
            // what the check at the top of the next one is entitled to look at.
            //
            // That check used to recompute them against a copy taken two lines
            // further up, with only residualFor(x, true) in between - and that
            // call re-evaluates the same x, so the "increments" were the
            // difference between two evaluations of one point, i.e. identically
            // zero. Four of the five tolerance terms therefore compared zero
            // against a positive number and could never fail: the test reduced
            // to `newtonResidual <= nonlinearTolerance` alone, and JFNK returned
            // the jfnkInitGuessOrder extrapolation without taking a single
            // Newton step whenever that extrapolation happened to land under the
            // coupled tolerance. Measured on M08 in 3-D: 24 of 50 timesteps
            // exited that way, 6.3e-08 from the solution the other two
            // nonlinear methods agree on, and no tolerance the user could set
            // made any difference - only nonlinearTolerance did.
            //
            // GREAT on entry means the first pass can never exit there, so at
            // least one Newton step is always taken. Verified by construction
            // against minNonlinearIterations = 2, which forces the same thing
            // from the outside and closes the gap to 2.0e-11.
            scalar VmIncrPrevIter = GREAT;
            scalar u1IncrPrevIter = GREAT;
            scalar u2IncrPrevIter = GREAT;
            scalar u3IncrPrevIter = GREAT;
            scalar IionIncrPrevIter = GREAT;

            scalar coupledResidualPrevIter = GREAT;

            for
            (
                label corr = 0;
                corr < max(implicitNonlinearIterations, label(1));
                ++corr
            )
            {
                const scalarField VmPrevious(VmGuess.primitiveField());
                const scalarField u1Previous(u1Guess.primitiveField());
                const scalarField u2Previous(u2Guess.primitiveField());
                const scalarField u3Previous(u3Guess.primitiveField());
                const scalarField u1IPPrevious(u1GuessIntegrationPoints);
                const scalarField u2IPPrevious(u2GuessIntegrationPoints);
                const scalarField u3IPPrevious(u3GuessIntegrationPoints);
                const scalarField IionPrevious(IionGuess.primitiveField());

                const EigVec R = residualFor(x, true);
                const EigVec rhsCurrent = rhsFromSource(sourceVm);
                const scalar newtonResidual = relativeL2Norm(R, rhsCurrent);

                // The increments that decide whether this iteration can be
                // skipped are the ones the PREVIOUS one measured; recomputing
                // them here would measure zero, see the note above the loop.
                const scalar VmResidualInitial = VmIncrPrevIter;
                const scalar u1ResidualInitial = u1IncrPrevIter;
                const scalar u2ResidualInitial = u2IncrPrevIter;
                const scalar u3ResidualInitial = u3IncrPrevIter;
                const scalar IionResidualInitial = IionIncrPrevIter;

                const bool residualConverged =
                    corr + 1 >= implicitMinNonlinearIterations
                 && newtonResidual <= nonlinearTolerance
                 && VmResidualInitial <= nonlinearVmTolerance
                 && IionResidualInitial <= nonlinearIionTolerance
                 && (
                        !nonlinearRequireStatesConvergence
                     || (
                            u1ResidualInitial <= nonlinearStatesTolerance
                         && u2ResidualInitial <= nonlinearStatesTolerance
                         && u3ResidualInitial <= nonlinearStatesTolerance
                        )
                    );

                if (residualConverged)
                {
                    writeAndPrintResiduals
                    (
                        corr + 1,
                        0,
                        0.0,
                        newtonResidual,
                        VmResidualInitial,
                        u1ResidualInitial,
                        u2ResidualInitial,
                        u3ResidualInitial,
                        IionResidualInitial,
                        true
                    );
                    nonlinearConverged = true;
                    nonlinearIters = corr + 1;
                    break;
                }

                // Globally reduced. An unreduced
                // norm here makes the finite-difference step, and therefore the
                // Jacobian approximation, depend on the rank count.
                const scalar epsScale = jfnkEpsilon*(1.0 + gNorm(x));

                auto matVec = [&](const EigVec& v) -> EigVec
                {
                    const scalar eps = epsScale/max(gNorm(v), SMALL);
                    return (residualFor(x + eps*v, false) - R)/eps;
                };

                label gmresIterations = 0;
                scalar gmresError = GREAT;
                EigVec delta(R.size());

                const auto gmresStart = std::chrono::steady_clock::now();

                if (jfnkPc.needsUpdate(corr))
                {
                    if (jfnkPc.usesSourceDerivative())
                    {
                        computeSourceDerivative();
                    }
                    jfnkPc.update
                    (
                        AImplicit,
                        jfnkPc.usesSourceDerivative()
                      ? &sourceDerivative.primitiveField()
                      : nullptr,
                        theta,
                        corr
                    );
                }

                delta = solveJfnkLinearSystem
                (
                    matVec,
                    R,
                    jfnkPc,
                    jfnkPetscShellSolver,
                    jfnkKrylov,
                    gmresIterations,
                    gmresError
                );

                if (profileTimings)
                {
                    fineTimings.preconditionerSetupWallTime = jfnkPc.setupWallTime();
                    fineTimings.preconditionerSetups = jfnkPc.setups();
                    fineTimings.preconditionerApplyWallTime = jfnkPc.applyWallTime();
                    fineTimings.preconditionerApplications = jfnkPc.applications();
                }

                if (profileTimings)
                {
                    const auto gmresEnd = std::chrono::steady_clock::now();
                    fineTimings.gmresWallTime +=
                        std::chrono::duration<scalar>
                        (
                            gmresEnd - gmresStart
                        ).count();
                    ++fineTimings.gmresCalls;
                }

                label lineSearchIters = 0;
                EigVec Rnew(R.size());
                scalar newtonResidualNew = GREAT;

                if (jfnkLineSearch)
                {
                    const EigVec xBackup = x;
                    const scalar Rold2 = gSquaredNorm(R);
                    scalar alpha = nonlinearRelaxation;

                    while (true)
                    {
                        x = xBackup + alpha*delta;
                        Rnew = residualFor(x, true);
                        const scalar Rnew2 = gSquaredNorm(Rnew);

                        if (Rnew2 <= (1.0 - 2.0*jfnkArmijoC*alpha)*Rold2)
                        {
                            break;
                        }

                        ++lineSearchIters;
                        if (lineSearchIters >= jfnkLineSearchMaxIter)
                        {
                            break;
                        }

                        alpha *= 0.5;
                        if (alpha < jfnkLineSearchAlphaMin)
                        {
                            alpha = jfnkLineSearchAlphaMin;
                            x = xBackup + alpha*delta;
                            Rnew = residualFor(x, true);
                            break;
                        }
                    }
                }
                else
                {
                    x += nonlinearRelaxation*delta;
                    Rnew = residualFor(x, true);
                }

                const EigVec rhsNew = rhsFromSource(sourceVm);
                newtonResidualNew = relativeL2Norm(Rnew, rhsNew);

                const scalar VmResidual =
                    relativeL2Difference(VmGuess.primitiveField(), VmPrevious);
                scalar u1Residual = GREAT;
                scalar u2Residual = GREAT;
                scalar u3Residual = GREAT;
                if (useHighOrder_Iion && dim > 1 && !reconstructStatesFromCellCentres)
                {
                    u1Residual =
                        relativeL2Difference(u1GuessIntegrationPoints, u1IPPrevious);
                    u2Residual =
                        relativeL2Difference(u2GuessIntegrationPoints, u2IPPrevious);
                    u3Residual =
                        relativeL2Difference(u3GuessIntegrationPoints, u3IPPrevious);
                }
                else
                {
                    // cellCentredReconstruct (and the low-order path) carry the
                    // authoritative states at cell centres -> measure there.
                    u1Residual =
                        relativeL2Difference(u1Guess.primitiveField(), u1Previous);
                    u2Residual =
                        relativeL2Difference(u2Guess.primitiveField(), u2Previous);
                    u3Residual =
                        relativeL2Difference(u3Guess.primitiveField(), u3Previous);
                }
                const scalar IionResidual =
                    relativeL2Difference(IionGuess.primitiveField(), IionPrevious);

                // Shared stopping test (monodomainBase.H); this used to be an
                // inline copy of it, without the residual-floor exit.
                const bool converged = nonlinearIterationStops
                (
                    nonlinearIterationExit
                    (
                        corr,
                        implicitMinNonlinearIterations,
                        newtonResidualNew,
                        coupledResidualPrevIter,
                        nonlinearTolerance,
                        VmResidual,
                        nonlinearVmTolerance,
                        IionResidual,
                        nonlinearIionTolerance,
                        nonlinearRequireStatesConvergence,
                        max(u1Residual, max(u2Residual, u3Residual)),
                        nonlinearStatesTolerance
                    ),
                    newtonResidualNew,
                    nonlinearTolerance,
                    nResidualFloorExits
                );
                coupledResidualPrevIter = newtonResidualNew;

                writeAndPrintResiduals
                (
                    corr + 1,
                    gmresIterations,
                    gmresError,
                    newtonResidualNew,
                    VmResidual,
                    u1Residual,
                    u2Residual,
                    u3Residual,
                    IionResidual,
                    converged,
                    lineSearchIters
                );

                // What the top-of-loop check of the next iteration may look at.
                VmIncrPrevIter = VmResidual;
                u1IncrPrevIter = u1Residual;
                u2IncrPrevIter = u2Residual;
                u3IncrPrevIter = u3Residual;
                IionIncrPrevIter = IionResidual;

                nonlinearIters = corr + 1;
                if (converged)
                {
                    nonlinearConverged = true;
                    break;
                }
            }
        }
        else if ((usePicard || useDiagonalIion) && !esdirkHandledStep)
        {
            // Persistent ACurrent buffer reused across the corr iterations
            // (and across timesteps). Only allocated when actually needed —
            // pure Picard with a cached PETSc solver never touches it. For
            // useDiagonalIion the value buffer is refreshed each iter from
            // AImplicit and then the diagonal source derivative correction
            // is applied in-place; sparsity is invariant.
            SpMat ACurrent;
            bool ACurrentInitialised = false;
            const bool needsACurrent = useDiagonalIion;

            scalar coupledResidualPrevIter = GREAT;

            for
            (
                label corr = 0;
                corr < max(implicitNonlinearIterations, label(1));
                ++corr
            )
            {
                const scalarField VmPrevious(VmGuess.primitiveField());
                const scalarField u1Previous(u1Guess.primitiveField());
                const scalarField u2Previous(u2Guess.primitiveField());
                const scalarField u3Previous(u3Guess.primitiveField());
                const scalarField u1IPPrevious(u1GuessIntegrationPoints);
                const scalarField u2IPPrevious(u2GuessIntegrationPoints);
                const scalarField u3IPPrevious(u3GuessIntegrationPoints);
                const scalarField IionPrevious(IionGuess.primitiveField());

                evaluateNonlinearFields
                (
                    VmGuess,
                    u1Guess,
                    u2Guess,
                    u3Guess,
                    VmGuessIntegrationPoints,
                    u1GuessIntegrationPoints,
                    u2GuessIntegrationPoints,
                    u3GuessIntegrationPoints,
                    sourceVm,
                    IionGuess
                );

                EigVec rhs = rhsFromSource(sourceVm);

                if (needsACurrent)
                {
                    if (!ACurrentInitialised)
                    {
                        ACurrent = AImplicit;
                        ACurrentInitialised = true;
                    }
                    else
                    {
                        // Refresh values only; sparsity pattern is invariant
                        // (M, K and AImplicit are time-constant and only the
                        // diagonal source term gets added below).
                        std::copy
                        (
                            AImplicit.valuePtr(),
                            AImplicit.valuePtr() + AImplicit.nonZeros(),
                            ACurrent.valuePtr()
                        );
                    }
                }

                if (useDiagonalIion)
                {
                    computeSourceDerivative();
                    const EigVec sourceDerivativeVec =
                        fieldToEigVec(sourceDerivative);
                    const EigVec VmLinearisationPoint =
                        fieldToEigVec(VmGuess);

                    addDiagonalToMatrix
                    (
                        sourceDerivative.primitiveField(),
                       -theta,
                        ACurrent
                    );

                    rhs -= theta
                       * sourceDerivativeVec.cwiseProduct(VmLinearisationPoint);
                }

                label linearIterations = 0;
                scalar linearError = GREAT;

                const auto sparseSolveStart =
                    std::chrono::steady_clock::now();
                EigVec Vsol(rhs.size());
                MatrixUpdate updatePolicy = MatrixUpdate::Reuse;
                if (linearSolverNeedsRebuild)
                {
                    updatePolicy = MatrixUpdate::Rebuild;
                    linearSolverNeedsRebuild = false;
                }
                else if (useDiagonalIion)
                {
                    updatePolicy = MatrixUpdate::ValuesOnly;
                }

                Vsol = solveSparseSystem
                (
                    useDiagonalIion ? ACurrent : AImplicit,
                    rhs,
                    linearSolverBackend,
                    implicitLinearSolver,
                    petscLinearKspType,
                    petscLinearPcType,
                    petscLinearRestart,
                    petscLinearOptionsPrefix,
                    petscUseOptions,
                    implicitTolerance,
                    implicitMaxIterations,
                    linearIterations,
                    linearError,
                    persistentPetscSolver,
                    persistentEigenSolvers,
                    updatePolicy,
                    petscReusePreconditioner
                );
                if (profileTimings)
                {
                    const auto sparseSolveEnd =
                        std::chrono::steady_clock::now();
                    fineTimings.sparseLinearSolveWallTime +=
                        std::chrono::duration<scalar>
                        (
                            sparseSolveEnd - sparseSolveStart
                        ).count();
                    ++fineTimings.sparseLinearSolveCalls;
                }

                const EigVec Vprev = fieldToEigVec(VmGuess);
                const EigVec Vrelaxed =
                    Vprev + nonlinearRelaxation*(Vsol - Vprev);

                eigVecToField(Vrelaxed, VmGuess);
                VmGuess.correctBoundaryConditions();

                if (useDiagonalIion)
                {
                    evaluateNonlinearFields
                    (
                        VmGuess,
                        u1Guess,
                        u2Guess,
                        u3Guess,
                        VmGuessIntegrationPoints,
                        u1GuessIntegrationPoints,
                        u2GuessIntegrationPoints,
                        u3GuessIntegrationPoints,
                        sourceVm,
                        IionGuess
                    );
                }
                else
                {
                    updateStateBoundaryValues(u1Guess, u2Guess, u3Guess, t + dt, dim);
                    u1Guess.correctBoundaryConditions();
                    u2Guess.correctBoundaryConditions();
                    u3Guess.correctBoundaryConditions();
                }

                const EigVec nonlinearResidualVec =
                    applyA(fieldToEigVec(VmGuess)) - rhsFromSource(sourceVm);
                const scalar newtonResidual =
                    relativeL2Norm(nonlinearResidualVec, rhsFromSource(sourceVm));

                const scalar VmResidual =
                    relativeL2Difference(VmGuess.primitiveField(), VmPrevious);
                scalar u1Residual = GREAT;
                scalar u2Residual = GREAT;
                scalar u3Residual = GREAT;
                if (useHighOrder_Iion && dim > 1 && !reconstructStatesFromCellCentres)
                {
                    u1Residual =
                        relativeL2Difference(u1GuessIntegrationPoints, u1IPPrevious);
                    u2Residual =
                        relativeL2Difference(u2GuessIntegrationPoints, u2IPPrevious);
                    u3Residual =
                        relativeL2Difference(u3GuessIntegrationPoints, u3IPPrevious);
                }
                else
                {
                    // cellCentredReconstruct (and the low-order path) carry the
                    // authoritative states at cell centres -> measure there.
                    u1Residual =
                        relativeL2Difference(u1Guess.primitiveField(), u1Previous);
                    u2Residual =
                        relativeL2Difference(u2Guess.primitiveField(), u2Previous);
                    u3Residual =
                        relativeL2Difference(u3Guess.primitiveField(), u3Previous);
                }
                const scalar IionResidual =
                    relativeL2Difference(IionGuess.primitiveField(), IionPrevious);

                // Shared stopping test (monodomainBase.H). Picard used to stop
                // here on the Vm and state increments alone.
                const bool converged = nonlinearIterationStops
                (
                    nonlinearIterationExit
                    (
                        corr,
                        implicitMinNonlinearIterations,
                        newtonResidual,
                        coupledResidualPrevIter,
                        nonlinearTolerance,
                        VmResidual,
                        nonlinearVmTolerance,
                        IionResidual,
                        nonlinearIionTolerance,
                        nonlinearRequireStatesConvergence,
                        max(u1Residual, max(u2Residual, u3Residual)),
                        nonlinearStatesTolerance
                    ),
                    newtonResidual,
                    nonlinearTolerance,
                    nResidualFloorExits
                );
                coupledResidualPrevIter = newtonResidual;

                writeAndPrintResiduals
                (
                    corr + 1,
                    linearIterations,
                    linearError,
                    newtonResidual,
                    VmResidual,
                    u1Residual,
                    u2Residual,
                    u3Residual,
                    IionResidual,
                    converged
                );

                nonlinearIters = corr + 1;
                if (converged)
                {
                    nonlinearConverged = true;
                    break;
                }
            }
        }

        bool nonlinearRolledBack = false;

        // Picard used to have its own branch here
        // that accepted the last iterate unconditionally, so the same solver
        // had two error semantics depending on the method. All three now obey
        // nonlinearAcceptUnconverged: false rolls the step back, true keeps the
        // iterate and says so.
        if (!nonlinearConverged && !nonlinearAcceptUnconverged)
        {
            WarningInFunction
                << "Nonlinear solver did not converge at t = " << t
                << " after " << nonlinearIters << " iterations."
                << " Rolling back Vm, Iion and states to the values"
                << " at the beginning of the time step." << nl
                << "Set 'nonlinearAcceptUnconverged true;' in"
                << " spatialIntegrationProperties to accept unconverged"
                << " iterates instead." << endl;

            VmGuess.primitiveFieldRef() = VmOld.primitiveField();
            VmGuess.correctBoundaryConditions();
            u1Guess.primitiveFieldRef() = u1Old.primitiveField();
            u2Guess.primitiveFieldRef() = u2Old.primitiveField();
            u3Guess.primitiveFieldRef() = u3Old.primitiveField();
            updateStateBoundaryValues(u1Guess, u2Guess, u3Guess, t, dim);
            u1Guess.correctBoundaryConditions();
            u2Guess.correctBoundaryConditions();
            u3Guess.correctBoundaryConditions();
            IionGuess.primitiveFieldRef() = Iion.primitiveField();
            IionGuess.correctBoundaryConditions();
            nonlinearRolledBack = true;

            if (useHighOrder_Iion && dim > 1)
            {
                u1GuessIntegrationPoints = u1OldIntegrationPoints;
                u2GuessIntegrationPoints = u2OldIntegrationPoints;
                u3GuessIntegrationPoints = u3OldIntegrationPoints;
            }

            if (++nConsecutiveRollbacks > nonlinearMaxConsecutiveRollbacks)
            {
                FatalErrorInFunction
                    << nConsecutiveRollbacks << " consecutive time steps"
                    << " rolled back (nonlinearMaxConsecutiveRollbacks = "
                    << nonlinearMaxConsecutiveRollbacks << "), the last at t = "
                    << t << ". The solution has been frozen since the first"
                    << " of them; tighten dt or loosen the nonlinear"
                    << " tolerances." << exit(FatalError);
            }
        }
        else if (!nonlinearConverged && nonlinearAcceptUnconverged)
        {
            WarningInFunction
                << "Nonlinear solver did not converge at t = " << t
                << " after " << nonlinearIters << " iterations."
                << " Accepting the unconverged iterate"
                << " (nonlinearAcceptUnconverged=true)." << endl;
        }

        if (!nonlinearRolledBack)
        {
            nConsecutiveRollbacks = 0;
        }

        NonlinearConvergenceRecord nonlinearRecord;
        nonlinearRecord.time = t + dt;
        nonlinearRecord.step = nSteps + 1;
        nonlinearRecord.iterations = nonlinearIters;
        nonlinearRecord.maxIterations =
            max(implicitNonlinearIterations, label(1));
        nonlinearRecord.converged = nonlinearConverged;
        nonlinearRecord.rolledBack = nonlinearRolledBack;
        nonlinearRecord.coupledResidual = finalCoupledResidual;
        nonlinearRecord.VmResidual = finalVmResidual;
        nonlinearRecord.u1Residual = finalU1Residual;
        nonlinearRecord.u2Residual = finalU2Residual;
        nonlinearRecord.u3Residual = finalU3Residual;
        nonlinearRecord.maxStateResidual = finalMaxStateResidual;
        nonlinearRecord.IionResidual = finalIionResidual;
        nonlinearRecord.linearIterations = finalLinearIterations;
        nonlinearRecord.linearError = finalLinearError;
        nonlinearRecord.lineSearchIterations = finalLineSearchIterations;
        nonlinearHistory.push_back(nonlinearRecord);

        Vm.primitiveFieldRef() = VmGuess.primitiveField();
        u1.primitiveFieldRef() = u1Guess.primitiveField();
        u2.primitiveFieldRef() = u2Guess.primitiveField();
        u3.primitiveFieldRef() = u3Guess.primitiveField();

        if (useHighOrder_Iion && dim > 1)
        {
            u1IntegrationPoints = u1GuessIntegrationPoints;
            u2IntegrationPoints = u2GuessIntegrationPoints;
            u3IntegrationPoints = u3GuessIntegrationPoints;
        }

        const scalar committedTime = nonlinearRolledBack ? t : t + dt;
        Vm.correctBoundaryConditions();
        updateStateBoundaryValues(u1, u2, u3, committedTime, dim);
        u1.correctBoundaryConditions();
        u2.correctBoundaryConditions();
        u3.correctBoundaryConditions();

        if (useHighOrder_Iion && dim > 1)
        {
            reconstructVmAtIionIntegrationPoints
            (
                Vm,
                useHighOrder_Vm,
                LREInterp_Vm,
                LREInterp_Iion,
                VmGuessIntegrationPoints,
                stateODEUseOpenMP,
                stateODEOpenMPThreshold,
                profileTimings,
                fineTimings
            );

            if (reconstructStatesFromCellCentres)
            {
                // Keep the Gauss-point states consistent with the committed
                // cell-centred states (matters after a rollback too).
                reconstructStatesAtIionIntegrationPoints
                (
                    u1,
                    u2,
                    u3,
                    LREInterp_statesPtr(),
                    LREInterp_Iion,
                    u1IntegrationPoints,
                    u2IntegrationPoints,
                    u3IntegrationPoints,
                    stateODEUseOpenMP,
                    stateODEOpenMPThreshold,
                    profileTimings,
                    fineTimings
                );
            }

            computeIionFromIntegrationPoints
            (
                VmGuessIntegrationPoints,
                u1IntegrationPoints,
                u2IntegrationPoints,
                u3IntegrationPoints,
                beta,
                chiVal,
                CmVal,
                IionIntegrationPoints,
                stateODEUseOpenMP,
                stateODEOpenMPThreshold
            );

            averageIntegrationPointFieldToCells
            (
                IionIntegrationPoints,
                LREInterp_Iion,
                Iion
            );
        }
        else
        {
            computeCellCentredIion
            (
                Vm,
                u1,
                u2,
                u3,
                beta,
                chiVal,
                CmVal,
                Iion,
                stateODEUseOpenMP,
                stateODEOpenMPThreshold
            );
        }

        // Rotate the BDF history. V^n becomes
        // V^{n-1} and the oldest level falls off the end. Done here, after the
        // step has been accepted, so that what enters the history is the value
        // the step actually kept.
        //
        // Modified for cardiacFoam: a ROLLED-BACK step (nonlinearAcceptUnconverged
        // false) advanced time but not the solution, so rotating would leave
        // two identical levels one dt apart - a history no solution has, which
        // the next k-1 steps would interpolate. It RESTARTS the history
        // instead: the next k-1 steps are taken by ESDIRK3 from the value kept.
        // Only with the ESDIRK3 startup; exact/constant seeding exist to
        // measure a bad startup and have no ESDIRK operator, so they keep the
        // old rotation. Same rule as the physiological solver.
        if (bdfOrder > 1)
        {
            const bool rolledBack =
                !nonlinearConverged && !nonlinearAcceptUnconverged;

            if (rolledBack && esdirkStartup)
            {
                bdfValidLevels = 0;
            }
            else
            {
                for (label j = VmBdfOlder.size() - 1; j >= 1; --j)
                {
                    VmBdfOlder[j] = VmBdfOlder[j - 1];
                }

                VmBdfOlder[0] = VmOld.primitiveField();

                bdfValidLevels = min(bdfValidLevels + 1, bdfOrder - 1);
            }
        }

        ++nSteps;
        ++runTime;

        const bool needsRhsVm =
            nSteps % stepLogInterval == 0 || nSteps <= 5 || runTime.outputTime();

        if (needsRhsVm)
        {
            updateNumericalLaplacian(t + dt, runTime.outputTime());
            rhsVm = lapVm/(chi*Cm) - Iion;
        }

        if (nSteps % stepLogInterval == 0 || nSteps <= 5)
        {
            Info<< "Time step " << nSteps
                << " : t = " << runTime.value()
                << " , RHS Linf = " << linfNorm(rhsVm)
                << endl;
        }

        if (runTime.outputTime())
        {
            fillExactFields(VmExact, u1Exact, u2Exact, runTime.value(), dim);
            VmError = Vm - VmExact;
            u1Error = u1 - u1Exact;
            u2Error = u2 - u2Exact;
            runTime.write();
        }
    }

    const auto tEndLoop = std::chrono::steady_clock::now();
    const auto tStartPost = tEndLoop;

    fillExactFields(VmExact, u1Exact, u2Exact, runTime.value(), dim);

    if (useHighOrder_Iion && dim > 1)
    {
        reconstructVmAtIionIntegrationPoints
        (
            Vm,
            useHighOrder_Vm,
            LREInterp_Vm,
            LREInterp_Iion,
            VmGuessIntegrationPoints,
            stateODEUseOpenMP,
            stateODEOpenMPThreshold,
            profileTimings,
            fineTimings
        );

        computeIionFromIntegrationPoints
        (
            VmGuessIntegrationPoints,
            u1IntegrationPoints,
            u2IntegrationPoints,
            u3IntegrationPoints,
            beta,
            chiVal,
            CmVal,
            IionIntegrationPoints,
            stateODEUseOpenMP,
            stateODEOpenMPThreshold
        );

        averageIntegrationPointFieldToCells
        (
            IionIntegrationPoints,
            LREInterp_Iion,
            Iion
        );
    }
    else
    {
        computeCellCentredIion
        (
            Vm,
            u1,
            u2,
            u3,
            beta,
            chiVal,
            CmVal,
            Iion,
            stateODEUseOpenMP,
            stateODEOpenMPThreshold
        );
    }

    updateNumericalLaplacian(runTime.value(), true);
    rhsVm = lapVm/(chi*Cm) - Iion;

    VmError = Vm - VmExact;
    u1Error = u1 - u1Exact;
    u2Error = u2 - u2Exact;

    const FieldErrorSummary VmCell = computeFieldErrorSummary(VmError);
    const FieldErrorSummary u1Cell = computeFieldErrorSummary(u1Error);
    const FieldErrorSummary u2Cell = computeFieldErrorSummary(u2Error);

    const ReconstructionErrorSummary VmHO =
            useHighOrder_Vm && dim > 1
        ? computeGaussReconstructedError(Vm, runTime.value(), mfVm, dim, LREInterp_Vm)
        : cellAsReconstructionError(VmCell, VmExact);

    const ReconstructionErrorSummary u1HO =
            useHighOrder_Iion && dim > 1
        ? computeIntegrationPointStateError
          (
              u1IntegrationPoints,
              runTime.value(),
              mfU1,
              dim,
              LREInterp_Iion,
              mesh
          )
        : cellAsReconstructionError(u1Cell, u1Exact);

    const ReconstructionErrorSummary u2HO =
            useHighOrder_Iion && dim > 1
        ? computeIntegrationPointStateError
          (
              u2IntegrationPoints,
              runTime.value(),
              mfU2,
              dim,
              LREInterp_Iion,
              mesh
          )
        : cellAsReconstructionError(u2Cell, u2Exact);

    const auto tEndPost = std::chrono::steady_clock::now();

    const scalar setupWallTime =
        std::chrono::duration_cast<std::chrono::duration<scalar>>
        (
            tEndSetup - tStartSetup
        ).count();

    const scalar timeLoopWallTime =
        std::chrono::duration_cast<std::chrono::duration<scalar>>
        (
            tEndLoop - tStartLoop
        ).count();

    const scalar postProcessWallTime =
        std::chrono::duration_cast<std::chrono::duration<scalar>>
        (
            tEndPost - tStartPost
        ).count();

    const scalar totalWallTime =
        std::chrono::duration_cast<std::chrono::duration<scalar>>
        (
            tEndPost - tStartTotal
        ).count();

    if (profileTimings)
    {
        Info<< "Fine-grained timing [s]:" << nl
            << "  nonlinear evaluations: calls = "
            << fineTimings.nonlinearEvalCalls
            << ", wall = " << fineTimings.nonlinearEvalWallTime << nl
            << "  JFNK KSP/GMRES: calls = " << fineTimings.gmresCalls
            << ", wall = " << fineTimings.gmresWallTime << nl
            << "  JFNK preconditioner setup: calls = "
            << fineTimings.preconditionerSetups
            << ", wall = " << fineTimings.preconditionerSetupWallTime << nl
            << "  JFNK preconditioner apply: calls = "
            << fineTimings.preconditionerApplications
            << ", wall = " << fineTimings.preconditionerApplyWallTime << nl
            << "  sparse PDE solves: calls = "
            << fineTimings.sparseLinearSolveCalls
            << ", wall = " << fineTimings.sparseLinearSolveWallTime << endl;
    }

    if (writeFinalFields)
    {
        // writeNow() fuerza la escritura: write() sola no hace nada si el
        // instante final no es un write time, que es el caso habitual.
        // Los campos exacto/error solo se refrescan dentro de la rama
        // outputTime, asi que sin esto quedarian con valores viejos.
        fillExactFields(VmExact, u1Exact, u2Exact, runTime.value(), dim);
        VmError = Vm - VmExact;
        u1Error = u1 - u1Exact;
        u2Error = u2 - u2Exact;
        runTime.writeNow();
    }
    else
    {
        runTime.write();
    }

    writeSummary
    (
        runTime,
        nSteps,
        dt,
        VmCell,
        u1Cell,
        u2Cell,
        VmHO,
        u1HO,
        u2HO,
        rhsVm,
        Iion,
        mesh,
        implicitScheme,
        bdfStartup,
        stateODEDriver,
        manufacturedSolution,
        massMatrixMode,
        useHighOrder_Vm,
        useHighOrder_Iion,
        stabilisationAlpha,
        memoryOptimization,
        memoryOptimizationEffective,
        memoryOptimizationCellThreshold,
        memoryOptimizationAdaptiveTripletReserve,
        memoryOptimizationFluxTripletsPerFace,
        memoryOptimizationStabilisedTripletsPerFace,
        memoryOptimizationTrimHeap,
        memoryOptimizationCompactMassAssembly,
        compactMassAssembly,
        stiffnessTripletsPerFaceReserve,
        allocatedIionIntegrationPoints,
        lreN,
        lreNn,
        lreK,
        lreMaxStencil,
        lreIionN,
        lreIionNn,
        lreIionK,
        lreIionMaxStencil,
        linearSolverBackend,
        implicitLinearSolver,
        petscLinearKspType,
        petscLinearPcType,
        implicitTolerance,
        implicitMaxIterations,
        max(implicitNonlinearIterations, label(1)),
        implicitMinNonlinearIterations,
        nonlinearTolerance,
        nonlinearVmTolerance,
        nonlinearStatesTolerance,
        nonlinearRequireStatesConvergence,
        nonlinearAcceptUnconverged,
        nonlinearIionTolerance,
        nonlinearRelaxation,
        jfnkMaxKrylovIterations,
        jfnkMaxRestarts,
        jfnkLinearTolerance,
        jfnkEpsilon,
        jfnkLinearSolverBackend,
        jfnkPetscKspType,
        jfnkPetscPcType,
        jfnkInitGuessOrder,
        jfnkInitGuessVmMin,
        jfnkInitGuessVmMax,
        jfnkClampODEInput,
        jfnkLineSearch,
        jfnkLineSearchMaxIter,
        jfnkLineSearchAlphaMin,
        jfnkArmijoC,
        jfnkPreconditioner,
        jfnkPreconditionerDropTolerance,
        jfnkPreconditionerFillFactor,
        jfnkPreconditionerUpdateFrequency,
        diagonalIionEpsilon,
        stateODESolver,
        stateODEInitialStep,
        stateODEAbsTol,
        stateODERelTol,
        stateODEMaxSteps,
        stateODEUseOpenMP,
        stateODENumThreads,
        nonlinearHistory,
        1024.0*currentPeakRSSMB(),
        setupWallTime,
        timeLoopWallTime,
        postProcessWallTime,
        totalWallTime,
        fineTimings,
        nonlinearMethod
    );

    Info<< "Nonlinear iterations stopped at the residual floor: "
        << nResidualFloorExits << nl;

    Info<< "End" << nl << endl;
    return 0;
}
