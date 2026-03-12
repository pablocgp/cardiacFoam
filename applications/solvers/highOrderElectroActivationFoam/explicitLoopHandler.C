/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "explicitLoopHandler.H"
#include "tmanufacturedFDA.H"
#include "fvm.H"
#include "fvc.H"

// HIGH ORDER //
#include "LRE.H"
#include "deltaVectors.H"
// HIGH ORDER //

explicitLoopHandler::explicitLoopHandler
(
    const fvMesh& mesh,
    ionicModel& model
)
:
    // mesh_(mesh),
    ionicModel_(model),
    dx_(0.0),
    dim_(0),
    Def_max_(0.0)
{}



//----------------------------------------------------------------------
// Compute stability limit dt <= CFL * dx^2 / (dim * Def_max)
//----------------------------------------------------------------------
scalar explicitLoopHandler::computeStableDt(const scalar CFL) const
{
    if (dim_ == 0 || Def_max_ <= SMALL)
        return GREAT;

    return CFL * dx_ * dx_ / (dim_ * Def_max_);
}


//----------------------------------------------------------------------
// Initialization: compute Def_max, stable dt
//----------------------------------------------------------------------
void explicitLoopHandler::initializeExplicit
(
    scalar& dt,
    int& nsteps,
    const scalar chiVal,
    const scalar CmVal,
    const volTensorField& conductivity,
    const scalar CFL,
    const scalar dx,
    const int dim
)
{
    // Store caller-provided dx and dim (includes MS corrections)
    dx_  = dx;
    dim_ = dim;

    // Extract max diffusion coefficient
    scalar Dxx = max(conductivity.component(tensor::XX)).value();
    scalar Dyy = max(conductivity.component(tensor::YY)).value();
    scalar Dzz = max(conductivity.component(tensor::ZZ)).value();

    scalar Dmax_raw = max(Dxx, max(Dyy, Dzz));

    // Effective diffusion coefficient: D / (chi * Cm)
    Def_max_ = Dmax_raw / (chiVal * CmVal);

    // Compute stable dt
    scalar dt_stable = computeStableDt(CFL);

    dt = min(dt, dt_stable);
    reduce(dt, maxOp<scalar>());

    Info << "ExplicitLoop: dx = " << dx_
         << ", dim = " << dim_
         << ", Def_max = " << Def_max_
         << ", stable dt = " << dt_stable
         << ", using dt = " << dt << nl;

    // nsteps recomputed in solver; we leave it untouched here
    (void)nsteps;
}

//----------------------------------------------------------------------
// Perform a single explicit step
//----------------------------------------------------------------------

// HIGH ORDER //

void explicitLoopHandler::highOrderExplicitLoop
(
    const scalar t0,
    const scalar dt,
    volScalarField& Vm,
    volScalarField& Iion,
    Field<Field<scalar>>& states,
    volScalarField& externalStimulusCurrent,
    const List<labelList>& stimulusCellIDsList,
    const List<scalar>& stimulusStartTimes,
    const scalar stimulusIntensity,
    const scalar stimulusDuration,
    const dimensionedScalar& chi,
    const dimensionedScalar& Cm,
    const volTensorField& conductivity,
    scalarField& VmIntegrationPoints,
    scalarField& IionIntegrationPoints,
    const LRE& LREInterp_Vm,
    const LRE& LREInterp_Iion,
    surfaceVectorField& surfaceGradVm_HO,
    const Switch useHighOrder_Iion,
    const Switch useHighOrder_Vm,
    volScalarField& lapVm
)
{
    const fvMesh& mesh = Vm.mesh();

    // 0) Manufactured-specific: store "old" internal states
    Info << "0) Manufactured-specific: store old internal states" << endl;
    if (ionicModel_.hasManufacturedSolution())
    {
        refCast<tmanufacturedFDA>(ionicModel_).updateStatesOld();
    }

    // 1) External stimulus (internal field only)
    Info << "1) External stimulus (internal field only)" << endl;
    scalarField& externalStimulusCurrentI = externalStimulusCurrent;
    externalStimulusCurrentI = 0.0;

    forAll(stimulusCellIDsList, bI)
    {
        const scalar tStart = stimulusStartTimes[bI];
        if (t0 < tStart || t0 > (tStart + stimulusDuration))
        {
            continue;
        }

        const labelList& stimulusCellIDs = stimulusCellIDsList[bI];
        forAll(stimulusCellIDs, cI)
        {
            const label id = stimulusCellIDs[cI];
            externalStimulusCurrentI[id] = stimulusIntensity;
        }
    }
    externalStimulusCurrent.correctBoundaryConditions();

    // 2) Ionic current using OLD Vm
    Info << "2) Ionic current using OLD Vm" << endl;
    if ( useHighOrder_Iion)
    {
        // Obtaining Vm at Gauss Points

        const vectorField& C = mesh.C();
        const CompactListList<scalar>& cellIionQuadW = LREInterp_Iion.cellQuadWeight();
        const CompactListList<point>& cellIionQuadP = LREInterp_Iion.cellQuadPoints();

        // Gradient
        tmp<volVectorField> tGradVm = LREInterp_Vm.grad(Vm);
        const vectorField& gradVm_HO = tGradVm->internalField();

        // Hessian
        tmp<volSymmTensorField> tHessVm = LREInterp_Vm.hessian(Vm);
        const symmTensorField& HessVm_HO = tHessVm->internalField();

        // Third derivative
        autoPtr<List<LRE::symmTensor3Order>> ThirdDerVm_HO = LREInterp_Vm.thirdDeriv(Vm);

        label integrationPointPos = 0;

        forAll(mesh.cells(), cellI)
        {
            const scalar Vc = Vm[cellI];
            const vector& gradVc = gradVm_HO[cellI];
            const vector& xc = C[cellI];
            forAll(cellIionQuadP[cellI], gI)
            {
                vector dx = cellIionQuadP[cellI][gI] - xc;
                // Linear part
                scalar Vg = Vc + (gradVc & dx) ;
                // Quadratic part
                if (LREInterp_Vm.order() >= 2)
                {
                    Vg += 0.5*(dx & (HessVm_HO[cellI] & dx));
                }
                // Cubic part
                if (LREInterp_Vm.order() >= 3)
                {
                    const bool twoD = (mesh.nGeometricD() == 2);
                    const scalar oneOverSix = 1.0/6.0;

                    Vg += oneOverSix *
                        LRE::cubicForm((*ThirdDerVm_HO)[cellI], dx, twoD);
                }
                VmIntegrationPoints[integrationPointPos] = Vg;
                integrationPointPos++;
            }
        }
        // Obtaining Iion at Gauss Points
        ionicModel_.calculateCurrent
        (
            t0,
            dt,
            VmIntegrationPoints,
            IionIntegrationPoints,
            states
        );
        // Obtaining Iion at Cell centers
        integrationPointPos = 0;
        forAll(mesh.cells(), cellI)
        {
            Iion[cellI] = 0.0;
            double totCellQuadW = 0.0;
            forAll(cellIionQuadP[cellI], gI)
            {
                Iion[cellI] += cellIionQuadW[cellI][gI]*IionIntegrationPoints[integrationPointPos];
                totCellQuadW += cellIionQuadW[cellI][gI];
                integrationPointPos++;
            }
            Iion[cellI] /= totCellQuadW; // this should be 1.. but just in case
        }
    }
    else
    {
        // Obtaining Iion at Cells centers
        ionicModel_.calculateCurrent
        (
            t0,
            dt,
            Vm.internalField(),   // Vm_old
            Iion,
            states
        );
    }
    Iion.correctBoundaryConditions();

    if (ionicModel_.hasManufacturedSolution())
    {
        Iion /= Cm.value();
    }

    // Discretization of the laplacian
    //
    if ( useHighOrder_Vm)
    {
        // Obtaining a HO laplacian aproximation
        autoPtr<List<List<vector>>> gradVmQuad_ptr = LREInterp_Vm.gradScalarFaceQuad(Vm);
        List<List<vector>>& gradVmQuad = gradVmQuad_ptr.ref();
        const CompactListList<scalar>& faceVmQuadW = LREInterp_Vm.faceQuadWeight();

        // Internal faces
        forAll(surfaceGradVm_HO, faceI)
        {
            const label ownerCell = mesh.owner()[faceI];
            surfaceGradVm_HO[faceI] = vector::zero;
            forAll(gradVmQuad[faceI], pI)
            {
                surfaceGradVm_HO[faceI] += (conductivity[ownerCell] & gradVmQuad[faceI][pI]) * faceVmQuadW[faceI][pI];
            }
        }

        // stabilization
        // const labelList& own = mesh.owner();
        // const labelList& nei = mesh.neighbour();

        // const surfaceVectorField n(mesh.Sf()/mesh.magSf());
        // const vectorField& nI = n.internalField();

        // const vectorField deltaI(deltaVectors(mesh));
        // const surfaceVectorField& Cf = mesh.Cf();
        // const vectorField& C = mesh.C();

        // const scalarField& VmI = Vm.internalField();

        // // High-order gradient
        // tmp<volVectorField> tGradVm = LREInterp_Vm.grad(Vm);
        // const vectorField& gradVmI = tGradVm->internalField();

        // // Hessian
        // tmp<volSymmTensorField> tHessVm = LREInterp_Vm.hessian(Vm);
        // const symmTensorField& hessVmI = tHessVm->internalField();

        // // Third derivative
        // autoPtr<List<LRE::symmTensor3Order>> tThirdDerVm;

        // if (LREInterp_Vm.order() >= 3)
        // {
        //     tThirdDerVm = LREInterp_Vm.thirdDeriv(Vm);
        // }

        // forAll(nei, faceI)
        // {
        //     const label ownCellID = own[faceI];
        //     const label neiCellID = nei[faceI];

        //     const scalar& ownVm = VmI[ownCellID];
        //     const scalar& neiVm = VmI[neiCellID];

        //     const vector& ownGradVm = gradVmI[ownCellID];
        //     const vector& neiGradVm = gradVmI[neiCellID];

        //     const vector& d = deltaI[faceI];

        //     const vector dOwn = Cf[faceI] - C[ownCellID];
        //     const vector dNei = Cf[faceI] - C[neiCellID];

        //     // Linear part
        //     scalar extrapOwnVm = ownVm + (ownGradVm & dOwn);
        //     scalar extrapNeiVm = neiVm + (neiGradVm & dNei);

        //     // Quadratic part
        //     if (LREInterp_Vm.order() >= 2)
        //     {
        //         extrapOwnVm += 0.5*(dOwn & (hessVmI[ownCellID] & dOwn));
        //         extrapNeiVm += 0.5*(dNei & (hessVmI[neiCellID] & dNei));
        //     }

        //     // Cubic part
        //     if (LREInterp_Vm.order() >= 3)
        //     {
        //         const bool twoD = (mesh.nGeometricD() == 2);
        //         const scalar oneOverSix = 1.0/6.0;

        //         extrapOwnVm += oneOverSix *
        //             LRE::cubicForm((*tThirdDerVm)[ownCellID], dOwn, twoD);

        //         extrapNeiVm += oneOverSix *
        //             LRE::cubicForm((*tThirdDerVm)[neiCellID], dNei, twoD);
        //     }

        //     const vector& n = nI[faceI];
        //     const scalar denom = max(mag(n & d), VSMALL);

        //     scalar h = mag(deltaI[faceI]);
        //     scalar alpha_val = 0.1;

        //     // alpha_val = scaleFactor*impKf_[faceI] this should be changed
        //     const scalar faceDamping =
        //         alpha_val*(extrapNeiVm - extrapOwnVm)/denom;

        //     surfaceGradVm_HO[faceI] += faceDamping * n;
        // }

        surfaceGradVm_HO.correctBoundaryConditions();


        lapVm = fvc::div(mesh.Sf() & surfaceGradVm_HO ) 
                // + fvc::laplacian(conductivity,Vm) 
                // - fvc::div(conductivity & fvc::grad(Vm))
                ;
    }
    else
    {
        // Obtaining a standar laplacian aproximation
        lapVm = fvc::laplacian(conductivity,Vm);
    }

    // 3) Solving equation
    Info << "// 3) Solving equation" << endl;
    solve
    (
        chi*Cm*fvm::ddt(Vm)
     == lapVm
      - chi*Cm*Iion
      + externalStimulusCurrent
    );

    Vm.correctBoundaryConditions();

    // 4) Manufactured-specific: reset internal states back to OLD
    Info << "4) Manufactured-specific: reset internal states back to OLD" << endl;
    if (ionicModel_.hasManufacturedSolution())
    {
        refCast<tmanufacturedFDA>(ionicModel_).resetStatesToStatesOld();
    }

    // 5) Advance ionic model in time (ODE solve) with NEW Vm, once per timestep
    Info << "5) Advance ionic model in time (ODE solve) with NEW Vm, once per timestep)" << endl;
    if ( useHighOrder_Iion)
    {

        const vectorField& C = mesh.C();
        const CompactListList<scalar>& cellIionQuadW = LREInterp_Iion.cellQuadWeight();
        const CompactListList<point>& cellIionQuadP = LREInterp_Iion.cellQuadPoints();

        // Gradient
        tmp<volVectorField> tGradVm = LREInterp_Vm.grad(Vm);
        const vectorField& gradVm_HO = tGradVm->internalField();
        tmp<volVectorField> tGradIion = LREInterp_Iion.grad(Iion);
        const vectorField& gradIion_HO = tGradIion->internalField();

        // Hessian
        tmp<volSymmTensorField> tHessVm = LREInterp_Vm.hessian(Vm);
        const symmTensorField& HessVm_HO = tHessVm->internalField();
        tmp<volSymmTensorField> tHessIion = LREInterp_Iion.hessian(Iion);
        const symmTensorField& HessIion_HO = tHessIion->internalField();

        // Third derivative
        autoPtr<List<LRE::symmTensor3Order>> ThirdDerVm_HO = LREInterp_Vm.thirdDeriv(Vm);
        autoPtr<List<LRE::symmTensor3Order>> ThirdDerIion_HO = LREInterp_Iion.thirdDeriv(Iion);

        label integrationPointPos = 0;

        forAll(mesh.cells(), cellI)
        {
            const scalar Vc = Vm[cellI];
            const scalar Iionc = Iion[cellI];
            const vector& gradVc = gradVm_HO[cellI];
            const vector& gradIionc = gradIion_HO[cellI];
            const vector& xc = C[cellI];

            forAll(cellIionQuadP[cellI], gI)
            {
                vector dx = cellIionQuadP[cellI][gI] - xc;

                // Linear part
                scalar Vg = Vc + (gradVc & dx);
                scalar Iiong = Iionc + (gradIionc & dx);
                // Quadratic part
                if (LREInterp_Vm.order() >= 2)
                {
                    Vg += 0.5*(dx & (HessVm_HO[cellI] & dx));
                }
                if (LREInterp_Iion.order() >= 2)
                {
                    Iiong += 0.5*(dx & (HessIion_HO[cellI] & dx));
                }
                // Cubic part
                if (LREInterp_Vm.order() >= 3)
                {
                    const bool twoD = (mesh.nGeometricD() == 2);
                    const scalar oneOverSix = 1.0/6.0;

                    Vg += oneOverSix *
                        LRE::cubicForm((*ThirdDerVm_HO)[cellI], dx, twoD);
                }
                if (LREInterp_Iion.order() >= 3)
                {
                    const bool twoD = (mesh.nGeometricD() == 2);
                    const scalar oneOverSix = 1.0/6.0;

                    Iiong += oneOverSix *
                        LRE::cubicForm((*ThirdDerIion_HO)[cellI], dx, twoD);
                }

                VmIntegrationPoints[integrationPointPos] = Vg;
                IionIntegrationPoints[integrationPointPos] = Iiong;
                integrationPointPos++;
            }
        }
        // Obtaining update of Iion at Gauss Points
        Info << "Ini solve ODE 5)" << endl;
        ionicModel_.solveODE
        (
            t0,
            dt,
            VmIntegrationPoints,    // Vm_new
            IionIntegrationPoints,
            states
        );
        Info << "End solve ODE 5)" << endl;
        // Obtaining update of Iion at cell centers
        integrationPointPos = 0;
        forAll(mesh.cells(), cellI)
        {
            Iion[cellI] = 0.0;
            double totCellQuadW = 0.0;
            forAll(cellIionQuadP[cellI], gI)
            {
                Iion[cellI] += cellIionQuadW[cellI][gI]*IionIntegrationPoints[integrationPointPos];
                totCellQuadW += cellIionQuadW[cellI][gI];
                integrationPointPos++;
            }
            Iion[cellI] /= totCellQuadW; // this should be 1.. but just in case
        }
        Iion.correctBoundaryConditions(); // correct????
    }
    else
    {
        // Obtaining update of Iion at cell centers
        ionicModel_.solveODE
        (
            t0,
            dt,
            Vm.internalField(),    // Vm_new
            Iion,
            states
        );
    }

}

// HIGH ORDER //

