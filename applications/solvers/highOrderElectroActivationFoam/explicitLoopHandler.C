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
void explicitLoopHandler::explicitLoop
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
    const volTensorField& conductivity
)
{
    // 0) Manufactured-specific: store "old" internal states
    if (ionicModel_.hasManufacturedSolution())
    {
        refCast<tmanufacturedFDA>(ionicModel_).updateStatesOld();
    }

    // 1) External stimulus (internal field only)
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
    ionicModel_.calculateCurrent
    (
        t0,
        dt,
        Vm.internalField(),   // Vm_old
        Iion,
        states
    );
    Iion.correctBoundaryConditions();

    if (ionicModel_.hasManufacturedSolution())
    {
        Iion /= Cm.value();
    }

    solve
    (
        chi*Cm*fvm::ddt(Vm)
     == fvc::laplacian(conductivity, Vm)
      - chi*Cm*Iion
      + externalStimulusCurrent
    );

    Vm.correctBoundaryConditions();

    // 4) Manufactured-specific: reset internal states back to OLD
    if (ionicModel_.hasManufacturedSolution())
    {
        refCast<tmanufacturedFDA>(ionicModel_).resetStatesToStatesOld();
    }

    // 5) Advance ionic model in time (ODE solve) with NEW Vm, once per timestep
    ionicModel_.solveODE
    (
        t0,
        dt,
        Vm.internalField(),    // Vm_new
        Iion,
        states
    );

}

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
    const label totalIntegrationPoints,
    const LRE& LREInterp,
    volVectorField& gradVm_HO,
    volScalarField& lapVm_HO,
    surfaceVectorField& surfaceGradVm_HO
)
{
    const fvMesh& mesh = Vm.mesh();

    // 0) Manufactured-specific: store "old" internal states
    if (ionicModel_.hasManufacturedSolution())
    {
        refCast<tmanufacturedFDA>(ionicModel_).updateStatesOld();
    }
    Info << "Fin 0)" << endl;
    // 1) External stimulus (internal field only)
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
    Info << "Fin 1)" << endl;
    // // 2) Ionic current using OLD Vm
    // ionicModel_.calculateCurrent
    // (
    //     t0,
    //     dt,
    //     Vm.internalField(),   // Vm_old
    //     Iion,
    //     states
    // );
    // Iion.correctBoundaryConditions();

    // if (ionicModel_.hasManufacturedSolution())
    // {
    //     Iion /= Cm.value();
    // }

    // gradVm_HO = LREInterp.grad(Vm);
    // gradVm_HO = fvc::grad(Vm);
    // gradVm_HO.correctBoundaryConditions();

    // surfaceVectorField gradVm_faces = fvc::interpolate(gradVm_HO);

    // surfaceScalarField flux = fvc::interpolate(conductivity & gradVm_HO) & mesh.Sf() ;// & gradVm_faces;

    // lapVm_HO = fvc::div(flux);

    const vectorField& C = mesh.C();  // centros de celda
    const CompactListList<scalar>& cellQuadW = LREInterp.cellQuadWeight();
    const CompactListList<point>& cellQuadP = LREInterp.cellQuadPoints();

    // SURFACE INTERPOLATION
    // Info << "surfaceGradVm_HO.dimensions() " << surfaceGradVm_HO.dimensions() << endl;

    scalarField VmIntegrationPoints( totalIntegrationPoints, 0.0);
    scalarField IionIntegrationPoints( totalIntegrationPoints, 0.0);

    gradVm_HO = LREInterp.grad(Vm);
    volVectorField gradIion_HO = LREInterp.grad(Iion);

    label integrationPointPos = 0;

    forAll(mesh.cells(), cellI)
    {
        // Info<< Iion[cellI] << " " << Vm[cellI] << " " << gradVm_HO[cellI] << endl;
        const scalar Vc = Vm[cellI];
        const scalar Iionc = Iion[cellI];
        const vector& gradVc = gradVm_HO[cellI];
        const vector& gradIionc = gradIion_HO[cellI];
        const vector& xc = C[cellI];
        // Field<Field<scalar>> statesGauss( cellQuadP[cellI].size(),Field<scalar>(nStates, 0.0));

        forAll(cellQuadP[cellI], gI)
        {
            vector dx = cellQuadP[cellI][gI] - xc;

            scalar Vg = Vc + (gradVc & dx);
            scalar Iiong = Iionc + (gradIionc & dx);

            VmIntegrationPoints[integrationPointPos] = Vg;
            IionIntegrationPoints[integrationPointPos] = Iiong;
            // Info << "    " << gI << "/" << cellQuadP[cellI].size() << " " << cellQuadP[cellI][gI] << " " << VmGauss[gI] << endl;
            // statesGauss[gI] = states[cellI];
            integrationPointPos++;
        }
    }
    // Info << "Llegue hasta aca" << endl;
    //         std::cin.get();
    ionicModel_.calculateCurrent
    (
        t0,
        dt,
        VmIntegrationPoints,
        IionIntegrationPoints,
        states
    );

    integrationPointPos = 0;
    forAll(mesh.cells(), cellI)
    {
        Iion[cellI] = 0.0;
        double totCellQuadW = 0.0;
        forAll(cellQuadP[cellI], gI)
        {
            Iion[cellI] += cellQuadW[cellI][gI]*IionIntegrationPoints[integrationPointPos];
            totCellQuadW += cellQuadW[cellI][gI];
            integrationPointPos++;
        }
        Iion[cellI] /= totCellQuadW; // this should be 1.. but just in case
        // Info << cellI << " " << totCellQuadW << " " << Iion[cellI] << endl;
    }
    Iion.correctBoundaryConditions();

    if (ionicModel_.hasManufacturedSolution())
    {
        Iion /= Cm.value();
    }
    Info << "Fin 2)" << endl;
    autoPtr<List<List<vector>>> gradQuadVm_ptr = LREInterp.gradScalarFaceQuad(Vm);
    List<List<vector>>& gradQuadVm = gradQuadVm_ptr.ref();
    const CompactListList<scalar>& faceQuadW = LREInterp.faceQuadWeight();

    // Info << "quadVm.size() = " << quadVm.size() << endl;
    // Info << "quadW.size() = " << quadW.size() << endl;
    // Info << "nInternalFaces = " << mesh.nInternalFaces() << endl;
    // Info << "nTotalFaces = " << mesh.nFaces() << endl;

    // Internal faces
    forAll(surfaceGradVm_HO, faceI)
    {
        const label ownerCell = mesh.owner()[faceI];
        surfaceGradVm_HO[faceI] = vector::zero;
        forAll(gradQuadVm[faceI], pI)
        {
            surfaceGradVm_HO[faceI] += (conductivity[ownerCell] & gradQuadVm[faceI][pI]) * faceQuadW[faceI][pI];
        }
    }

    surfaceGradVm_HO.correctBoundaryConditions();

    // volScalarField lapVm_HO2 = fvc::div(surfaceGradVm_HO & mesh.Sf());

    // Info << "Vm.dimensions()               " << Vm.dimensions() << endl;
    // Info << "conductivity.dimensions()     " << conductivity.dimensions() << endl;
    // Info << "surfaceGradVm_HO.dimensions() " << surfaceGradVm_HO.dimensions() << endl;
    // Info << "lapVm_HO2.dimensions()        " << lapVm_HO2.dimensions() << endl;
    // Info << "externalStimulusCurrent.dimensions() " << externalStimulusCurrent.dimensions() << endl;

    lapVm_HO = fvc::div(mesh.Sf() & surfaceGradVm_HO );

    solve
    (
        chi*Cm*fvm::ddt(Vm)
     == lapVm_HO
      - chi*Cm*Iion
      + externalStimulusCurrent
    );

    Vm.correctBoundaryConditions();
    Info << "Fin 3)" << endl;
    // 4) Manufactured-specific: reset internal states back to OLD
    if (ionicModel_.hasManufacturedSolution())
    {
        refCast<tmanufacturedFDA>(ionicModel_).resetStatesToStatesOld();
    }

    VmIntegrationPoints = 0.0;
    IionIntegrationPoints = 0.0;

    gradVm_HO = LREInterp.grad(Vm);
    gradIion_HO = LREInterp.grad(Iion);
    integrationPointPos = 0;

    forAll(mesh.cells(), cellI)
    {
        const scalar Vc = Vm[cellI];
        const scalar Iionc = Iion[cellI];
        const vector& gradVc = gradVm_HO[cellI];
        const vector& gradIionc = gradIion_HO[cellI];
        const vector& xc = C[cellI];
        
        // Field<Field<scalar>> statesGauss( cellQuadP[cellI].size(),Field<scalar>(nStates, 0.0));

        forAll(cellQuadP[cellI], gI)
        {
            vector dx = cellQuadP[cellI][gI] - xc;

            scalar Vg = Vc + (gradVc & dx);
            scalar Iiong = Iionc + (gradIionc & dx);

            VmIntegrationPoints[integrationPointPos] = Vg;
            IionIntegrationPoints[integrationPointPos] = Iiong;
            // Info << "    " << gI << "/" << cellQuadP[cellI].size() << " " << cellQuadP[cellI][gI] << " " << VmGauss[gI] << endl;
            // statesGauss[gI] = states[cellI];
            integrationPointPos++;
        }
    }
    Info << "Fin 4)" << endl;
    // 5) Advance ionic model in time (ODE solve) with NEW Vm, once per timestep
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

    integrationPointPos = 0;
    forAll(mesh.cells(), cellI)
    {
        Iion[cellI] = 0.0;
        double totCellQuadW = 0.0;
        forAll(cellQuadP[cellI], gI)
        {
            Iion[cellI] += cellQuadW[cellI][gI]*IionIntegrationPoints[integrationPointPos];
            totCellQuadW += cellQuadW[cellI][gI];
            integrationPointPos++;
        }
        Iion[cellI] /= totCellQuadW; // this should be 1.. but just in case
    }
    Iion.correctBoundaryConditions(); // correct????

    // if (ionicModel_.hasManufacturedSolution())
    // {
    //     Iion /= Cm.value();
    // }
    Info << "Fin 5)" << endl;
    // Info << "Llegue hasta aca" << endl;
    //         std::cin.get();
}

// HIGH ORDER //

