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
    const volTensorField& conductivity,
    volVectorField& gradVm_HO,
    volScalarField& lapVm_HO,
    volScalarField& lapVm_standar
)
{
    const fvMesh& mesh = Vm.mesh();
    
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

    // gradVm_HO = fvc::grad(Vm);
    // gradVm_HO.correctBoundaryConditions();

    // surfaceScalarField flux = fvc::interpolate(conductivity & gradVm_HO) & mesh.Sf() ;// & gradVm_faces;

    // lapVm_HO = fvc::div(flux);

    // lapVm_standar = fvc::laplacian(conductivity, Vm);

    // fvc::div( conductivity & HOgrad(Vm) )

    // fvc::laplacian(conductivity, Vm) = fvc::div(conductivity & fvc::grad(Vm))

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

    // gradVm_HO = LREInterp.grad(Vm);
    // gradVm_HO = fvc::grad(Vm);
    // gradVm_HO.correctBoundaryConditions();

    // surfaceVectorField gradVm_faces = fvc::interpolate(gradVm_HO);

    // surfaceScalarField flux = fvc::interpolate(conductivity & gradVm_HO) & mesh.Sf() ;// & gradVm_faces;

    // lapVm_HO = fvc::div(flux);

    // SURFACE INTERPOLATION
    // Info << "surfaceGradVm_HO.dimensions() " << surfaceGradVm_HO.dimensions() << endl;

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

    // Boundary faces
    // forAll(surfaceGradVm_HO.boundaryField(), patchI)
    // {
    //     vectorField& tfPatch = surfaceGradVm_HO.boundaryFieldRef()[patchI];

    //     forAll(tfPatch, faceI)
    //     {
    //         const label globalFaceID = mesh.boundaryMesh()[patchI].start() + faceI;
    //         const label ownerCell = mesh.owner()[globalFaceID];
    //         const List<vector>& faceQuadVm = quadVm[globalFaceID];

    //         tfPatch[faceI] = vector::zero;
    //         forAll(faceQuadVm, pI)
    //         {
    //             tfPatch[faceI] += (conductivity[ownerCell] & faceQuadVm[pI]) * quadW[globalFaceID][pI];
    //         }
    //     }
    // }

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


void explicitLoopHandler::setSpatialIntegration
(
    const Switch useHighOrder,
    const dictionary& solutionVariablesMemory
)
{
    solutionVariablesMemory_ = solutionVariablesMemory;
    useHighOrder_ = useHighOrder;

    if (useHighOrder_)
    {
        highOrderCoeffs_ = solutionVariablesMemory.subDict("highOrderCoeffs");
        LRECoeffs_       = highOrderCoeffs_.subDict("LRECoeffs");
        
        Info<< "Spatial integration: HIGH ORDER\n";

        // Info<< "  Jacobian  : "
        //     << highOrderCoeffs_.lookup("highOrderJacobian") << nl
        //     << "  Residual  : "
        //     << highOrderCoeffs_.lookup("highOrderResidual") << nl
        //     << "  LRE order : "
        //     << LRECoeffs_.lookup("N") << nl;
    }
    else
    {
        Info<< "Spatial integration: STANDARD (OpenFOAM)\n";
    }
}

// HIGH ORDER //

