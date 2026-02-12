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

Solver
    highOrderElectroActivationFoam

Description
    Generic model-agnostic solver for cardiac electrophysiology
    based on the monodomain reaction–diffusion equation. The solver:

    - Delegates all ionic-model state indexing and ODE evaluation to the
      run-time selectable ionicModel.

    - Supports multiple time-integration strategies (explicit, implicit)
      via dedicated loop-handler classes.

    - Provides infrastructure for manufactured-solution verification through
      model-supplied export functions, avoiding any solver-side indexing of
      ionic states.

    - Stores ionic state vectors externally (one N-state vector per cell),
      ensuring complete separation between solver logic and model detail.

Authors
   Simão Nieto de Castro, UCD.
   Philip Cardiff, UCD.
   Pablo Castrillo, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "manufacturedSolutionHandler.H" // helper class for manufactured solution
#include "explicitLoopHandler.H"         // helper class for explicit Loop
#include "implicitLoopHandler.H"         // helper class for implicit Loop
#include "ionicModel.H"
#include "pimpleControl.H"
#include "Field.H"
#include "volFields.H"

// HIGH ORDER //
#include "hofvc.H"
#include "LRE.H"
// HIGH ORDER //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // Initialisation of handler functions
    manufacturedSolutionHandler msHandler(mesh, ionicModel());
    explicitLoopHandler explicitHandler(mesh, ionicModel());
    implicitLoopHandler implicitHandler(mesh, ionicModel());

    // External state storage: one N-state vector per cell.
    const label nStates = ionicModel->nEqns();
    Field<Field<scalar>> states
    (
        mesh.nCells(),Field<scalar>(nStates, 0.0)
    );

    // Initialisation
    pimpleControl pimple(mesh);

    // Structured mesh initialization
    scalar dx  = Foam::cbrt(mesh.V().average().value());
    int dim = mesh.nGeometricD();

    scalar dt = runTime.deltaTValue();
    int nsteps = int(ceil(runTime.endTime().value()/dt));

    if (ionicModel->hasManufacturedSolution())
    {
        msHandler.initializeManufactured(Vm, outFields, dx ,dim);
    }

    // Solution methodology flag
    const Switch solveExplicit(solutionVariablesMemory.lookup("solveExplicit"));

    // Extract the names of the fields to be exported
    const wordList exportNames = ionicModel->exportedFieldNames();
    if (!exportNames.empty())
    {
        Info<< "Exporting fields: " << exportNames << nl;
    }

    // =============== CASE 1: EXPLICIT SOLVER ====================
    if (solveExplicit)
    {
        const scalar CFL = readScalar(solutionVariablesMemory.lookup("CFL"));

        explicitHandler.initializeExplicit
        (
            dt,
            nsteps,
            chi.value(),
            Cm.value(),
            conductivity,
            CFL,
            dx,
            dim
        );

        runTime.setDeltaT(dt);
        nsteps = int(std::ceil(runTime.endTime().value()/dt));

        while (runTime.loop())
        {
            const scalar t0 = runTime.value() - dt;

            if (useHighOrder)
            {   

                // gradVm_HO = LREInterp.grad(Vm);
                // gradVm_HO.correctBoundaryConditions();

                // surfaceVectorField gradVm_faces = fvc::interpolate(gradVm_HO);

                // surfaceVectorField flux = fvc::interpolate(conductivity) & gradVm_faces;

                // volScalarField lapVm_HO = fvc::div(flux & mesh.Sf());






                // SURFACE INTERPOLATION
                // autoPtr<List<List<vector>>> quadVm_ptr = LREInterp.gradScalarQuad(Vm);
                // List<List<vector>>& quadVm = quadVm_ptr.ref();
                // const CompactListList<scalar>& quadW = LREInterp.faceQuadWeight();

                // // Info << "quadVm.size() = " << quadVm.size() << endl;
                // // Info << "quadW.size() = " << quadW.size() << endl;
                // // Info << "nInternalFaces = " << mesh.nInternalFaces() << endl;
                // // Info << "nTotalFaces = " << mesh.nFaces() << endl;

                // // Internal faces
                // forAll(surfaceGradVm_HO, faceI)
                // {
                //     const List<vector>& faceQuadVm = quadVm[faceI];
                //     const label ownerCell = mesh.owner()[faceI];

                //     forAll(faceQuadVm, pI)
                //     {
                //         surfaceGradVm_HO[faceI] += (conductivity[ownerCell] & faceQuadVm[pI]) * quadW[faceI][pI];
                //     }
                // }

                // // Boundary faces
                // forAll(surfaceGradVm_HO.boundaryField(), patchI)
                // {
                //     vectorField& tfPatch = surfaceGradVm_HO.boundaryFieldRef()[patchI];

                //     forAll(tfPatch, faceI)
                //     {
                //         const label globalFaceID = mesh.boundaryMesh()[patchI].start() + faceI;
                //         const label ownerCell = mesh.owner()[globalFaceID];
                //         const List<vector>& faceQuadVm = quadVm[globalFaceID];

                //         forAll(faceQuadVm, pI)
                //         {
                //             tfPatch[faceI] += (conductivity[ownerCell] & faceQuadVm[pI]) * quadW[globalFaceID][pI];
                //         }
                //     }
                // }

                // surfaceGradVm_HO.correctBoundaryConditions();

                // surfaceVectorField flux = fvc::interpolate(conductivity) & surfaceGradVm_HO;

                // volScalarField lapVm_HO = fvc::div(flux & mesh.Sf());

                // Info << "Llegue hasta aca" << endl;
                // std::cin.get();

                explicitHandler.highOrderExplicitLoop
                (
                    t0,
                    dt,
                    Vm,
                    Iion,
                    states,
                    externalStimulusCurrent,
                    stimulusCellIDsList,
                    stimulusStartTimes,
                    stimulusIntensity.value(),
                    stimulusDuration.value(),
                    chi,
                    Cm,
                    conductivity,
                    totalIntegrationPoints,
                    LREInterp,
                    gradVm_HO,
                    lapVm_HO,
                    surfaceGradVm_HO
                );
                // Info << "Llegue hasta aca" << endl;
                // std::cin.get();
                
            }
            else
            {
                explicitHandler.explicitLoop
                (
                    t0,
                    dt,
                    Vm,
                    Iion,
                    states,
                    externalStimulusCurrent,
                    stimulusCellIDsList,
                    stimulusStartTimes,
                    stimulusIntensity.value(),
                    stimulusDuration.value(),
                    chi,
                    Cm,
                    conductivity
                );
            }

            ionicModel->exportStates(states, outFields);

            #include "updateActivationTimes.H"

            runTime.write();
        }
    }
    // =============== CASE 2: IMPLICIT SOLVER ====================
    else
    {
        Info<< "\nUsing implicit solver\n" << endl;

        while (runTime.loop())
        {
            const scalar currentTime = runTime.value();
            const scalar deltaT      = runTime.deltaTValue();
            const scalar t0          = currentTime - deltaT;

            implicitHandler.implicitLoop
            (
                t0,
                deltaT,
                Vm,
                Iion,
                states,
                externalStimulusCurrent,
                stimulusCellIDsList,
                stimulusStartTimes,
                stimulusIntensity.value(),
                stimulusDuration.value(),
                chi,
                Cm,
                conductivity,
                pimple
            );

            ionicModel->exportStates(states, outFields);

            #include "updateActivationTimes.H"

            runTime.write();
        }
    }

    // Manufactured-solution post-processing
    if (ionicModel->hasManufacturedSolution())
    {
        msHandler.postProcess(Vm, outFields, dt, nsteps, solveExplicit);
    }

    runTime.printExecutionTime(Info);

    Info<< "End" << nl << endl;

    return 0;
}
