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

#include "manufacturedSolutionHandler.H"
#include "tmanufacturedFields.H"


// Constructor
manufacturedSolutionHandler::manufacturedSolutionHandler
(
    const fvMesh& mesh,
    ionicModel& model
)
:
    mesh_(mesh),
    ionicModel_(model),
    dx_(0.0),
    dim_(0),
    N_(1)
{}


// 1. Manufactured-solution INITIALIZATION
// ===============================================================
void manufacturedSolutionHandler::initializeManufactured
(
    volScalarField& Vm,
    PtrList<volScalarField>& outFields,
    scalar dxstructured,
    int dim,
    const vectorField& pos
)
{
    Info << "Manufactured solution needs to export u1,u2,u3 always for error computation" << nl;

    // Get exported variable names from the ionic model
    Foam::wordList names = ionicModel_.exportedFieldNames();

    const label iu1 = names.find("u1");
    const label iu2 = names.find("u2");
    const label iu3 = names.find("u3");

    volScalarField& u1 = outFields[iu1];
    volScalarField& u2 = outFields[iu2];
    volScalarField& u3 = outFields[iu3];

    // Analytic initialization from tmanufacturedFields via the ionic model
    ionicModel_.initializeFields(Vm, u1, u2, u3, pos);
    dim_ = dim;     // from the mesh-based general dimension
    dx_  = dxstructured;      // the generalized dx you computed in main()
    const double totalCells =
         returnReduce(mesh_.nCells(), sumOp<int>());

    N_ = std::round(Foam::pow(totalCells, 1.0/dim_));
}

void manufacturedSolutionHandler::initializeManufacturedIntegrationPoints
(
    scalarField& VmIntegrationPoints,
    PtrList<scalarField>& outFieldsIntegrationPoints,
    scalar dxstructured,
    int dim,
    const vectorField& posIntegrationPoints
)
{
    Info << "Manufactured solution needs to export u1,u2,u3 always for error computation" << nl;

    // Get exported variable names from the ionic model
    Foam::wordList names = ionicModel_.exportedFieldNames();

    const label iu1 = names.find("u1");
    const label iu2 = names.find("u2");
    const label iu3 = names.find("u3");

    scalarField& u1IntegrationPoints = outFieldsIntegrationPoints[iu1];
    scalarField& u2IntegrationPoints = outFieldsIntegrationPoints[iu2];
    scalarField& u3IntegrationPoints = outFieldsIntegrationPoints[iu3];

    // Analytic initialization from tmanufacturedFields via the ionic model
    ionicModel_.initializeFieldsIntegrationPoints(VmIntegrationPoints, u1IntegrationPoints, 
        u2IntegrationPoints, u3IntegrationPoints, posIntegrationPoints);
    dim_ = dim;     // from the mesh-based general dimension
    dx_  = dxstructured;      // the generalized dx you computed in main()
    const double totalCells =
         returnReduce(mesh_.nCells(), sumOp<int>());

    N_ = std::round(Foam::pow(totalCells, 1.0/dim_));
}


// 3. MS POST-PROCESSING
// ===============================================================
void manufacturedSolutionHandler::postProcess
(
    const volScalarField& Vm,
    const PtrList<volScalarField>& outFields,
    const scalar dt,
    const int nsteps,
    const bool solveExplicit,
    const vectorField& pos,
    const scalar Tfinal
) const
{
    Info << "\nCalculating manufactured-solution errors..." << nl;

    // Lookup u1 and u2 from exported fields (numerical solution)
    Foam::wordList names = ionicModel_.exportedFieldNames();

    const label iu1 = names.find("u1");
    const label iu2 = names.find("u2");

    const volScalarField& u1 = outFields[iu1];
    const volScalarField& u2 = outFields[iu2];

    const scalarField x = pos.component(vector::X);
    const scalarField y = pos.component(vector::Y);
    const scalarField z = pos.component(vector::Z);
    
    // Compare numerical Vm,u1,u2 with analytic manufactured solution at volume cells
    computeAndPrintErrors
    (
        Vm,
        u1,
        u2,
        x, y, z,
        Tfinal,
        ionicModel_.tissue(),
        N_, dx_, dt, nsteps,
        solveExplicit
    );
}

void manufacturedSolutionHandler::postProcessIntegrationPoints
(
    const scalarField& VmIntegrationPoints,
    const PtrList<scalarField>& outFieldsIntegrationPoints,
    const scalar dt,
    const int nsteps,
    const bool solveExplicit,
    const vectorField& posIntegrationPoints,
    const scalar Tfinal
) const
{
    Info << "\nCalculating manufactured-solution errors..." << nl;

    // Lookup u1 and u2 from exported fields (numerical solution)
    Foam::wordList names = ionicModel_.exportedFieldNames();

    const label iu1 = names.find("u1");
    const label iu2 = names.find("u2");

    const scalarField& u1IntegrationPoints = outFieldsIntegrationPoints[iu1];
    const scalarField& u2IntegrationPoints = outFieldsIntegrationPoints[iu2];

    const scalarField x = posIntegrationPoints.component(vector::X);
    const scalarField y = posIntegrationPoints.component(vector::Y);
    const scalarField z = posIntegrationPoints.component(vector::Z);
    
    // Compare numerical Vm,u1,u2 with analytic manufactured solution at Integration Points
    computeAndPrintErrors
    (
        VmIntegrationPoints,
        u1IntegrationPoints,
        u2IntegrationPoints,
        x, y, z,
        Tfinal,
        ionicModel_.tissue(),
        N_, dx_, dt, nsteps,
        solveExplicit
    );
}
