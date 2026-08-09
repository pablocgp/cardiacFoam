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
    testHighOrderScalarGrad

Description
    Test the high order scalar gradient

Authors
   Pablo Castrillo, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Field.H"
#include "volFields.H"

// HIGH ORDER //
// Modified for cardiacFoam: LRE replaced by highOrderInterp, the adapter over
// solids4foam's parallel movingLeastSquares. hofvc.H dropped: it was included
// but never used, and it now requires a registered solidModel.
#include "highOrderInterp.H"
// HIGH ORDER //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "analyticalFunctions.H"

    #include "gradFromFvSchemes.H"

    // HIGH ORDER //
    Info << "High Order Computing numerical gradient" << endl;
    boolList includePatchInStencils(mesh.boundaryMesh().size(), false);
    forAll(includePatchInStencils, patchI)
    {
        if
        (
            isA<fixedValueFvPatchScalarField>
            (
                volTestFun.boundaryField()[patchI]
            )
        )
        {
            includePatchInStencils[patchI] = true;
        }
    }

    #include "highOrderInterpP1.H"

    // Added for cardiacFoam: guard on the quadrature weight convention.
    //
    // LRE normalised its quadrature weights, so each face summed to 1 and each
    // cell summed to 1, and callers multiplied by |Sf| or by the cell volume
    // themselves. fvMeshQuadrature returns PHYSICAL weights that already carry
    // the area/volume. Getting this wrong rescales a flux by |Sf| without
    // changing its convergence rate, which is exactly the kind of error that
    // survives a convergence study, so assert the convention here rather than
    // trusting it.
    {
        const CompactListList<scalar>& fw =
            LREInterp_p1.faceQuadWeightPhysical();
        const CompactListList<scalar>& cw =
            LREInterp_p1.cellQuadWeightPhysical();

        // magSf() is a surfaceScalarField: its internal field only covers the
        // internal faces, so boundary faces must go through the boundary field.
        // faceQuadPoints/faceQuadWeights, by contrast, are indexed by global
        // face ID over all mesh.nFaces().
        scalar maxFaceErr = 0.0;

        const surfaceScalarField& magSf = mesh.magSf();

        auto checkFace = [&](const label faceI, const scalar area)
        {
            if (fw[faceI].empty() || area < SMALL)
            {
                return;   // empty patches carry no quadrature points
            }
            scalar sum = 0.0;
            forAll(fw[faceI], qpI)
            {
                sum += fw[faceI][qpI];
            }
            maxFaceErr = max(maxFaceErr, mag(sum - area)/area);
        };

        for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
        {
            checkFace(faceI, magSf[faceI]);
        }

        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchI];
            const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchI];

            forAll(pp, i)
            {
                checkFace(pp.start() + i, pMagSf[i]);
            }
        }

        scalar maxCellErr = 0.0;
        forAll(cw, cellI)
        {
            scalar sum = 0.0;
            forAll(cw[cellI], qpI)
            {
                sum += cw[cellI][qpI];
            }
            maxCellErr = max(maxCellErr, mag(sum - mesh.V()[cellI])/mesh.V()[cellI]);
        }

        reduce(maxFaceErr, maxOp<scalar>());
        reduce(maxCellErr, maxOp<scalar>());

        Info<< "Quadrature weight convention (physical, not normalised):" << nl
            << "    max |sum(w_face) - magSf|/magSf = " << maxFaceErr << nl
            << "    max |sum(w_cell) - V|/V         = " << maxCellErr << endl;

        if (maxFaceErr > 1e-10 || maxCellErr > 1e-10)
        {
            FatalErrorInFunction
                << "Quadrature weights are not the physical weights this code"
                << " assumes." << nl
                << "    face error = " << maxFaceErr
                << ", cell error = " << maxCellErr
                << abort(FatalError);
        }
    }

    #include "highOrderInterpP2.H"
    #include "highOrderInterpP3.H"

    Info << "End\n" << endl;
    // HIGH ORDER //

    runTime++;
    runTime.write();

    return 0;
}
