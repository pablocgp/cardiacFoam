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
#include "hofvc.H"
#include "LRE.H"
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
                TestFun.boundaryField()[patchI]
            )
        )
        {
            includePatchInStencils[patchI] = true;
        }
    }

    #include "highOrderInterpP1.H"
    #include "highOrderInterpP2.H"
    #include "highOrderInterpP3.H"

    Info << "End\n" << endl;
    // HIGH ORDER //

    runTime++;
    runTime.write();

    return 0;
}
