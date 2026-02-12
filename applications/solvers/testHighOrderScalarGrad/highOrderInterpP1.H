dictionary LRECoeffs_p1;

LRECoeffs_p1.add("N", 1);
LRECoeffs_p1.add("Nn", 40);
LRECoeffs_p1.add("weightFunction", word("radiallySymmetricExponential"));
LRECoeffs_p1.add("k", 6);
LRECoeffs_p1.add("maxStencilSize", 80);
LRECoeffs_p1.add("useQRDecomposition", true);
LRECoeffs_p1.add("useGlobalStencils", true);
LRECoeffs_p1.add("calcConditionNumber", false);

LRE LREInterp_p1(
    mesh,
    includePatchInStencils,
    LRECoeffs_p1
);

volVectorField gradTestFun_num_HO_p1
(
    IOobject
    (
        "gradTestFun_num_HO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    LREInterp_p1.grad(TestFun)
);

volVectorField error_HO_p1
(
    IOobject
    (
        "gradError_HO_p1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gradTestFun_num_HO_p1 - gradTestFun_exact
);

// Normas del error
double L2_HO_p1 =
    Foam::sqrt
    (
        gSum(magSqr(error_HO_p1) * mesh.V())
        / gSum(mesh.V())
    );

volScalarField magError_HO_p1
(
    IOobject
    (
        "magError_HO_p1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(error_HO_p1)
);

double Linf_HO_p1 = gMax(magError_HO_p1);

autoPtr<List<List<vector>>> gradQuadPtr_p1 = LREInterp_p1.gradScalarFaceQuad(TestFun);
const List<List<vector>>& gradQuad_p1 = gradQuadPtr_p1();
const CompactListList<point>& faceQuadPts_p1 = LREInterp_p1.faceQuadPoints();
const CompactListList<scalar>& faceQuadW_p1 =  LREInterp_p1.faceQuadWeight();

// Info << "\nfaceQuadPts_p1:" << endl;
// Info << faceQuadPts_p1 << endl;

// forAll(faceQuadPts_p1, faceI)
// {
//     Info << "Face " << faceI
//          << " has " << faceQuadPts_p1[faceI].size()
//          << " quadrature points" << nl;
//     Info << "Face " << faceI
//          << " has " << faceQuadW_p1[faceI].size()
//          << " quadrature weights" << nl;
//     Info << "Face " << faceI
//          << " has " << gradQuad_p1[faceI].size()
//          << " gradient points" << nl;

//     forAll(faceQuadPts_p1[faceI], pI)
//     {
//         const point& p = faceQuadPts_p1[faceI][pI];
//         Info << "    p[" << pI << "] = "
//              << p << nl;
//         Info << "    w[" << pI << "] = "
//              << faceQuadW_p1[faceI][pI] << nl;
//         Info << "    grad[" << pI << "] = "
//              << gradQuad_p1[faceI][pI] << nl;
//     }
// }

// Info << "Number of faces with (1) points, (2) weights and (3) gradQuad: "
//      << faceQuadPts_p1.size() << " " << faceQuadW_p1.size() << " " << gradQuad_p1.size() << nl << endl;

// const List<CompactListList<vector>>& pointQRGradCoeffs =
//         LREInterp_p1.QRGradFaceGPCoeffs();

// const List<labelList>& stencils = LREInterp_p1.globalFaceStencils();

// =====================
// Interior faces
// =====================
// forAll(mesh.owner(), faceI)
// {
//     const labelList& curStencil = stencils[faceI];
//     const List<point>& fGP = quadPts[faceI];
//     Info << "faceI: " << faceI << " quadPts[faceI].size() " << quadPts[faceI].size() << endl;
//     forAll(quadPts[faceI], pointI)
//     {
//         Info << "faceI: " << faceI << " pointI: " << pointI << endl;
//         forAll(curStencil, cI)
//         {
//             Info << "faceI: " << faceI << " pointI: " << pointI << " cI: " << cI << endl;
//             const label neiGlobalCellI = curStencil[cI];
            
//             Info << "VI[neiGlobalCellI]: " << pointQRGradCoeffs[faceI][pointI][cI] << " " << TestFun[neiGlobalCellI] << endl;
//         }
//     }
// }

scalar L2_srf_HO_p1 = 0.0;
scalar Linf_srf_HO_p1 = 0.0;
scalar wSum_p1 = 0.0;

forAll(mesh.owner(), faceI)
{
    forAll(faceQuadPts_p1[faceI], pI)
    {
        const point& p = faceQuadPts_p1[faceI][pI];

        const scalar x = p.x();
        const scalar y = p.y();
        const scalar z = p.z();

        // gradiente analítico exacto
        vector gradEx
        (
           -coefPI
            * Foam::sin(coefPI*x)
            * Foam::cos(coefPI*y)
            * Foam::cos(coefPI*z),

           -coefPI
            * Foam::cos(coefPI*x)
            * Foam::sin(coefPI*y)
            * Foam::cos(coefPI*z),

           -coefPI
            * Foam::cos(coefPI*x)
            * Foam::cos(coefPI*y)
            * Foam::sin(coefPI*z)
        );

        const vector& gradNum = gradQuad_p1[faceI][pI];
        const vector err = gradNum - gradEx;
        const scalar errMag = mag(err);

        // normas
        L2_srf_HO_p1   += faceQuadW_p1[faceI][pI] * Foam::sqr(errMag);
        Linf_srf_HO_p1  = max(Linf_srf_HO_p1, errMag);
        wSum_p1    += faceQuadW_p1[faceI][pI];
    }
}

L2_srf_HO_p1 = Foam::sqrt(L2_srf_HO_p1 / wSum_p1);
    
Info << "HO p1: "
        << "N = " << numCells
        << "   h = " << h
        << "   L2_HO_p1 = " << L2_HO_p1
        << "   Linf_HO_p1 = " << Linf_HO_p1
        << "   L2_srf_HO_p1 = " << L2_srf_HO_p1
        << "   Linf_srf_HO_p1 = " << Linf_srf_HO_p1
        << endl;
        