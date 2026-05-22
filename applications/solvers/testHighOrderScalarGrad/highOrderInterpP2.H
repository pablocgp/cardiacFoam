dictionary LRECoeffs_p2;

LRECoeffs_p2.add("N", 2);
LRECoeffs_p2.add("Nn", 50);
LRECoeffs_p2.add("weightFunction", word("radiallySymmetricExponential"));
LRECoeffs_p2.add("k", 6);
LRECoeffs_p2.add("maxStencilSize", 80);
LRECoeffs_p2.add("useQRDecomposition", true);
LRECoeffs_p2.add("useGlobalStencils", true);
LRECoeffs_p2.add("calcConditionNumber", false);
// Declare zeroGradient boundaries as mirror patches for full high-order
// convergence at the boundary (replaces the O(1)/clipping fallback).
LRECoeffs_p2.add("mirrorPatches", wordList(1, word("patch0")));

LRE LREInterp_p2(
    mesh,
    includePatchInStencils,
    LRECoeffs_p2
);

volVectorField gradVolTestFun_num_HO_p2
(
    IOobject
    (
        "gradVolTestFun_num_HO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    LREInterp_p2.grad(volTestFun)
);

volVectorField error_HO_p2_vol
(
    IOobject
    (
        "error_HO_p2_vol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gradVolTestFun_num_HO_p2 - gradVolTestFun_exact
);

// Normas del error
double L2_HO_p2 =
    Foam::sqrt
    (
        gSum(magSqr(error_HO_p2_vol) * mesh.V())
        / gSum(mesh.V())
    );

volScalarField magError_HO_p2_vol
(
    IOobject
    (
        "magError_HO_p2_vol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(error_HO_p2_vol)
);

double Linf_HO_p2 = gMax(magError_HO_p2_vol);

autoPtr<List<List<vector>>> gradQuadPtr_p2 = LREInterp_p2.gradScalarFaceQuad(volTestFun);
const List<List<vector>>& gradQuad_p2 = gradQuadPtr_p2();
const CompactListList<point>& faceQuadPts_p2 = LREInterp_p2.faceQuadPoints();
const CompactListList<scalar>& faceQuadW_p2 =  LREInterp_p2.faceQuadWeight();

scalar L2_srf_HO_p2 = 0.0;
scalar Linf_srf_HO_p2 = 0.0;
scalar wSum_p2 = 0.0;

forAll(mesh.owner(), faceI)
{
    forAll(faceQuadPts_p2[faceI], pI)
    {
        const point& p = faceQuadPts_p2[faceI][pI];

        const scalar x = p.x();
        const scalar y = p.y();
        const scalar z = p.z();

        // gradiente analítico exacto en punto de gauss pI de la cara faceI
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

        const vector& gradNum = gradQuad_p2[faceI][pI];
        const vector err = gradNum - gradEx;
        const scalar errMag = mag(err);

        // normas
        L2_srf_HO_p2   += faceQuadW_p2[faceI][pI] * Foam::sqr(errMag);
        Linf_srf_HO_p2  = max(Linf_srf_HO_p2, errMag);
        wSum_p2    += faceQuadW_p2[faceI][pI];
    }
}

L2_srf_HO_p2 = Foam::sqrt(L2_srf_HO_p2 / wSum_p2);

Info << "HO p2: "
        << "N = " << numCells
        << "   h = " << h
        << "   L2_HO_p2 = " << L2_HO_p2
        << "   Linf_HO_p2 = " << Linf_HO_p2
        << "   L2_srf_HO_p2 = " << L2_srf_HO_p2
        << "   Linf_srf_HO_p2 = " << Linf_srf_HO_p2
        << endl;