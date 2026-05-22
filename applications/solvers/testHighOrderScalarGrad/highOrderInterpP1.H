dictionary LRECoeffs_p1;

LRECoeffs_p1.add("N", 1);
LRECoeffs_p1.add("Nn", 40);
LRECoeffs_p1.add("weightFunction", word("radiallySymmetricExponential"));
LRECoeffs_p1.add("k", 6);
LRECoeffs_p1.add("maxStencilSize", 80);
LRECoeffs_p1.add("useQRDecomposition", true);
LRECoeffs_p1.add("useGlobalStencils", true);
LRECoeffs_p1.add("calcConditionNumber", false);
// Declare zeroGradient boundaries as mirror patches for full high-order
// convergence at the boundary (replaces the O(1)/clipping fallback).
LRECoeffs_p1.add("mirrorPatches", wordList(1, word("patch0")));

LRE LREInterp_p1(
    mesh,
    includePatchInStencils,
    LRECoeffs_p1
);

volVectorField gradVolTestFun_num_HO_p1
(
    IOobject
    (
        "gradVolTestFun_num_HO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    LREInterp_p1.grad(volTestFun)
);

volVectorField error_HO_p1_vol
(
    IOobject
    (
        "error_HO_p1_vol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gradVolTestFun_num_HO_p1 - gradVolTestFun_exact
);

// Normas del error
double L2_HO_p1 =
    Foam::sqrt
    (
        gSum(magSqr(error_HO_p1_vol) * mesh.V())
        / gSum(mesh.V())
    );

volScalarField magError_HO_p1_vol
(
    IOobject
    (
        "magError_HO_p1_vol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(error_HO_p1_vol)
);

double Linf_HO_p1 = gMax(magError_HO_p1_vol);

autoPtr<List<List<vector>>> gradQuadPtr_p1 = LREInterp_p1.gradScalarFaceQuad(volTestFun);
const List<List<vector>>& gradQuad_p1 = gradQuadPtr_p1();
const CompactListList<point>& faceQuadPts_p1 = LREInterp_p1.faceQuadPoints();
const CompactListList<scalar>& faceQuadW_p1 =  LREInterp_p1.faceQuadWeight();

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