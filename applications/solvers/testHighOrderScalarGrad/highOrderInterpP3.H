dictionary LRECoeffs_p3;

LRECoeffs_p3.add("N", 3);
LRECoeffs_p3.add("Nn", 60);
LRECoeffs_p3.add("weightFunction", word("radiallySymmetricExponential"));
LRECoeffs_p3.add("k", 6);
LRECoeffs_p3.add("maxStencilSize", 80);
LRECoeffs_p3.add("useQRDecomposition", true);
LRECoeffs_p3.add("useGlobalStencils", true);
LRECoeffs_p3.add("calcConditionNumber", false);

LRE LREInterp_p3(
    mesh,
    includePatchInStencils,
    LRECoeffs_p3
);

volVectorField gradTestFun_num_HO_p3
(
    IOobject
    (
        "gradTestFun_num_HO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    LREInterp_p3.grad(TestFun)
);

volVectorField error_HO_p3
(
    IOobject
    (
        "gradError_HO_p3",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gradTestFun_num_HO_p3 - gradTestFun_exact
);

// Normas del error
double L2_HO_p3 =
    Foam::sqrt
    (
        gSum(magSqr(error_HO_p3) * mesh.V())
        / gSum(mesh.V())
    );

volScalarField magError_HO_p3
(
    IOobject
    (
        "magError_HO_p3",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(error_HO_p3)
);

double Linf_HO_p3 = gMax(magError_HO_p3);

autoPtr<List<List<vector>>> gradQuadPtr_p3 = LREInterp_p3.gradScalarFaceQuad(TestFun);
const List<List<vector>>& gradQuad_p3 = gradQuadPtr_p3();
const CompactListList<point>& faceQuadPts_p3 = LREInterp_p3.faceQuadPoints();
const CompactListList<scalar>& faceQuadW_p3 =  LREInterp_p3.faceQuadWeight();

scalar L2_srf_HO_p3 = 0.0;
scalar Linf_srf_HO_p3 = 0.0;
scalar wSum_p3 = 0.0;

forAll(mesh.owner(), faceI)
{
    forAll(faceQuadPts_p3[faceI], pI)
    {
        const point& p = faceQuadPts_p3[faceI][pI];

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

        const vector& gradNum = gradQuad_p3[faceI][pI];
        const vector err = gradNum - gradEx;
        const scalar errMag = mag(err);

        // normas
        L2_srf_HO_p3   += faceQuadW_p3[faceI][pI] * Foam::sqr(errMag);
        Linf_srf_HO_p3  = max(Linf_srf_HO_p3, errMag);
        wSum_p3    += faceQuadW_p3[faceI][pI];
    }
}

L2_srf_HO_p3 = Foam::sqrt(L2_srf_HO_p3 / wSum_p3);

Info << "HO p3: "
        << "N = " << numCells
        << "   h = " << h
        << "   L2_HO_p3 = " << L2_HO_p3
        << "   Linf_HO_p3 = " << Linf_HO_p3
        << "   L2_srf_HO_p3 = " << L2_srf_HO_p3
        << "   Linf_srf_HO_p3 = " << Linf_srf_HO_p3
        << endl;