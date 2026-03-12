Info << "Computing numerical gradient from fvSchemes" << endl;

int numCells = mesh.nCells();

volVectorField gradVolTestFun_num
(
    IOobject
    (
        "gradVolTestFun_num",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(volTestFun)
);

volVectorField gradVolError
(
    IOobject
    (
        "gradVolError",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gradVolTestFun_num - gradVolTestFun_exact
);

// Normas del error
double L2 =
    Foam::sqrt
    (
        gSum(magSqr(gradVolError) * mesh.V())
        / gSum(mesh.V())
    );

volScalarField magGradVolError
(
    IOobject
    (
        "magGradVolError",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(gradVolError)
);

double Linf = gMax(magGradVolError);

double h = Foam::pow
(
    gSum(mesh.V()) / numCells,
    1.0/3.0
);

// gradiente estándar según fvSchemes
volVectorField gradV = fvc::grad(volTestFun);

// interpolación a caras
surfaceVectorField gradF = fvc::interpolate(gradV);

const vectorField& Cf = mesh.Cf();

scalar L2_srf = 0.0;
scalar Linf_srf = 0.0;
scalar wSum = 0.0;

forAll(mesh.owner(), faceI)
{
    const point& p = Cf[faceI];

    const scalar x = p.x();
    const scalar y = p.y();
    const scalar z = p.z();

    // gradiente analítico exacto en el centro de cara
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

    const vector& gradNum = gradF[faceI];

    const vector err = gradNum - gradEx;
    const scalar errMag = mag(err);

    // peso = área de la cara
    const scalar w = mesh.magSf()[faceI];

    L2_srf   += w * Foam::sqr(errMag);
    Linf_srf  = max(Linf_srf, errMag);
    wSum     += w;
}

L2_srf = Foam::sqrt(L2_srf / wSum);

Info << "fvSchemes grad: "
        << "N = " << numCells
        << "   h = " << h
        << "   L2 = " << L2
        << "   Linf = " << Linf
        << "   L2_srf = " << L2_srf
        << "   Linf_srf = " << Linf_srf
        << endl;