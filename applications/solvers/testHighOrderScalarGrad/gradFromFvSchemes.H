Info << "Computing numerical gradient from fvSchemes" << endl;

int numCells = mesh.nCells();

volVectorField gradTestFun_num
(
    IOobject
    (
        "gradTestFun_num",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(TestFun)
);

volVectorField error
(
    IOobject
    (
        "gradError",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gradTestFun_num - gradTestFun_exact
);

// Normas del error
double L2 =
    Foam::sqrt
    (
        gSum(magSqr(error) * mesh.V())
        / gSum(mesh.V())
    );

volScalarField magError
(
    IOobject
    (
        "magError",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(error)
);

double Linf = gMax(magError);

double h = Foam::pow
(
    gSum(mesh.V()) / numCells,
    1.0/3.0
);

Info << "fvSchemes grad: "
        << "N = " << numCells
        << "   h = " << h
        << "   L2 = " << L2
        << "   Linf = " << Linf
        << endl;