Info << "Creating analytic scalar field TestFun" << endl;

double coefPI = Foam::constant::mathematical::pi;

volScalarField TestFun
(
    IOobject
    (
        "TestFun",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("TestFun", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

const vectorField& C = mesh.C();

forAll(TestFun, cellI)
{
    const double x = C[cellI].x();
    const double y = C[cellI].y();
    const double z = C[cellI].z();

    TestFun[cellI] = Foam::cos(coefPI*x)*Foam::cos(coefPI*y)*Foam::cos(coefPI*z);
}

TestFun.correctBoundaryConditions();

Info << "Computing exact gradient of TestFun" << endl;

volVectorField gradTestFun_exact
(
    IOobject
    (
        "gradTestFun_exact",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "zero",
        dimless/dimLength,
        vector::zero
    ),
    zeroGradientFvPatchVectorField::typeName
);

forAll(gradTestFun_exact, cellI)
{
    const double x = C[cellI].x();
    const double y = C[cellI].y();
    const double z = C[cellI].z();

    gradTestFun_exact[cellI] =
        vector
        (
            -coefPI*Foam::sin(coefPI*x)*Foam::cos(coefPI*y)*Foam::cos(coefPI*z),
            -coefPI*Foam::cos(coefPI*x)*Foam::sin(coefPI*y)*Foam::cos(coefPI*z),
            -coefPI*Foam::cos(coefPI*x)*Foam::cos(coefPI*y)*Foam::sin(coefPI*z)
        );
}

gradTestFun_exact.correctBoundaryConditions();

