Info << "Creating analytic cell center scalar field volTestFun" << endl;

double coefPI = Foam::constant::mathematical::pi;

volScalarField volTestFun
(
    IOobject
    (
        "volTestFun",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("volTestFun", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

const vectorField& C = mesh.C();

forAll(volTestFun, cellI)
{
    const double x = C[cellI].x();
    const double y = C[cellI].y();
    const double z = C[cellI].z();

    volTestFun[cellI] = Foam::cos(coefPI*x)*Foam::cos(coefPI*y)*Foam::cos(coefPI*z);
}

volTestFun.correctBoundaryConditions();

Info << "Computing exact gradient of volTestFun" << endl;

volVectorField gradVolTestFun_exact
(
    IOobject
    (
        "gradVolTestFun_exact",
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

forAll(gradVolTestFun_exact, cellI)
{
    const double x = C[cellI].x();
    const double y = C[cellI].y();
    const double z = C[cellI].z();

    gradVolTestFun_exact[cellI] =
        vector
        (
            -coefPI*Foam::sin(coefPI*x)*Foam::cos(coefPI*y)*Foam::cos(coefPI*z),
            -coefPI*Foam::cos(coefPI*x)*Foam::sin(coefPI*y)*Foam::cos(coefPI*z),
            -coefPI*Foam::cos(coefPI*x)*Foam::cos(coefPI*y)*Foam::sin(coefPI*z)
        );
}

gradVolTestFun_exact.correctBoundaryConditions();

