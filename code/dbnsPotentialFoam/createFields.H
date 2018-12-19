    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> thermo
    (
        basicPsiThermo::New(mesh)
    );

    // Primitive variables

    volScalarField& h = thermo->h();
    volScalarField& p = thermo->p();
    const volScalarField& T = thermo->T();

    Info<< "\nReading field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo->rho()
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading potential field\n" << endl;
    volScalarField potential
    (
        IOobject
        (
            "potential",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Conservative variables

    volVectorField rhoU
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*U
    );

    volScalarField rhoE
    (
        IOobject
        (
            "rhoE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*(h + 0.5*magSqr(U)) - p
    );

    const bool isBalancedFlux = readBool
    (
        mesh.schemesDict().subDict("divSchemes").subDict("dbns")
            .lookup("balanced")
    );

    autoPtr<basicNumericFlux> dbnsFluxPtr;
    basicBalancedNumericFlux* balancedDbnsFluxPtr = NULL;

    if ( isBalancedFlux )
    {
        // Create balanced numeric flux
        balancedDbnsFluxPtr =
        (
            basicBalancedNumericFlux::New
            (
                potential,
                U,
                rho,
                thermo()
            )
        ).ptr();
        dbnsFluxPtr.set(dynamic_cast<basicNumericFlux*>(balancedDbnsFluxPtr));
    }
    else
    {
        // Create standard numeric flux
        dbnsFluxPtr = basicNumericFlux::New
        (
            p,
            U,
            T,
            thermo()
        );
    }

    basicNumericFlux& dbnsFlux = dbnsFluxPtr();

    // Create mass flux alias for easier coupling with other code components
    const surfaceScalarField& phi = dbnsFlux.rhoFlux();
