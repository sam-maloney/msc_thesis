    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> thermo
    (
        basicPsiThermo::New(mesh)
    );

    // Primitive variables

    volScalarField& h = thermo->h();
    volScalarField& p = thermo->p();
    const volScalarField& T = thermo->T();

    Info<< "Reading field rho\n" << endl;
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

    // Store limiter type
    const word limiterType = word(
      mesh.schemesDict().subDict("divSchemes").subDict("dbns").lookup("limiter")
    );
    const bool isBalancedFlux = ( limiterType.substr(0,8) == "balanced" );

    autoPtr<basicNumericFlux> dbnsFluxPtr;

    if ( isBalancedFlux )
    {
        // Create balanced numeric flux
        dbnsFluxPtr = basicNumericFlux::New
        (
            potential, // p will be extracted from thermo
            U,
            T,
            thermo()
        );
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
//    const bool isBalancedFlux = ( dbnsFlux.type() == "balancedNumericFlux" );

    // Create mass flux alias for easier coupling with other code components
    const surfaceScalarField& phi = dbnsFlux.rhoFlux();
