/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "numericFlux.H"
#include "balancedPotentialLimiter.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Flux>
Foam::numericFlux<Flux, Foam::balancedPotentialLimiter>::numericFlux
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    const volScalarField& potential,
    basicThermo& thermo
)
:
    numericFluxBase<Flux>(p.mesh()),
    p_(p),
    U_(U),
    T_(T),
    potential_(potential),
    thermo_(thermo),
    rhoFlux_
    (
        IOobject
        (
            "phi",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(thermo_.rho()*U_) & this->mesh().Sf())
    ),
    rhoUFlux_
    (
        IOobject
        (
            "rhoUFlux",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux_*linearInterpolate(U_)
    ),
    rhoEFlux_
    (
        IOobject
        (
            "rhoEFlux",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux_*linearInterpolate(thermo.Cv()*T_ + 0.5*magSqr(U_))
    ),
    rhoUSource_
    (
        IOobject
        (
            "rhoUSource",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        -thermo_.rho()*fvc::grad(potential)
    ),
    volumeInverse_( 1.0/this->mesh().V() ),
    tol_
    (
        readScalar
        (
            p.mesh().solutionDict().subDict("solvers").subDict("rho").lookup("tolerance")
        )
    )
{
    Info << "Constructing balancedNumericFlux" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Flux>
void Foam::numericFlux<Flux, Foam::balancedPotentialLimiter>::computeFlux()
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const unallocLabelList& owner = this->mesh().owner();
    const unallocLabelList& neighbour = this->mesh().neighbour();

    // Get the face area vector and cell volumes
    const surfaceVectorField& Sf = this->mesh().Sf();
    const surfaceScalarField& magSf = this->mesh().magSf();

    // Thermodynamics
    const volScalarField Cv  = thermo_.Cv();
    const volScalarField R   = thermo_.Cp() - Cv;
    const volScalarField rho = thermo_.rho();

    // Interpolate the potential to the cell faces
    const surfaceScalarField potentialFace = linearInterpolate(potential_);

    // Compute adiabatic index, momentum squared, Bernoulli, and EOS constants
    const volScalarField gamma = R/Cv + 1.0;
    const volScalarField m2 = magSqr(rho*U_);
    const volScalarField B = 0.5*magSqr(U_) + gamma/(gamma-1.0)*R*T_ + potential_;
    const volScalarField K = p_ / pow(rho, gamma);

    // Reset momentum source back to zero
    rhoUSource_ *= scalar(0.0);

    // Calculate fluxes at internal faces
    forAll (owner, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        scalar pLeft, TLeft, pRight, TRight, rhoUSourceLeft, rhoUSourceRight;
        vector ULeft, URight;

        computePrimitives
        (
            pLeft,
            ULeft,
            TLeft,
            rhoUSourceLeft,
            U_[own],
            R[own],
            B[own],
            K[own],
            m2[own],
            rho[own],
            gamma[own],
            potential_[own],
            potentialFace[faceI]
        );

        computePrimitives
        (
            pRight,
            URight,
            TRight,
            rhoUSourceRight,
            U_[nei],
            R[nei],
            B[nei],
            K[nei],
            m2[nei],
            rho[nei],
            gamma[nei],
            potential_[nei],
            potentialFace[faceI]
        );

        // Calculate local momemtum source contributions
        rhoUSource_[own] += Sf[faceI]*rhoUSourceLeft;
        rhoUSource_[nei] -= Sf[faceI]*rhoUSourceRight;

        // Calculate fluxes with reconstructed primitive variables at faces
        Flux::evaluateFlux
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
            pLeft,
            pRight,
            ULeft,
            URight,
            TLeft,
            TRight,
            R[own],
            R[nei],
            Cv[own],
            Cv[nei],
            Sf[faceI],
            magSf[faceI]
        );
    }

    // Update boundary field and values
    forAll (rhoFlux_.boundaryField(), patchi)
    {
        // Fluxes
        fvsPatchScalarField& pRhoFlux  = rhoFlux_.boundaryField()[patchi];
        fvsPatchVectorField& pRhoUFlux = rhoUFlux_.boundaryField()[patchi];
        fvsPatchScalarField& pRhoEFlux = rhoEFlux_.boundaryField()[patchi];

        // Patch fields
        const fvPatchScalarField& pp = p_.boundaryField()[patchi];
        const vectorField& pU = U_.boundaryField()[patchi];
        const scalarField& pT = T_.boundaryField()[patchi];

        const scalarField& pCv = Cv.boundaryField()[patchi];
        const scalarField& pR = R.boundaryField()[patchi];

        // Face areas
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        const unallocLabelList& pFaceCells =
            this->mesh().boundary()[patchi].faceCells();

        // Cell face potential
        const fvsPatchScalarField& ppotentialFace =
            potentialFace.boundaryField()[patchi];

        if (pp.coupled())
        {
            // Coupled patch
            const scalarField ppLeft  =
                p_.boundaryField()[patchi].patchInternalField();
            const scalarField ppRight =
                p_.boundaryField()[patchi].patchNeighbourField();
            const vectorField pULeft  =
                U_.boundaryField()[patchi].patchInternalField();
            const vectorField pURight =
                U_.boundaryField()[patchi].patchNeighbourField();
            const scalarField pTLeft  =
                T_.boundaryField()[patchi].patchInternalField();
            const scalarField pTRight =
                T_.boundaryField()[patchi].patchNeighbourField();
            const scalarField prhoLeft  =
                rho.boundaryField()[patchi].patchInternalField();
            const scalarField prhoRight =
                rho.boundaryField()[patchi].patchNeighbourField();
            const scalarField pgammaLeft  =
                gamma.boundaryField()[patchi].patchInternalField();
            const scalarField pgammaRight =
                gamma.boundaryField()[patchi].patchNeighbourField();
            const scalarField pBLeft  =
                B.boundaryField()[patchi].patchInternalField();
            const scalarField pBRight =
                B.boundaryField()[patchi].patchNeighbourField();
            const scalarField pKLeft  =
                K.boundaryField()[patchi].patchInternalField();
            const scalarField pKRight =
                K.boundaryField()[patchi].patchNeighbourField();
            const scalarField pm2Left  =
                m2.boundaryField()[patchi].patchInternalField();
            const scalarField pm2Right =
                m2.boundaryField()[patchi].patchNeighbourField();
            const scalarField ppotentialLeft  =
                potential_.boundaryField()[patchi].patchInternalField();
            const scalarField ppotentialRight =
                potential_.boundaryField()[patchi].patchNeighbourField();

            forAll (pp, facei)
            {
                scalar pLeft, TLeft, pRight, TRight, rhoUSourceLocal;
                vector ULeft, URight;

                computePrimitives
                (
                    pLeft,
                    ULeft,
                    TLeft,
                    rhoUSourceLocal,
                    pULeft[facei],
                    pR[facei],
                    pBLeft[facei],
                    pKLeft[facei],
                    pm2Left[facei],
                    prhoLeft[facei],
                    pgammaLeft[facei],
                    ppotentialLeft[facei],
                    ppotentialFace[facei]
                );

                rhoUSource_[pFaceCells[facei]] += pSf[facei]*rhoUSourceLocal;

                computePrimitives
                (
                    pRight,
                    URight,
                    TRight,
                    rhoUSourceLocal,
                    pURight[facei],
                    pR[facei],
                    pBRight[facei],
                    pKRight[facei],
                    pm2Right[facei],
                    prhoRight[facei],
                    pgammaRight[facei],
                    ppotentialRight[facei],
                    ppotentialFace[facei]
                );

                Flux::evaluateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    pLeft,
                    pRight,
                    ULeft,
                    URight,
                    TLeft,
                    TRight,
                    pR[facei],
                    pR[facei],
                    pCv[facei],
                    pCv[facei],
                    pSf[facei],
                    pMagSf[facei]
                );
            }
        }
        else // if (!pp.coupled())
        {
            forAll (pp, facei)
            {
                // Calculate fluxes
                Flux::evaluateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    pp[facei],
                    pp[facei],
                    pU[facei],
                    pU[facei],
                    pT[facei],
                    pT[facei],
                    pR[facei],
                    pR[facei],
                    pCv[facei],
                    pCv[facei],
                    pSf[facei],
                    pMagSf[facei]
                );

                scalar rhoLocal = pp[facei]/(pR[facei]*pT[facei]);
                rhoUSource_[pFaceCells[facei]] +=
                    pSf[facei] * (rhoLocal*magSqr(pU[facei]) + pp[facei]);
            }
        }
    }

    Field<vector>& irhoUSource = rhoUSource_;
    irhoUSource *= volumeInverse_;
    rhoUSource_.correctBoundaryConditions();
}


template<class Flux>
void Foam::numericFlux<Flux, Foam::balancedPotentialLimiter>::computePrimitives
(
    scalar& pFace,
    vector& UFace,
    scalar& TFace,
    scalar& rhoUSource,
    const vector& U,
    const scalar& R,
    const scalar& B,
    const scalar& K,
    const scalar& m2,
    const scalar& rho,
    const scalar& gamma,
    const scalar& potential,
    const scalar& potentialFace
)
{
    // Initalize rhoFace with value at the cell centre
    scalar rhoPrev, rhoFace = rho;

    // Use Newton's method to solve for rhoFace at the cell face
    unsigned int iter = 0;
    const scalar tolerance = rho*tol_;
    do
    {
        rhoPrev = rhoFace;
        rhoFace -= ( 0.5*m2/sqr(rhoFace) + gamma/(gamma-1.0)*K*pow(rhoFace, gamma-1.0) +
                 potentialFace - B )
             / ( -m2/pow3(rhoFace) + gamma*K*pow(rhoFace, gamma-2.0) );
        iter++;
    }
    while ( ( mag(rhoFace-rhoPrev) > tolerance ) && ( iter < 10000 ) );

    if ( iter == 10000 )
    {
        Info << "maxIter (10000) reached computing face primitives.\n"
             << "rhoFace = " << rhoFace << ",    rhoPrev = " << rhoPrev
             << ",    residual = " << mag(rhoFace-rhoPrev) << endl;
    }

    // Use new value of rhoFace to compute other primitives at the cell face
    pFace = K*pow(rhoFace, gamma);
    UFace = U*rho/rhoFace;
    TFace = pFace/(rhoFace*R);
    rhoUSource = m2/rhoFace + pFace;
}


// ************************************************************************ //
