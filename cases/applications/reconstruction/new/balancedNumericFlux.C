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

    // Get the face area vector
    const surfaceVectorField& Sf = this->mesh().Sf();
    const surfaceScalarField& magSf = this->mesh().magSf();

    // Thermodynamics
    const volScalarField Cv = thermo_.Cv();
    const volScalarField R  = thermo_.Cp() - Cv;

    // Interpolate the potential to the cell faces
    const surfaceScalarField potFace = linearInterpolate(potential_);

    // Calculate fluxes at internal faces
    forAll (owner, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        scalar pLeft, TLeft, pRight, TRight;
        vector ULeft, URight;

        computePrimitives
        (
            pLeft,
            ULeft,
            TLeft,
            p_[own],
            U_[own],
            T_[own],
            R[own],
            Cv[own],
            potential_[own],
            potFace[faceI]
        );

        computePrimitives
        (
            pRight,
            URight,
            TRight,
            p_[nei],
            U_[nei],
            T_[nei],
            R[nei],
            Cv[nei],
            potential_[nei],
            potFace[faceI]
        );

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

        // Cell face potential
        const fvsPatchScalarField& ppotFace = potFace.boundaryField()[patchi];

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
            const scalarField ppotentialLeft  =
                potential_.boundaryField()[patchi].patchInternalField();
            const scalarField ppotentialRight =
                potential_.boundaryField()[patchi].patchNeighbourField();

            forAll (pp, facei)
            {
                scalar pLeft, TLeft, pRight, TRight;
                vector ULeft, URight;
                
                computePrimitives
                (
                    pLeft,
                    ULeft,
                    TLeft,
                    ppLeft[facei],
                    pULeft[facei],
                    pTLeft[facei],
                    pR[facei],
                    pCv[facei],
                    ppotentialLeft[facei],
                    ppotFace[facei]
                );
                
                computePrimitives
                (
                    pRight,
                    URight,
                    TRight,
                    ppRight[facei],
                    pURight[facei],
                    pTRight[facei],
                    pR[facei],
                    pCv[facei],
                    ppotentialRight[facei],
                    ppotFace[facei]
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
            }
        }
    }
}


template<class Flux>
void Foam::numericFlux<Flux, Foam::balancedPotentialLimiter>::computePrimitives
(
    scalar& pFace,
    vector& UFace,
    scalar& TFace,
    const scalar p,
    const vector U,
    const scalar T,
    const scalar R,
    const scalar Cv,
    const scalar potential,
    const scalar potentialFace
)
{
    // Initalize rho with value at the cell centre
    scalar rhoPrev, rho = p / (R*T);
    
    // Compute adiabatic index and equation of state constant
    const scalar gamma = R/Cv + 1.0;
    const scalar K = p / pow(rho, gamma);

    // Compute Bernoulli constant at cell centre
    const scalar B = 0.5*magSqr(U) + gamma/(gamma-1.0)*R*T + potential;

    // Compute momentum and its square at cell centre (assumed constant over cell)
    const vector m  = rho*U;
    const scalar m2 = magSqr(m);
    
    // Use Newton's method to solve for rho at the cell face
    do
    {
        rhoPrev = rho;
        rho -= ( 0.5*m2/sqr(rho) + gamma/(gamma-1.0)*K*pow(rho, gamma-1.0) +
                 potentialFace - B )
             / ( -m2/pow3(rho) + gamma*K*pow(rho, gamma-2.0) );
    } 
    while ( abs(rho-rhoPrev) > tol_ );

    // Use the new value of rho to compute other primitives at the cell face
    pFace = K*pow(rho, gamma);
    UFace = m/rho;
    TFace = pFace/(rho*R);
}


// ************************************************************************ //
