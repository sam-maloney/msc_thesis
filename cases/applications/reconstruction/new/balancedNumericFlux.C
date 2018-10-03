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
#include "MDLimiter.H"
#include "balancedPotentialLimiter.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Flux>
Foam::numericFlux<Flux, Foam::balancedPotentialLimiter>::numericFlux
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    const volScalarField& rho,
    const volScalarField& potential,
    basicThermo& thermo
)
:
    numericFluxBase<Flux>(p.mesh()),
    p_(p),
    U_(U),
    T_(T),
    rho_(rho),
    potential_(potential),
    thermo_(thermo),
    pLeft_
    (
        IOobject
        (
            "pLeft",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(p_)
    ),
    pRight_
    (
        IOobject
        (
            "pRight",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(p_)
    ),
    ULeft_
    (
        IOobject
        (
            "ULeft",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(U_)
    ),
    URight_
    (
        IOobject
        (
            "URight",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(U_)
    ),
    TLeft_
    (
        IOobject
        (
            "TLeft",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(T_)
    ),
    TRight_
    (
        IOobject
        (
            "TRight",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(T_)
    ),
    gradP_
    (
        IOobject
        (
            "gradP",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(p_)
    ),
    gradU_
    (
        IOobject
        (
            "gradU",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(U_)
    ),
    gradT_
    (
        IOobject
        (
            "gradT",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(T_)
    ),
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

    const volVectorField& cellCentre = this->mesh().C();
    const surfaceVectorField& faceCentre = this->mesh().Cf();

    // Thermodynamics
    const volScalarField Cv  = thermo_.Cv();
    const volScalarField R   = thermo_.Cp() - Cv;

    // Interpolate the potential to the cell faces
    const surfaceScalarField potentialFace = linearInterpolate(potential_);

    // Compute adiabatic index, momentum squared, Bernoulli, and EOS constants
    const volScalarField gamma = R/Cv + 1.0;
    const volScalarField m2 = magSqr(rho_*U_);
    const volScalarField B = 0.5*magSqr(U_) + gamma/(gamma-1.0)*R*T_ + potential_;
    const volScalarField K = p_ / pow(rho_, gamma);

    // Reset momentum source and gradients back to zero
    rhoUSource_ *= scalar(0.0);
    gradP_ *= scalar(0.0);
    gradU_ *= scalar(0.0);
    gradT_ *= scalar(0.0);

    // Create a data structure for storing local contributions
    dataStruct localData;

    // Store min and max values from neighbourhood
    volScalarField pMinValue("pMinValue", p_*scalar(0.0));
    volScalarField pMaxValue("pMaxValue", p_*scalar(0.0));
    volVectorField UMinValue("UMinValue", U_*scalar(0.0));
    volVectorField UMaxValue("UMaxValue", U_*scalar(0.0));
    volScalarField TMinValue("TMinValue", T_*scalar(0.0));
    volScalarField TMaxValue("TMaxValue", T_*scalar(0.0));

    Field<scalar>& pMinIn = pMinValue.internalField();
    Field<scalar>& pMaxIn = pMaxValue.internalField();
    Field<vector>& UMinIn = UMinValue.internalField();
    Field<vector>& UMaxIn = UMaxValue.internalField();
    Field<scalar>& TMinIn = TMinValue.internalField();
    Field<scalar>& TMaxIn = TMaxValue.internalField();

    // Calculate fluxes at internal faces
    forAll (owner, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        computePrimitives
        (
            pLeft_[faceI],
            ULeft_[faceI],
            TLeft_[faceI],
            localData,
            U_[own],
            R[own],
            B[own],
            K[own],
            m2[own],
            rho_[own],
            gamma[own],
            potentialFace[faceI],
            potential_[nei]
        );

        // Calculate local momemtum source contribution
        rhoUSource_[own] += Sf[faceI]*localData.rhoUSource;
        gradP_[own] += Sf[faceI]*(p_[nei] - localData.p);
        gradU_[own] += Sf[faceI]*(U_[nei] - localData.U);
        gradT_[own] += Sf[faceI]*(T_[nei] - localData.T);

        pMinIn[own] = min(pMinIn[own], p_[nei] - localData.p);
        pMaxIn[own] = max(pMaxIn[own], p_[nei] - localData.p);
        UMinIn[own] = min(UMinIn[own], U_[nei] - localData.U);
        UMaxIn[own] = max(UMaxIn[own], U_[nei] - localData.U);
        TMinIn[own] = min(TMinIn[own], T_[nei] - localData.T);
        TMaxIn[own] = max(TMaxIn[own], T_[nei] - localData.T);

        computePrimitives
        (
            pRight_[faceI],
            URight_[faceI],
            TRight_[faceI],
            localData,
            U_[nei],
            R[nei],
            B[nei],
            K[nei],
            m2[nei],
            rho_[nei],
            gamma[nei],
            potentialFace[faceI],
            potential_[own]
        );

        rhoUSource_[nei] -= Sf[faceI]*localData.rhoUSource;
        gradP_[nei] -= Sf[faceI]*(p_[own] - localData.p);
        gradU_[nei] -= Sf[faceI]*(U_[own] - localData.U);
        gradT_[nei] -= Sf[faceI]*(T_[own] - localData.T);

        pMinIn[nei] = min(pMinIn[nei], p_[own] - localData.p);
        pMaxIn[nei] = max(pMaxIn[nei], p_[own] - localData.p);
        UMinIn[nei] = min(UMinIn[nei], U_[own] - localData.U);
        UMaxIn[nei] = max(UMaxIn[nei], U_[own] - localData.U);
        TMinIn[nei] = min(TMinIn[nei], T_[own] - localData.T);
        TMaxIn[nei] = max(TMaxIn[nei], T_[own] - localData.T);
    }

    // Update boundary field and values
    forAll (rhoFlux_.boundaryField(), patchi)
    {
        // Reconstructed values
        fvsPatchScalarField& ppLeft_  = pLeft_.boundaryField()[patchi];
        fvsPatchScalarField& ppRight_ = pRight_.boundaryField()[patchi];
        fvsPatchVectorField& pULeft_  = ULeft_.boundaryField()[patchi];
        fvsPatchVectorField& pURight_ = URight_.boundaryField()[patchi];
        fvsPatchScalarField& pTLeft_  = TLeft_.boundaryField()[patchi];
        fvsPatchScalarField& pTRight_ = TRight_.boundaryField()[patchi];

        // Patch fields
        const fvPatchScalarField& pp = p_.boundaryField()[patchi];
        const vectorField& pU = U_.boundaryField()[patchi];
        const scalarField& pT = T_.boundaryField()[patchi];
        const scalarField& pR = R.boundaryField()[patchi];

        // Face areas
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
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
                rho_.boundaryField()[patchi].patchInternalField();
            const scalarField prhoRight =
                rho_.boundaryField()[patchi].patchNeighbourField();
            const scalarField pgammaLeft  =
                gamma.boundaryField()[patchi].patchInternalField();
            const scalarField pgammaRight =
                gamma.boundaryField()[patchi].patchNeighbourField();
            const scalarField pRLeft  =
                R.boundaryField()[patchi].patchInternalField();
            const scalarField pRRight =
                R.boundaryField()[patchi].patchNeighbourField();
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
                computePrimitives
                (
                    ppLeft_[facei],
                    pULeft_[facei],
                    pTLeft_[facei],
                    localData,
                    pULeft[facei],
                    pRLeft[facei],
                    pBLeft[facei],
                    pKLeft[facei],
                    pm2Left[facei],
                    prhoLeft[facei],
                    pgammaLeft[facei],
                    ppotentialFace[facei],
                    ppotentialRight[facei]
                );

                const label own = pFaceCells[facei];

                rhoUSource_[own] += pSf[facei]*localData.rhoUSource;
                gradP_[own] += pSf[facei]*(ppRight[facei] - localData.p);
                gradU_[own] += pSf[facei]*(pURight[facei] - localData.U);
                gradT_[own] += pSf[facei]*(pTRight[facei] - localData.T);

                pMinIn[own] = min(pMinIn[own], ppRight[facei] - localData.p);
                pMaxIn[own] = max(pMaxIn[own], ppRight[facei] - localData.p);
                UMinIn[own] = min(UMinIn[own], pURight[facei] - localData.U);
                UMaxIn[own] = max(UMaxIn[own], pURight[facei] - localData.U);
                TMinIn[own] = min(TMinIn[own], pTRight[facei] - localData.T);
                TMaxIn[own] = max(TMaxIn[own], pTRight[facei] - localData.T);

                computePrimitives
                (
                    ppRight_[facei],
                    pURight_[facei],
                    pTRight_[facei],
                    localData,
                    pURight[facei],
                    pRRight[facei],
                    pBRight[facei],
                    pKRight[facei],
                    pm2Right[facei],
                    prhoRight[facei],
                    pgammaRight[facei],
                    ppotentialFace[facei],
                    ppotentialLeft[facei]
                );
            }
        }
        else // if (!pp.coupled())
        {
            forAll (pp, facei)
            {
                scalar rhoLocal = pp[facei]/(pR[facei]*pT[facei]);
                rhoUSource_[pFaceCells[facei]] +=
                    pSf[facei] * (rhoLocal*magSqr(pU[facei]) + pp[facei]);
            }
        }
    }

    Field<vector>& irhoUSource = rhoUSource_;
    Field<vector>& igradP = gradP_;
    Field<tensor>& igradU = gradU_;
    Field<vector>& igradT = gradT_;

    irhoUSource *= volumeInverse_;
    igradP *= 0.5*volumeInverse_;
    igradU *= 0.5*volumeInverse_;
    igradT *= 0.5*volumeInverse_;

    /// TODO: replace 0.5 factor with actual interpolation to faces

    rhoUSource_.correctBoundaryConditions();
    gradP_.correctBoundaryConditions();
    gradU_.correctBoundaryConditions();
    gradT_.correctBoundaryConditions();

    // Get limiters; alternative: VenkatakrishnanLimiter
    MDLimiter<scalar, Foam::BarthJespersenLimiter> scalarPLimiter
    (
        gradP_,
        pMaxValue,
        pMinValue
    );
    MDLimiter<vector, Foam::BarthJespersenLimiter> vectorULimiter
    (
        gradU_,
        UMaxValue,
        UMinValue
    );
    MDLimiter<scalar, Foam::BarthJespersenLimiter> scalarTLimiter
    (
        gradT_,
        TMaxValue,
        TMinValue
    );

    const volScalarField& pLimiter = scalarPLimiter.phiLimiter();
    const volVectorField& ULimiter = vectorULimiter.phiLimiter();
    const volScalarField& TLimiter = scalarTLimiter.phiLimiter();


    // Calculate fluxes at internal faces
    forAll (owner, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        const vector deltaRLeft = faceCentre[faceI] - cellCentre[own];
        const vector deltaRRight = faceCentre[faceI] - cellCentre[nei];

        // Calculate fluxes with reconstructed primitive variables at faces
        Flux::evaluateFlux
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
            pLeft_ [faceI] + pLimiter[own]*(deltaRLeft  & gradP_[own]),
            pRight_[faceI] + pLimiter[nei]*(deltaRRight & gradP_[nei]),
            ULeft_ [faceI] + cmptMultiply(ULimiter[own], (deltaRLeft  & gradU_[own])),
            URight_[faceI] + cmptMultiply(ULimiter[nei], (deltaRRight & gradU_[nei])),
            TLeft_ [faceI] + TLimiter[own]*(deltaRLeft  & gradT_[own]),
            TRight_[faceI] + TLimiter[nei]*(deltaRRight & gradT_[nei]),
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
        const fvPatch& curPatch = p_.boundaryField()[patchi].patch();

        // Fluxes
        fvsPatchScalarField& pRhoFlux  = rhoFlux_.boundaryField()[patchi];
        fvsPatchVectorField& pRhoUFlux = rhoUFlux_.boundaryField()[patchi];
        fvsPatchScalarField& pRhoEFlux = rhoEFlux_.boundaryField()[patchi];

        // Reconstructed values
        fvsPatchScalarField& ppLeft  = pLeft_.boundaryField()[patchi];
        fvsPatchScalarField& ppRight = pRight_.boundaryField()[patchi];
        fvsPatchVectorField& pULeft  = ULeft_.boundaryField()[patchi];
        fvsPatchVectorField& pURight = URight_.boundaryField()[patchi];
        fvsPatchScalarField& pTLeft  = TLeft_.boundaryField()[patchi];
        fvsPatchScalarField& pTRight = TRight_.boundaryField()[patchi];

        // Patch fields
        const fvPatchScalarField& pp = p_.boundaryField()[patchi];
        const vectorField& pU = U_.boundaryField()[patchi];
        const scalarField& pT = T_.boundaryField()[patchi];

        const scalarField& pCv = Cv.boundaryField()[patchi];
        const scalarField& pR = R.boundaryField()[patchi];

        // Gradients
        const fvPatchVectorField& pGradP = gradP_.boundaryField()[patchi];
        const fvPatchTensorField& pGradU = gradU_.boundaryField()[patchi];
        const fvPatchVectorField& pGradT = gradT_.boundaryField()[patchi];

        // Limiters
        const fvPatchScalarField& pPatchLim = pLimiter.boundaryField()[patchi];
        const fvPatchVectorField& UPatchLim = ULimiter.boundaryField()[patchi];
        const fvPatchScalarField& TPatchLim = TLimiter.boundaryField()[patchi];

        // Face areas
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        if (pp.coupled())
        {
            // Gradients
            const vectorField pgradPLeft = pGradP.patchInternalField();
            const vectorField pgradPRight = pGradP.patchNeighbourField();

            const tensorField pgradULeft = pGradU.patchInternalField();
            const tensorField pgradURight = pGradU.patchNeighbourField();

            const vectorField pgradTLeft = pGradT.patchInternalField();
            const vectorField pgradTRight = pGradT.patchNeighbourField();

            // Geometry: call the raw cell-to-face vector by calling
            // the base patch (cell-to-face) delta coefficient
            // Work out the right delta from the cell-to-cell delta
            // across the coupled patch and left delta
            vectorField pDeltaRLeft = curPatch.fvPatch::delta();
            vectorField pDdeltaRRight = pDeltaRLeft - curPatch.delta();

            // Limiters
            const scalarField ppLimiterLeft = pPatchLim.patchInternalField();
            const scalarField ppLimiterRight = pPatchLim.patchNeighbourField();

            const vectorField pULimiterLeft = UPatchLim.patchInternalField();
            const vectorField pULimiterRight = UPatchLim.patchNeighbourField();

            const scalarField pTLimiterLeft = TPatchLim.patchInternalField();
            const scalarField pTLimiterRight = TPatchLim.patchNeighbourField();

            forAll (pp, facei)
            {
                Flux::evaluateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],

                    ppLeft[facei]
                  + ppLimiterLeft[facei]*
                    (pDeltaRLeft[facei] & pgradPLeft[facei]),

                    ppRight[facei]
                  + ppLimiterRight[facei]*
                    (pDdeltaRRight[facei] & pgradPRight[facei]),

                    pULeft[facei]
                  + cmptMultiply
                    (
                        pULimiterLeft[facei],
                        pDeltaRLeft[facei] & pgradULeft[facei]
                    ),

                    pURight[facei]
                  + cmptMultiply
                    (
                        pULimiterRight[facei],
                        pDdeltaRRight[facei] & pgradURight[facei]
                    ),

                    pTLeft[facei]
                  + pTLimiterLeft[facei]*
                    (pDeltaRLeft[facei] & pgradTLeft[facei]),

                    pTRight[facei]
                  + pTLimiterRight[facei]*
                    (pDdeltaRRight[facei] & pgradTRight[facei]),

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


//template<class Flux>
//void Foam::numericFlux<Flux, Foam::balancedPotentialLimiter>::computePrimitives
//(
//    scalar& pFace,
//    vector& UFace,
//    scalar& TFace,
//    scalar& rhoUSource,
//    const vector& U,
//    const scalar& R,
//    const scalar& B,
//    const scalar& K,
//    const scalar& m2,
//    const scalar& rho,
//    const scalar& gamma,
//    const scalar& potential,
//    const scalar& potentialFace
//)
//{
//    // Initalize rhoFace with value at the cell centre
//    scalar rhoPrev, rhoFace = rho;
//
//    // Use Newton's method to solve for rhoFace at the cell face
//    unsigned int iter = 0;
//    const scalar tolerance = rho*tol_;
//    do
//    {
//        rhoPrev = rhoFace;
//        rhoFace -= ( 0.5*m2/sqr(rhoFace) + gamma/(gamma-1.0)*K*pow(rhoFace, gamma-1.0) +
//                 potentialFace - B )
//             / ( -m2/pow3(rhoFace) + gamma*K*pow(rhoFace, gamma-2.0) );
//        iter++;
//    }
//    while ( ( mag(rhoFace-rhoPrev) > tolerance ) && ( iter < 10000 ) );
//
//    if ( iter == 10000 )
//    {
//        Info << "maxIter (10000) reached computing face primitives.\n"
//             << "rhoFace = " << rhoFace << ",    rhoPrev = " << rhoPrev
//             << ",    residual = " << mag(rhoFace-rhoPrev) << endl;
//    }
//
//    // Use new value of rhoFace to compute other primitives at the cell face
//    pFace = K*pow(rhoFace, gamma);
//    UFace = U*rho/rhoFace;
//    TFace = pFace/(rhoFace*R);
//    rhoUSource = m2/rhoFace + pFace;
//}

template<class Flux>
void Foam::numericFlux<Flux, Foam::balancedPotentialLimiter>::computePrimitives
(
    scalar& pFace,
    vector& UFace,
    scalar& TFace,
    dataStruct& localData,
    const vector& U,
    const scalar& R,
    const scalar& B,
    const scalar& K,
    const scalar& m2,
    const scalar& rho,
    const scalar& gamma,
    const scalar& potentialFace,
    const scalar& potentialNei
)
{
    // Initalize rhoNext with value at the cell centre
    scalar rhoPrev, rhoNext = rho;

    // Use Newton's method to solve for rhoNext at the cell face
    unsigned int iter = 0;
    const scalar tolerance = rho*tol_;
    do
    {
        rhoPrev = rhoNext;
        rhoNext -= ( 0.5*m2/sqr(rhoNext) + gamma/(gamma-1.0)*K*pow(rhoNext, gamma-1.0) +
                 potentialFace - B )
             / ( -m2/pow3(rhoNext) + gamma*K*pow(rhoNext, gamma-2.0) );
        iter++;
    }
    while ( ( mag(rhoNext-rhoPrev) > tolerance ) && ( iter < 10000 ) );

    if ( iter == 10000 )
    {
        Info << "maxIter (10000) reached computing face primitives.\n"
             << "rhoNext = " << rhoNext << ",    rhoPrev = " << rhoPrev
             << ",    residual = " << mag(rhoNext-rhoPrev) << endl;
    }

    // Use new value of rhoNext to compute other primitives at the cell face
    pFace = K*pow(rhoNext, gamma);
    UFace = U*rho/rhoNext;
    TFace = pFace/(rhoNext*R);
    localData.rhoUSource = m2/rhoNext + pFace;

    // Use Newton's method to solve for rhoNext in the neighbour cell
    iter = 0;
    do
    {
        rhoPrev = rhoNext;
        rhoNext -= ( 0.5*m2/sqr(rhoNext) + gamma/(gamma-1.0)*K*pow(rhoNext, gamma-1.0) +
                 potentialNei - B )
             / ( -m2/pow3(rhoNext) + gamma*K*pow(rhoNext, gamma-2.0) );
        iter++;
    }
    while ( ( mag(rhoNext-rhoPrev) > tolerance ) && ( iter < 10000 ) );

    if ( iter == 10000 )
    {
        Info << "maxIter (10000) reached computing neighbour primitives.\n"
             << "rhoNext = " << rhoNext << ",    rhoPrev = " << rhoPrev
             << ",    residual = " << mag(rhoNext-rhoPrev) << endl;
    }

    // Use new value of rhoNext to compute other primitives in neighbour cell
    localData.p = K*pow(rhoNext, gamma);
    localData.U = U*rho/rhoNext;
    localData.T = pFace/(rhoNext*R);
}


// ************************************************************************ //
