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

#include "balancedNumericFlux.H"
#include "balancedMDLimiter.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Flux, class Limiter>
Foam::balancedNumericFlux<Flux, Limiter>::balancedNumericFlux
(
    const volScalarField& potential,
    const volVectorField& U,
    const volScalarField& rho,
    basicThermo& thermo
)
:
    balancedNumericFluxBase<Flux>(U.mesh()),
    p_(thermo.p()),
    U_(U),
    T_(thermo.T()),
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
        -fvc::grad(p_)
    ),
    volumeInverse_( 1.0/this->mesh().V() ),
    tol_
    (
        readScalar
        (
            this->mesh().solutionDict().subDict("solvers")
                        .subDict("rho").lookup("tolerance")
        )
    ),
    maxIter_
    (
        readLabel
        (
            this->mesh().solutionDict().subDict("solvers")
                        .subDict("rho").lookup("maxIter")
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Flux, class Limiter>
void Foam::balancedNumericFlux<Flux, Limiter>::computeFlux()
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

    // Compute adiabatic index, momentum squared, Bernoulli, and EOS constants
    const volScalarField gamma = R/Cv + 1.0;
    const volScalarField m2 = magSqr(rho_*U_);
    const volScalarField B = 0.5*magSqr(U_) + gamma/(gamma-1.0)*R*T_ + potential_;
    const volScalarField K = p_ / pow(rho_, gamma);

    // Reset momentum source and gradients back to zero
    Field<vector>& irhoUSource = rhoUSource_;
    Field<vector>& igradP = gradP_;
    Field<tensor>& igradU = gradU_;
    Field<vector>& igradT = gradT_;

    rhoUSource_.Field<vector>::operator=(vector::zero);
    gradP_.Field<vector>::operator=(vector::zero);
    gradU_.Field<tensor>::operator=(tensor::zero);
    gradT_.Field<vector>::operator=(vector::zero);

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

        // Note: mag in the dot-product.
        // For all valid meshes, the non-orthogonality will be less than
        // 90 deg and the dot-product will be positive.  For invalid
        // meshes (d & s <= 0), this will stabilise the calculation
        // but the result will be poor.
        const scalar SfdOwn = mag(Sf[faceI]&(faceCentre[faceI] - cellCentre[own]));
        const scalar SfdNei = mag(Sf[faceI]&(cellCentre[nei] - faceCentre[faceI]));
        const scalar weightLeft = SfdOwn/(SfdOwn + SfdNei);
        const scalar weightRight = SfdNei/(SfdOwn + SfdNei);

        computePrimitives
        (
            localData,
            U_[own],
            R[own],
            B[own],
            K[own],
            m2[own],
            rho_[own],
            gamma[own],
            potential_[nei]
        );
//        scalar rhoFace = weightLeft*localData.rho + weightRight*rho_[own];
        scalar rhoFace = (localData.rho + rho_[own])/2.0;
        pLeft_[faceI] = K[own]*pow(rhoFace, gamma[own]);
        ULeft_[faceI] = U_[own]*rho_[own]/rhoFace;
        TLeft_[faceI] = pLeft_[faceI]/(rhoFace*R[own]);

        // Calculate local momemtum source contribution
        rhoUSource_[own] += Sf[faceI]*(m2[own]/rhoFace + pLeft_[faceI]);
        gradP_[own] += Sf[faceI]*(p_[nei] - localData.p)*weightLeft;
        gradU_[own] += Sf[faceI]*(U_[nei] - localData.U)*weightLeft;
        gradT_[own] += Sf[faceI]*(T_[nei] - localData.T)*weightLeft;

        pMinIn[own] = min(pMinIn[own], p_[nei] - localData.p);
        pMaxIn[own] = max(pMaxIn[own], p_[nei] - localData.p);
        UMinIn[own] = min(UMinIn[own], U_[nei] - localData.U);
        UMaxIn[own] = max(UMaxIn[own], U_[nei] - localData.U);
        TMinIn[own] = min(TMinIn[own], T_[nei] - localData.T);
        TMaxIn[own] = max(TMaxIn[own], T_[nei] - localData.T);

        computePrimitives
        (
            localData,
            U_[nei],
            R[nei],
            B[nei],
            K[nei],
            m2[nei],
            rho_[nei],
            gamma[nei],
            potential_[own]
        );

//        rhoFace = weightRight*localData.rho + weightLeft*rho_[nei];
        rhoFace = (localData.rho + rho_[nei])/2.0;
        pRight_[faceI] = K[nei]*pow(rhoFace, gamma[nei]);
        URight_[faceI] = U_[nei]*rho_[nei]/rhoFace;
        TRight_[faceI] = pRight_[faceI]/(rhoFace*R[nei]);

        rhoUSource_[nei] -= Sf[faceI]*(m2[nei]/rhoFace + pRight_[faceI]);
        gradP_[nei] -= Sf[faceI]*(p_[own] - localData.p)*weightRight;
        gradU_[nei] -= Sf[faceI]*(U_[own] - localData.U)*weightRight;
        gradT_[nei] -= Sf[faceI]*(T_[own] - localData.T)*weightRight;

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
        const fvsPatchVectorField& pCf = faceCentre.boundaryField()[patchi];
        const unallocLabelList& pFaceCells =
            this->mesh().boundary()[patchi].faceCells();

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
            const vectorField pcellCentreLeft  =
                cellCentre.boundaryField()[patchi].patchInternalField();
            const vectorField pcellCentreRight =
                cellCentre.boundaryField()[patchi].patchNeighbourField();


            forAll (pp, facei)
            {
                computePrimitives
                (
                    localData,
                    pULeft[facei],
                    pRLeft[facei],
                    pBLeft[facei],
                    pKLeft[facei],
                    pm2Left[facei],
                    prhoLeft[facei],
                    pgammaLeft[facei],
                    ppotentialRight[facei]
                );

                const scalar SfdOwn = mag(pSf[facei]&(pCf[facei] - pcellCentreLeft[facei]));
                const scalar SfdNei = mag(pSf[facei]&(pcellCentreRight[facei] - pCf[facei]));
                const scalar weightLeft = SfdOwn/(SfdOwn + SfdNei);
                const scalar weightRight = SfdNei/(SfdOwn + SfdNei);

                const label own = pFaceCells[facei];

//                scalar rhoFace = weightLeft*localData.rho + weightRight*prhoLeft[facei];
                scalar rhoFace = (localData.rho + prhoLeft[facei])/2.0;
                ppLeft_[facei] = pKLeft[facei]*pow(rhoFace, pgammaLeft[facei]);
                pULeft_[facei] = pULeft[facei]*prhoLeft[facei]/rhoFace;
                pTLeft_[facei] = ppLeft_[facei]/(rhoFace*pRLeft[facei]);

                rhoUSource_[own] += pSf[facei]*(pm2Left[facei]/rhoFace + ppLeft_[facei]);
                gradP_[own] += pSf[facei]*(ppRight[facei] - localData.p)*weightLeft;
                gradU_[own] += pSf[facei]*(pURight[facei] - localData.U)*weightLeft;
                gradT_[own] += pSf[facei]*(pTRight[facei] - localData.T)*weightLeft;

                pMinIn[own] = min(pMinIn[own], ppRight[facei] - localData.p);
                pMaxIn[own] = max(pMaxIn[own], ppRight[facei] - localData.p);
                UMinIn[own] = min(UMinIn[own], pURight[facei] - localData.U);
                UMaxIn[own] = max(UMaxIn[own], pURight[facei] - localData.U);
                TMinIn[own] = min(TMinIn[own], pTRight[facei] - localData.T);
                TMaxIn[own] = max(TMaxIn[own], pTRight[facei] - localData.T);

                computePrimitives
                (
                    localData,
                    pURight[facei],
                    pRRight[facei],
                    pBRight[facei],
                    pKRight[facei],
                    pm2Right[facei],
                    prhoRight[facei],
                    pgammaRight[facei],
                    ppotentialLeft[facei]
                );

//                rhoFace = weightRight*localData.rho + weightLeft*prhoRight[facei];
                rhoFace = (localData.rho + prhoRight[facei])/2.0;
                ppRight_[facei] = pKRight[facei]*pow(rhoFace, pgammaRight[facei]);
                pURight_[facei] = pURight[facei]*prhoRight[facei]/rhoFace;
                pTRight_[facei] = ppRight_[facei]/(rhoFace*pRRight[facei]);
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

    irhoUSource *= volumeInverse_;
    igradP *= volumeInverse_;
    igradU *= volumeInverse_;
    igradT *= volumeInverse_;

    /// TODO: replace 0.5 factor with actual interpolation to faces

    rhoUSource_.correctBoundaryConditions();
    gradP_.correctBoundaryConditions();
    gradU_.correctBoundaryConditions();
    gradT_.correctBoundaryConditions();

    balancedMDLimiter<scalar, Limiter> scalarPLimiter
    (
        gradP_,
        pMaxValue,
        pMinValue
    );
    balancedMDLimiter<vector, Limiter> vectorULimiter
    (
        gradU_,
        UMaxValue,
        UMinValue
    );
    balancedMDLimiter<scalar, Limiter> scalarTLimiter
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


template<class Flux, class Limiter>
void Foam::balancedNumericFlux<Flux, Limiter>::computePrimitives
(
    dataStruct& localData,
    const vector& U,
    const scalar& R,
    const scalar& B,
    const scalar& K,
    const scalar& m2,
    const scalar& rho,
    const scalar& gamma,
    const scalar& potentialNei
)
{
    // Initalize rhoNext with value at the cell centre
    scalar rhoPrev, rhoNext = rho;

    // Use Newton's method to solve for rhoNext in the neighbour cell
    label iter = 0;
    const scalar tolerance = rho*tol_;
    do
    {
        rhoPrev = rhoNext;
        rhoNext -= ( 0.5*m2/sqr(rhoNext) + gamma/(gamma-1.0)*K*pow(rhoNext, gamma-1.0) +
                 potentialNei - B )
             / ( -m2/pow3(rhoNext) + gamma*K*pow(rhoNext, gamma-2.0) );
        iter++;
    }
    while ( ( mag(rhoNext-rhoPrev) > tolerance ) && ( iter < maxIter_ ) );

    if ( iter == maxIter_ )
    {
        Info << "Maximum number of iterations exceeded.\n"
             << "rhoNext = " << rhoNext << ",    rhoPrev = " << rhoPrev
             << ",    residual = " << mag(rhoNext-rhoPrev) << endl;
    }

    // Use new value of rhoNext to compute other primitives in neighbour cell
    localData.rho = rhoNext;
    localData.p = K*pow(rhoNext, gamma);
    localData.U = U*rho/rhoNext;
    localData.T = localData.p/(rhoNext*R);
}


// ************************************************************************ //
