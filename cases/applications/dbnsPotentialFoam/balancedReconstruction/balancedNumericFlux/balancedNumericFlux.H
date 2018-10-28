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

Class
    balancedNumericFlux

Description
    Single level numeric flux class for density-based solvers.
    Uses well-balanced reconstruction for external potential.
    Reference: http://dx.doi.org/10.1016/j.jcp.2013.11.028

Author
    Samuel Maloney
    Modified from numericFlux by Aleksandar Jemcov
      which was rewritten by Hrvoje Jasak

SourceFiles
    balancedNumericFlux.H
    balancedNumericFlux.C

\*---------------------------------------------------------------------------*/

#ifndef balancedNumericFlux_H
#define balancedNumericFlux_H

#include "balancedNumericFluxBase.H"
#include "basicThermo.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class balancedNumericFlux Declaration
\*---------------------------------------------------------------------------*/

template<class Flux, class Limiter>
class balancedNumericFlux
:
    public balancedNumericFluxBase<Flux>
{
    // Private data

        // Reference to primitive fields

            //- Static pressure
            const volScalarField& p_;

            //- Velocity
            const volVectorField& U_;

            //- Static temperature
            const volScalarField& T_;

            //- Mass density
            const volScalarField rho_;

            //- External potential
            const volScalarField& potential_;

            //- Reference to the thermophysicalModel
            basicThermo& thermo_;


        // Reconstructed Primitives at cell faces

            //- Static pressure
            surfaceScalarField pLeft_;
            surfaceScalarField pRight_;

            //- Velocity
            surfaceVectorField ULeft_;
            surfaceVectorField URight_;

            //- Static temperature
            surfaceScalarField TLeft_;
            surfaceScalarField TRight_;


        // Gradients of reconstruction

            //- Static pressure
            volVectorField gradP_;

            //- Velocity
            volTensorField gradU_;

            //- Static temperature
            volVectorField gradT_;


        // Fluxes

            //- Density flux
            surfaceScalarField rhoFlux_;

            //- Velocity flux
            surfaceVectorField rhoUFlux_;

            //- Energy flux
            surfaceScalarField rhoEFlux_;


        // Sources

            //- Momentum Source
            volVectorField rhoUSource_;

            //- Inverse volume for source computation
            const DimensionedField<scalar, volMesh> volumeInverse_;


        // Other method parameters

            // Newton's method tolerance for rho
            const scalar tol_;

            // Newton's method maximum number of iterations
            const label maxIter_;


        // Gradients: no longer stored.  HJ, 22/Apr/2016


    // Private Member Functions

        //- Disallow default bitwise copy construct
        balancedNumericFlux(const balancedNumericFlux&);

        //- Disallow default bitwise assignment
        void operator=(const balancedNumericFlux&);


        //- Return internal field of mass flux
        const scalarField& rhoFluxI() const
        {
            return rhoFlux_.internalField();
        }

        //- Return access to internal field of mass flux
        scalarField& rhoFluxI()
        {
            return rhoFlux_.internalField();
        }

        //- Return internal field of momentum flux
        const vectorField& rhoUFluxI() const
        {
            return rhoUFlux_.internalField();
        }

        //- Return access to internal field of momentum flux
        vectorField& rhoUFluxI()
        {
            return rhoUFlux_.internalField();
        }

         //- Return access to internal field of energy flux
        const scalarField& rhoEFluxI() const
        {
            return rhoEFlux_.internalField();
        }

        //- Return access to internal field of energy flux
        scalarField& rhoEFluxI()
        {
            return rhoEFlux_.internalField();
        }


public:

    //- Runtime type information
    TypeName("balancedNumericFlux");

    // Constructors

        //- Construct from components
        balancedNumericFlux
        (
            const volScalarField& potential,
            const volVectorField& U,
            const volScalarField& rho,
            basicThermo& thermo
        );


    //- Destructor
    virtual ~balancedNumericFlux()
    {}


    // Typedef of a structure for passing local data
    struct dataStruct
    {
        scalar rhoUSource;
        scalar p;
        vector U;
        scalar T;
    };


    // Member Functions

        // Return fluxes

            //- Return density flux
            virtual const surfaceScalarField& rhoFlux() const
            {
                return rhoFlux_;
            }

            //- Return velocity flux
            virtual const surfaceVectorField& rhoUFlux() const
            {
                return rhoUFlux_;
            }

            //- Return energy flux
            virtual const surfaceScalarField& rhoEFlux() const
            {
                return rhoEFlux_;
            }


       // Return residuals

            //- Return density equation residual
            virtual tmp<scalarField> rhoResidual() const
            {
                return fvc::div(rhoFlux_)().internalField();
            }

            //- Return momentum equation flux
            virtual tmp<vectorField> rhoUResidual() const
            {
                return fvc::div(rhoUFlux_)().internalField();
            }

            //- Return energy equation flux
            virtual tmp<scalarField> rhoEResidual() const
            {
                return fvc::div(rhoEFlux_)().internalField();
            }


        // Return source terms

            //- Return well-balanced momentum source term
            virtual const volVectorField& rhoUSource() const
            {
                return rhoUSource_;
            }


        // Update fluxes based on current state

            //- Compute flux
            virtual void computeFlux();


        // Computes primitive values for the given cell face potential

            //- Compute primitives at cell face
            void computePrimitives
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
                const scalar& potential,
                const scalar& potentialFace
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "balancedNumericFlux.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // balancedNumericFlux_H

// ************************************************************************* //
