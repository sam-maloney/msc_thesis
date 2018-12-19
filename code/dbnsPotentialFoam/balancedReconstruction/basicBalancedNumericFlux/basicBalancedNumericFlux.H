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
    Foam::basicBalancedNumericFlux

Description
    Base class for runtime-selectable balanced flux and limiters

Author
    Henrik Rusche
    Modified by Samuel Maloney

SourceFiles
    basicBalancedNumericFlux.C
    newBasicBalancedNumericFlux.C

\*---------------------------------------------------------------------------*/

#ifndef basicBalancedNumericFlux_H
#define basicBalancedNumericFlux_H

#include "basicNumericFlux.H"
#include "fvMesh.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class basicBalancedNumericFlux Declaration
\*---------------------------------------------------------------------------*/

class basicBalancedNumericFlux
:
    public basicNumericFlux
{
protected:

    // Protected member functions

        //- Construct as copy (not implemented)
        basicBalancedNumericFlux(const basicBalancedNumericFlux&);


public:

    //- Runtime type information
    TypeName("basicBalancedNumericFlux");


    //- Declare run-time constructor selection table

#ifndef SWIG
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicBalancedNumericFlux,
        state,
        (
            const volScalarField& potential,
            const volVectorField& U,
            const volScalarField& rho,
            basicThermo& thermo
        ),
        (potential, U, rho, thermo)
    );
#endif


    // Constructors

        //- Construct from components
        basicBalancedNumericFlux
        (
            const fvMesh& mesh
        );


    //- Selector
    static autoPtr<basicBalancedNumericFlux> New
    (
        const volScalarField& potential,
        const volVectorField& U,
        const volScalarField& rho,
        basicThermo& thermo
    );


    //- Destructor
    virtual ~basicBalancedNumericFlux();


    // Member functions

        //- Return well-balanced momentum source term
        virtual const volVectorField& rhoUSource() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // basicBalancedNumericFlux_H

// ************************************************************************* //
