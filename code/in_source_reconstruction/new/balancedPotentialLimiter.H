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
    balancedPotentialLimiter

Description
    Dummy limiter class to force specialization of the numericFlux class.
    Specialization implemented in .../src/dbns/numericFlux/balancedNumericFlux
    files, where the actual balanced reconstruction is performed

Author
    Samuel Maloney
    Based on class structure from
        Aleksandar Jemcov
        Rewritten by Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef balancedPotentialLimiter_H
#define balancedPotentialLimiter_H

#include "firstOrderLimiter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class balancedPotentialLimiter Declaration
\*---------------------------------------------------------------------------*/

class balancedPotentialLimiter
:
    public firstOrderLimiter
{
public:

    // Constructor

        //- Construct null
        balancedPotentialLimiter()
        {}


    // Destructor - default


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
