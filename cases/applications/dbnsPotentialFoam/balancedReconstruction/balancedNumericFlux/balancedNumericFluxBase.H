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
    balancedNumericFluxBase

Description
    Base balanced numeric flux class for density-based solvers

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak
    Modified by Samuel Maloney

SourceFiles
    balancedNumericFluxBase.H
    balancedNumericFluxBase.C

\*---------------------------------------------------------------------------*/

#ifndef balancedNumericFluxBase_H
#define balancedNumericFluxBase_H

#include "basicBalancedNumericFlux.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class balancedNumericFluxBase Declaration
\*---------------------------------------------------------------------------*/

template<class Flux>
class balancedNumericFluxBase
:
    public basicBalancedNumericFlux,
    public Flux
{
public:

    // Constructors

        //- Construct from mesh
        balancedNumericFluxBase(const fvMesh& mesh)
        :
            basicBalancedNumericFlux(mesh)
        {}


    //- Destructor
    virtual ~balancedNumericFluxBase()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // balancedNumericFluxBase_H

// ************************************************************************* //
