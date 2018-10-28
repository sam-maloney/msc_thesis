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

#include "basicBalancedNumericFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicBalancedNumericFlux> Foam::basicBalancedNumericFlux::New
(
    const volScalarField& potential,
    const volVectorField& U,
    const volScalarField& rho,
    basicThermo& thermo
)
{
    const dictionary& subDict =
        U.mesh().schemesDict().subDict("divSchemes").subDict("dbns");

    word name = word(subDict.lookup("flux")) + "Flux"
        + word(subDict.lookup("limiter")) + "Limiter";

    Info<< "Selecting balancedNumericFlux " << name << endl;

    stateConstructorTable::iterator cstrIter =
        stateConstructorTablePtr_->find(name);

    if (cstrIter == stateConstructorTablePtr_->end())
    {
        FatalErrorIn("basicBalancedNumericFlux::New(const fvMesh&)")
            << "Unknown basicBalancedNumericFlux type " << name << nl << nl
            << "Valid basicBalancedNumericFlux types are:" << nl
            << stateConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<basicBalancedNumericFlux>
    (
        cstrIter()
        (
            potential,
            U,
            rho,
            thermo
        )
    );
}


// ************************************************************************* //
