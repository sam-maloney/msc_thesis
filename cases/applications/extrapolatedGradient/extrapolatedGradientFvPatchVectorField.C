/*---------------------------------------------------------------------------*\
| File modified by Engys Ltd 2010                                             |
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "extrapolatedGradientFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

extrapolatedGradientFvPatchVectorField::
extrapolatedGradientFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF)
{}


extrapolatedGradientFvPatchVectorField::
extrapolatedGradientFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF)
{
    fvPatchField<vector>::operator=(patchInternalField());
    gradient() = pTraits<vector>::zero;
}


extrapolatedGradientFvPatchVectorField::
extrapolatedGradientFvPatchVectorField
(
    const extrapolatedGradientFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(ptf, p, iF, mapper)
{}


extrapolatedGradientFvPatchVectorField::
extrapolatedGradientFvPatchVectorField
(
    const extrapolatedGradientFvPatchVectorField& ptf
)
:
    fixedGradientFvPatchVectorField(ptf)
{}


extrapolatedGradientFvPatchVectorField::
extrapolatedGradientFvPatchVectorField
(
    const extrapolatedGradientFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void extrapolatedGradientFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    word field = this->dimensionedInternalField().name();

    vectorField nf = patch().nf();

    volTensorField gradField
        = fvc::grad(db().lookupObject<volVectorField>(field));

    gradient()
        = (nf & gradField.boundaryField()[patch().index()]
        .patchInternalField()());

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void extrapolatedGradientFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    extrapolatedGradientFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
