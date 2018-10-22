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
    balancedMDLimiter

Description
    Multidimensional limiter for density based solver with balancedNumericFlux

Author
    Samuel Maloney
    Modified from MDLimiter by Aleksandar Jemcov
      which was updated by Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef balancedMDLimiter_H
#define balancedMDLimiter_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "BarthJespersenLimiter.H"
#include "VenkatakrishnanLimiter.H"
#include "MichalakGoochLimiter.H"
#include "WangLimiter.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class MDLimiter Declaration
\*---------------------------------------------------------------------------*/

template
<
    class Type,
    typename LimiterFunctionType
>
class balancedMDLimiter
{
public:

    // Public typedefs

    typedef Field<Type> FieldType;
    typedef GeometricField<Type, fvPatchField, volMesh> GeoFieldType;

    typedef Field<typename outerProduct<vector, Type>::type> GradFieldType;
    typedef GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    > GeoGradFieldType;


private:

    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Grad field
        const GeoGradFieldType& gradPhi_;

        //- Limiter
        GeoFieldType phiLimiter_;


public:

    // Constructors

        //- Construct from gradient field and min/max fields
        balancedMDLimiter
        (
            const GeoGradFieldType& gradPhi,
            const GeoFieldType& phiMaxValue,
            const GeoFieldType& phiMinValue
        )
        :
            mesh_(gradPhi.mesh()),
            gradPhi_(gradPhi),
            phiLimiter_
            (
                IOobject
                (
                    "phiLimiter",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensioned<Type>("pLimiter", dimless, pTraits<Type>::one)
            )
        {
            // Calculate limiter

            // Get geometrical information
            const DimensionedField<scalar, volMesh>& cellVolume = mesh_.V();
            const volVectorField& cellCentre = mesh_.C();
            const surfaceVectorField& faceCentre = mesh_.Cf();

            // Create limiter function
            LimiterFunctionType limitFunction;

            FieldType& phiLimiterIn = phiLimiter_.internalField();

            const GradFieldType& gradPhiIn = gradPhi.internalField();

            const unallocLabelList& owner = mesh_.owner();
            const unallocLabelList& neighbour = mesh_.neighbour();

            const Field<Type>& phiMinIn = phiMinValue.internalField();
            const Field<Type>& phiMaxIn = phiMaxValue.internalField();

            // Compute limiter values, internal faces
            forAll (owner, faceI)
            {
                const label& own = owner[faceI];
                const label& nei = neighbour[faceI];

                vector deltaRLeft = faceCentre[faceI] - cellCentre[own];
                vector deltaRRight = faceCentre[faceI] - cellCentre[nei];

                // Find minimal limiter value in each cell

                // Owner side
                limitFunction.limiter
                (
                    phiLimiterIn[own],
                    cellVolume[own],
                    phiMaxIn[own],
                    phiMinIn[own],
                    (deltaRLeft & gradPhiIn[own])
                );

                // Neighbour side
                limitFunction.limiter
                (
                    phiLimiterIn[nei],
                    cellVolume[nei],
                    phiMaxIn[nei],
                    phiMinIn[nei],
                    (deltaRRight & gradPhiIn[nei])
                );
            }

            // Coupled boundaries
            forAll (gradPhi.boundaryField(), patchI)
            {
                if (gradPhi.boundaryField()[patchI].coupled())
                {
                    // Get patch
                    const fvPatch& p = gradPhi.boundaryField()[patchI].patch();

                    const labelList& fc = p.faceCells();

                    const vectorField deltaR = p.Cf() - p.Cn();

                    // Get gradients
                    const GradFieldType gradPhiLeft =
                        gradPhi.boundaryField()[patchI].patchInternalField();

                    const GradFieldType gradPhiRight =
                        gradPhi.boundaryField()[patchI].patchNeighbourField();

                    // Find minimal limiter value in each cell
                    forAll (fc, faceI)
                    {
                        const label& curFC = fc[faceI];

                        limitFunction.limiter
                        (
                            phiLimiterIn[curFC],
                            cellVolume[curFC],
                            phiMaxIn[curFC],
                            phiMinIn[curFC],
                            (deltaR[faceI] & gradPhiLeft[faceI])
                        );
                    }
                }
            }

            // Do parallel communication to correct limiter on
            // coupled boundaries
            phiLimiter_.correctBoundaryConditions();
        }


    // Destructor - default


    // Member functions

        //- Return limiter
        const GeoFieldType& phiLimiter() const
        {
            return phiLimiter_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif // balancedMDLimiter_H

// ************************************************************************* //
