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

Application
    dbnsPotentialFoam

Description
    Density-based compressible explicit time-marching flow solver
    using enthalpy based thermo packages and acting under an
    externally applied potential field

    Reference for time-marching scheme: https://hal.inria.fr/inria-00075385

Author
    Hrvoje Jasak
    Modified by Samuel Maloney

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "bound.H"
#include "numericFlux.H"
#include "balancedNumericFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "createTimeControls.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // Runge-Kutta coefficient
    scalarList beta(4);
    beta[0] = 0.1100;
    beta[1] = 0.2766;
    beta[2] = 0.5000;
    beta[3] = 1.0000;

    // Switch off solver messages
    lduMatrix::debug = 0;

    while (runTime.run())
    {
#       include "readTimeControls.H"
#       include "readFieldBounds.H"
#       include "compressibleCourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "\n Time = " << runTime.value() << endl;

        // Low storage Runge-Kutta time integration
        forAll (beta, i)
        {
            // Solve the approximate Riemann problem for this time step
            dbnsFlux.computeFlux();

            // Time integration
            solve
            (
                1.0/beta[i]*fvm::ddt(rho)
              + fvc::div(dbnsFlux.rhoFlux())
            );

            if ( isBalancedFlux )
            {
                solve
                (
                    1.0/beta[i]*fvm::ddt(rhoU)
                  + fvc::div(dbnsFlux.rhoUFlux())
                  - balancedDbnsFluxPtr->rhoUSource()
                );
            }
            else
            {
                solve
                (
                    1.0/beta[i]*fvm::ddt(rhoU)
                  + fvc::div(dbnsFlux.rhoUFlux())
                  + rho*fvc::grad(potential)
                );
            }

            solve
            (
                1.0/beta[i]*fvm::ddt(rhoE)
              + fvc::div(dbnsFlux.rhoEFlux())
              + ( rhoU & fvc::grad(potential) )
            );

#           include "updateFields.H"
        }

        runTime.write();

        Info<< "    ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl;
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
