/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::MultiphaseCavitations::Kunz

Description
    Kunz cavitation model slightly modified so that the condensation term
    is switched off when the pressure is less than the saturation vapour
    pressure.  This change allows the condensation term to be formulated as
    a coefficient multiplying (p - p_sat) so that it can be included as an
    implicit term in the pressure equation.

    Reference:
    \verbatim
        Kunz, R.F., Boger, D.A., Stinebring, D.R., Chyczewski, Lindau. J.W.,
        Gibeling, H.J., Venkateswaran, S., Govindan, T.R.,
        "A Preconditioned Implicit Method for Two-Phase Flows with Application
         to Cavitation Prediction,"
        Computers and Fluids,
        29(8):849-875, 2000.
    \endverbatim

SourceFiles
    Kunz.C

\*--------------------------------------------------------------------*/

#ifndef Kunz_H
#define Kunz_H

#include "MultiphaseCavitation.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace MultiphaseCavitations
{

/*--------------------------------------------------------------------*\
                              Class Kunz
\*--------------------------------------------------------------------*/

class Kunz
:
public MultiphaseCavitation
{
    // Private data

        dimensionedScalar UInf_;
        dimensionedScalar tInf_;
        dimensionedScalar Cc_;
        dimensionedScalar Cv_;

        dimensionedScalar p0_;


public:

    //- Runtime type information
    TypeName("Kunz");


    // Constructors

        //- Construct from components
        Kunz
        (
                const volVectorField& U,
                const surfaceScalarField& phi,
    			const volScalarField& rho1,
    			const volScalarField& rho2,
    			const phaseModel& alpha1,
    			const phaseModel& alpha2
        );


    //- Destructor
    virtual ~Kunz()
    {}


    // Member Functions

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField>> mDotAlphaW() const;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        virtual Pair<tmp<volScalarField>> mDotP() const;

        //- Correct the Kunz phaseChange model
        virtual void correct();

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace MultiphaseCavitations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
