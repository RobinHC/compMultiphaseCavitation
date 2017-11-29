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

\*---------------------------------------------------------------------------*/

#include "MultiphaseCavitation.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MultiphaseCavitation, 0);
    defineRunTimeSelectionTable(MultiphaseCavitation, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::MultiphaseCavitation::MultiphaseCavitation
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi,
	const volScalarField& rho1,
	const volScalarField& rho2,
	const phaseModel& alpha1,
	const phaseModel& alpha2
)
:
	IOdictionary
    (
        IOobject
        (
            "thermophysicalProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    MultiphaseCavitationCoeffs_(subDict(type + "Coeffs")),
    pSat_("pSat", dimPressure, lookup("pSat")),
	rho1_(rho1),
	rho2_(rho2),
	alpha1_(alpha1),
	alpha2_(alpha2)

{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::MultiphaseCavitation::vDotAlphaW() const
{
	volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    volScalarField alphaWCoeff(1.0/rho1_ - limitedAlpha1*(1.0/rho1_ - 1.0/rho2_));
    Pair<tmp<volScalarField>> mDotAlphaW = this->mDotAlphaW();

    //- Pair[0]: production term
	//  Pair[1}: destruction term
    return Pair<tmp<volScalarField>>
    (
        alphaWCoeff*mDotAlphaW[0],
        alphaWCoeff*mDotAlphaW[1]
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::MultiphaseCavitation::vDotAlphaV() const
{
	volScalarField limitedAlpha2(min(max(alpha2_, scalar(0)), scalar(1)));
    volScalarField alphaVCoeff(1.0/rho2_ + limitedAlpha2*(1.0/rho1_ - 1.0/rho2_));

    Pair<tmp<volScalarField>> mDotAlphaW = this->mDotAlphaW();
    // Destruction term of the vapor mass transfer rate
    volScalarField mDotAlphaVdes(-1.0 * mDotAlphaW[0]);
    // Production term of the vapor mass transfer rate
    volScalarField mDotAlphaVprod(-1.0 * mDotAlphaW[1]);

    //- Pair[0]: production term
	//  Pair[1}: destruction term
    return Pair<tmp<volScalarField>>
    (
        alphaVCoeff*mDotAlphaVprod,
        alphaVCoeff*mDotAlphaVdes
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::MultiphaseCavitation::vDotP() const
{
    volScalarField pCoeff(1.0/rho1_ - 1.0/rho2_);
    Pair<tmp<volScalarField>> mDotP = this->mDotP();

    return Pair<tmp<volScalarField>>(pCoeff*mDotP[0], pCoeff*mDotP[1]);
}


bool Foam::MultiphaseCavitation::read()
{
		if (regIOobject::read())
		{
	        MultiphaseCavitationCoeffs_ = subDict(type() + "Coeffs");
	        pSat_ = lookup("pSat");

	        return true;
		}
		else
		{
			return false;
		}




}


// ************************************************************************* //
