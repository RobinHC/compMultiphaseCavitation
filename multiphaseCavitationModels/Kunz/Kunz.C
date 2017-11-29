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

#include "Kunz.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace MultiphaseCavitations
{
    defineTypeNameAndDebug(Kunz, 0);
    addToRunTimeSelectionTable(MultiphaseCavitation, Kunz, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::MultiphaseCavitations::Kunz::Kunz
(
        const volVectorField& U,
        const surfaceScalarField& phi,
		const volScalarField& rho1,
		const volScalarField& rho2,
		const phaseModel& alpha1,
		const phaseModel& alpha2
)
:
	MultiphaseCavitation(typeName, U, phi,
			  rho1,
			  rho2,
			  alpha1,
			  alpha2),

    UInf_(MultiphaseCavitationCoeffs_.lookup("UInf")),
    tInf_(MultiphaseCavitationCoeffs_.lookup("tInf")),
    Cc_(MultiphaseCavitationCoeffs_.lookup("Cc")),
    Cv_(MultiphaseCavitationCoeffs_.lookup("Cv")),
    p0_("0", pSat().dimensions(), 0.0)
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::MultiphaseCavitations::Kunz::mDotAlphaW() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));

    // The coefficients must be calculated for each time step, because the density changes
    volScalarField mcCoeff_(Cc_*rho2_/tInf_);
    volScalarField mvCoeff_(Cv_*rho2_/(0.5*rho1_*sqr(UInf_)*tInf_));

    return Pair<tmp<volScalarField>>
    (
        mcCoeff_*sqr(limitedAlpha1)
       *max(p - pSat(), p0_)/max(p - pSat(), 0.01*pSat()),

        mvCoeff_*min(p - pSat(), p0_)
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::MultiphaseCavitations::Kunz::mDotP() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    volScalarField limitedAlpha2(min(max(alpha2_, scalar(0)), scalar(1)));
    // The coefficients must be calculated for each time step, because the density changes
    volScalarField mcCoeff_(Cc_*rho2_/tInf_);
    volScalarField mvCoeff_(Cv_*rho2_/(0.5*rho1_*sqr(UInf_)*tInf_));


    return Pair<tmp<volScalarField>>
    (
        mcCoeff_*sqr(limitedAlpha1)*limitedAlpha2
       *pos(p - pSat())/max(p - pSat(), 0.01*pSat()),

        (-mvCoeff_)*limitedAlpha1*neg(p - pSat())
    );
}


void Foam::MultiphaseCavitations::Kunz::correct()
{}


bool Foam::MultiphaseCavitations::Kunz::read()
{
	if (MultiphaseCavitation::read())
	{
	    MultiphaseCavitationCoeffs_ = this->subDict(type() + "Coeffs");

	    MultiphaseCavitationCoeffs_.lookup("UInf") >> UInf_;
	    MultiphaseCavitationCoeffs_.lookup("tInf") >> tInf_;
	    MultiphaseCavitationCoeffs_.lookup("Cc") >> Cc_;
	    MultiphaseCavitationCoeffs_.lookup("Cv") >> Cv_;
	    return true;
	}
	else
	{
		return false;
	}

}


// ************************************************************************* //
