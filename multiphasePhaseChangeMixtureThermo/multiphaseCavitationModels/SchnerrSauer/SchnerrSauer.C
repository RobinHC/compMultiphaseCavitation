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

#include "SchnerrSauer.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace MultiphaseCavitations
{
    defineTypeNameAndDebug(SchnerrSauer, 0);
    addToRunTimeSelectionTable
    (
    	MultiphaseCavitation,
        SchnerrSauer,
        components
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MultiphaseCavitations::SchnerrSauer::SchnerrSauer
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

    n_(MultiphaseCavitationCoeffs_.lookup("n")),
    dNuc_(MultiphaseCavitationCoeffs_.lookup("dNuc")),
    Cc_(MultiphaseCavitationCoeffs_.lookup("Cc")),
    Cv_(MultiphaseCavitationCoeffs_.lookup("Cv")),
    p0_("0", pSat().dimensions(), 0.0)
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::MultiphaseCavitations::SchnerrSauer::rRb
(
    const volScalarField& limitedAlpha1
) const
{
    return pow
    (
        ((4*constant::mathematical::pi*n_)/3)
       *limitedAlpha1/(1.0 + alphaNuc() - limitedAlpha1),
        1.0/3.0
    );
}


Foam::dimensionedScalar
Foam::MultiphaseCavitations::SchnerrSauer::alphaNuc() const
{
    dimensionedScalar Vnuc = n_*constant::mathematical::pi*pow3(dNuc_)/6;
    return Vnuc/(1 + Vnuc);
}


Foam::tmp<Foam::volScalarField>
Foam::MultiphaseCavitations::SchnerrSauer::pCoeff
(
    const volScalarField& p
) const
{
    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    // Using the rho field of the mixture instead of the weighted rho field of the bulk densities for
    // each phase.
    // Doing this captures the incompressbile effects calculated by the solver
    const volScalarField& rho = alpha1_.db().lookupObject<volScalarField>("rho");

    return
        (3*rho1_*rho2_)*sqrt(2/(3*rho1_))
       *rRb(limitedAlpha1)/(rho*sqrt(mag(p - pSat()) + 0.01*pSat()));
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::MultiphaseCavitations::SchnerrSauer::mDotAlphaW() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeff(this->pCoeff(p));

    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));

    //- Pair[0]: production term
	//  Pair[1}: destruction term
    return Pair<tmp<volScalarField>>
    (
        Cc_*limitedAlpha1*pCoeff*max(p - pSat(), p0_),

        Cv_*(1.0 + alphaNuc() - limitedAlpha1)*pCoeff*min(p - pSat(), p0_)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::MultiphaseCavitations::SchnerrSauer::mDotP() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeff(this->pCoeff(p));

    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    volScalarField apCoeff(limitedAlpha1*pCoeff);

    //- Pair[0]: production term
	//  Pair[1}: destruction term
    return Pair<tmp<volScalarField>>
    (
        Cc_*(1.0 - limitedAlpha1)*pos(p - pSat())*apCoeff,

        (-Cv_)*(1.0 + alphaNuc() - limitedAlpha1)*neg(p - pSat())*apCoeff
    );
}


void Foam::MultiphaseCavitations::SchnerrSauer::correct()
{}


bool Foam::MultiphaseCavitations::SchnerrSauer::read()
{
	if (MultiphaseCavitation::read())
	{
		MultiphaseCavitationCoeffs_ = subDict(type() + "Coeffs");

		MultiphaseCavitationCoeffs_.lookup("n") >> n_;
        MultiphaseCavitationCoeffs_.lookup("dNuc") >> dNuc_;
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
