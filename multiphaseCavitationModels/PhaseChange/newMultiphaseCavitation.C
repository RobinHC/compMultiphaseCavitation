/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::MultiphaseCavitation>
Foam::MultiphaseCavitation::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
	const volScalarField& rho1,
	const volScalarField& rho2,
	const phaseModel& alpha1,
	const phaseModel& alpha2
)
{
    IOdictionary transportPropertiesDict
    (
        IOobject
        (
            "thermophysicalProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word MultiphaseCavitationTypeName
    (
        transportPropertiesDict.lookup("MultiphaseCavitation")
    );

    Info<< "Selecting phaseChange model "
        << MultiphaseCavitationTypeName << endl;

    componentsConstructorTable::iterator cstrIter =
        componentsConstructorTablePtr_
            ->find(MultiphaseCavitationTypeName);

    if (cstrIter == componentsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown MultiphaseCavitation type "
            << MultiphaseCavitationTypeName << endl << endl
            << "Valid  MultiphaseCavitations are : " << endl
            << componentsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<MultiphaseCavitation>(cstrIter()(U, phi,rho1,rho2,alpha1,alpha2));
}


// ************************************************************************* //
