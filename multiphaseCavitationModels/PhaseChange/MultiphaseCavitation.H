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
    Foam::MultiphaseCavitation

Description

SourceFiles
    MultiphaseCavitation.C
    newPhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef MultiphaseCavitation_H
#define MultiphaseCavitation_H


#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

#include "regIOobject.H"
#include "dictionary.H"

#include "phaseModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class MultiphaseCavitation Declaration
\*---------------------------------------------------------------------------*/

class phaseModel;

class MultiphaseCavitation
:
		public IOdictionary
{

protected:

    // Protected data

        dictionary MultiphaseCavitationCoeffs_;

        //- Saturation vapour pressure
        dimensionedScalar pSat_;

        // Fields
        const volScalarField& rho1_;
        const volScalarField& rho2_;

        //Using phaseModel class instead of volScalarField
        //- water field
        const phaseModel& alpha1_;
        //- water vapor feild
        const phaseModel& alpha2_;


    // Private Member Functions

        //- Disallow copy construct
        MultiphaseCavitation(const MultiphaseCavitation&);

        //- Disallow default bitwise assignment
        void operator=(const MultiphaseCavitation&);


public:

    //- Runtime type information
    TypeName("MultiphaseCavitation");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            MultiphaseCavitation,
            components,
            (
                const volVectorField& U,
                const surfaceScalarField& phi,
				const volScalarField& rho1,
				const volScalarField& rho2,
				const phaseModel& alpha1,
				const phaseModel& alpha2
            ),
            (U, phi,
    			  rho1,
				  rho2,
    			  alpha1,
    			  alpha2)
        );

    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<MultiphaseCavitation> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
			const volScalarField& rho1,
			const volScalarField& rho2,
			const phaseModel& alpha1,
			const phaseModel& alpha2
        );


    // Constructors

        //- Construct from components
        MultiphaseCavitation
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi,
			const volScalarField& rho1,
			const volScalarField& rho2,
			const phaseModel& alpha1,
			const phaseModel& alpha2
        );


    //- Destructor
    virtual ~MultiphaseCavitation()
    {}


    // Member Functions

        //- Return const-access to the saturation vapour pressure
        const dimensionedScalar& pSat() const
        {
            return pSat_;
        }

		void testFunction();

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
	    //- Pair[0]: production term
		//  Pair[1}: destruction term
        virtual Pair<tmp<volScalarField>> mDotAlphaW() const = 0;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        //- Pair[0]: production term
    	//  Pair[1}: destruction term
        virtual Pair<tmp<volScalarField>> mDotP() const = 0;

        //- Return the volumetric condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        //- Pair[0]: production term
    	//  Pair[1}: destruction term
        Pair<tmp<volScalarField>> vDotAlphaW() const;

        //- Returns the volumetric condensation and vaporization rates of vapor base.
        //  Calculate based on the water phase, therefore using the complementary values for
        //  mass transfer

        Pair<tmp<volScalarField>> vDotAlphaV() const;

        //- Return the volumetric condensation and vaporisation rates as
        //  coefficients to multiply (p - pSat)
        //- Pair[0]: production term
    	//  Pair[1}: destruction term
        Pair<tmp<volScalarField>> vDotP() const;

        //- Correct the phaseChange model
        virtual void correct() = 0;

        //- Read the transportProperties dictionary and update
        virtual bool read() = 0;



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
