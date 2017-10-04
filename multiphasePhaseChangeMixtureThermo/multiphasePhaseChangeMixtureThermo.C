/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "alphaContactAngleFvPatchScalarField.H"
#include "Time.H"
#include "subCycle.H"
#include "MULES.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcMeshPhi.H"
#include "multiphasePhaseChangeMixtureThermo.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphasePhaseChangeMixtureThermo, 0);
}


const Foam::scalar Foam::multiphasePhaseChangeMixtureThermo::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiphasePhaseChangeMixtureThermo::calcAlphas()
{
    scalar level = 0.0;
    alphas_ == 0.0;

    forAllIter(PtrDictionary<phaseModel>, phases_, phase)
    {
        alphas_ += level*phase();
        level += 1.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphasePhaseChangeMixtureThermo::multiphasePhaseChangeMixtureThermo
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    psiThermo(U.mesh(), word::null),

    phases_(lookup("phases"), phaseModel::iNew(p_, T_)),

    mesh_(U.mesh()),
    U_(U),
    phi_(phi),

    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("rhoPhi", dimMass/dimTime, 0.0)
    ),

    alphas_
    (
        IOobject
        (
            "alphas",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphas", dimless, 0.0)
    ),
    alphaSum_
    (
        IOobject
        (
            "alphaSum",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphaSum_", dimless, 0.0)
    ),

    sigmas_(lookup("sigmas")),
    dimSigma_(1, 0, -2, 0, 0),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    )
{

	//Initializing the cavitation model
	forAllIter(PtrDictionary<phaseModel>, phases_, phasei)
	{
		if (phasei().name() == "water")
		{
			//Info << "phasei name: " << phasei().name() << "\n";
			const phaseModel& alpha1 = phasei();
			tmp<volScalarField> TMPrho1(phasei().thermo().rho());
			const volScalarField& rho1(TMPrho1());

			forAllIter(PtrDictionary<phaseModel>, phases_, phasej)
			{
				if (phasej().name() == "vapor")
				{
					//Info << "phasej name: " << phasej().name() << "\n";
					const phaseModel& alpha2 = phasej();
					tmp<volScalarField> TMPrho2(phasej().thermo().rho());
					const volScalarField& rho2(TMPrho2());

					cavitationModel_ = MultiphaseCavitation::New(U, 	phi,
																		rho1,
																		rho2,
																		alpha1,
																		alpha2);
				}
			}
		}
	}

    calcAlphas();
    alphas_.write();
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::multiphasePhaseChangeMixtureThermo::correct()
{
    forAllIter(PtrDictionary<phaseModel>, phases_, phasei)
    {
        phasei().correct();
    }

    PtrDictionary<phaseModel>::iterator phasei = phases_.begin();

    psi_ = phasei()*phasei().thermo().psi();
    mu_ = phasei()*phasei().thermo().mu();
    alpha_ = phasei()*phasei().thermo().alpha();

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        psi_ += phasei()*phasei().thermo().psi();
        mu_ += phasei()*phasei().thermo().mu();
        alpha_ += phasei()*phasei().thermo().alpha();
    }
}


void Foam::multiphasePhaseChangeMixtureThermo::correctRho(const volScalarField& dp)
{
	// Implementing minimum rho value. Otherwise, negative densities are allowed
	dimensionedScalar minDensity("0.01", dimDensity, 0.01);

    forAllIter(PtrDictionary<phaseModel>, phases_, phasei)
    {

    	forAll(phasei(), celli)
		{
    		phasei().thermo().rho()[celli] +=  phasei().thermo().psi()[celli]*dp[celli];
    		if (phasei().thermo().rho()[celli] < minDensity.value())
    		{
    			phasei().thermo().rho()[celli] = minDensity.value();
    		}
		}

    }
}


bool Foam::multiphasePhaseChangeMixtureThermo::incompressible() const
{
    bool ico = true;

    forAllConstIter(PtrDictionary<phaseModel>, phases_, phase)
    {
        ico &= phase().thermo().incompressible();
    }

    return ico;
}


bool Foam::multiphasePhaseChangeMixtureThermo::isochoric() const
{
    bool iso = true;

    forAllConstIter(PtrDictionary<phaseModel>, phases_, phase)
    {
        iso &= phase().thermo().incompressible();
    }

    return iso;
}


Foam::tmp<Foam::volScalarField> Foam::multiphasePhaseChangeMixtureThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<volScalarField> the(phasei()*phasei().thermo().he(p, T));

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        the.ref() += phasei()*phasei().thermo().he(p, T);
    }

    return the;
}


Foam::tmp<Foam::scalarField> Foam::multiphasePhaseChangeMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<scalarField> the
    (
        scalarField(phasei(), cells)*phasei().thermo().he(p, T, cells)
    );

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        the.ref() +=
            scalarField(phasei(), cells)*phasei().thermo().he(p, T, cells);
    }

    return the;
}


Foam::tmp<Foam::scalarField> Foam::multiphasePhaseChangeMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<scalarField> the
    (
        phasei().boundaryField()[patchi]*phasei().thermo().he(p, T, patchi)
    );

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        the.ref() +=
            phasei().boundaryField()[patchi]*phasei().thermo().he(p, T, patchi);
    }

    return the;
}


Foam::tmp<Foam::volScalarField> Foam::multiphasePhaseChangeMixtureThermo::hc() const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<volScalarField> thc(phasei()*phasei().thermo().hc());

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        thc.ref() += phasei()*phasei().thermo().hc();
    }

    return thc;
}


Foam::tmp<Foam::scalarField> Foam::multiphasePhaseChangeMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::scalarField> Foam::multiphasePhaseChangeMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::volScalarField> Foam::multiphasePhaseChangeMixtureThermo::rho() const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<volScalarField> trho(phasei()*phasei().thermo().rho());

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        trho.ref() += phasei()*phasei().thermo().rho();
    }

    return trho;
}


Foam::tmp<Foam::scalarField> Foam::multiphasePhaseChangeMixtureThermo::rho
(
    const label patchi
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<scalarField> trho
    (
        phasei().boundaryField()[patchi]*phasei().thermo().rho(patchi)
    );

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        trho.ref() +=
            phasei().boundaryField()[patchi]*phasei().thermo().rho(patchi);
    }

    return trho;
}


Foam::tmp<Foam::volScalarField> Foam::multiphasePhaseChangeMixtureThermo::Cp() const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<volScalarField> tCp(phasei()*phasei().thermo().Cp());

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        tCp.ref() += phasei()*phasei().thermo().Cp();
    }

    return tCp;
}


Foam::tmp<Foam::scalarField> Foam::multiphasePhaseChangeMixtureThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<scalarField> tCp
    (
        phasei().boundaryField()[patchi]*phasei().thermo().Cp(p, T, patchi)
    );

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        tCp.ref() +=
            phasei().boundaryField()[patchi]*phasei().thermo().Cp(p, T, patchi);
    }

    return tCp;
}


Foam::tmp<Foam::volScalarField> Foam::multiphasePhaseChangeMixtureThermo::Cv() const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<volScalarField> tCv(phasei()*phasei().thermo().Cv());

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        tCv.ref() += phasei()*phasei().thermo().Cv();
    }

    return tCv;
}


Foam::tmp<Foam::scalarField> Foam::multiphasePhaseChangeMixtureThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<scalarField> tCv
    (
        phasei().boundaryField()[patchi]*phasei().thermo().Cv(p, T, patchi)
    );

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        tCv.ref() +=
            phasei().boundaryField()[patchi]*phasei().thermo().Cv(p, T, patchi);
    }

    return tCv;
}


Foam::tmp<Foam::volScalarField> Foam::multiphasePhaseChangeMixtureThermo::gamma() const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<volScalarField> tgamma(phasei()*phasei().thermo().gamma());

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        tgamma.ref() += phasei()*phasei().thermo().gamma();
    }

    return tgamma;
}


Foam::tmp<Foam::scalarField> Foam::multiphasePhaseChangeMixtureThermo::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<scalarField> tgamma
    (
        phasei().boundaryField()[patchi]*phasei().thermo().gamma(p, T, patchi)
    );

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        tgamma.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().gamma(p, T, patchi);
    }

    return tgamma;
}


Foam::tmp<Foam::volScalarField> Foam::multiphasePhaseChangeMixtureThermo::Cpv() const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<volScalarField> tCpv(phasei()*phasei().thermo().Cpv());

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        tCpv.ref() += phasei()*phasei().thermo().Cpv();
    }

    return tCpv;
}


Foam::tmp<Foam::scalarField> Foam::multiphasePhaseChangeMixtureThermo::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<scalarField> tCpv
    (
        phasei().boundaryField()[patchi]*phasei().thermo().Cpv(p, T, patchi)
    );

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        tCpv.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().Cpv(p, T, patchi);
    }

    return tCpv;
}


Foam::tmp<Foam::volScalarField> Foam::multiphasePhaseChangeMixtureThermo::CpByCpv() const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<volScalarField> tCpByCpv(phasei()*phasei().thermo().CpByCpv());

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        tCpByCpv.ref() += phasei()*phasei().thermo().CpByCpv();
    }

    return tCpByCpv;
}


Foam::tmp<Foam::scalarField> Foam::multiphasePhaseChangeMixtureThermo::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<scalarField> tCpByCpv
    (
        phasei().boundaryField()[patchi]*phasei().thermo().CpByCpv(p, T, patchi)
    );

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        tCpByCpv.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().CpByCpv(p, T, patchi);
    }

    return tCpByCpv;
}


Foam::tmp<Foam::volScalarField> Foam::multiphasePhaseChangeMixtureThermo::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::scalarField> Foam::multiphasePhaseChangeMixtureThermo::nu
(
    const label patchi
) const
{
    return mu(patchi)/rho(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::multiphasePhaseChangeMixtureThermo::kappa() const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<volScalarField> tkappa(phasei()*phasei().thermo().kappa());

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        tkappa.ref() += phasei()*phasei().thermo().kappa();
    }

    return tkappa;
}


Foam::tmp<Foam::scalarField> Foam::multiphasePhaseChangeMixtureThermo::kappa
(
    const label patchi
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<scalarField> tkappa
    (
        phasei().boundaryField()[patchi]*phasei().thermo().kappa(patchi)
    );

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        tkappa.ref() +=
            phasei().boundaryField()[patchi]*phasei().thermo().kappa(patchi);
    }

    return tkappa;
}


Foam::tmp<Foam::volScalarField> Foam::multiphasePhaseChangeMixtureThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<volScalarField> tkappaEff(phasei()*phasei().thermo().kappaEff(alphat));

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        tkappaEff.ref() += phasei()*phasei().thermo().kappaEff(alphat);
    }

    return tkappaEff;
}


Foam::tmp<Foam::scalarField> Foam::multiphasePhaseChangeMixtureThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<scalarField> tkappaEff
    (
        phasei().boundaryField()[patchi]
       *phasei().thermo().kappaEff(alphat, patchi)
    );

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        tkappaEff.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().kappaEff(alphat, patchi);
    }

    return tkappaEff;
}


Foam::tmp<Foam::volScalarField> Foam::multiphasePhaseChangeMixtureThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<volScalarField> talphaEff(phasei()*phasei().thermo().alphaEff(alphat));

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        talphaEff.ref() += phasei()*phasei().thermo().alphaEff(alphat);
    }

    return talphaEff;
}


Foam::tmp<Foam::scalarField> Foam::multiphasePhaseChangeMixtureThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<scalarField> talphaEff
    (
        phasei().boundaryField()[patchi]
       *phasei().thermo().alphaEff(alphat, patchi)
    );

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        talphaEff.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().alphaEff(alphat, patchi);
    }

    return talphaEff;
}


Foam::tmp<Foam::volScalarField> Foam::multiphasePhaseChangeMixtureThermo::rCv() const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<volScalarField> trCv(phasei()/phasei().thermo().Cv());

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        trCv.ref() += phasei()/phasei().thermo().Cv();
    }

    return trCv;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphasePhaseChangeMixtureThermo::surfaceTensionForce() const
{
    tmp<surfaceScalarField> tstf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "surfaceTensionForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "surfaceTensionForce",
                dimensionSet(1, -2, -2, 0, 0),
                0.0
            )
        )
    );

    surfaceScalarField& stf = tstf.ref();

    forAllConstIter(PtrDictionary<phaseModel>, phases_, phase1)
    {
        const phaseModel& alpha1 = phase1();

        PtrDictionary<phaseModel>::const_iterator phase2 = phase1;
        ++phase2;

        for (; phase2 != phases_.end(); ++phase2)
        {
            const phaseModel& alpha2 = phase2();

            sigmaTable::const_iterator sigma =
                sigmas_.find(interfacePair(alpha1, alpha2));

            if (sigma == sigmas_.end())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << " in list of sigma values"
                    << exit(FatalError);
            }

            stf += dimensionedScalar("sigma", dimSigma_, sigma())
               *fvc::interpolate(K(alpha1, alpha2))*
                (
                    fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                  - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                );
        }
    }

    return tstf;
}


void Foam::multiphasePhaseChangeMixtureThermo::solve()
{
    const Time& runTime = mesh_.time();

    const dictionary& alphaControls = mesh_.solverDict("alpha");
    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
    scalar cAlpha(readScalar(alphaControls.lookup("cAlpha")));

    volScalarField& alpha = phases_.first();

    if (nAlphaSubCycles > 1)
    {
        surfaceScalarField rhoPhiSum(0.0*rhoPhi_);
        dimensionedScalar totalDeltaT = runTime.deltaT();

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            solveAlphas(cAlpha);
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi_;
        }

        rhoPhi_ = rhoPhiSum;
    }
    else
    {
        solveAlphas(cAlpha);
    }
}


Foam::tmp<Foam::surfaceVectorField> Foam::multiphasePhaseChangeMixtureThermo::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    /*
    // Cell gradient of alpha
    volVectorField gradAlpha =
        alpha2*fvc::grad(alpha1) - alpha1*fvc::grad(alpha2);

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf = fvc::interpolate(gradAlpha);
    */

    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN_);
}


Foam::tmp<Foam::surfaceScalarField> Foam::multiphasePhaseChangeMixtureThermo::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::multiphasePhaseChangeMixtureThermo::correctContactAngle
(
    const phaseModel& alpha1,
    const phaseModel& alpha2,
    surfaceVectorField::Boundary& nHatb
) const
{
    const volScalarField::Boundary& gbf
        = alpha1.boundaryField();

    const fvBoundaryMesh& boundary = mesh_.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(gbf[patchi]))
        {
            const alphaContactAngleFvPatchScalarField& acap =
                refCast<const alphaContactAngleFvPatchScalarField>(gbf[patchi]);

            vectorField& nHatPatch = nHatb[patchi];

            vectorField AfHatPatch
            (
                mesh_.Sf().boundaryField()[patchi]
               /mesh_.magSf().boundaryField()[patchi]
            );

            alphaContactAngleFvPatchScalarField::thetaPropsTable::
                const_iterator tp =
                acap.thetaProps().find(interfacePair(alpha1, alpha2));

            if (tp == acap.thetaProps().end())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << "\n    in table of theta properties for patch "
                    << acap.patch().name()
                    << exit(FatalError);
            }

            bool matched = (tp.key().first() == alpha1.name());

            scalar theta0 = convertToRad*tp().theta0(matched);
            scalarField theta(boundary[patchi].size(), theta0);

            scalar uTheta = tp().uTheta();

            // Calculate the dynamic contact angle if required
            if (uTheta > SMALL)
            {
                scalar thetaA = convertToRad*tp().thetaA(matched);
                scalar thetaR = convertToRad*tp().thetaR(matched);

                // Calculated the component of the velocity parallel to the wall
                vectorField Uwall
                (
                    U_.boundaryField()[patchi].patchInternalField()
                  - U_.boundaryField()[patchi]
                );
                Uwall -= (AfHatPatch & Uwall)*AfHatPatch;

                // Find the direction of the interface parallel to the wall
                vectorField nWall
                (
                    nHatPatch - (AfHatPatch & nHatPatch)*AfHatPatch
                );

                // Normalise nWall
                nWall /= (mag(nWall) + SMALL);

                // Calculate Uwall resolved normal to the interface parallel to
                // the interface
                scalarField uwall(nWall & Uwall);

                theta += (thetaA - thetaR)*tanh(uwall/uTheta);
            }


            // Reset nHatPatch to correspond to the contact angle

            scalarField a12(nHatPatch & AfHatPatch);

            scalarField b1(cos(theta));

            scalarField b2(nHatPatch.size());

            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatPatch = a*AfHatPatch + b*nHatPatch;

            nHatPatch /= (mag(nHatPatch) + deltaN_.value());
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::multiphasePhaseChangeMixtureThermo::K
(
    const phaseModel& alpha1,
    const phaseModel& alpha2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(alpha1, alpha2);

    correctContactAngle(alpha1, alpha2, tnHatfv.ref().boundaryFieldRef());

    // Simple expression for curvature
    return -fvc::div(tnHatfv & mesh_.Sf());
}


Foam::tmp<Foam::volScalarField>
Foam::multiphasePhaseChangeMixtureThermo::nearInterface() const
{
    tmp<volScalarField> tnearInt
    (
        new volScalarField
        (
            IOobject
            (
                "nearInterface",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("nearInterface", dimless, 0.0)
        )
    );

    forAllConstIter(PtrDictionary<phaseModel>, phases_, phase)
    {
        tnearInt.ref() =
            max(tnearInt(), pos(phase() - 0.01)*pos(0.99 - phase()));
    }

    return tnearInt;
}


void Foam::multiphasePhaseChangeMixtureThermo::solveAlphas
(
    const scalar cAlpha
)
{
    static label nSolves=-1;
    nSolves++;

    // - Defining the schemes to be used for the bounded flux (alphaScheme) and the unbounded flux (alpharScheme)
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");


    surfaceScalarField phic(mag(phi_/mesh_.magSf()));				// creating compression flux
    phic = min(cAlpha*phic, max(phic));								// setting maximum values of the compression flux

    PtrList<surfaceScalarField> alphaPhiCorrs(phases_.size());
    int phasei = 0;


    //- /////////////////////////////////////////////////////////////////
	//  Building the mass transfer coefficients for the current iteration
	//- /////////////////////////////////////////////////////////////////


    //- Creating the mass transfer values for the water phase
	Pair<tmp<volScalarField>> vDotAlphaW =
			cavitationModel_->vDotAlphaW();
	// volume condesation rate
	// water volume transfer rate production
	const volScalarField& vDotProdAlphaW = vDotAlphaW[0]();
	// volume evaporation rate
	// water volume transfer rate desctruction
	const volScalarField& vDotDestAlphaW = vDotAlphaW[1]();
	//- volume evaporation rate minus the volume condensation rate
	//  water volume transfer rate, destruction minus production term
	//  used for implicit source
	const volScalarField vDotDmPAlphaW(vDotDestAlphaW - vDotProdAlphaW);

	//- Creating the mass transfer values for the water phase
	Pair<tmp<volScalarField>> vDotAlphaV =
			cavitationModel_->vDotAlphaV();
	//- vapor volume transfer rate production
	const volScalarField& vDotProdAlphaV = vDotAlphaV[0]();
	//- vapor volume transfer rate desctruction
	const volScalarField& vDotDestAlphaV = vDotAlphaV[1]();
	//- vapor volume transfer rate, destruction minus prodution term
	//  used for implicit source
	const volScalarField vDotDmPAlphaV(vDotDestAlphaV - vDotProdAlphaV);


    //- Forward declaration of variables to store the mass transfer rates calculated from the water phase
    //- For re-use of the vapor phase
    tmp<volScalarField> tmpWaterCRate;
    tmp<volScalarField> tmpWaterVRate;

    //- Field to store the implicit and explicit source terms due to the mass transfer for phase alpha
    tmp<volScalarField> implWaterSource;
    tmp<volScalarField> explWaterSource;
    tmp<volScalarField> implVaporSource;
    tmp<volScalarField> explVaporSource;

    //- /////////////////////////////////////////////////////////////////
    //  Source term validation fields
    //- /////////////////////////////////////////////////////////////////

    //- Calculating the source terms for the phases
    	//- water
    	//- vapor
    forAllIter(PtrDictionary<phaseModel>, phases_, phase)
    {
    	phaseModel& alpha = phase();
    	if (alpha.name() == "water")
    	{
    		// Initializing fields
    		implWaterSource = vDotProdAlphaW * 0;
    		explWaterSource = vDotProdAlphaW * 0;

    		volScalarField& Sp = implWaterSource.ref();
			volScalarField& Su = explWaterSource.ref();

			volScalarField limitedAlpha(min(max(alpha, scalar(0)), scalar(1)));

        	// Implicit source term
            forAll(vDotDmPAlphaW, celli)
            {

            		Sp[celli] += vDotDmPAlphaW[celli]*limitedAlpha[celli];
                	if (vDotDmPAlphaW[celli] > 0.0)
                    {
					Info << "Error: The implicit source term of the water phase is positive: " << vDotDmPAlphaW[celli] << "\n";
					Info << "Error: Water phase fraction: " << limitedAlpha[celli] << "\n";
					Info << "Error: vDotDestAlphaW: " << vDotDestAlphaW[celli] << "\n";
					Info << "Error: vDotDestAlphaV: " << vDotDestAlphaV[celli] << "\n";
					Info << "Error: vDotProdAlphaW: " << vDotProdAlphaW[celli] << "\n";
					Info << "Error: vDotProdAlphaV: " << vDotProdAlphaV[celli] << "\n";
				}
            }

            //- Explicit source term
            forAll(vDotProdAlphaW, celli)
            {

            		Su[celli] += vDotProdAlphaW[celli];
        	if (vDotProdAlphaW[celli] < 0.0)
			{
					Info << "Error: The explicit source term of the water phase is negative: " << vDotProdAlphaW[celli] << "\n";
					Info << "Error: Water phase fraction: " << limitedAlpha[celli] << "\n";
					Info << "Error: vDotDestAlphaW: " << vDotDestAlphaW[celli] << "\n";
					Info << "Error: vDotDestAlphaV: " << vDotDestAlphaV[celli] << "\n";
					Info << "Error: vDotProdAlphaW: " << vDotProdAlphaW[celli] << "\n";
					Info << "Error: vDotProdAlphaV: " << vDotProdAlphaV[celli] << "\n";
				}
            }
    	}

    	if (alpha.name() == "vapor")
    	{
    		// Initializing fields
    		implVaporSource = vDotProdAlphaV * 0;
    		explVaporSource = vDotProdAlphaV * 0;

    		volScalarField& Sp = implVaporSource.ref();
			volScalarField& Su = explVaporSource.ref();

			volScalarField limitedAlpha(min(max(alpha, scalar(0)), scalar(1)));

        	//- Implicit source term
            forAll(vDotDmPAlphaV, celli)
            {

            			Sp[celli] += vDotDmPAlphaV[celli]*limitedAlpha[celli];
            		if (vDotDmPAlphaV[celli] > 0.0) {
    					Info << "Error: The implicit source term of the vapor phase is positive: " << vDotProdAlphaV[celli] << "\n";
    					Info << "Error: Vapor phase fraction: " << limitedAlpha[celli] << "\n";
    					Info << "Error: vDotDestAlphaW: " << vDotDestAlphaW[celli] << "\n";
    					Info << "Error: vDotDestAlphaV: " << vDotDestAlphaV[celli] << "\n";
    					Info << "Error: vDotProdAlphaW: " << vDotProdAlphaW[celli] << "\n";
    					Info << "Error: vDotProdAlphaV: " << vDotProdAlphaV[celli] << "\n";
            		}
            }

            //- Explicit source term
            forAll(vDotProdAlphaV, celli)
            {

    				Su[celli] += vDotProdAlphaV[celli];
            	if (vDotProdAlphaV[celli] < 0.0)
            	{
					Info << "Error: The expicit source term of the vapor phase is negative: " << vDotProdAlphaV[celli] << "\n";
					Info << "Error: Vapor phase fraction: " << limitedAlpha[celli] << "\n";
					Info << "Error: vDotDestAlphaW: " << vDotDestAlphaW[celli] << "\n";
					Info << "Error: vDotDestAlphaV: " << vDotDestAlphaV[celli] << "\n";
					Info << "Error: vDotProdAlphaW: " << vDotProdAlphaW[celli] << "\n";
					Info << "Error: vDotProdAlphaV: " << vDotProdAlphaV[celli] << "\n";
            	}
            }
    	}
    }

    //- /////////////////////////////////////////////////////////////////
	//  Build the limited fluxes of each phase
	//- /////////////////////////////////////////////////////////////////


    forAllIter(PtrDictionary<phaseModel>, phases_, phase)
    {
        phaseModel& alpha = phase();

        // bounded upwind flux of phase alpha
        alphaPhiCorrs.set
        (
            phasei,
            new surfaceScalarField
            (
                phi_.name() + alpha.name(),
                fvc::flux
                (
                    phi_,
                    alpha,
                    alphaScheme
                )
            )
        );

        // Reference to the corrector flux of phase alpha
        surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phasei];


        forAllIter(PtrDictionary<phaseModel>, phases_, phase2)
        {
            phaseModel& alpha2 = phase2();
            if (&alpha2 == &alpha) continue;	// Skip the currently calculated phase alpha for this iterator

            surfaceScalarField phir(phic*nHatf(alpha, alpha2));		// compression flux for the phases alpha and alpha2

            // substract the compression flux of alpha,alpha2 from the corrector flux of alpha
            alphaPhiCorr += fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha,
                alpharScheme
            );
        }

        if (alpha.name() == "water")
        {

        	volScalarField& implSource = implWaterSource.ref();
        	volScalarField& explSource = explWaterSource.ref();

            MULES::limit
            (
                1.0/mesh_.time().deltaT().value(),
                geometricOneField(),
                alpha,
                phi_,
                alphaPhiCorr,
				implSource,
				explSource,
                1,
                0,
                true
            );

        }

        if (alpha.name() == "vapor")
        {
        	volScalarField& implSource = implVaporSource.ref();
			volScalarField& explSource = explVaporSource.ref();

            MULES::limit
            (
                1.0/mesh_.time().deltaT().value(),
                geometricOneField(),
                alpha,
                phi_,
                alphaPhiCorr,
				implSource,
				explSource,
                1,
                0,
                true
            );
        }


//            MULES::limit
//            (
//                1.0/mesh_.time().deltaT().value(),
//                geometricOneField(),
//                alpha,
//                phi_,
//                alphaPhiCorr,
//    			zeroField(),
//    			zeroField(),
//                1,
//                0,
//                true
//            );

        //- /////////////////////////////////////////////////////////////////
    	//- END - Build the limited fluxes of each phase
    	//- /////////////////////////////////////////////////////////////////
        phasei++;
    }


	//- /////////////////////////////////////////////////////////////////
	// Limits the fluxes between the phases
	// http://cpp.openfoam.org/v4/a05954_source.html
    //- /////////////////////////////////////////////////////////////////
    MULES::limitSum(alphaPhiCorrs);


    rhoPhi_ = dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0);	// initializing the fluxes

    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("sumAlpha", dimless, 0)
    );


    volScalarField divU(fvc::div(fvc::absolute(phi_, U_)));				// Calcs the divergence of the velocity field
    phasei = 0;


	//- ///////////////////////////////////////////////////////////////////////////////
	//- Solve-loop for all phases by explicitly iterating over the phases_ dictionary
	//- ///////////////////////////////////////////////////////////////////////////////


    forAllIter(PtrDictionary<phaseModel>, phases_, phase)
    {
        phaseModel& alpha = phase();									// Rename the phase given by the iterator to alpha
        surfaceScalarField& alphaPhi = alphaPhiCorrs[phasei];			// Rename the corrector flux of alpha to alphaPhi
        alphaPhi += upwind<scalar>(mesh_, phi_).flux(alpha);			// Add the bounded upwind fluxes to the corrected flux of each phase


    	//- ///////////////////////////////////////////////////////////////////////////////
    	//  Define the implicit and explicit source on the internal field for each phase
    	//- ///////////////////////////////////////////////////////////////////////////////

        //- Implicit source term
        volScalarField Sp
        (
            IOobject
            (
                "Sp",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("Sp", alpha.dgdt().dimensions(), 0.0)
        );

        //- Explicit source term
        volScalarField Su
        (
            IOobject
            (
                "Su",
                mesh_.time().timeName(),
                mesh_
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            divU*min(alpha, scalar(1))
        );

        //- Adds the divergence of the phase dgdt either to the implicit or the explicit source,
        //  depending on the sign
        //  Initially,
        	// Sp (impl) is zero and
        	// Su (expl) is divU*min(alpha, 1)

        {
            const volScalarField::Internal& dgdt = alpha.dgdt().internalField();
            volScalarField limitedAlpha(min(max(alpha, scalar(0)), scalar(1)));

            forAll(dgdt, celli)
            {
                if (dgdt[celli] < 0.0 && limitedAlpha[celli] > 0.0)
                {
                    Sp[celli] += dgdt[celli]*limitedAlpha[celli];
                    Su[celli] -= dgdt[celli]*limitedAlpha[celli];
                }
                else if (dgdt[celli] > 0.0 && limitedAlpha[celli] < 1.0)
                {
                    Sp[celli] -= dgdt[celli]*(1.0 - limitedAlpha[celli]);
                    // In total, this yields:
                    // Su + Sp 	= divU*alpha - divU * (1 - alpha)
                    //			= positive total div of alpha is in source - pos div of all other phases is
                }
            }
        }


        //- Adds or substracts the divergence of all the other phases dgdt either to the implicit or the explicit source,
        //  depending on the sign
        forAllConstIter(PtrDictionary<phaseModel>, phases_, phase2)
        {
            const phaseModel& alpha2 = phase2();
            volScalarField limitedAlpha2(min(max(alpha2, scalar(0)), scalar(1)));

            if (&alpha2 == &alpha) continue;

            const scalarField& dgdt2 = alpha2.dgdt();

            forAll(dgdt2, celli)
            {
                if (dgdt2[celli] > 0.0 && limitedAlpha2[celli] < 1.0)
                {
                    Sp[celli] -= dgdt2[celli]*(1.0 - limitedAlpha2[celli]);
                    Su[celli] += dgdt2[celli]*alpha[celli];
                }
                else if (dgdt2[celli] < 0.0 && limitedAlpha2[celli] > 0.0)
                {
                    Sp[celli] += dgdt2[celli]*limitedAlpha2[celli];
                }
            }
        }

        //- /////////////////////////////////////////
		//- Adding the volumetric mass transfer rates to the source terms
		//- /////////////////////////////////////////


        if (alpha.name() == "water")
        {

        	volScalarField& implSource = implWaterSource.ref();
        	volScalarField& explSource = explWaterSource.ref();


        	forAll(implSource, celli)
			{
        		Sp[celli] += implSource[celli];
			}

        	forAll(explSource, celli)
        	{
        		Su[celli] += explSource[celli];
        	}
        }

        if (alpha.name() == "vapor")
        {
        	volScalarField& implSource = implVaporSource.ref();
			volScalarField& explSource = explVaporSource.ref();

        	forAll(implSource, celli)
			{
        		Sp[celli] += implSource[celli];
			}

        	forAll(explSource, celli)
        	{
        		Su[celli] += explSource[celli];
        	}
        }


        //- /////////////////////////////////////////
		//- Solving the transport equation for phase alpha explicitly
		//- /////////////////////////////////////////


        MULES::explicitSolve
        (
            geometricOneField(),
            alpha,
            alphaPhi,
            Sp,
            Su
        );

        // Further limiting the phase values to increase stability
        forAllIter(PtrDictionary<phaseModel>, phases_, phase)
        {
        	forAll(phase(), celli)
			{
        		if (phase()[celli] < 0.0)
        		{
        			phase()[celli] = 0.0;
        		}
        		if (phase()[celli] > 1.0)
        		{
        			phase()[celli] = 1.0;
        		}
			}
        }

        // Adds the flux of the current phase to the total flux
        rhoPhi_ += fvc::interpolate(alpha.thermo().rho())*alphaPhi;

        Info<< alpha.name() << " volume fraction, min, max = "
            << alpha.weightedAverage(mesh_.V()).value()
            << ' ' << min(alpha).value()
            << ' ' << max(alpha).value()
            << endl;

        // Adds the alpha values of all phases together
        sumAlpha += alpha;

        phasei++;
    }

    Info<< "Phase-sum volume fraction, min, max = "
        << sumAlpha.weightedAverage(mesh_.V()).value()
        << ' ' << min(sumAlpha).value()
        << ' ' << max(sumAlpha).value()
        << endl;

    calcAlphas();


    alphaSum_ *= 0.0;
    forAllIter(PtrDictionary<phaseModel>, phases_, phasei)
    {
    	alphaSum_ += phasei();
    }
    alphaSum_ -= 1;
}


// ************************************************************************* //
