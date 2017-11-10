{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    PtrList<fvScalarMatrix> p_rghEqnComps(mixture.phases().size());

    label phasei = 0;
    forAllConstIter
    (
        PtrDictionary<phaseModel>,
        mixture.phases(),
        phase
    )
    {
        const rhoThermo& thermo = phase().thermo();
        const volScalarField& rho = thermo.rho()();

        p_rghEqnComps.set
        (
            phasei,
            (
                fvc::ddt(rho) + thermo.psi()*correction(fvm::ddt(p_rgh))
              + fvc::div(phi, rho) - fvc::Sp(fvc::div(phi), rho)
            ).ptr()
        );

        phasei++;
    }

    //Coefficients of the mass transfer
    //Info << "calling mixture.cavitationModel()->pSat() now \n";
    const dimensionedScalar pSat = mixture.cavitationModel()->pSat();
    volScalarField pSatField(
							IOobject
							(
								"pSatField",
								mesh.time().timeName(),
								mesh,
								IOobject::NO_READ,
								IOobject::NO_WRITE
							),
							mesh,
							dimensionedScalar("pSatField", pSat.dimensions(), 0.0));
    pSatField.dimensions().reset(pSat.dimensions());
    pSatField = mixture.cavitationModel()->pSat();


    //Info << "calling mixture.cavitationModel()->vDotP() now \n";
    Pair<tmp<volScalarField>> vDotP = mixture.cavitationModel()->vDotP();
    //Info << "vDotP[0]() = " << vDotP[0]() << "\n";
    //Info << "vDotP[1]() = " << vDotP[1]() << "\n";
    const volScalarField& vDotcP = vDotP[0]();
    const volScalarField& vDotvP = vDotP[1]();


    // Cache p_rgh prior to solve for density update
    volScalarField p_rgh_0(p_rgh);
//
    while (pimple.correctNonOrthogonal())
    {

 		fvScalarMatrix p_rghEqnIncomp
		(
			fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh)
		  - (vDotvP - vDotcP)*(pSat- rho*gh)
		  + fvm::Sp(vDotvP - vDotcP, p_rgh)
		);

        tmp<fvScalarMatrix> p_rghEqnComp;

        phasei = 0;
        forAllConstIter
        (
            PtrDictionary<phaseModel>,
            mixture.phases(),
            phase
        )
        {
            tmp<fvScalarMatrix> hmm
            (
                (max(phase(), scalar(0))/phase().thermo().rho())
               *p_rghEqnComps[phasei]
            );

            if (phasei == 0)
            {
                p_rghEqnComp = hmm;
            }
            else
            {
                p_rghEqnComp.ref() += hmm;
            }

            phasei++;
        }

        solve
        (
            p_rghEqnComp
          + p_rghEqnIncomp,
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phasei = 0;
            forAllIter
            (
                PtrDictionary<phaseModel>,
                mixture.phases(),
                phase
            )
            {
				phase().dgdt() =
						pos(phase())*
					  (p_rghEqnComps[phasei] & p_rgh)/phase().thermo().rho();
            }

            phi = phiHbyA + p_rghEqnIncomp.flux();

            U = HbyA
              + rAU*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
            U.correctBoundaryConditions();
        }
    }

    p = max(p_rgh + mixture.rho()*gh, pMin);

    mixture.correctRho(p_rgh - p_rgh_0);
    rho = mixture.rho();

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();

    K = 0.5*magSqr(U);

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(p_rgh) " << min(p_rgh).value() << endl;
}
