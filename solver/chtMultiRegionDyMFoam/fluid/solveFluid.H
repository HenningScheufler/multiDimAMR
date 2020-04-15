if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

if (frozenFlow)
{
    #include "EEqn.H"
}
else
{
    bool correctPhi = pimple.lookupOrDefault
    (
        "correctPhi",
        correctPhi
    );

    // Store divrhoU from the previous mesh so that it can be mapped
    // and used in correctPhi to ensure the corrected phi has the
    // same divergence
    volScalarField divrhoU("divrhoU",fvc::div(fvc::absolute(phi, rho, U)));

    if (oCorr == 0)
    {
        // Store momentum to set rhoUf for introduced faces.
        volVectorField rhoU("rhoU", rho*U);

        // Do any mesh changes
        mesh.update();

        if (mesh.changing())
        {
            MRF.update();

            if (correctPhi)
            {
                // Calculate absolute flux
                // from the mapped surface velocity
                phi = mesh.Sf() & rhoUf;

                #include "correctPhi.H"

                // Make the fluxes relative to the mesh-motion
                fvc::makeRelative(phi, rho, U);
            }

            // if (checkMeshCourantNo)
            // {
            //     #include "meshCourantNo.H"
            // }
        }
    }

    if (oCorr == 0)
    {
        #include "rhoEqn.H"
    }

    #include "UEqn.H"
    #include "YEqn.H"
    #include "EEqn.H"

    // --- PISO loop
    for (int corr=0; corr<nCorr; corr++)
    {
        #include "pEqn.H"
    }

    turbulence.correct();

    rho = thermo.rho();
}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
