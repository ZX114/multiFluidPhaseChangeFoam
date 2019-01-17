/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "multiFluidSystem.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "fixedValueFvsPatchFields.H"
#include "Time.H"
#include "subCycle.H"
#include "MULES.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcAverage.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::multiFluidSystem::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiFluidSystem::calcAlphas()
{
    scalar level = 0.0;
    alphas_ == 0.0;

    forAllIter(PtrDictionary<phaseModel>, phases_, iter)
    {
        alphas_ += level*iter();
        level += 1.0;
    }
}


void Foam::multiFluidSystem::calcCAlpha()
{
    cAlpha_ == 0.0;

    forAll(cAlpha_.primitiveField(),cellI)
    {
        cAlpha_.primitiveFieldRef()[cellI] = 0.0;
    }
}


void Foam::multiFluidSystem::solveAlphas()
{
    PtrList<surfaceScalarField> alphaPhiCorrs(phases_.size());
    int phasei = 0;

    forAllIter(PtrDictionary<phaseModel>, phases_, iter)
    {
        //if (phasei == (phases_.size()-1)) break;
        
        phaseModel& phase = iter();
        volScalarField& alpha1 = phase;

        alphaPhiCorrs.set
        (
            phasei,
            new surfaceScalarField
            (
                "phi" + alpha1.name() + "Corr",
                fvc::flux
                (
                    phi_,
                    phase,
                    "div(phi," + alpha1.name() + ')'
                )
            )
        );

        surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phasei];

        forAllIter(PtrDictionary<phaseModel>, phases_, iter2)
        {
           
            phaseModel& phase2 = iter2();
            volScalarField& alpha2 = phase2;

            if (&phase2 == &phase) continue;

            surfaceScalarField phir(phase.phi() - phase2.phi());

            /*
            scalarCoeffSymmTable::const_iterator cAlpha
            (
                cAlphas_.find(interfacePair(phase, phase2))
            );

            if (cAlpha != cAlphas_.end())
            {
                surfaceScalarField phic
                (
                    (mag(phi_) + mag(phir))/mesh_.magSf()
                );

                phir += min(cAlpha()*phic, max(phic))*nHatf(phase, phase2);
            }
            */
            if (phasei != (phases_.size()-1))
            {
                surfaceScalarField phic
                (
                    (mag(phi_) + mag(phir))/mesh_.magSf()
                );
                            
                surfaceScalarField gamma = cAlphaGamma(phase, phase2);   

                //cAlphaSwitch() returns a binary value, which is then multiplied by
                //compressionCoeff to give the right interface compression value
                
                phir += min(compressionCoeff*cAlphaSwitch(gamma)*phic, max(phic))
                        *nHatf(phase, phase2);

                volScalarField gammav = cAlphaGammav(phase, phase2);
                volScalarField cAlphav = cAlphaSwitchv(gammav);

                //multiplied compressionCoeff
                forAll(cAlpha_.primitiveField(),cellI)
                {
                    cAlpha_.primitiveFieldRef()[cellI]
                    = compressionCoeff*cAlphav.primitiveField()[cellI];
                }
            }
            
            word phirScheme
            (
                "div(phir," + alpha2.name() + ',' + alpha1.name() + ')'
            );

            //current alphaPhiCorr = alphaPhiCorrs[phasei]
            alphaPhiCorr += fvc::flux
            (
                -fvc::flux(-phir, phase2, phirScheme),
                phase,
                phirScheme
            );
        }

        const volScalarField Sp(phaseChangeSp(phase, cCoeffs(), vCoeffs()));
        const volScalarField Su(phaseChangeSu(phase, cCoeffs(), vCoeffs()));
        
        phase.correctInflowOutflow(alphaPhiCorr);
        
        MULES::limit
        (
            1.0/mesh_.time().deltaT().value(),
            geometricOneField(),
            phase,
            phi_,
            alphaPhiCorr,
            Sp(),
            Su(),
            //zeroField(),
            //zeroField(),
            oneField(),
            zeroField(),
            true
        );

        phasei++;
    }

    MULES::limitSum(alphaPhiCorrs);

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

    phasei = 0;

    forAllIter(PtrDictionary<phaseModel>, phases_, iter)
    {
    	/*
        if (phasei == (phases_.size()-1))
        {
            phaseModel& lastPhase = iter();
            
            // forAll(lastPhase.primitiveField(),cellI)
            // {
            //     lastPhase.primitiveFieldRef()[cellI] = 
            //     1.0 - sumAlpha.primitiveField()[cellI];
            // }
            // forAll(lastPhase.boundaryField(),patchI)
            // {
            //     lastPhase.boundaryFieldRef()[patchI] = 
            //     1.0 - sumAlpha.boundaryField()[patchI];
            // }

            
            scalarField& lastPRef = lastPhase;
            lastPRef = (oneField() - sumAlpha);
            lastPhase.correctBoundaryConditions();
            
            sumAlpha += lastPhase;
            
            Info<< lastPhase.name() << " volume fraction, min, max = "
            << lastPhase.weightedAverage(mesh_.V()).value()
            << ' ' << min(lastPhase).value()
            << ' ' << max(lastPhase).value()
            << endl;
            break;
        }
        */

        phaseModel& phase = iter();

        surfaceScalarField& alphaPhi = alphaPhiCorrs[phasei];
        alphaPhi += upwind<scalar>(mesh_, phi_).flux(phase);
        phase.correctInflowOutflow(alphaPhi);
        
        const volScalarField Sp(phaseChangeSp(phase, cCoeffs(), vCoeffs()));
        const volScalarField Su(phaseChangeSu(phase, cCoeffs(), vCoeffs()));

        MULES::explicitSolve
        (
            geometricOneField(),
            phase,
            alphaPhi,
            Sp(),
            Su()
        );

        phase.alphaPhi() = alphaPhi;

        Info<< phase.name() << " volume fraction, min, max = "
            << phase.weightedAverage(mesh_.V()).value()
            << ' ' << min(phase).value()
            << ' ' << max(phase).value()
            << endl;

        sumAlpha += phase;

        phasei++;
    }
    
    
    Info<< "Phase-sum volume fraction, min, max = "
        << sumAlpha.weightedAverage(mesh_.V()).value()
        << ' ' << min(sumAlpha).value()
        << ' ' << max(sumAlpha).value()
        << endl;

    // Correct the sum of the phase-fractions to avoid 'drift'
    
    volScalarField sumCorr(1.0 - sumAlpha);
    forAllIter(PtrDictionary<phaseModel>, phases_, iter)
    {
        phaseModel& phase = iter();
        volScalarField& alpha = phase;
        alpha += alpha*sumCorr;
    }
    
    calcAlphas();
}


Foam::tmp<Foam::surfaceVectorField> Foam::multiFluidSystem::nHatfv
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


Foam::tmp<Foam::surfaceScalarField> Foam::multiFluidSystem::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


Foam::tmp<Foam::surfaceScalarField> Foam::multiFluidSystem::cAlphaGamma
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );
    
    return mag(gradAlphaf)/(max(mag(gradAlphaf))+deltaN_);
}



// binary value to determine whether or not
// interface capturing will be used (based on gamma field)
Foam::tmp<Foam::surfaceScalarField> Foam::multiFluidSystem::cAlphaSwitch
(
    surfaceScalarField& gamma
) const
{   
    forAll(gamma.internalField(),faceI)
    {
        if (gamma.internalField()[faceI] < gammaStar)
        {
            gamma.primitiveFieldRef()[faceI] = 0.0;
        }
        else
        {
            gamma.primitiveFieldRef()[faceI] = 1.0;
        }
    }
    
    return gamma;
}


// volScalarField
Foam::tmp<Foam::volScalarField> Foam::multiFluidSystem::cAlphaGammav
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    volVectorField gradAlpha
    (
        alpha2*fvc::grad(alpha1)
      - alpha1*fvc::grad(alpha2)
    );

    return mag(gradAlpha)/(max(mag(gradAlpha))+deltaN_);
}


Foam::tmp<Foam::volScalarField> Foam::multiFluidSystem::cAlphaSwitchv
(
    volScalarField& gammav
) const
{    
    forAll(gammav.internalField(),cellI)
    {
        if (gammav.internalField()[cellI] < gammaStar)
        {
            gammav.primitiveFieldRef()[cellI] = 0.0;
        }
        else
        {
            gammav.primitiveFieldRef()[cellI] = 1.0;
        }
    }

    Info << "compressionCoeff: " << compressionCoeff << endl;
    
    return gammav;
}


// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::multiFluidSystem::correctContactAngle
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    surfaceVectorField::Boundary& nHatb
) const
{
    const volScalarField::Boundary& gbf
        = phase1.boundaryField();

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
                acap.thetaProps().find(interfacePair(phase1, phase2));

            if (tp == acap.thetaProps().end())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(phase1, phase2)
                    << "\n    in table of theta properties for patch "
                    << acap.patch().name()
                    << exit(FatalError);
            }

            bool matched = (tp.key().first() == phase1.name());

            scalar theta0 = convertToRad*tp().theta0(matched);
            scalarField theta(boundary[patchi].size(), theta0);

            scalar uTheta = tp().uTheta();

            // Calculate the dynamic contact angle if required
            if (uTheta > small)
            {
                scalar thetaA = convertToRad*tp().thetaA(matched);
                scalar thetaR = convertToRad*tp().thetaR(matched);

                // Calculated the component of the velocity parallel to the wall
                vectorField Uwall
                (
                    phase1.U().boundaryField()[patchi].patchInternalField()
                  - phase1.U().boundaryField()[patchi]
                );
                Uwall -= (AfHatPatch & Uwall)*AfHatPatch;

                // Find the direction of the interface parallel to the wall
                vectorField nWall
                (
                    nHatPatch - (AfHatPatch & nHatPatch)*AfHatPatch
                );

                // Normalise nWall
                nWall /= (mag(nWall) + small);

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


Foam::tmp<Foam::volScalarField> Foam::multiFluidSystem::K
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(phase1, phase2);

    correctContactAngle(phase1, phase2, tnHatfv.ref().boundaryFieldRef());

    // Simple expression for curvature
    return -fvc::div(tnHatfv & mesh_.Sf());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiFluidSystem::multiFluidSystem
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    phases_(lookup("phases"), phaseModel::iNew(U.mesh())),

    mesh_(U.mesh()),
    phi_(phi),

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
    
    cAlpha_
    (
        IOobject
        (
            "cAlpha",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("cAlpha", dimless, 0.0)
    ),

    // lookup interfaceCompressionCoeff
    compressionCoeff(lookupOrDefault<scalar>("interfaceCompressionCoeff",1.0)),

    // lookup gammaStar
    gammaStar(lookupOrDefault<scalar>("gammaStar", 0.4)),
    
    sigmas_(lookup("sigmas")),
    dimSigma_(1, 0, -2, 0, 0),
    // cAlphas_(lookup("interfaceCompression")),
    Cvms_(lookup("virtualMass")),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    )
{
    calcAlphas();
    alphas_.write();
    
    // initialize and write cAlpha_
    calcCAlpha();
    cAlpha_.write();
    
    interfaceDictTable dragModelsDict(lookup("drag"));

    forAllConstIter(interfaceDictTable, dragModelsDict, iter)
    {
        dragModels_.insert
        (
            iter.key(),
            dragModel::New
            (
                iter(),
                *phases_.lookup(iter.key().first()),
                *phases_.lookup(iter.key().second())
            ).ptr()
        );
    }

    interfaceDictTable phaseChangeModelsDict(lookup("phaseChange"));
    
    forAllConstIter(interfaceDictTable, phaseChangeModelsDict, iter)
    {
        phaseChangeModels_.insert
        (
            iter.key(),
            phaseChangeModel::New
            (
                iter(),
                *phases_.lookup(iter.key().first()),
                *phases_.lookup(iter.key().second())
            ).ptr()
        );
    }

    forAllConstIter(PtrDictionary<phaseModel>, phases_, iter1)
    {
        const phaseModel& phase1 = iter1();

        forAllConstIter(PtrDictionary<phaseModel>, phases_, iter2)
        {
            const phaseModel& phase2 = iter2();

            if (&phase2 != &phase1)
            {
                scalarCoeffSymmTable::const_iterator sigma
                (
                    sigmas_.find(interfacePair(phase1, phase2))
                );

                /*
                if (sigma != sigmas_.end())
                {
                    scalarCoeffSymmTable::const_iterator cAlpha
                    (
                        cAlphas_.find(interfacePair(phase1, phase2))
                    );

                    if (cAlpha == cAlphas_.end())
                    {
                        WarningInFunction
                          << "Compression coefficient not specified for "
                             "phase pair ("
                          << phase1.name() << ' ' << phase2.name()
                          << ") for which a surface tension "
                             "coefficient is specified"
                          << endl;
                    }
                }
                */
            }
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::multiFluidSystem::rho() const
{
    PtrDictionary<phaseModel>::const_iterator iter = phases_.begin();

    tmp<volScalarField> trho = iter()*iter().rho();
    volScalarField& rho = trho.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        rho += iter()*iter().rho();
    }

    return trho;
}


Foam::tmp<Foam::scalarField>
Foam::multiFluidSystem::rho(const label patchi) const
{
    PtrDictionary<phaseModel>::const_iterator iter = phases_.begin();

    tmp<scalarField> trho = iter().boundaryField()[patchi]*iter().rho().value();
    scalarField& rho = trho.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        rho += iter().boundaryField()[patchi]*iter().rho().value();
    }

    return trho;
}


Foam::tmp<Foam::volScalarField> Foam::multiFluidSystem::nu() const
{
    PtrDictionary<phaseModel>::const_iterator iter = phases_.begin();

    tmp<volScalarField> tmu = iter()*(iter().rho()*iter().nu());
    volScalarField& mu = tmu.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        mu += iter()*(iter().rho()*iter().nu());
    }

    return tmu/rho();
}


Foam::tmp<Foam::scalarField>
Foam::multiFluidSystem::nu(const label patchi) const
{
    PtrDictionary<phaseModel>::const_iterator iter = phases_.begin();

    tmp<scalarField> tmu =
        iter().boundaryField()[patchi]
       *(iter().rho().value()*iter().nu().value());
    scalarField& mu = tmu.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        mu +=
            iter().boundaryField()[patchi]
           *(iter().rho().value()*iter().nu().value());
    }

    return tmu/rho(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::multiFluidSystem::Cvm
(
    const phaseModel& phase
) const
{
    tmp<volScalarField> tCvm
    (
        new volScalarField
        (
            IOobject
            (
                "Cvm",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "Cvm",
                dimensionSet(1, -3, 0, 0, 0),
                0
            )
        )
    );

    forAllConstIter(PtrDictionary<phaseModel>, phases_, iter)
    {
        const phaseModel& phase2 = iter();

        if (&phase2 != &phase)
        {
            scalarCoeffTable::const_iterator Cvm
            (
                Cvms_.find(interfacePair(phase, phase2))
            );

            if (Cvm != Cvms_.end())
            {
                tCvm.ref() += Cvm()*phase2.rho()*phase2;
            }
            else
            {
                Cvm = Cvms_.find(interfacePair(phase2, phase));

                if (Cvm != Cvms_.end())
                {
                    tCvm.ref() += Cvm()*phase.rho()*phase2;
                }
            }
        }
    }

    return tCvm;
}


Foam::tmp<Foam::volVectorField> Foam::multiFluidSystem::Svm
(
    const phaseModel& phase
) const
{
    tmp<volVectorField> tSvm
    (
        new volVectorField
        (
            IOobject
            (
                "Svm",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector
            (
                "Svm",
                dimensionSet(1, -2, -2, 0, 0),
                Zero
            )
        )
    );

    forAllConstIter(PtrDictionary<phaseModel>, phases_, iter)
    {
        const phaseModel& phase2 = iter();

        if (&phase2 != &phase)
        {
            scalarCoeffTable::const_iterator Cvm
            (
                Cvms_.find(interfacePair(phase, phase2))
            );

            if (Cvm != Cvms_.end())
            {
                tSvm.ref() += Cvm()*phase2.rho()*phase2*phase2.DDtU();
            }
            else
            {
                Cvm = Cvms_.find(interfacePair(phase2, phase));

                if (Cvm != Cvms_.end())
                {
                    tSvm.ref() += Cvm()*phase.rho()*phase2*phase2.DDtU();
                }
            }
        }
    }

    volVectorField::Boundary& SvmBf =
        tSvm.ref().boundaryFieldRef();

    // Remove virtual mass at fixed-flux boundaries
    forAll(phase.phi().boundaryField(), patchi)
    {
        if
        (
            isA<fixedValueFvsPatchScalarField>
            (
                phase.phi().boundaryField()[patchi]
            )
        )
        {
            SvmBf[patchi] = Zero;
        }
    }

    return tSvm;
}


Foam::autoPtr<Foam::multiFluidSystem::dragCoeffFields>
Foam::multiFluidSystem::dragCoeffs() const
{
    autoPtr<dragCoeffFields> dragCoeffsPtr(new dragCoeffFields);

    forAllConstIter(dragModelTable, dragModels_, iter)
    {
        const dragModel& dm = *iter();

        volScalarField* Kptr =
            (
                max
                (
                    // fvc::average(dm.phase1()*dm.phase2()),
                    // fvc::average(dm.phase1())*fvc::average(dm.phase2()),
                    max(
                        dm.phase1()*dm.phase2(),

                        // artificial drag
                       (cAlpha_ / compressionCoeff) /
                        (dm.residualPhaseFraction())
                       ),

                    // stabalization
                    dm.residualPhaseFraction()
                    
                    // cAlpha_ / dm.residualPhaseFraction()
                )
               *dm.K
                (
                    max
                    (
                        mag(dm.phase1().U() - dm.phase2().U()),
                        dm.residualSlip()
                    )
                )
            ).ptr();
        

        Info<< "max cAlpha/residual: " << max(cAlpha_ / 
                                            dm.residualPhaseFraction())
            << endl;
            
        Info<< "min cAlpha/residual: " << min(cAlpha_ / 
                                            dm.residualPhaseFraction())
            << endl;
        
        volScalarField::Boundary& Kbf = Kptr->boundaryFieldRef();

        // Remove drag at fixed-flux boundaries
        forAll(dm.phase1().phi().boundaryField(), patchi)
        {
            if
            (
                isA<fixedValueFvsPatchScalarField>
                (
                    dm.phase1().phi().boundaryField()[patchi]
                )
            )
            {
                Kbf[patchi] = 0.0;
            }
        }

        dragCoeffsPtr().insert(iter.key(), Kptr);
    }

    return dragCoeffsPtr;
}


Foam::autoPtr<Foam::multiFluidSystem::cFields>
Foam::multiFluidSystem::cCoeffs() const
{
    autoPtr<cFields> cCoeffsPtr(new cFields);

    forAllConstIter(phaseChangeTable, phaseChangeModels_, iter)
    {
        const phaseChangeModel& pc = *iter();

        volScalarField* Cptr =
            (
                pc.mDotAlphal()[0]
            ).ptr();
        
        cCoeffsPtr().insert(iter.key(), Cptr);
    }

    return cCoeffsPtr;
}


Foam::autoPtr<Foam::multiFluidSystem::vFields>
Foam::multiFluidSystem::vCoeffs() const
{
    autoPtr<vFields> vCoeffsPtr(new vFields);

    forAllConstIter(phaseChangeTable, phaseChangeModels_, iter)
    {
        const phaseChangeModel& pc = *iter();

        volScalarField* Vptr =
            (
                pc.mDotAlphal()[1]
            ).ptr();
        
        vCoeffsPtr().insert(iter.key(), Vptr);
    }

    return vCoeffsPtr;
}


Foam::tmp<Foam::volScalarField> Foam::multiFluidSystem::dragCoeff
(
    const phaseModel& phase,
    const dragCoeffFields& dragCoeffs
) const
{
    tmp<volScalarField> tdragCoeff
    (
        new volScalarField
        (
            IOobject
            (
                "dragCoeff",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "dragCoeff",
                dimensionSet(1, -3, -1, 0, 0),
                0
            )
        )
    );

    dragModelTable::const_iterator dmIter = dragModels_.begin();
    dragCoeffFields::const_iterator dcIter = dragCoeffs.begin();
    for
    (
        ;
        dmIter != dragModels_.end() && dcIter != dragCoeffs.end();
        ++dmIter, ++dcIter
    )
    {
        if
        (
            &phase == &dmIter()->phase1()
         || &phase == &dmIter()->phase2()
        )
        {
            tdragCoeff.ref() += *dcIter();
        }
    }

    return tdragCoeff;
}


Foam::tmp<Foam::volScalarField> Foam::multiFluidSystem::phaseChangeSp
(
    const phaseModel& phase,
    const cFields& cCoeffs,
    const vFields& vCoeffs
) const
{
    tmp<volScalarField> tpcSp
    (
        new volScalarField
        (
            IOobject
            (
                "tpcSp",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "tpcSp",
                dimensionSet(0, 0, -1, 0, 0),
                0
            )
        )
    );    
    
    phaseChangeTable::const_iterator pcIter = phaseChangeModels_.begin();
    cFields::const_iterator cIter = cCoeffs.begin();
    vFields::const_iterator vIter = vCoeffs.begin();
    for
    (
        ;
        pcIter != phaseChangeModels_.end() && cIter != cCoeffs.end() && vIter != vCoeffs.end();
        ++pcIter, ++cIter, ++vIter
    )
    {
        if
        (
            &phase == &pcIter()->phase1()
            || &phase == &pcIter()->phase2()
        )
        {
            // vapor phase
            if(&phase == &pcIter()->phase1())
            {
                const volScalarField Sp = *cIter();
                // condensation is an implicit source term for vapor phase
                tpcSp.ref() = -Sp / phase.rho();
                // tpcSp.ref() = -Sp;                
            }
            // continuous liquid phase
            if(&phase == &pcIter()->phase2())
            {
                const volScalarField Sp = *vIter();
                // vaporization is an implicit source term for liquid phase
                tpcSp.ref()= Sp / phase.rho();
                // tpcSp.ref()= Sp;
            }
        }
    }
    
    return tpcSp;
}


Foam::tmp<Foam::volScalarField> Foam::multiFluidSystem::phaseChangeSu
(
    const phaseModel& phase,
    const cFields& cCoeffs,
    const vFields& vCoeffs
) const
{
    tmp<volScalarField> tpcSu
    (
        new volScalarField
        (
            IOobject
            (
                "tpcSu",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "tpcSu",
                dimensionSet(0, 0, -1, 0, 0),
                0
            )
        )
    );    
    
    phaseChangeTable::const_iterator pcIter = phaseChangeModels_.begin();
    cFields::const_iterator cIter = cCoeffs.begin();
    vFields::const_iterator vIter = vCoeffs.begin();
    for
    (
        ;
        pcIter != phaseChangeModels_.end() && cIter != cCoeffs.end() && vIter != vCoeffs.end();
        ++pcIter, ++cIter, ++vIter
    )
    {
        const phaseModel* phasePtr = nullptr;
        if
        (
            &phase == &pcIter()->phase1()
            || &phase == &pcIter()->phase2()
        )
        {
            // vapor phase
            if(&phase == &pcIter()->phase1())
            {
                // point to liquid phase
                phasePtr = &pcIter()->phase2();
                const volScalarField alpha = *phasePtr;
                const volScalarField Su = *vIter();
                // alphal * vCoeffs
                // equivalent to: alphal * mDotAlphal()[1]
                // vaporization is a explicit source term for vapor phase
                tpcSu.ref() = -Su*alpha / phase.rho();
                // tpcSu.ref() = -Su*alpha;
            }
            
            // continuous liquid phase
            if(&phase == &pcIter()->phase2())
            {
                // point to vapor phase
                phasePtr = &pcIter()->phase1();
                const volScalarField alpha = *phasePtr;
                const volScalarField Su = *cIter();
                // alphav * cCoeffs
                // equivalent to: alphav * mDotAlphal()[0]
                // condensation is a explicit source term for liquid phase
                tpcSu.ref()= Su*alpha / phase.rho();
                // tpcSu.ref()= Su*alpha;
            }
        }
    }
        
    return tpcSu;
}


Foam::tmp<Foam::surfaceScalarField> Foam::multiFluidSystem::surfaceTension
(
    const phaseModel& phase1
) const
{
    tmp<surfaceScalarField> tSurfaceTension
    (
        new surfaceScalarField
        (
            IOobject
            (
                "surfaceTension",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "surfaceTension",
                dimensionSet(1, -2, -2, 0, 0),
                0
            )
        )
    );

    // read-only volScalarField alpha1
    // const volScalarField& alpha1 = phase1;
    
    forAllConstIter(PtrDictionary<phaseModel>, phases_, iter)
    {
        const phaseModel& phase2 = iter();
        
        // read-only volScalarField alpha2
        // const volScalarField& alpha2 = phase2;

        if (&phase2 != &phase1)
        {
            scalarCoeffSymmTable::const_iterator sigma
            (
                sigmas_.find(interfacePair(phase1, phase2))
            );

            if (sigma != sigmas_.end())
            {
                // surface tension force only exists where cAlphaSwitch=1
                surfaceScalarField gamma = cAlphaGamma(phase1, phase2);                

                tSurfaceTension.ref() +=
                    dimensionedScalar("sigma", dimSigma_, sigma())
                    // insert cAlphaSwitch() here!!!!!!!!!!!!!!!!!!!!!!
                   *cAlphaSwitch(gamma)
                   *fvc::interpolate(K(phase1, phase2))*
                    (
                        fvc::interpolate(phase2)*fvc::snGrad(phase1)
                      - fvc::interpolate(phase1)*fvc::snGrad(phase2)
                    );
            }
        }
    }

    return tSurfaceTension;
}


Foam::tmp<Foam::volScalarField>
Foam::multiFluidSystem::nearInterface() const
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

    forAllConstIter(PtrDictionary<phaseModel>, phases_, iter)
    {
        tnearInt.ref() =
            max(tnearInt(), pos0(iter() - 0.01)*pos0(0.99 - iter()));
    }

    return tnearInt;
}


void Foam::multiFluidSystem::solve()
{
    forAllIter(PtrDictionary<phaseModel>, phases_, iter)
    {
        iter().correct();
    }

    const Time& runTime = mesh_.time();

    const dictionary& alphaControls = mesh_.solverDict("alpha");
    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();

        PtrList<volScalarField> alpha0s(phases_.size());
        PtrList<surfaceScalarField> alphaPhiSums(phases_.size());

        int phasei = 0;
        forAllIter(PtrDictionary<phaseModel>, phases_, iter)
        {
            phaseModel& phase = iter();
            volScalarField& alpha = phase;

            alpha0s.set
            (
                phasei,
                new volScalarField(alpha.oldTime())
            );

            alphaPhiSums.set
            (
                phasei,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "phiSum" + alpha.name(),
                        runTime.timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
                )
            );

            phasei++;
        }

        for
        (
            subCycleTime alphaSubCycle
            (
                const_cast<Time&>(runTime),
                nAlphaSubCycles
            );
            !(++alphaSubCycle).end();
        )
        {
            solveAlphas();

            int phasei = 0;
            forAllIter(PtrDictionary<phaseModel>, phases_, iter)
            {
                alphaPhiSums[phasei] += iter().alphaPhi()/nAlphaSubCycles;
                phasei++;
            }
        }

        phasei = 0;
        forAllIter(PtrDictionary<phaseModel>, phases_, iter)
        {
            phaseModel& phase = iter();
            volScalarField& alpha = phase;

            phase.alphaPhi() = alphaPhiSums[phasei];

            // Correct the time index of the field
            // to correspond to the global time
            alpha.timeIndex() = runTime.timeIndex();

            // Reset the old-time field value
            alpha.oldTime() = alpha0s[phasei];
            alpha.oldTime().timeIndex() = runTime.timeIndex();

            phasei++;
        }
    }
    else
    {
        solveAlphas();
    }
}


bool Foam::multiFluidSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        PtrList<entry> phaseData(lookup("phases"));
        label phasei = 0;

        forAllIter(PtrDictionary<phaseModel>, phases_, iter)
        {
            readOK &= iter().read(phaseData[phasei++].dict());
        }

        lookup("sigmas") >> sigmas_;
        //lookup("interfaceCompression") >> cAlphas_;
        lookup("virtualMass") >> Cvms_;

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
