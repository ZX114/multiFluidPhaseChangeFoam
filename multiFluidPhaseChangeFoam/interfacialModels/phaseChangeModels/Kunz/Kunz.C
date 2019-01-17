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

#include "Kunz.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeModels
{
    defineTypeNameAndDebug(Kunz, 0);
    addToRunTimeSelectionTable(phaseChangeModel, Kunz, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModels::Kunz::Kunz
(
    const dictionary& interfaceDict,
    const phaseModel& phase1,
    const phaseModel& phase2
)
:
    phaseChangeModel(interfaceDict, phase1, phase2),

    UInf_("UInf", dimVelocity, phaseChangeModelCoeffs_.lookup("UInf") ),
    tInf_("tInf", dimTime, phaseChangeModelCoeffs_.lookup("tInf") ),
    Cc_("Cc", dimless, phaseChangeModelCoeffs_.lookup("Cc") ),
    Cv_("Cv", dimless, phaseChangeModelCoeffs_.lookup("Cv") ),

    p0_("0", pSat().dimensions(), 0.0),

    mcCoeff_(Cc_*phase1.rho()/tInf_),
    mvCoeff_(Cv_*phase1.rho()/(0.5*phase2.rho()*sqr(UInf_)*tInf_))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::Kunz::mDotAlphal() const
{
    const volScalarField& p = phase2_.db().lookupObject<volScalarField>("p_rgh");
    volScalarField limitedAlpha2(min(max(phase2_, scalar(0)), scalar(1)));

    return Pair<tmp<volScalarField>>
    (
        mcCoeff_*sqr(limitedAlpha2)
       *max(p - pSat(), p0_)/max(p - pSat(), 0.01*pSat()),

        mvCoeff_*min(p - pSat(), p0_)
    );
}

/*
Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::Kunz::mDotP() const
{
    const volScalarField& p = phase2_.db().lookupObject<volScalarField>("p");
    volScalarField limitedAlpha2(min(max(phase2_, scalar(0)), scalar(1)));

    return Pair<tmp<volScalarField>>
    (
        mcCoeff_*sqr(limitedAlpha2)*(1.0 - limitedAlpha2)
       *pos0(p - pSat())/max(p - pSat(), 0.01*pSat()),

        (-mvCoeff_)*limitedAlpha2*neg(p - pSat())
    );
}


bool Foam::phaseChangeModels::Kunz::read()
{
    if (phaseChangeModels::read())
    {
        phaseChangeModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        phaseChangeModelCoeffs_.lookup("UInf") >> UInf_;
        phaseChangeModelCoeffs_.lookup("tInf") >> tInf_;
        phaseChangeModelCoeffs_.lookup("Cc") >> Cc_;
        phaseChangeModelCoeffs_.lookup("Cv") >> Cv_;

        mcCoeff_ = Cc_*phase1.rho()/tInf_;
        mvCoeff_ = Cv_*phase1.rho()/(0.5*phase2.rho()*sqr(UInf_)*tInf_);

        return true;
    }
    else
    {
        return false;
    }
}
*/

// ************************************************************************* //
