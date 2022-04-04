/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "adaptCriteria.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(adaptCriteria, 0);
defineRunTimeSelectionTable(adaptCriteria, dictionary);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adaptCriteria::adaptCriteria
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    coeffDict_(dict),
    maxCellLevel_(dict.lookupOrDefault<label>("maxCellLevel",labelMax)),
    minCellLevel_(dict.lookupOrDefault<label>("minCellLevel",0)),
    negate_(dict.lookupOrDefault<bool>("negate",false))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::adaptCriteria> Foam::adaptCriteria::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    // Get the name of the desired refinement selection algorithm
    const word adaptCriteriaTypeName(dict.lookup("type"));
    Info<< "Creating adaptCriteria " << adaptCriteriaTypeName << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(adaptCriteriaTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
            (
            "adaptCriteria::adaptCriteria::New\n"
            "(\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
            )   << "Unknown adaptCriteria type "
                << adaptCriteriaTypeName << endl << endl
                << "Valid adaptCriteria types are :" << endl
                << dictionaryConstructorTablePtr_->toc()
                << exit(FatalError);
    }

    return autoPtr<adaptCriteria>(cstrIter()(mesh, dict));
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::adaptCriteria::~adaptCriteria()
{}


// ************************************************************************* //
