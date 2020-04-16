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

#include "fieldBounds.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "labelIOField.H"
#include "syncTools.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(fieldBounds, 0);
addToRunTimeSelectionTable
(
    adaptCriteria,
    fieldBounds,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldBounds::fieldBounds
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    adaptCriteria(mesh, dict),
    fieldName_(coeffDict().lookup("fieldName")),
    lowerBound_(coeffDict().get<scalar>("lowerBound")),
    upperBound_(coeffDict().get<scalar>("upperBound")),
    nLayer_(coeffDict().lookupOrDefault<scalar>("nLayer",0))
{}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::fieldBounds::~fieldBounds()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

void Foam::fieldBounds::extendMarkedCells
(
    bitSet& markedCell
) const
{
    // Mark faces using any marked cell
    bitSet markedFace(mesh_.nFaces());

    for (const label celli : markedCell)
    {
        markedFace.set(mesh_.cells()[celli]);  // set multiple faces
    }

    syncTools::syncFaceList(mesh_, markedFace, orEqOp<unsigned int>());

    // Update cells using any markedFace
    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        if (markedFace.test(facei))
        {
            markedCell.set(mesh_.faceOwner()[facei]);
            markedCell.set(mesh_.faceNeighbour()[facei]);
        }
    }
    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); ++facei)
    {
        if (markedFace.test(facei))
        {
            markedCell.set(mesh_.faceOwner()[facei]);
        }
    }
}

Foam::bitSet
Foam::fieldBounds::refinementCellCandidates() const
{
    bitSet refineCells(mesh().nCells(),false);

    // Get the field
    const volScalarField& vField =
        mesh().lookupObject<volScalarField>(fieldName_);
    // get cellLevels
    const labelIOList& curCellLevel =
        mesh().lookupObject<labelIOList>("cellLevel");

    // Loop through internal field and collect cells to refine
    const scalarField& vfIn = vField.internalField();

    forAll (vfIn, cellI)
    {
        // Get current cell value
        const scalar& cellValue = vfIn[cellI];

        if
        (
            cellValue > lowerBound_
            && cellValue < upperBound_
            && curCellLevel[cellI] < maxCellLevel_
        )
        {
            // Cell value is within the bounds, append cell for potential
            // refinement
            refineCells.set(cellI);
        }
    }

    // Extend with a buffer layer to prevent neighbouring points
    // being unrefined.
    for (label i = 0; i < nLayer_; ++i)
    {
        extendMarkedCells(refineCells);
    }

    // Print out some information
    // Info<< "Selection algorithm " << type() << " selected "
    //     << returnReduce(refineCells.toc().size(), sumOp<label>())
    //     << " cells as refinement candidates."
    //     << endl;

    // Return the list in the Xfer container to prevent copying
    if(negate_)
    {
        refineCells = ~refineCells;
    }
    return refineCells;
}


Foam::bitSet
Foam::fieldBounds::unrefinementPointCandidates() const
{
    bitSet unrefinePoints(mesh().nPoints(),false);

    // Get the field
    const volScalarField& vField =
        mesh().lookupObject<volScalarField>(fieldName_);

    // get cellLevels
    const labelIOList& curCellLevel =
        mesh().lookupObject<labelIOList>("cellLevel");

    // Loop through internal field and collect cells to refine
    const scalarField& vfIn = vField.internalField();

    forAll (mesh().points(), pI)
    {
        // Get point value
        for (const auto celli : mesh().pointCells()[pI])
        {
            bool unRefPoint = true;
            const scalar pointValue = vfIn[celli];
            // all cells need to fullfill this criteria
            if
            (
                !(pointValue > upperBound_
                || pointValue < lowerBound_
                || curCellLevel[celli] > maxCellLevel_)
            )
            {
                unRefPoint = false;
                break;
            }
            
            if(unRefPoint)
            {
                unrefinePoints.set(pI);
            }
        }
    }

    syncTools::syncPointList(mesh(), unrefinePoints, minEqOp<unsigned int>(), 0);

    // Print out some information
    // Info<< "Selection algorithm " << type() << " selected "
    //     << returnReduce(unrefinePoints.toc().size(), sumOp<label>())
    //     << " points as unrefinement candidates."
    //     << endl;

    if(negate_)
    {
        unrefinePoints = ~unrefinePoints;
    }

    return unrefinePoints;
}


// ************************************************************************* //
