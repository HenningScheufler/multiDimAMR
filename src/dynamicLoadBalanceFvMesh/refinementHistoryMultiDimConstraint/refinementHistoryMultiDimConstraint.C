/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

#include "refinementHistoryMultiDimConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"
#include "refinementHistoryMultiDim.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{
    defineTypeName(refinementHistoryMultiDim);

    addToRunTimeSelectionTable
    (
        decompositionConstraint,
        refinementHistoryMultiDim,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionConstraints::refinementHistoryMultiDim::refinementHistoryMultiDim
(
    const dictionary& dict
)
:
    decompositionConstraint(dict, typeName)
{
    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : setting constraints to refinement history" << endl;
    }
}


Foam::decompositionConstraints::refinementHistoryMultiDim::refinementHistoryMultiDim()
:
    decompositionConstraint(dictionary(), typeName)
{
    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : setting constraints to refinement history" << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::decompositionConstraints::refinementHistoryMultiDim::add
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    // The refinement history type
    typedef ::Foam::refinementHistoryMultiDim HistoryType;

    // Local storage if read from file
    autoPtr<const HistoryType> readFromFile;

    const HistoryType* historyPtr =
        mesh.findObject<HistoryType>("refinementHistoryMultiDim");

    if (historyPtr)
    {
        if (decompositionConstraint::debug)
        {
            Info<< type() << " : found refinementHistoryMultiDim" << endl;
        }
    }
    else
    {
        if (decompositionConstraint::debug)
        {
            Info<< type() << " : reading refinementHistoryMultiDim from time "
                << mesh.facesInstance() << endl;
        }

        readFromFile.reset
        (
            new HistoryType
            (
                IOobject
                (
                    "refinementHistoryMultiDim",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh.nCells()
            )
        );

        historyPtr = readFromFile.get();  // get(), not release()
    }

    const auto& history = *historyPtr;

    if (history.active())
    {
        // refinementHistoryMultiDim itself implements decompositionConstraint
        history.add
        (
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            explicitConnections
        );
    }
}


void Foam::decompositionConstraints::refinementHistoryMultiDim::apply
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    // The refinement history type
    typedef ::Foam::refinementHistoryMultiDim HistoryType;

    // Local storage if read from file
    autoPtr<const HistoryType> readFromFile;

    const HistoryType* historyPtr =
        mesh.findObject<HistoryType>("refinementHistoryMultiDim");

    if (!historyPtr)
    {
        readFromFile.reset
        (
            new HistoryType
            (
                IOobject
                (
                    "refinementHistoryMultiDim",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh.nCells()
            )
        );

        historyPtr = readFromFile.get();  // get(), not release()
    }

    const auto& history = *historyPtr;

    if (history.active())
    {
        // refinementHistoryMultiDim itself implements decompositionConstraint
        history.apply
        (
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            explicitConnections,
            decomposition
        );
    }
}


// ************************************************************************* //
