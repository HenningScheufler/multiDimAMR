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

Class
    Foam::decompositionConstraints::refinementHistoryMultiDim

Description
    Constraint to keep all cells originating from refining the same cell
    onto the same processor. Reads polyMesh/refinementHistoryMultiDim.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required  | Default
        type        | refinementHistoryMultiDim                 | yes   |
    \endtable

SourceFiles
    refinementHistoryMultiDimConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef refinementHistoryMultiDimConstraint_H
#define refinementHistoryMultiDimConstraint_H

#include "decompositionConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{

/*---------------------------------------------------------------------------*\
                      Class refinementHistoryMultiDim Declaration
\*---------------------------------------------------------------------------*/

class refinementHistoryMultiDim
:
    public decompositionConstraint
{
public:

    //- Runtime type information
    TypeName("refinementHistoryMultiDim");


    // Constructors

        //- Construct with constraint dictionary
        explicit refinementHistoryMultiDim(const dictionary& dict);

        //- Construct from components
        refinementHistoryMultiDim();


    //- Destructor
    virtual ~refinementHistoryMultiDim() = default;


    // Member Functions

        //- Add this constraint to list of constraints
        virtual void add
        (
            const polyMesh& mesh,
            boolList& blockedFace,
            PtrList<labelList>& specifiedProcessorFaces,
            labelList& specifiedProcessor,
            List<labelPair>& explicitConnections
        ) const;

        //- Apply this constraint post-decomposition
        virtual void apply
        (
            const polyMesh& mesh,
            const boolList& blockedFace,
            const PtrList<labelList>& specifiedProcessorFaces,
            const labelList& specifiedProcessor,
            const List<labelPair>& explicitConnections,
            labelList& decomposition
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decompositionConstraints
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
