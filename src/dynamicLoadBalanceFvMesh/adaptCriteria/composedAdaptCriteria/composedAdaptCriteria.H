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

Class
    Foam::composedAdaptCriteria

Description
    Selection of refinement cells based on an arbitrary number of combined
    "basic" selection algorithms.

SourceFiles
    composedAdaptCriteria.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.
    modified by Henning Scheufler, DLR (2020)

\*---------------------------------------------------------------------------*/

#ifndef composedAdaptCriteria_H
#define composedAdaptCriteria_H

#include "adaptCriteria.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class composedAdaptCriteria Declaration
\*---------------------------------------------------------------------------*/

class composedAdaptCriteria
:
    public adaptCriteria
{
    // Private data

        enum operationType
        {
            opOr,
            opAnd,
            opXor,
        };

        //- Operation type names
        static const Enum<operationType> operationTypeNames_;

        //- Operation to apply to values
        operationType operation_;

        //- List of basic refinement selection algorithms
        PtrList<adaptCriteria> criteriaSelections_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        composedAdaptCriteria(const composedAdaptCriteria&);

        //- Disallow default bitwise assignment
        void operator=(const composedAdaptCriteria&);


public:

    //- Runtime type information
    TypeName("composedAdaptCriteria");


    // Constructors

        //- Construct from components
        composedAdaptCriteria
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~composedAdaptCriteria();


    // Member Functions

        // Selection of refinement/unrefinement candidates

            //- Return bitSet of refinement cell candidates
            virtual bitSet refinementCellCandidates() const;

            //- Return bitSet of protectedCells
            virtual bitSet protectedCells() const
            {
                bitSet protectedCells(mesh().nCells(),false);
                return protectedCells;
            }

            //- Return bitSet of unrefinement split point candidates
            virtual bitSet unrefinementPointCandidates() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
