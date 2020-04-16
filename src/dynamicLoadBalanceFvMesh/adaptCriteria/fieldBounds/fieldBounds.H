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
    Foam::fieldBounds

Description
    Selection of refinement cells based on a given scalar field bounds:
    if a given (scalar) field is between lowerBound and upperBound, cell will
    be marked for refinement.

SourceFiles
    fieldBounds.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.
    modified by Henning Scheufler, DLR (2020)

\*---------------------------------------------------------------------------*/

#ifndef fieldBounds_H
#define fieldBounds_H

#include "adaptCriteria.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fieldBounds Declaration
\*---------------------------------------------------------------------------*/

class fieldBounds
:
    public adaptCriteria
{
    // Private data

        //- Name of the field to refine
        word fieldName_;

        //- Lower bound for the field
        scalar lowerBound_;

        //- Upper bound for the field
        scalar upperBound_;

        //- layers to extend the refinementregion
        scalar nLayer_;


    // Private Member Functions

        //- Extend markedCell with cell-face-cell.
        void extendMarkedCells(bitSet& markedCell) const;

        //- Disallow default bitwise copy construct
        fieldBounds(const fieldBounds&);

        //- Disallow default bitwise assignment
        void operator=(const fieldBounds&);


public:

    //- Runtime type information
    TypeName("fieldBounds");


    // Constructors

        //- Construct from components
        fieldBounds(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~fieldBounds();


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
