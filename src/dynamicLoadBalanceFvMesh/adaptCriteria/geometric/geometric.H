/*---------------------------------------------------------------------------*\
            Copyright (c) 2020, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
	unofficial extension to OpenFOAM.
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
    Foam::geometric

Description
    Selection of refinement cells based on geometric constrains:
    They are selected if the are inside the geometry

SourceFiles
    geometric.C

Author
    Henning Scheufler, DLR  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef geometric_H
#define geometric_H

#include "adaptCriteria.H"
#include "Switch.H"
#include "searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class geometric Declaration
\*---------------------------------------------------------------------------*/

class geometric
:
    public adaptCriteria
{
    // Private data
        
    //- surface pointer
    autoPtr<searchableSurface> surfPtr_;


    // Private Member Functions

        //- Extend markedCell with cell-face-cell.
        void extendMarkedCells(bitSet& markedCell) const;

        //- Disallow default bitwise copy construct
        geometric(const geometric&);

        //- Disallow default bitwise assignment
        void operator=(const geometric&);


public:

    //- Runtime type information
    TypeName("geometric");


    // Constructors

        //- Construct from components
        geometric(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~geometric();


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
