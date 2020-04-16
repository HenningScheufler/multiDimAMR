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
    Foam::adaptCriteria

Description
    An abstract base class for providing an interface for refinement selection
    algorithms. The interface provides two functionalities:

    1. Selects candidate cells to refine based on chosen criteria through
       refinementCellCandidates() member function. Returns a list of all cells
       that satisfy given criteria,
    2. Selects candidate split points to unrefine based on chosen criteria
       through unrefinementPointCandidates() member function. Returns a list of
       all points (not just split points!) that satisfy the criteria.

    Note: Here, we do not check that the unrefinement point candidates do not
    clash the points of cells marked as candidates for refinement. This is taken
    care of in polyhedralRefinement class.

SourceFiles
    adaptCriteria.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.
    modified by Henning Scheufler, DLR (2020)

\*---------------------------------------------------------------------------*/

#ifndef adaptCriteria_H
#define adaptCriteria_H

#include "dictionary.H"
#include "labelList.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "bitSet.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class adaptCriteria Declaration
\*---------------------------------------------------------------------------*/

class adaptCriteria
{
    // Protected data
    protected:
        //- Const reference to fvMesh
        const fvMesh& mesh_;

        //- Refinement selection dictionary (subdictionary of
        //  dynamicPolyRefinementFvMeshCoeffs dictionary)
        const dictionary coeffDict_;

    private:
    // Private Member Functions

        //- Disallow default bitwise copy construct
        adaptCriteria(const adaptCriteria&);

        //- Disallow default bitwise assignment
        void operator=(const adaptCriteria&);


protected:

    // protected data

    // optional maxLevel entry
    const label maxCellLevel_;

    // optional minLevel entry
    const label minCellLevel_;

    // optional minLevel entry
    const bool negate_;

    // Protected member functions

        // Access functions for derived classes

            //- Const access to the fvMesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Const access to the coefficient dictionary
            const dictionary& coeffDict() const
            {
                return coeffDict_;
            }




public:

    //- Runtime type information
    TypeName("adaptCriteria");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            adaptCriteria,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct from components
        adaptCriteria(const fvMesh& mesh, const dictionary& dict);


    // Selectors

        //- Return an autoPtr to the selected adaptCriteria
        static autoPtr<adaptCriteria> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~adaptCriteria();


    // Member Functions

        // Selection of refinement/unrefinement candidates

            //- Return bitSet of refinement cell candidates
            virtual bitSet refinementCellCandidates() const = 0;

            //- Return bitSet of protectedCells
            virtual bitSet protectedCells() const = 0;

            //- Return bitSet of unrefinement split point candidates
            virtual bitSet unrefinementPointCandidates() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
