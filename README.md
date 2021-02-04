# dynamicAMR

## Descrition

dynamic meshing with load balancing for hexahedral meshes in 3D and 2D

Library is based on:

Rettenmaier, Daniel, et al. "Load balanced 2D and 3D adaptive mesh refinement in OpenFOAM." SoftwareX 10 (2019): 100317.

link:
https://www.sciencedirect.com/science/article/pii/S2352711018301699

port to the OpenFOAM+ version v1812 and v2006,v2012

refinement selection algoritm is based on foam extended 4.1

## Getting Started

install OpenFOAM v1812 and v2006

compile the library
```
./Allwmake
```
in case of v2012:
```
git checkout v2012
./Allwmake
```
### Prerequisites

Requires OpenFOAM v1812 or v2006,v2012:

```
https://www.openfoam.com/download/release-history.php
```

### Installing

```
 git clone https://github.com/HenningScheufler/multiDimAMR
 cd multiDimAMR
 ./Allwmake
```
### Usage

add following lines to the contolDict:
```
libs
(
   "libdynamicLoadBalanceFvMesh.so"
);
or depending on the openfoam version
libs
(
   dynamicLoadBalanceFvMesh
);
```

cat decomposeParDict
```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                   |
|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains 4;

method          scotch;

constraints
{
    refinementHistoryMultiDim
    {
        //- Decompose cells such that all cell originating from single cell
        //  end up on same processor
        type    refinementHistoryMultiDim;
    }
}
```

cat balanceParDict:
```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      balanceParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains 4;

method          ptscotch; //clsutered //scotch

constraints
{
    refinementHistoryMultiDim
    {
        //- Decompose cells such that all cell originating from single cell
        //  end up on same processor
        type    refinementHistoryMultiDim;
    }
}
// ************************************************************************* //
```
### refinement selection

the cells to refine are specified in the adaptCriteria dictionary
```
adaptCriteria
{
    type fieldBounds; // options fieldBounds geometric
    fieldName alpha.water;
    lowerBound 0.001;
    upperBound 0.999;
    nLayer 2; // extends refinement zone by two layers
    maxCellLevel 2; // default value very high number limits the maxium refinement level
    minCellLevel 0; // default value 0 specify minimum refinement level
    negate false; // default false // negates the selection
}
```
the composedAdaptCriteria enables us to combine the simple functions above by logial operator:
```
adaptCriteria
{

    type composedAdaptCriteria;
    operation or; // and or xor
    criteria
    (
        interface // refine to the maxRefinement
        {
            type fieldBounds;
            fieldName alpha.water;
            lowerBound 0.01;
            upperBound 0.99;
            nLayer     2;
        }
        fluid // refLvl 2 in fluid
        {
            type fieldBounds;
            fieldName alpha.water;
            lowerBound 0.01;
            upperBound 2;
            maxCellLevel 2;
        }
    );
}
```
This way we can specify a higher refinement level at the interface than in the fluid or in combination with the geometric refinement option chose only to refine in a given region. Note that the composedAdaptCriteria can also consists of multiple composedAdaptCriteria


### Known Issues

*   boundary condition that use function e.g. uniformFixedValue
    crash if the table function is used see issue #2

    Solution: the issue can be resolved by switching to v2012

*   the load balancing crashes if a domain with no cells is created

    Solution: use more cells or try a different method in balanceParDict kahip seems to work better than ptscotch
I

## License

This project is licensed under the GPL v3 License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments



