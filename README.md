# dynamicAMR

## Descrition

dynamic meshing with load balancing for hexahedral meshes in 3D and 2D

Library is based on:

Rettenmaier, Daniel, et al. "Load balanced 2D and 3D adaptive mesh refinement in OpenFOAM." SoftwareX 10 (2019): 100317.

link:
https://www.sciencedirect.com/science/article/pii/S2352711018301699

port to the OpenFOAM+ version v1812 and v2006

refinement selection algoritm is based on foam extended 4.1

## Getting Started

install OpenFOAM v1812

compile the library

./Allwmake

### Prerequisites

Requires OpenFOAM v1812:

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

add to contolDict:
```
libs
(
   "libdynamicLoadBalanceFvMesh.so"	
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

## License

This project is licensed under the GPL v3 License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments



