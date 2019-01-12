# PHtools

PHtools provides a set of Matlab functions for analyzing the hypsometry of drainage basins in a provided set of digital elevation models. It performs the [Progressive Hypsometry algorithm (PH)](https://mcunningham917.github.io/PHdoc/Method/) described in Cunningham et al. (2019). The aim is specifically to demonstrate how hypsometric maxima are distributed in the landscape. 

Example output is available at [PHdoc](https://mcunningham917.github.io/PHdoc/). Users can create PH plots by supplying catchment DEMs.

## Requirements

PHtools is implemented in Matlab (R2017b) and requires Topotoolbox (v. 2.1).

## Structure

Code is separated into two folders: [`Core/`](https://github.com/mcunningham917/PHtools/tree/master/Core), which contains the scripts needed to run PH, and [`Example/`](https://github.com/mcunningham917/PHtools/tree/master/Example), which contains a demo, [`CostaRica.m`](https://github.com/mcunningham917/PHtools/blob/master/Example/CostaRica.m), and a [`Defaults.m`](https://github.com/mcunningham917/PHtools/blob/master/Example/Defaults.m) script that be shared by multiple driver scripts.

## Basic Operation

The PH algorithm involves two routines: 

### Step 1: Hypsometry of progressive (nested) subcatchments 

Traverse upstream from base level to main drainage divide along chains (flowpaths) in supplied DEM.

Record the modal elevation of the catchment draining to the (progressively higher) position on each chain

### Step 2: Progressive Hypsometric Bench (PHB) Identification

Identify nested subcatchments with similar modal elevation, i.e., PHBs.

## Output

PHtools generates and writes to a data repository, [`PHanalysis`](https://github.com/mcunningham917/PHanalysis). 

**Step 1:** Folder containing a txt file for each flow path: [`PHanalysis/CostaRica/Subcatchments/25mStep`](https://github.com/mcunningham917/PHanalysis/tree/master/CostaRica/Subcatchments/25mStep)

**Step 2:** Two folders:

One contains a .txt file for each supercatchment and one .txt for the entire ROI: [`PHanalysis/CostaRica/PHBs/Cusum02_BenchLength3Steps/AllSupercatchmentsTxt`](https://github.com/mcunningham917/PHanalysis/tree/master/CostaRica/PHBs/Cusum02_BenchLength3Steps/AllSupercatchmentsTxt) 

The second contains maps of PHBs (geotiffs) for each supercatchment:[`PHanalysis/CostaRica/PHBs/Cusum02_BenchLength3Steps/AllSupercatchmentTiffs`](https://github.com/mcunningham917/PHanalysis/tree/master/CostaRica/PHBs/Cusum02_BenchLength3Steps/AllSupercatchmentsTxt).

## Documentation

[PHdoc](https://mcunningham917.github.io/PHdoc/)

   - Core documentation of PH method and standard output
   - Links to all code and data repositories
