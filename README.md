# PHtools

PHtools provides a set of Matlab functions for analyzing the hypsometry of drainage basins in a provided set of digital elevation models. It performs the Progressive Hypsometry Algorithm (PH) described in Cunningham et al. (2019). The aim is specifically to demonstrate how hypsometric maxima are distributed in the landscape. 

Example output is available at [PHdoc](https://mcunningham917.github.io/PHdoc/). Users can create PH plots by supplying catchment DEMs.

## Requirements

PHtools is implemented in Matlab (R2017b) and requires Topotoolbox (v. 2.1).

## Structure

Code is separated into two folders: `Core/`, which contains the scripts need to run PH, and `Example`, which contains a demo, `CostaRica.m` and a `Defaults.m` script that be shared by multiple driver scripts.

## Basic Operation

The PH algorithm involves two routines: 

### Step 1:Hypsometry of progressive (nested) subcatchments 

Traverse upstream along chains (flowpaths) in supplied DEM and record the modal elevation of the catchment draining to (progressively higher) position on stream network.

### Step 2:Progressive Hypomsetric Bench (PHB) Identification

Indentify groups of nested subcatchments with similar modal elevation, PHBs.

### Output
PHtools writes the data repository, PHanalysis. 

**Step 1:** Folder containing a txt file for each flow path: `PHanalysis/ROI/Subcatchments/25mStep`
**Step 2:** Folder containing a txt file for each PHB: `PHanalysis/ROI/PHBs/AllSupercatchments`

## Documentation

[PHdoc](https://mcunningham917.github.io/PHdoc/)

   - Core documentation of PH method and standard output
   - Links to all code and data repositories
