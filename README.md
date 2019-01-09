# PHtools
Matlab toolbox for identifying patterns in catchment hypsometry 

PHtools provides a set of Matlab functions for analyzing the hypsometry of drainage basins in a provided set of digital elevation models. The perform the Progressive Hypsometry Algorithm (PH) described in Cunningham et al. (2019). The aim is specifically to demonstrate how hypsometric maxima are distributed in the landscape. Users can create PH plots by supplying catchment DEMs.

## Requirements

PHtools is implemented in Matlab (R2017b) and requires Topotoolbox (v. 2.1).

## Structure

Code is separated into two folders: `Core/`, which contains the scripts need to run PH, and an example folder, which contains a demo, `CostaRica.m` and a `Defaults` script that be shared by multiple driver scripts.
