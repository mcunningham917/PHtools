# PHtools

PHtools provides a set of Matlab functions for analyzing the hypsometry of drainage 
basins in a provided set of digital elevation models. 
It performs the 
[Progressive Hypsometry algorithm (PH)](https://mcunningham917.github.io/PHdoc/Method) 
described in 
[Cunningham et al. (2019)](https://mcunningham917.github.io/PHdoc/Publications). 
The aim is specifically to demonstrate how 
hypsometric maxima are distributed in the landscape. 


## Requirements

PHtools is implemented in [`Matlab`](https://www.mathworks.com/products/matlab.html)
 (R2017b) and requires 
[`TopoToolbox`](https://topotoolbox.wordpress.com/) (v. 2.1).

## Structure

Code is separated into two folders: 
[`Core/`](https://github.com/mcunningham917/PHtools/tree/master/Core), which contains
 the scripts needed to run PH, and 
 [`Example/`](https://github.com/mcunningham917/PHtools/tree/master/Example), 
 which contains a demo, 
 [`CostaRica.m`](https://github.com/mcunningham917/PHtools/blob/master/Example/CostaRica.m), 
 and a 
 [`Defaults.m`](https://github.com/mcunningham917/PHtools/blob/master/Example/Defaults.m) 
 script that be shared by multiple driver scripts.

## Basic Operation

The PH algorithm involves two routines: 

 <dl>
  <dt>Step 1: Hypsometry of progressive (nested) subcatchments </dt>
  <dd> 
	Traverse upstream from base level to main drainage divide along chains (flowpaths) in 
	supplied DEM.
  <br> 
	Record the modal elevation of the catchment draining to the (progressively higher) 
	position on each chain.
  </dd>
  <dt>Step 2: Progressive Hypsometric Bench (PHB) Identification</dt>
  <dd>
  	Identify nested subcatchments with similar modal elevation, i.e., PHBs.
  </dd>
</dl> 


## Output

PHtools generates and writes to a results folder (archived example at 
[`PHanalysis`](https://github.com/mcunningham917/PHanalysis)) with two subfolders:

 <dl>
  <dt><em>Subcatchments</em> subfolder:</dt>
  <dd> 
	 which contains a text file for each flow path
  </dd>
  <dt><em>PHBs</em> subfolder:</dt>
  <dd> 
  
  which contains two subfolders:
  
  <ul>
	<li>
		contains a .txt file for each supercatchment and one .txt for the entire ROI: 
		<a href="https://github.com/mcunningham917/PHanalysis/tree/master/CostaRica/PHBs/Cusum02_BenchLength3Steps/AllSupercatchmentsTxt">
		PHanalysis/CostaRica/PHBs/Cusum02_BenchLength3Steps/AllSupercatchmentsTxt</a> 
	</li>
	<li>
		contains maps of PHBs (geotiffs) for each supercatchment:
		<a href="https://github.com/mcunningham917/PHanalysis/tree/master/CostaRica/PHBs/Cusum02_BenchLength3Steps/AllSupercatchmentsTxt">
		PHanalysis/CostaRica/PHBs/Cusum02_BenchLength3Steps/AllSupercatchmentTiffs</a> 
	</li>
</ul>
  </dd>
</dl> 


## Documentation

[PHdoc](https://mcunningham917.github.io/PHdoc/) provides summary documentation, 
and includes:

   - core description of PH method and standard output
   - links to all code and data repositories
