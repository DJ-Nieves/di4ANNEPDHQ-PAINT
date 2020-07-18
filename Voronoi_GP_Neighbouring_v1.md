Voronoi_GP_Neighbouring v1.0

Author: Daniel J. Nieves, University of Birmingham
E-mail: d.j.nieves@bham.ac.uk
 
This code is used to perform Voronoi tesselation analysis of ”marked" single molecule localisation data, e.g., here each point is associated with a value for membrane polarity. This code was implemented in "Quantitative mapping of nanoenvironments through single-molecule imaging of solvatochromic probes", Nieves et al., (2020), bioarxiv, doi: .

Here, the Voronoi tesselation of the marked data (in the form [x,y,GP value]) is calculated, and then neighbours within a certain GP range (rg) are found. If the tiles neighbour each other, and are within the GP range, they are combined. This proceeds until the neighbouring converges and no new neighbours within the GP range can be made. This can also be compared to a random neighbouring situation, which is generate by setting randGP equal to 1.
 
%%%%%%%%% NOTE: this code requires the "uniquecell" function to run completely %%%%%%%%%%%%%%
a copy will be provided at the repository with this code, or it can be downloaded from here: 

https://uk.mathworks.com/matlabcentral/fileexchange/31718-unique-elements-in-cell-array

make sure it is added to the path, or in the directory. 

Inputs and Outputs
 
Data to provide: 
- a matrix of the form (x,y, GP value), called "GP_data".
 
Inputs:
 
 - randGP:         Set to either 0 or 1. Default is 0. Choose 1 to run
                   	analysis again, with randomised neighbouring for
                   	comparison
 
 - rg:            	This sets the +/- range of GP values to include in
                   	the neighbouring. E.g. a value of 0.1 would mean values
                   	0.1 either side of a GP value will be included in the
                   	neighbouring group.
 
 - area_co:        This sets the maximum area (nm^2) to allow into the Voronoi
                   	neighbouring, as larger tiles maybe from areas of low data 
                   	coverage.
 
 - pos:            	Coordinates of the ROI from which to take [x,y,GP]
                  	data.
 
 Outputs:
 
 - GP_v:           Output from Voronoi tesselation of GP_data.
                   	- GP_v{:,1} - coordinates of tile vertices.
                  	- GP_v{:,2} - GP value.
                  	- GP_v{:,3} - Tile Area.
                  	- GP_v{:,4} - Indices of tiles within +/- rg of GP value.
                  	- GP_v{:,5} - Indices of the initial neighbouring result.
 
 - compGPv:      Final indices from complete neighbouring of Voronoi tiles by
                   	GP value.
 
 - compMeanGP:     Final mean GP value for combined tiles.
 
 - compGPareas:     Final combined tile areas.
 



