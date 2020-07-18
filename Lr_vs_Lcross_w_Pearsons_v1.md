Lr_vs_Lcross_w_Pearsons v1.0
 
This is used to calculate the L(r) vs Lcross(r) (according to Rossy et al., (2014) Histochem Cell Biol 141, 605-612) from "marked" single molecule localisation data, e.g., here each point is associated with a value for membrane polarity. This code was implemented in  "Quantitative mapping of nanoenvironments through single-molecule imaging of solvatochromic probes", Nieves et al., (2020), bioarxiv, doi:
 
The L(r) vs Lcross(r) is calculated for the GP segregated channels, and also for points randomly assigned into the channels as a control. Finally, the Pearson's correlation coefficient for the L(r) vs Lcross(r) data is calculated for each channel.
 
Author: Daniel J. Nieves, University of Birmingham
email: d.j.nieves@bham.ac.uk
 
Data to provide: 
- a matrix of the form (x,y, GP value), called "GP_data".

Inputs and Outputs:

Inputs:

- r :	         Radius of capture for Lr and Lcross (nm)

- GP_co:       GP value to use split the data. A vector using a series of cut-off values may also be used.

- pos:            Coordinates of the ROI from which to take [x,y,GP] data.

 
Outputs:
 
- Lr_plus :       		list of L(r) values for each point in the above GP cut-off
                  		data.
- Lr_minus :      		list of L(r) values for each point in the above GP cut-off
                  		data.
- Lr_rand_plus :  	list of L(r) values for each point assigned
                  		randomly in the above GP cut-off data.
- Lr_rand_minus : 	list of L(r) values for each point assigned
                  		randomly in the below GP cut-off data.
 
- Lc_plus :      		list of Lcross(r) values for each point in the above GP cut-off
                  		data.
- Lc_minus :      		list of Lcross(r) values for each point in the below GP cut-off
                  		data.
- Lc_rand_plus :  	list of Lcross(r) values for each point assigned
                  		randomly in the above GP cut-off data.
- Lc_rand_minus : 	list of Lcross(r) values for each point assigned
                  		randomly in the below GP cut-off data.
 
- pLrLc_plus :            	Pearson's coefficent of L(r) vs Lcross(r) for above GP cut-off data
- pLrLc_minus :           	Pearson's coefficent of L(r) vs Lcross(r) for below GP cut-off data
- pLrLc_rand_plus :      Pearson's coefficent of L(r) vs Lcross(r) for randomly assigned above GP cut-off data
- pLrLc_rand_minus :   Pearson's coefficent of L(r) vs Lcross(r) for randomly assigned below GP cut-off data
 

