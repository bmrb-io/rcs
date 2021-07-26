VDW APPROACH CODE

####INTRODUCTION####

This section of the code is concerned with calculating the van der Waal's
interaction between the amide proton of an ALA approaching a PHE aromatic ring.
Each of the .pdb files contains the atomic coordinates for such an assembly with
the ALA N-H aligned with the PHE ring normal. MoSART  
(Hoch, J. C. and Stern, A. S.: MoSART [code], 2003.) has been used to perform the
calculations. 



####RUNNING####
To generate the plot found in Figure 9 of the manuscript, simply run vdw_plot.py .
It has been tested on Python 3.8.0

In order to generate the energies anew, run MoSART using the input file 
energy-amide-aromatic.inp. This can be done on NMRbox by running:
/usr/software/mosart/mosart energy-amide-aromatic.inp
If not on NMRbox, mosart can be found at https://simtk.org/projects/mosart
To recreate the plot, the van der Waals terms for each assembly would then have
to be copied into energies_fine.csv. 