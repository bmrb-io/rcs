# Anomalous Amide Proton Chemical Shifts as Signatures of Hydrogen Bonding to Aromatic Sidechains

Find the manuscript at : \
https://doi.org/10.5194/mr-2021-53 \
\
The code found in this repository can be used to generate Figures 3, 4, 5, 6, 9, and 11 in the manuscript. It is broken up into 4 sections: \
\
  -chemical_shift_analysis (Figures 3 and 4): code for the chemical shift/geometric analysis portion of the paper \
  -survey (Figures 5 and 6): code for the restraint analysis portion of the paper (a survey of NOE restraints between amide protons and aromatic ring protons in paired BMRB and PDB entries) \
  -mosart (Figure 9): code for analysis of the van der Waals interaction between an ALA amide proton and a PHE aromatic ring as the proton approaches the ring \
  -dep (Figure 11): code for analysis of BMRB deposition trends, specifically structure depositions, structure depositions using CS Rosetta, and runs executed using the BMRB CS-Rosetta server by year 
  
README files can be found in each of the above subdirectories (except images) with more specific descriptions and instructions.

## Installation

All code was built and run in Python 3.8.0 \
Installed python packages and their versions can be found in requirements.txt


  
## Running:

All figures can be found already made in images. To recreate them:  


  -survey (Figures 5 and 6): (within survey subdir) python noe_analysis.py (also outputs to terminal exceptions listed in Supp. Table 1). Figures 5 and 6 saved to images/combo_plot.pdf and images/noes_by_num.pdf , respectively \
  -mosart (Figure 9): (within mosart subdir) python vdw_plot.py . Figure 9 saved to images/plot_vdw.pdf \
  -dep (Figure 11): (within dep subdir) python dep_trends.py . Figure 11 saved to images/dep_plot.pdf 
  
## Demo:
As an example of data federation, demo.ipynb contains an interactive, pared-down version of some of the survey code. It parses a BMRB entry for chemical shifts, a restraint file for amide-aromatic restraints, and combines the two to fine the chemical shifts of amides restrained to aromatic ring protons. Excluded from the demo is much of the code that deals with errors that arise from surveying/federating two large databases. 

To run the demo, only PyNMRSTAR needs to be installed. It was last run with PyNMRSTAR version 3.2.1.

# Chemical shift analysis
The procedure to calculate Z-score and azimuthal angle has been explained in the Jupyter notbook in chemical_shift_analysis folder. The underlying data for Figure 3 & 4 in the manuscript can be found in chemical_shift_analysis/data

## Order parameter
Order parameter can be fetched via API for a given brmb ID. Once we have the combined information from PDB and BMRB 
as a CSV file, order parameter can be appended using the script order_parameter.py. This data is used to generate 
Fig 10 in the manuscript

## Software information
The list of sfotware used to generate BMRB can also be extracted via API. Simalr to order parameter the combined CSV file 
can be used as input to generate Fig 12 in the manuscript. 