# Anomalous Amide Proton Chemical Shifts as Signatures of Hydrogen Bonding to Aromatic Sidechains

Find the manuscript at : \
https://doi.org/10.5194/mr-2021-53 \
\
The code found in this repository can be used to generate Figures 5, 6, 9, and 11 in the manuscript. It is broken up into 4 sections: \
\
  -survey (Figures 5 and 6): code for the restraint analysis portion of the paper (a survey of NOE restraints between amide protons and aromatic ring protons in paired BMRB and PDB entries) \
  -mosart (Figure 9): code for analysis of the van der Waals interaction between an ALA amide proton and a PHE aromatic ring as the proton approaches the ring \
  -dep (Figure 11): code for analysis of BMRB deposition trends, specifically structure depositions, structure depositions using CS Rosetta, and runs executed using the BMRB CS-Rosetta server by year \
  -images: stores images produced by other code\


## Installation

All code was built and run in Python 3.8.0 \
Installed python packages and their versions can be found in requirements.txt


  
## Running:

All figures can be found already made in images. To recreate them: \
\

  -survey (Figures 5 and 6): (within survey subdir) python noe_analysis.py (also outputs to terminal exceptions listed in Supp. Table 1). Figures 5 and 6 saved to images/combo_plot.pdf and images/noes_by_num.pdf , respectively \
  -mosart (Figure 9): (within mosart subdir) python vdw_plot.py . Figure 9 saved to images/plot_vdw.pdf \
  -dep (Figure 11): (within dep subdir) python dep_trends.py . Figure 11 saved to images/dep_plot.pdf 
