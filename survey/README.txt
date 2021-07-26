NOE SURVEY CODE

####INTRODUCTION####

This section of the code is concerned with finding pairs of aromatic rings
and amide protons that:

    -have coordinates listed in the PDB and chemical shifts in the BMRB 
    ("paired entries")
    -have a NOE restraint between the amide proton and at least one atom
    in the ring

Further analysis indicates how many restraints exist and whether or not the
chemical shift of the amide proton is an "outlier," which can be defined in
noe_analysis.py (default 2 sigma).
All files (BMRB entry, PDB entry, restraint file, and a list of paired entries)
are extracted from the June 2021 version of reboxitory, found in 
/reboxitory/2021/06
If not on NMRbox, this can be changed to download from PDB/BMRB with a small
amount of effort.

####RUNNING####

To run the default analysis, run noe_analysis.py . It has been tested on 
Python 3.8.0

The output will list the various exceptions raised during the analysis. These
exceptions occur for any paired entry that could not be accurately surveyed for
amide-aromatic NOEs. It will also generate Figures 5 & 6 from the manuscript, which
will be saved in ../images/combo_plot.pdf and ../images/noes_by_num.pdf , 
respectively.

Running the script for the first time will take a time on the order of hours as the
protein dump files in proteins/ and the data files in output/ are created. 
Subsequent runs should be quite fast.

