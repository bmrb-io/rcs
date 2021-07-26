DEPOSITION TRENDS CODE

####INTRODUCTION####

This section of the code is concerned with analyzing deposition trends to BMRB,
particularly:

    (a) structural depositions to BMRB by year
    (b) structural depositions to BMRB citing CS Rosetta by year
    (c) runs executed using the BMRB CS-Rosetta server

(a) was found at 
https://api.bmrb.io/v2/meta/release_statistics

(b) can be found in CS-Rosetta_Entries.csv

(c) was found at 
https://csrosetta.bmrb.io/

####RUNNING####

To generate Figure 11, run dep_trends.py . It will be saved to ../images/dep_plot.pdf
