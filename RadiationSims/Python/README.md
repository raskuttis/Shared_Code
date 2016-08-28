Package containing all the Python code used to analyse the outputs of Hyperion simulations. Code is stored in the following directories

    1) Hyperion: Most of the code required to read in Hyperion outputs (particularly hyp_read) and analyze them. Output plot directories, the clusters and folders in which Hyperion data are located and the names of .hst and .out files are all in hyp_models along with the lists of models used in our analysis. All of these can either be changed in hyp_models or in the individual python functions.
    2) Experimental: A number of scripts used to plot hst outputs, PDFs and other utilities
    3) Paper_I, PI_Referee, Paper_II, Thesis: Figures from Papers and Thesis, all labeled according to their figure numbers
    4) Tables: Scripts to output tables that can be input to Latex for Papers and Thesis

    Because the subdirectories are set up as a package, scripts need to be run from the directory above the Python directory using the command
        python -m Python.subdir.script
        e.g. to plot Figure 1 from Paper 1 you'd run
        python -m Python.Paper_I.f1

    Modules that you'll need to have installed to make everything run smoothly are:
        numpy, scipy, subprocess, re, copy, sys, matplotlib