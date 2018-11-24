# Quantitative analysis of the yeast pheromone response pathway
Code for analysis of experiments that quantify the yeast pheromone response

## A Note for MBoC reviewers:

We are currently compiling the code used by the different authors for their analyses and commenting it appropriately. 
We will include: 

* Code for Bem1-GFP kymographs, rose diagrams, and single cell traces
* Code for the Fus3 translocation analysis
* Code for the flow cytometry analysis using 'FlowCytometryTools' package
* Code for the imaging cytometry analysis 

## Gradient Tracking

Scripts used for the analysis of gradient tracking experiments done in microfluidic chambers

### Kymographs

Cells were segmented manually in ImageJ based on the GFP images. The kymograph function requires a *.tif file for the mask of the cells and a *.tif file for the GFP channel. 

Run the scripts in the following order:

1. Wholecell_cap_v4.m
2. Wholecell_anglefix.m
3. Combokeeper.m

These three scripts require 

* jtrackv3.m
* meannan.m
* stderrornan.m
* jthresh.m
* lineprof.m

Example data is included to show the functionality of the scripts. 

These scripts were developed by Joshua B. Kelley and have been modified by Amy E. Pomeroy. 

### Single cell polar cap traces and polar histograms 

This notebook requires the following packages: numpy, pandas, math, and matplotlib

Single cell traces were generated using the manual tracking Plug-in in ImageJ

These notebooks were developed by Amy E. Pomeroy. 
