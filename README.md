# Quantitative analysis of the yeast pheromone response pathway
Code for analysis of experiments that quantify the yeast pheromone response

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

## Translocation reporter

1.To run the code, open the start_file_cyto_new.m and put in details in the scripts, like directory of the data, channel information. 

e.g. 
Folders{i} = '~/Desktop/Data/'
BaseFileNameExample{i} = 'alphaFactor1um'

2.For each data set, the code will generate several csv files contains the fluorescense information. The fluorescence intensity of nuclear and cyto has been calculated with  different method, as it shown in the csv file name ( e.g. mean.csv, median.csv stand for nuclear fluorescence intensity, cyto_mean.csv, cyto_median.csv stand for the cytoplasmic fluorescence intensity). 

Code from Yang Li. 
