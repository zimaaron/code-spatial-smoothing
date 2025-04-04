# History
This is work started during Sam's DVP visit to Bucknell during Nov 2023

# Motivation

The initial ideas was to work on spatially smoothing of spatial ligand and receptor gene expression to:
1) come up with "better" images of spatial cell-cell signaling
2) use the smoothing as a crude means of imputation

# Method
Use GP kriging (with fixed/specified parameters) to interpolate the spatial omics data to higher resolution. Initially we are running this independently per gene expression and without fitting parameters so we can try running NICHES on a stack of higher resolution spatially smoothed inputs.

# File Descriptions
* ExtractionScriptDemo - runs LRextract
* LRextract - function to pull out ligand receptor spatial data and xy coords
* gp-no-estimation - runs spatial smoothing on data from LRextract, but user supplies cov params
# spatial-transcriptomics
