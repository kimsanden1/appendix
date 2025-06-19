# Sensitivity analysis of missing data approaches
Code for sensitivity analysis used in my internship report appendix.

This repository contains the code and output related to the sensitivity analysis performed for the appendix of my internship report. The analysis checks what happens to the networks visually under different combinations of missing data strategies and participant inclusion thresholds.

## Overview
Three approaches to handling missing data were compared:
1. No imputation
2. Kalman filtering
3. multilevel MICE

Each method was applied to datasets filtered at four participant compliance thresholds:
- 0% (i.e., all participants are included)
- 25%
- 50%
- 75% (i.e., high compliance)

For each resulting dataset, mlVAR models were estimated to extract temporal and contemporaneous network structures. The main objective was to examine the visual consistency of networks across the methods and thresholds.

## Contents
- Output of the 12 approaches (3 methods x 4 thresholds) for both network types (temporal and contemporaneous). See `output.pdf` for all network plots.
- The full R code (`analysis.R`) used for the sensitivity analysis pipeline.
- Interactive plots are available for centrality indices (strength, in-strength, and out-strength) under `contemp_cent_plot.html` and `temp_cent_plot.html` (download files to open).

### How to use
The R script (`analysis.R`) contains the full pipeline. 

### Citation
If you use this code or analysis in your own work, please cite this repository or reference it as part of the appendix of the internship report or, if available, the published version.
