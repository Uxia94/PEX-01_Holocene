# PEX-01_Holocene
Version archived for publication in Paleoceanography and Paleoclimatology. Includes data processing scripts, figure generation code, and analysis workflows associated with the manuscript.

## Scripts and Files
- **XRF_analysis_Peixao.R**: processes XRF elemental data (clr transformation), performs PCA, broken stick analysis and ANOSIM, and generates PC1 and Kclr time-series plots with 200-year smoothing. 
- **XRF_Peixao_clean.xlsx**: input file for `XRF_analysis_Peixao.R`
- **Plots_figure_5.R**: processes paleoclimate data. The script computes z-score normalization for δ²Hterr (dD) data and visualizes temporal variability using conditional areas fills (positive vs negative anomalies). Additionally, it bins NAO-like index data into 500-year intervals, calculates median values, and displays them as bar plots to emphasize centennial-scale variability. 
- **Data_figure_5.xlsx**: input file for `Plots_figure_5.R` containing: (1) sheet 1, δ²H (dD) data used for z-score standardization and visualization; and (2) sheet 2, NAO-like index data used for binning and median aggregation.
- **Statistical_assessment_deuterio_carbon_covariation.R**: this script evaluates the covariation between δ2Hterr (dD) and δ13Cterr time series using multiple statistical approaches: (1) Pearson correlation, (2) linear detrending and correlattion of residuals, and (3) first-differencing to assss high-frequency variability. 

## Output Files

### From `XRF_analysis_Peixao.R`:
- `XRF_clr.csv` (clr-transformed dataset containing elemental compositions used for PCA and subsequent analyses); biplot of PC1 vs PC2 including loadings and confidence ellipses; broken-stick scree plot; ANOSIM output (plot and R statistic and significance values); time-series plot of PC1 vs age (raw and smoothed); and time-series plot of clr-transformed K vs age (smoothed).

### From `Plots_figure_5.R`: 
- Generates two base R plots: a continuous z-score time series with filled anomalies and a bar-style representation of binned NAO variability.

### From `Statistical_assessment_deuterio_carbon_covariation.R`: 
- This script does not generate external files by default. Results are returned as: (a) Console output (Pearson correlation coefficints and associated statistical tests, and correlation results for detrended (residual) series and for first_differenced series), (b) Graphical output (scatter plot with linear regression line), and (c) Diagnostic plots: autocorrelation functions (ACF) fot both time series. 

## Usage

### General Setup:
- Ensure R and packages are installed.
- Place input files in yout working directory or update paths in scripts.
- Run scripts via `source("script_name.R")` in R.

### Run `XRF_analysis_Peixao.R`: 
- To export the clr-transformed dataset, add the following line to the script: `write.csv(xrf_clr, "XRF_clr.csv`, row.names = FALSE)

### Run `Statistical_assessment_deuterio_carbon_covariation.R`: 
- Uses the input file `Data_figure_5.xlsx` (sheet 1), which contains δ²H and δ¹³C isotope data from the terrestrial endmember derived from the NMF analysis.

## Notes
- **Input Data**: Contact the corresponding author (m.eugenia.fernandezp@udc.es) for supporting information or full datasets.
- **Customization**: Adjust parameters for other datasets.

## Citation
If you use this repository, please cite the associated manuscript (in preparation): Fernández-Pérez, U., et al. (in prep.). North Atlantic atmospheric circulation controlled the hydroclimate in western Iberia during the Holocene. 
A full reference will be updated upon publication.

## Contact
For questions, contact U. Fernández-Pérez (m.eugenia.fernandezp@udc.es).

