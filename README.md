# PEX-01_Holocene
Version archived for publication in Paleoceanography and Paleoclimatology. Includes data processing scripts, figure generation code, and analysis workflows associated with the manuscript.

## Scripts and Files
- **XRF_analysis_Peixao.R**: processes XRF elemental data (clr transformation), performs PCA, broken stick analysis and ANOSIM, and generates PC1 and Kclr time-series plots with 200-year smoothing. 
- **XRF_Peixao_clean.xlsx**: input file for `XRF_analysis_Peixao.R`

## Output Files

### From `XRF_analysis_Peixao.R`:
- `d2H_precipitation_simulations.csv`: XRF_clr.csv (clr-transformed dataset containing elemental compositions used for PCA and subsequent analyses); biplot of PC1 vs PC2 including loadings and confidence ellipses; broken-stick scree plot; ANOSIM output (plot and R statistic and significance values); time-series plot of PC1 vs age (raw and smoothed); and time-series plot of clr-transformed K vs age (smoothed).


## Usage

### General Setup:
- Ensure R and packages are installed.
- Place input files in yout working directory or update paths in scripts.
- Run scripts via `source("script_name.R")` in R.

### Run `XRF_analysis_Peixao.R`: 
- To export the clr-transformed dataset, add the following line to the script: `write.csv(xrf_clr, "XRF_clr.csv`, row.names = FALSE)

## Notes
- **Input Data**: Contact the corresponding author (m.eugenia.fernandezp@udc.es) for supporting information or full datasets.


## Contact
For questions, contact U. Fernández-Pérez (m.eugenia.fernandezp@udc.es).

