# PEX-01_Holocene
Version archived for publication in Paleoceanography and Paleoclimatology. Includes data processing scripts, figure generation code, and analysis workflows associated with the manuscript.

## Scripts and Files
-**Plots_alkanes.R**:generates three figures from n-alkanes and NMF-derived isotope data: (1) modern plant n-alkane distributions by vegetation endmember, (2) reconstructed percipitation δD (δD_prc) from terrestrial endmembers (EM1+EM2), and (3) NMF-derived terrestrial isotope curves (δ²Hterr and δ13Cterr).
- **XRF_analysis_Peixao.R**: processes XRF elemental data (clr transformation), performs PCA, broken stick analysis and ANOSIM, and generates PC1 and Kclr time-series plots with 200-year smoothing. 
- **XRF_Peixao_clean.xlsx**: input file for `XRF_analysis_Peixao.R`, containing selected elemental data (cps) used in the study. 
- **Plots_figure_5.R**: processes paleoclimate data. The script computes z-score normalization for δ²Hterr (dD) data and visualizes temporal variability using conditional areas fills (positive vs negative anomalies). Additionally, it bins NAO-like index data into 500-year intervals, calculates median values, and displays them as bar plots to emphasize centennial-scale variability. 
- **Data_figure_5.xlsx**: input file for `Plots_figure_5.R` containing: (1) sheet 1, δ²H (dD) data used for z-score standardization and visualization; and (2) sheet 2, NAO-like index data used for binning and median aggregation.
- **Statistical_assessment_deuterio_carbon_covariation.R**: this script evaluates the covariation between δ2Hterr (dD) and δ13Cterr time series using multiple statistical approaches: (1) Pearson correlation, (2) linear detrending and correlation of residuals, and (3) first-differencing to assess high-frequency variability.
- **NMF_Peixao_2ka.m**: performs non-negative matrix factorization (NMF) on n-alkane distributions. The script requires an external NMF function (see below).
- **Peixao_2000_alkanes.xlsx**: input file for `NMF_Peixao_2ka.m`, containing n-alkane concentration (ng/g).
- **NMF_Peixao_11ka.m**: this script performs Non-negative Matrix Factorization (NMF) on n-alkane distributions from the Peixao sedimentary record and integrates complementary analyses to reconstruct paleoclimate signals and interpret organic matter sources. 
- **Peixao_11ka_alkanes.xlsx**: input file for `NMF_Peixao_11ka.m`. This file contains all datasets required to reproduce the NMF-based analyses of n-alkane distributions and associated isotope and proxy data from the Peixao sedimentary record (~11 ka). The excel file is organized into multiple sheets: (1) n-alkane distributions, (2) age, (3) n-alkane indices, (4) n-alkane hydrogen isotopes, (5) n-alkane carbon isotopes, and (6) distributions of the modern vegetation.

## Output Files

### From `Plots_alkanes.R`:
- Generates three ggplot2 figures: bar plot of modern n-alkane distributions by endmember group, time-series plot of reconstructed precipitation δ²H with uncertainty band, and two independent time-series plots of terrestrial δ²H and δ13C.
  
### From `XRF_analysis_Peixao.R`:
- `XRF_clr.csv` (clr-transformed dataset containing elemental compositions used for PCA and subsequent analyses); biplot of PC1 vs PC2 including loadings and confidence ellipses; broken-stick scree plot; ANOSIM output (plot and R statistic and significance values); time-series plot of PC1 vs age (raw and smoothed); and time-series plot of clr-transformed K vs age (smoothed).

### From `Plots_figure_5.R`: 
- Generates two base R plots: a continuous z-score time series with filled anomalies and a bar-style representation of binned NAO variability.

### From `Statistical_assessment_deuterio_carbon_covariation.R`: 
- This script does not generate external files by default. Results are returned as: (a) Console output (Pearson correlation coefficints and associated statistical tests, and correlation results for detrended (residual) series and for first_differenced series), (b) Graphical output (scatter plot with linear regression line), and (c) Diagnostic plots: autocorrelation functions (ACF) for both time series.

### From `NMF_Peixao_2ka.m`:
- NMF-derived endmember distributions from the Peixão dataset (Peixao_2000_alkanes). 
- Evaluates models with 1-4 endmembers and compares residual variance with PCA.
- Figures: constraints, scree plot, endmember distributions.

### From `NMF_Peixao_11ka.m`:
- NMF-derived endmember distributions from the Peixão dataset (Peixao_11ka_alkanes). 
- Evaluates models with 1-4 endmembers and compares residual variance with PCA.
- Figures: constraints, scree plot, endmember distributions, temporal evolution of endmember contributions.
- Projects NMF results into PCA space for comparison, which incorporates independent environmental indices (ACL, CPI and Paq). 
- Model evaluation: sum of squared errors (SSE), explained variance (R2) and relative improvement between models.
- Reconstructs endmember-specific isotope signals (carbon and hydrogen from n-alkanes) using NMF weights.
- Computes terrestrial isotope signal by combining endmembers 1 and 2 (EM1+EM2). 
- Applies moving-average smoothing and generates time-series plots.
- Modern plant n-alkane distributions comparison with NMF-derived endmembers using: Pearson and Spearman correlations, Non-negative least squares (NNLS), and Bray-Curtis similarity. 

## Usage

### General Setup (R):
- Ensure R and packages are installed.
- Place input files in yout working directory or update paths in scripts.
- Run scripts via `source("script_name.R")` in R.
- Outputs are generated as figures displayed in the R environment unless otherwise specified.

### General Setup (MATLAB):
- Ensure MATLAB (R2023b or later) is installed with required toolboxes (if applicable).
- Clone or download the repository and set the working directory to the project folder.
- Place input files (e.g., `.mat` or `.xlsx`) in the specified directories or update file paths within the scripts.
- Run scripts from the MATLAB Command Window or Editor (e.g., `script_name`).
- Outputs are generated as figures displayed in the MATLAB environment unless otherwise specified.
- Set the random seed for reproducibility: `rng(1)`.

### Run `XRF_analysis_Peixao.R`: 
- To export the clr-transformed dataset, add the following line to the script: `write.csv(xrf_clr, "XRF_clr.csv`, row.names = FALSE)

### Run `Statistical_assessment_deuterio_carbon_covariation.R`: 
- Uses the input file `Data_figure_5.xlsx` (sheet 1), which contains δ²H and δ¹³C isotope data from the terrestrial endmember derived from the NMF analysis.

### Run `NMF_Peixao_alkanes.m`: 
- Requires the input file `Peixao_2000_alkanes.xlsx`.
- Requires the external function `NMF.m` from Polissar: Polissar, P; Karp, A. Tyler; D'Andrea, William (2025), “Mixed messages: Unmixing sedimentary molecular distributions reveals source contributions and isotopic values”, Mendeley Data, V2, doi: 10.17632/3hymgv47jv.2. The `NMF.m` file must be downloaded separately and placed in the same directory as the script or added to the MATLAB path before execution.
- Run the script in MATLAB; figures will be generated in the MATLAB environment.

### Run `NMF_Peixao_11ka.m`: 
- Requires the input file `Peixao_11ka_alkanes.xlsx`.
- Requires the external function `NMF.m` from Polissar et al. (2025), cited above. 

### Run `Plots_alkanes.m`: 
- Uses three independent input files (set paths in the USER INPUT section).
- Data associated with figures are linked to a manuscript currently under review and will be accessible in the PANGAEA repository upon publication. 

## Notes
- **Input Data**: Contact the corresponding author (m.eugenia.fernandezp@udc.es) for supporting information or full datasets.
- **Customization**: Adjust parameters for other datasets.
- **Reproducibility**: All scripts are designed to be fully reproducible using the provided datasets and specified external dependencies.

## References
The MATLAB implementation of non-negative matrix factorization (NMF) used in this repository is based on the code provided by Polissar, P.J., Karp, A.T., D’Andrea, W.J., 2025. Mixed messages: Unmixing sedimentary molecular distributions reveals source contributions and isotopic values. Geochimica et Cosmochimica Acta 396, 122–134. https://doi.org/10.1016/j.gca.2025.03.001. The corresponding `NMF.m` function must be obtained from the original source: Polissar, P; Karp, A. Tyler; D'Andrea, William (2025), “Mixed messages: Unmixing sedimentary molecular distributions reveals source contributions and isotopic values”, Mendeley Data, V2, doi: 10.17632/3hymgv47jv.2 [dataset]

## Citation
If you use this repository, please cite the associated manuscript (in preparation): Fernández-Pérez, U., et al. (in preparation). North Atlantic atmospheric circulation controlled the hydroclimate in western Iberia during the Holocene. 
A full reference will be updated upon publication.

## Contact
For questions, contact U. Fernández-Pérez (m.eugenia.fernandezp@udc.es).

