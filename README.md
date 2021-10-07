<!-- badges: starts -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

# Omics Untargeted Key Script *(OUKS)* <img src="GH logo .png" align="right" height="250" width="300"/> 
## Brief Description :old_key:
R based open-source collection of scripts called *OUKS* (*Omics Untargeted Key Script*) providing comprehensive nine step LC-MS untargeted metabolomic profiling data processing toolbox :toolbox:

Script | Purpose
------------ | -------------
[1. Randomization.R](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/1.%20Randomization.R) | experimental design and sample randomization
[2. Integration.R](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/2.%20Integration.R) | peaks integration and time alignment
[3. Imputation.R](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/3.%20Imputation.R) | missing value imputation (MVI) and artifacts removal
[4. Correction.R](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/4.%20Correction.R) | signal drift correction and batch effect removal
[5. Annotation.R](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/5.%20Annotation.R) | feature annotation and tentative identification by database search
[6. Filtering.R](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/6.%20Filtering.R) | peaks filtering for quality checking and accounting of technical variation
[7. Normalization.R](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/7.%20Normalization.R) | data normalization and adjusting of biological variation
[8. Grouping.R](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/8.%20Grouping.R) | peaks grouping and molecular features clustering
[9. Statistics.R](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/9.%20Statistics.R) | statistical analysis and hypothesis testing

## Table of Contents :clipboard:
- Instruction and introduction into the *OUKS* toolbox is provided by [Basic tutorial](https://github.com/plyush1993/OUKS/blob/main/Basic%20tutorial.pdf) file. [Session info](https://github.com/plyush1993/OUKS/blob/main/Session%20Info.txt) and [installed packages](https://github.com/plyush1993/OUKS/blob/main/Used%20packages.pdf) are listed in corresponding files.
- [Scripts](https://github.com/plyush1993/OUKS/tree/main/Scripts%20(R)) with comments, notes and references are stored in Scripts folder at a previously defined order along with [code](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/Figures%20for%20OUKS.R) for plotting figures associated with [article](https://doi.org/10.1021/acs.jproteome.1c00392).
- [MS2 spectra](https://github.com/plyush1993/OUKS/tree/main/MS2%20spectra%20(mzXML)) for selected potential biomarkers of bladder cancer are stored in mzXML format at corresponding folder.
- [Datasets](https://github.com/plyush1993/OUKS/tree/main/Datasets%20(csv)) in .csv and [other files](https://github.com/plyush1993/OUKS/tree/main/Auxiliary%20files%20(RData)) (.RData, .R) are available for reproducibility from corresponding folders. Files descriptions are provived by [Roadmap file](https://github.com/plyush1993/OUKS/blob/main/Roadmap.pdf).
- [Report](https://github.com/plyush1993/OUKS/tree/main/Report%20(Rmd)) in [.Rmd](https://github.com/plyush1993/OUKS/blob/main/Report%20(Rmd)/Report%20example%20OUKS.Rmd), [.pdf](https://github.com/plyush1993/OUKS/blob/main/Report%20(Rmd)/Report-example-OUKS.pdf) and [.docx](https://github.com/plyush1993/OUKS/blob/main/Report%20(Rmd)/Report-example-OUKS.docx) formats were provided as an example to reproduce the *OUKS* code script.

## Requirements :building_construction:
The only requirements are to be familiar with the basic syntax of the R language, PC with Internet connection and Windows OS (desirable), [RStudio](https://www.rstudio.com/products/rstudio/download/) and [R](https://cloud.r-project.org/) (≥ 4.0.0).

## Release notes :new:
### **V. 0.1** :one: 2021.06.23
* Freely available at [link](https://doi.org/10.1021/acs.jproteome.1c00392) (Supporting Information File 2).
### **V. 0.2** :two: 2021.07.12
* "9. Statistics": Outlier detection method implementation (by Mahalanobis distance) via ClassDiscovery package (3.3.13, CRAN) was added. OutlierDetection package require spatstat package version 1.64-1 (CRAN).
* "9. Statistics": Add adjusted p-value for multiple comparisons in all cases.
* "9. Statistics": Multigroup Fold Change (structToolbox package) was replaced by base packages implementation.
* "7. Normalization": Add adjusted p-value for multiple comparisons in all cases.
* "4. Correction": Add PCA with gradient color.
* [“Deleted functionality.R”](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/Deleted%20functionality.R) was created for storing deleted code strings.
### **V. 0.3** :three: 2021.07.17
* "7. Normalization": GAM (mgcv, 1.8-32, CRAN) and GAMM (gamm4, 0.2-6, CRAN) were added as new biological factor adjustment algorithms.
* "9. Statistics": Section “Signal Modeling” was added for LM, LMM, GAM (mgcv, 1.8-32, CRAN), GAMM (gamm4, 0.2-6, CRAN) and some other nonlinear functions for Dose-Response curve analysis (drc, 3.0-1, CRAN) modeling.
* "9. Statistics": In section “Time series” Dose-Response curve analysis and modeling was added (DRomics, 2.2-0, CRAN).
### **V. 0.4** :four: 2021.08.04
* "5. Annotation": mWISE (0.1.0, GH, forked from b2slab/mWISE to plyush1993/mWISE and depends were manually changed to R (>= 4.0)). 
* "9. Statistics": Add tdfdr (0.1, GH) for two-dimensional false discovery rate control in filtration and multigroup analysis.
* All scripts (from 5. Annotation to 9. Statistics) and files were updated.
### **V. 0.5** :five: 2021.08.18
* "4. Correction": Add PC-PR2 for correction evaluation (pcpr2, 0.0.0.1, GH).
* "9. Statistics": Add PC-PR2 and PVCA for multigroup analysis.
### **V. 0.6** :six: 2021.08.31
* "4. Correction": Add box-plot, mean Silhouette Score and One-Sample Test metric.
* "7. Normalization": Box-plot construction updated.
### **V. 0.7** :seven: 2021.09.10
* "5. Annotation": metID (1.1.0, GH) for database identification from peak table.
* "7. Normalization": GPBoost (gpboost, 0.6.7, CRAN) with boosting and mixed effects boosting were added as new biological factor adjustment algorithms.
* "9. Statistics": Slightly changed Fold Change calculations and canonical limma implementation was added.
### **V. 0.8** :eight: 2021.10.07
* "9. Statistics": In section “Time series” TOXcms (1.0.3, GH) and timeOmics (1.0.1, BC) were added, also DRomics part was updated.
* Add reports (by R Markdown, folder [Report (Rmd)](https://github.com/plyush1993/OUKS/tree/main/Report%20(Rmd))).

## Citation :link:
*OUKS* has been published in the [Journal of Proteome Research](https://pubs.acs.org/journal/jprobs). If you use this software to analyze your own data, please cite it as below, thanks:

Ivan V. Plyushchenko, Elizaveta S. Fedorova, Natalia V. Potoldykova, Konstantin A. Polyakovskiy, Alexander I. Glukhov, Igor A. Rodin. *Omics Untargeted Key Script*: R‑Based Software Toolbox for Untargeted Metabolomics with Bladder Cancer Biomarkers Discovery Case Study, Journal of Proteome Research, 2021, [https://doi.org/10.1021/acs.jproteome.1c00392](https://doi.org/10.1021/acs.jproteome.1c00392).

## Contact :memo:
Please send any comment, suggestion or question you may have to the author (Mr. Ivan Plyushchenko :man_scientist:): :e-mail: plyushchenko.ivan@gmail.com, <img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png"> [0000-0003-3883-4695](https://orcid.org/0000-0003-3883-4695).
