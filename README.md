<!-- badges: starts -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/plyush1993/OUKS/blob/main/CHANGELOG.md)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![](https://img.shields.io/badge/devel%20version-1.8-yellow.svg)](https://github.com/plyush1993/OUKS/)
[![](https://img.shields.io/badge/doi-10.1021/acs.jproteome.1c00392-blueviolet.svg)](https://doi.org/10.1021/acs.jproteome.1c00392)
![](https://img.shields.io/github/languages/code-size/plyush1993/OUKS)
![](https://img.shields.io/github/repo-size/plyush1993/OUKS)
<!-- badges: end -->

# Omics Untargeted Key Script *(OUKS)* <img src="GH logo .png" align="right" height="250" width="300"/> 
## Brief Description :old_key:
R based open-source collection of scripts called :red_circle:*OUKS*:large_blue_circle: (*Omics Untargeted Key Script*) providing comprehensive nine step LC-MS untargeted metabolomic profiling data processing toolbox :toolbox:

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
- Instruction and introduction into the :red_circle:*OUKS* toolbox:large_blue_circle: is provided by [Basic tutorial](https://github.com/plyush1993/OUKS/blob/main/Basic%20tutorial.pdf) file. [Session info](https://github.com/plyush1993/OUKS/blob/main/Session%20Info.txt) and [installed packages](https://github.com/plyush1993/OUKS/blob/main/Used%20packages.pdf) are listed in corresponding files.
- [Scripts](https://github.com/plyush1993/OUKS/tree/main/Scripts%20(R)) with comments, notes and references are stored in Scripts folder at a previously defined order along with [code](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/Figures%20for%20OUKS.R) for plotting figures associated with [article](https://doi.org/10.1021/acs.jproteome.1c00392).
- [MS2 spectra](https://github.com/plyush1993/OUKS/tree/main/MS2%20spectra%20(mzXML)) for selected potential biomarkers of bladder cancer are stored in mzXML format at corresponding folder.
- [Datasets](https://github.com/plyush1993/OUKS/tree/main/Datasets%20(csv)) in .csv and [other files](https://github.com/plyush1993/OUKS/tree/main/Auxiliary%20files%20(RData)) (.RData, .R) are available for reproducibility from corresponding folders. Files descriptions are provived by [Roadmap file](https://github.com/plyush1993/OUKS/blob/main/Roadmap.pdf). Raw data (.CDF format) are available from [Metabolomics Workbench Repository](https://www.metabolomicsworkbench.org/), study ID: [ST001682](http://doi.org/10.21228/M8ZT4C). [Metadata](https://github.com/plyush1993/OUKS/blob/main/Datasets%20(csv)/metadata.csv) table is also provided.
- [Report](https://github.com/plyush1993/OUKS/tree/main/Report%20(Rmd)) in [.Rmd](https://github.com/plyush1993/OUKS/blob/main/Report%20(Rmd)/Report%20example%20OUKS.Rmd), [.pdf](https://github.com/plyush1993/OUKS/blob/main/Report%20(Rmd)/Report-example-OUKS.pdf) and [.docx](https://github.com/plyush1993/OUKS/blob/main/Report%20(Rmd)/Report-example-OUKS.docx) formats were provided as an example to reproduce the *OUKS* code script.
- [Changelog](https://github.com/plyush1993/OUKS/blob/main/CHANGELOG.md) is provided.

## Requirements :building_construction:
The only requirements are to be familiar with the basic syntax of the R language, PC with Internet connection and Windows OS (desirable), [RStudio](https://www.rstudio.com/products/rstudio/download/) and [R](https://cloud.r-project.org/) (≥ 4.0.0).

## Citation :link:
*OUKS* has been published in the [Journal of Proteome Research](https://pubs.acs.org/journal/jprobs). If you use this software to analyze your own data, please cite it as below, thanks:

Ivan V. Plyushchenko, Elizaveta S. Fedorova, Natalia V. Potoldykova, Konstantin A. Polyakovskiy, Alexander I. Glukhov, Igor A. Rodin. *Omics Untargeted Key Script*: R‑Based Software Toolbox for Untargeted Metabolomics with Bladder Cancer Biomarkers Discovery Case Study, Journal of Proteome Research, 2021, [https://doi.org/10.1021/acs.jproteome.1c00392](https://doi.org/10.1021/acs.jproteome.1c00392).

## Contact :memo:
Please send any comment, suggestion or question you may have to the author (Mr. Ivan Plyushchenko :man_scientist:): :e-mail: plyushchenko.ivan@gmail.com, <img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png"> [0000-0003-3883-4695](https://orcid.org/0000-0003-3883-4695).
