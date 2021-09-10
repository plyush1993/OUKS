# *OUKS* <img src="Spider 2.jpg" align="right" height="304" width="280"/> 
# Brief Description
R based open-source collection of scripts called *OUKS* (*Omics Untargeted Key Script*) providing comprehensive nine step LC-MS untargeted metabolomic profiling data processing:

1. experimental design and sample randomization
2. peaks integration and time alignment
3. missing value imputation (MVI) and artifacts removal
4. signal drift correction and batch effect removal
5. feature annotation and tentative identification by database search
6. peaks filtering for quality checking and accounting of technical variation
7. data normalization and adjusting of biological variation
8. peaks grouping and molecular features clustering
9. statistical analysis and hypothesis testing

# Table of Contents
- Instruction and introduction into the *OUKS* toolbox is provided by [Basic tutorial](https://github.com/plyush1993/OUKS/blob/main/Basic%20tutorial.pdf) file. [Session info](https://github.com/plyush1993/OUKS/blob/main/Session%20Info.txt) and [installed packages](https://github.com/plyush1993/OUKS/blob/main/Used%20packages.pdf) are listed in corresponding files.
- [Scripts](https://github.com/plyush1993/OUKS/tree/main/Scripts%20(R)) with comments, notes and references are stored in Scripts folder at a previously defined order along with code for figures construction.
- [MS2 spectra](https://github.com/plyush1993/OUKS/tree/main/MS2%20spectra%20(mzXML)) for selected potential biomrkers of bladder cancer are stored in mzXML format at corresponding folder.
- [Datasets](https://github.com/plyush1993/OUKS/tree/main/Datasets%20(csv)) in .csv and [other files](https://github.com/plyush1993/OUKS/tree/main/Auxiliary%20files%20(RData)) (.RData, .R) are available for reproducibility from corresponding folders. Files descriptions are provived by [Roadmap file](https://github.com/plyush1993/OUKS/blob/main/Roadmap.pdf).

# Requirements
The only requirements are to be familiar with the basic syntax of the R language, PC with Internet connection and Windows OS (desirable), [RStudio](https://www.rstudio.com/products/rstudio/download/) and [R](https://cloud.r-project.org/) (≥ 4.0.0).

# Release notes
## **V. 0.1** 
* freely available at [https://doi.org/10.1021/acs.jproteome.1c00392](https://doi.org/10.1021/acs.jproteome.1c00392) (Supporting Information File 2)
## **V. 0.2**  
* "9. Statistics": Outlier detection method implementation (by Mahalanobis distance) via ClassDiscovery package (3.3.13, CRAN) was added. OutlierDetection package require spatstat package version 1.64-1 (CRAN).
* "9. Statistics": Add adjusted p-value for multiple comparisons in all cases.
* "9. Statistics": Multigroup Fold Change (structToolbox package) was replaced by base packages implementation.
* "7. Normalization": Add adjusted p-value for multiple comparisons in all cases.
* "4. Correction": Add PCA with gradient color.
* “Deleted functionality.R” was created for storing deleted code strings.
## **V. 0.3** 
* "7. Normalization": GAM (mgcv, 1.8-32, CRAN) and GAMM (gamm4, 0.2-6, CRAN) were added as new biological factor adjustment algorithms.
* "9. Statistics": Section “Signal Modeling” was added for LM, LMM, GAM (mgcv, 1.8-32, CRAN), GAMM (gamm4, 0.2-6, CRAN) and some other nonlinear functions for Dose-Response curve analysis (drc, 3.0-1, CRAN) modeling.
* "9. Statistics": In section “Time series” Dose-Response curve analysis and modeling was added (DRomics, 2.2-0, CRAN).
## **V. 0.4** 
* "5. Annotation": mWISE (0.1.0, GitHub, forked from b2slab/mWISE to plyush1993/mWISE and depends were manually changed to R (>= 4.0). 
* "9. Statistics": Add tdfdr (0.1, GitHub) for two-dimensional false discovery rate control in filtration and multigroup analysis.
* All scripts (from 5. Annotation to 9. Statistics) and files were updated.
## **V. 0.5** 
* "4. Correction": Add PC-PR2 for correction evaluation (pcpr2, 0.0.0.1, GitHub).
* "9. Statistics": Add PC-PR2 and PVCA for multigroup analysis.
## **V. 0.6** 
* "4. Correction": Add box-plot, mean Silhouette Score and One-Sample Test metric.
* "7. Normalization": box plots updated.
## **V. 0.7**
* "5. Annotation": metID (1.1.0, GitHub) for database identification from peak table.
* "7. Normalization": GPBoost (gpboost, 0.6.7, CRAN) with boosting and mixed effects boosting were added as new biological factor adjustment algorithms.
* "9. Statistics": GPBoost (gpboost, 0.6.7, CRAN) with boosting and mixed effects boosting were added as new algorithms in Section “Signal Modeling”.
     
## Citation
*OUKS* has been published in the [Journal of Proteome Research](https://pubs.acs.org/journal/jprobs). If you use this software to analyze your own data, please cite it as below, thanks:

Ivan V. Plyushchenko, Elizaveta S. Fedorova, Natalia V. Potoldykova, Konstantin A. Polyakovskiy, Alexander I. Glukhov, Igor A. Rodin. *Omics Untargeted Key Script*: R‑Based Software Toolbox for Untargeted Metabolomics with Bladder Cancer Biomarkers Discovery Case Study, Journal of Proteome Research, 2021, [https://doi.org/10.1021/acs.jproteome.1c00392](https://doi.org/10.1021/acs.jproteome.1c00392).

## Contact
Please send any comment, suggestion or question you may have to the author (Mr. Ivan Plyushchenko), email: plyushchenko.ivan@gmail.com.

