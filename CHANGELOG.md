# Changelog :new:
## **V. 0.1** 
* :spiral_calendar: 2021.04.29
* Creation date :clapper:
## **V. 1.0** 
* :spiral_calendar: 2021.06.23
* Freely available at [link](https://doi.org/10.1021/acs.jproteome.1c00392) (Supporting Information File 2).
## **V. 1.1** 
* :spiral_calendar: 2021.07.12
### Added
* "9. Statistics": Outlier detection method implementation (by Mahalanobis distance) via ClassDiscovery package (3.3.13, CRAN) was added. OutlierDetection package require spatstat package version 1.64-1 (CRAN).
* "9. Statistics": Add adjusted p-value for multiple comparisons in all cases.
* "7. Normalization": Add adjusted p-value for multiple comparisons in all cases.
* "4. Correction": Add PCA with gradient color.
* [“Deleted functionality.R”](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/Deleted%20functionality.R) was created for storing deleted code strings.
### Changed
* "9. Statistics": Multigroup Fold Change (structToolbox package) was replaced by base packages implementation.
## **V. 1.2**
* :spiral_calendar: 2021.07.17
### Added
* "7. Normalization": GAM (mgcv, 1.8-32, CRAN) and GAMM (gamm4, 0.2-6, CRAN) were added as new biological factor adjustment algorithms.
* "9. Statistics": Section “Signal Modeling” was added for LM, LMM, GAM (mgcv, 1.8-32, CRAN), GAMM (gamm4, 0.2-6, CRAN) and some other nonlinear functions for Dose-Response curve analysis (drc, 3.0-1, CRAN) modeling.
* "9. Statistics": In section “Time series” Dose-Response curve analysis and modeling was added (DRomics, 2.2-0, CRAN).
## **V. 1.3**
* :spiral_calendar: 2021.08.04
### Added
* "5. Annotation": mWISE (0.1.0, GH, forked from b2slab/mWISE to plyush1993/mWISE and depends were manually changed to R (>= 4.0)). 
* "9. Statistics": Add tdfdr (0.1, GH) for two-dimensional false discovery rate control in filtration and multigroup analysis.
### Fixed
* All scripts (from 5. Annotation to 9. Statistics) and files were updated.
## **V. 1.4**
* :spiral_calendar: 2021.08.18
### Added
* "4. Correction": Add PC-PR2 for correction evaluation (pcpr2, 0.0.0.1, GH).
* "9. Statistics": Add PC-PR2 and PVCA for multigroup analysis.
### **V. 1.5** 
* :spiral_calendar: 2021.08.31
### Added
* "4. Correction": Add box-plot, mean Silhouette Score and One-Sample Test metric.
### Fixed
* "7. Normalization": Box-plot construction updated.
## **V. 1.6** 
* :spiral_calendar: 2021.09.10
### Added
* "5. Annotation": metID (1.1.0, GH) for database identification from peak table.
* "7. Normalization": GPBoost (gpboost, 0.6.7, CRAN) with boosting and mixed effects boosting were added as new biological factor adjustment algorithms.
### Changed
* "9. Statistics": Slightly changed Fold Change calculations and canonical limma implementation was added.
## **V. 1.7**
* :spiral_calendar: 2021.10.07
### Added
* "9. Statistics": In section “Time series” TOXcms (1.0.3, GH) and timeOmics (1.0.1, BC) were added, also DRomics part was updated.
* "4. Correction": Add WaveICA2.0 (0.1.0, GH) correction method.
* Add reports (by R Markdown, folder [Report (Rmd)](https://github.com/plyush1993/OUKS/tree/main/Report%20(Rmd))).
## **V. 1.8**
* :spiral_calendar: 2021.10.25
### Added
* "4. Correction": Two features (tf) and subtraction (s) modes were implemented for QC-BT/DT/KNN/GBM in all combinations.
### Changed
* "4. Correction": QC-norm was modified and division-/subtraction-based versions are provided (original version of the function is available at [“Deleted functionality.R”](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/Deleted%20functionality.R)).
### Fixed
* "6. Filtering": Small changes in Descriptive Statistics, Missing Value and Blank filtering approaches.
* "4. Correction": Fixed dependencies for QC-SVR, QC-LOESS and QC-RF.
