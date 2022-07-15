# Changelog :new:
## **V. 0.1** 
:spiral_calendar: 2021.04.29
* Creation date :clapper:
## **V. 1.0** 
:spiral_calendar: 2021.06.23
* Freely available at [link](https://doi.org/10.1021/acs.jproteome.1c00392) (Supporting Information File 2).
## **V. 1.1** 
:spiral_calendar: 2021.07.12
### Added
* "9. Statistics": Outlier detection method implementation (by Mahalanobis distance) via ClassDiscovery package (3.3.13, CRAN) was added. OutlierDetection package require spatstat package version 1.64-1 (CRAN). Adjusted p-value for multiple comparisons in all cases were added.
* "7. Normalization": Add adjusted p-value for multiple comparisons in all cases.
* "4. Correction": Add PCA with gradient color.
* [“Deprecated functionality.R”](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/Deprecated%20functionality.R) was created for storing deleted/deprecated code strings.
### Changed
* "9. Statistics": Multigroup Fold Change (structToolbox package) was replaced by base packages implementation.
## **V. 1.2**
:spiral_calendar: 2021.07.17
### Added
* "7. Normalization": GAM (mgcv, 1.8-32, CRAN) and GAMM (gamm4, 0.2-6, CRAN) were added as new biological factor adjustment algorithms.
* "9. Statistics": Section “Signal Modeling” was added for LM, LMM, GAM (mgcv, 1.8-32, CRAN), GAMM (gamm4, 0.2-6, CRAN) and some other nonlinear functions for Dose-Response curve analysis (drc, 3.0-1, CRAN) modeling. In section “Time series” Dose-Response curve analysis and modeling was added (DRomics, 2.2-0, CRAN).
## **V. 1.3**
:spiral_calendar: 2021.08.04
### Added
* "5. Annotation": mWISE (0.1.0, GH, forked from b2slab/mWISE to plyush1993/mWISE and depends were manually changed to R (>= 4.0)). 
* "9. Statistics": Add tdfdr (0.1, GH) for two-dimensional false discovery rate control in filtration and multigroup analysis.
### Fixed
* All scripts (from 5. Annotation to 9. Statistics) and files were updated.
## **V. 1.4**
:spiral_calendar: 2021.08.18
### Added
* "4. Correction": Add PC-PR2 for correction evaluation (pcpr2, 0.0.0.1, GH).
* "9. Statistics": Add PC-PR2 and PVCA for multigroup analysis.
## **V. 1.5** 
:spiral_calendar: 2021.08.31
### Added
* "4. Correction": Add box-plot, mean Silhouette Score and One-Sample Test metric.
### Fixed
* "7. Normalization": Box-plot construction updated.
## **V. 1.6** 
:spiral_calendar: 2021.09.10
### Added
* "5. Annotation": metID (1.1.0, GH) for database identification from peak table.
* "7. Normalization": GPBoost (gpboost, 0.6.7, CRAN) with boosting and mixed effects boosting were added as new biological factor adjustment algorithms.
### Changed
* "9. Statistics": Slightly changed Fold Change calculations and canonical limma implementation was added.
## **V. 1.7**
:spiral_calendar: 2021.10.07
### Added
* "9. Statistics": In section “Time series” TOXcms (1.0.3, GH) and timeOmics (1.0.1, BC) were added, also DRomics part was updated.
* "4. Correction": Add WaveICA2.0 (0.1.0, GH) correction method.
* Add reports (by R Markdown, folder [Report (Rmd)](https://github.com/plyush1993/OUKS/tree/main/Report%20(Rmd))).
## **V. 1.8**
:spiral_calendar: 2021.10.26
### Added
* "4. Correction": Two features (tf) and subtraction (s) modes were implemented for QC-BT/DT/KNN/GBM in all combinations.
* [Releases](https://github.com/plyush1993/OUKS/releases) page and [tags](https://github.com/plyush1993/OUKS/tags).
### Changed
* "4. Correction": QC-norm was modified and division-/subtraction-based versions are provided (original version of the function is available at [“Deprecated functionality.R”](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/Deprecated%20functionality.R)).
### Fixed
* "6. Filtering": Small changes in Descriptive Statistics, Missing Value and Blank filtering approaches.
* "4. Correction": Fixed dependencies for QC-SVR, QC-LOESS, QC-RF and updated one sample test quality metric.
## **V. 1.9**
:spiral_calendar: 2021.11.19
### Added
* "9. Statistics": LM, LMM, GLM, GLMM and correlation analysis (MWASTools, 1.12.0, BC; glmmsr, 0.2.3, CRAN) are added in "MWAS/ANCOVA" section (with categorical feature as dependent variable). Add PLS and sPLS for multigroup analysis.
## **V. 1.10**
:spiral_calendar: 2022.05.10
### Fixed
* Updated to R 4.1.2.
### Changed
* "3. Imputation": Class-specific NRMSE calculation is moved to [“Deprecated functionality.R”](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/Deprecated%20functionality.R)).
* "6. Filtering": Minor changes to reduce dependencies.
* "4. Correction": Script was rewritten to reduce dependencies. 
### Added
* "2. Integration": cpc for filtering low quality peaks, Paramounter optimization algorithm, “groupFeatures” function from xcms, data smoothing and refinement of the centroid’s m/z values (by MSnbase package) were added. 
* "3. Imputation": tWLSA multivariate MVI algorithm.
* "4. Correction": SERRF, TIGERr, QC-LOESS and QC-RF (in subtraction mode) algorithms were added. Manual ANOVA version.
* "7. Normalization": Manual LM modeling version.
* "9. Statistics": Nested feature selection and Boruta ("Multiple statistical filtration" section), Nested cross-validation ("Classification Machine Learning task" and "Regression Machine Learning task" sections), and DBSCAN, HDBSCAN, Spectral Clustering, UMAP, MCLUST, MDS, LLE, IsoMap, Laplacian Score, Diffusion Maps, kernel PCA, sparse PCA, ICA, FA, NMF, PAM, CLARA, Fuzzy Clustering ("Unsupervised Data Projection" section) were added. Manual LM, GLM modeling versions.
## **V. 1.10.1**
:spiral_calendar: 2022.06.29
### Added
CompDb.Hsapiens.HMDB.5.0 (from [jorainer/MetaboAnnotationTutorials](https://github.com/jorainer/MetaboAnnotationTutorials/releases/tag/2021-11-02))
## **V. 1.11**
:spiral_calendar: 2022.07.12
### Fixed
* Some packages were updated.
### Changed
* "7. Normalization": GBM, GBMM modeling were moved to [“Deprecated functionality.R”](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/Deprecated%20functionality.R)).
### Added
* "2. Integration": Target peak integration in EIC chromatograms was implemented.
* "3. Imputation": Automatic missingness type determination and imputation were implemented by MAI package.
* "4. Correction": D-Ratio metric, RLA-plot, correlogram, 2-factors PCA are added.
* "5. Annotation": MetaboAnnotation for database-based annotation was added.
* "6. Filtering": D-Ratio filtering was added.
* "7. Normalization": Adaptive Box-Cox Transformation normalization was implemented.
* "9. Statistics": Hotelling Ellipse with T-squared statistic and DModX metric were added for outlier detection. polyPK for pharmacokinetic analysis was implemented. Some minor changes.
