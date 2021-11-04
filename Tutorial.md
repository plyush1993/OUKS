# Omics Untargeted Key Scripts
![logo](/vignette/logo.png) 
## Preface
R based open-source collection of scripts called OUKS (Omics Untargeted Key Script) providing comprehensive nine step LC-MS untargeted metabolomic profiling data processing toolbox.

The session info snapshot, information about used packages with comments and scripts are available from https://github.com/plyush1993/OUKS. Each script consists of contents, comments, references and notes in places where functions should be adjusted for experimental parameters (usually “# adjust to your data”). Working directory was set in each script by setwd function, all data tables were loaded into environment in csv format.

Description and instruction for every function can be obtained by the code: ?function. Self-written functions are commented in the script. We utilized a special format for the name of every experimental injection file: “150. qc30 b6 QC30.CDF” for QC sample and “167. s64_2 b7 MS 45.CDF” for study file. Thus, the first value is a run order, then the index of sample with repeat number, batch value and clinical/experimental ID (or QC ID). This form allows to obtaining key information directly from the file name and automatically constructs data frame in the R environment. You can load a similar table to your environment manually if you prefer a different name format.

Raw data (.CDF format) and table with metadata (“metadata.csv”) are the only requirements to reproduce the entire code and are available from Metabolomics Workbench, study ID: ST001682 (http://doi.org/10.21228/M8ZT4C). Other input files (.csv data tables and .R/.RData files) can be generated in the appropriate script or downloaded from the GitHub repository. The entire processing scheme can be easily adapted to the researcher’s needs and their raw data or data tables. The information about used packages and some troubleshooting procedures are described in “installed packages.docx” file. Other useful information, references and significant additions are listed below by each script file name. 

To reproduce any script code, download all files to the working directory folder, set the working directory and load (install) the corresponding packages. The R Markdown document and reports files were provided as an example to reproduce the OUKS code script (available in “Report (Rmd)” folder).

The only requirements are to be familiar with the basic syntax of the R language, PC with Internet connection and Windows OS (desirable), RStudio and R (≥ 4.0.0).

OUKS has been published in the Journal of Proteome Research. If you use this software to analyze your own data, please cite it as below, thanks:

Ivan V. Plyushchenko, Elizaveta S. Fedorova, Natalia V. Potoldykova, Konstantin A. Polyakovskiy, Alexander I. Glukhov, Igor A. Rodin. Omics Untargeted Key Script: R‑Based Software Toolbox for Untargeted Metabolomics with Bladder Cancer Biomarkers Discovery Case Study, Journal of Proteome Research, 2021, https://doi.org/10.1021/acs.jproteome.1c00392.

Please send any comment, suggestion or question you may have to the author (Mr. Ivan Plyushchenko), email: plyushchenko.ivan@gmail.com.

## 1. Randomization 
All samples were analyzed at random order to prevent systematic bias [1,2]. Each analytical batch consisted of ten samples in two technical repeats (repeat samples were acquired after last tenth sample and repeats were analyzed at the same order as first repeat). The QC samples were acquired at the beginning of the batch and after every five injections (overall five QC samples for each batch). The code generates a random sequence of samples in accordance with the user's conditions: the number of samples, technical and biological repeats and batch size.

## 2. Integration
In our case 311 injects in 13 batches were acquired. All raw data were stored in a specific folder. Seven injects were selected for integration and alignment parameters optimization via IPO [3] (1 (batch 1, 1st QC from 5), 75 (3, 5), 101 (5, 1), 163 (7, 3), 219 (9, 4), 257 (11, 2), 311 (13, 3 of 3)) in order to decrease time of computing. QC samples spanning all study were randomly selected (at the beginning, end and middle of batch).

Best parameters were selected for XCMS processing [4] according to IPO optimization. Only bw and minfrac parameters were manually checked according to [5,6]. The “bw” parameter did not introduce any benefits in terms of the number of peaks. The decreasing of “minfrac” parameter sequentially increase the number of peaks and “minfrac” was finally set equal to minimum proportion of the smallest experimental group (0.2 for QC group). Peaks integration, grouping and retention time alignment were performed for all 311 sample injections (without blank injections). The final peak table was described by total number of peaks and fraction of missing values.

Warpgroup algorithm [7] for increasing of precision XCMS results integration was also tested. However, time of computing was approximately 6 days (test was performed via 5 random peak groups). The algorithm has not been implemented for this reason.

ncGTW algorithm [8] for increasing of precision XCMS alignment results was also implemented. Alternative way is checking different “bw “values.

Other approaches for XCMS parameters selection (Autotuner [9] in semi-automatic mode and MetaboAnalystR [10] in automatic) were also implemented as alternative options. Also, old version of XCMS functions (in consistent with xcmsSet objects instead of XCMSnExp objects in recent version) and new capabilities, that are provided by “refineChromPeaks” function and subset-based alignment, were integrated into the script.

## 3. Imputation
Artifacts checking was performed via MetProc [11] which is based on the comparison of missing rates (all peaks were retained). Nine methods of missing value imputation (MVI) were tested for previously obtained peak table, closely to [12]. The half minimum, KNN (k-nearest neighbors), PLS (partial least squares), PPCA (probabilistic principal component analysis, or other PCA-based), CART (Classification and Regression Tree), three types of RF (random forest) and QRILC (quantile regression for the imputation of left-censored missing data) algorithms were implemented for data table without any NA values and with randomly introducing NA values (proportion of introducing was equal to proportion of NA in the original dataset) [13]. The best algorithm was selected by the normalized root-mean-square error (NRMSE, value between 0 and 1) value close to zero (the RF implementation from missForest package). Assumption about decreasing of mean weighted NRMSE value after dividing the procedure for each group was sequentially tested and rejected (no statistically significant increment). Sum of squared error in Procrustes analysis on principal component scores (minimal value required) and correlation coefficient (maximum value required) were implemented as other quality metrics for MVI.

Moreover, other univariate MVI methods are implemented: converting NA values into the one single number (for example, 0), replacing by mean or median or minimum values and by random generated numbers. Values for random generation are produced by normal distribution with mean equal to the noise level (which was determined by IPO optimization or set manually) and standard deviation which was calculated from the noise (for example, 0.3*noise). 

## 4.	Correction
31 methods for signal drift correction and batch effects removal were implemented in this work, which include: model-based (WaveICA [14], WaveICA 2.0 [15], EigenMS [16], RUVSeqs family [17,18] (RUVSeqs, RUVSeqg, RUVSeqr), Parametric/Non-Parametric Combat, Ber, Ber-Bagging [19]), internal standards based [20,21] (SIS, NOMIS, CCMN, BMIS), QC metabolites based [21] (ruvrand, ruvrandclust), QC samples based regression algorithms are: cubic smoothing splines [22,23] (in two different implementation via pmp and notame packages respectively), local polynomial regression fitting (LOESS) and RF [24], support vector machine (SVM) [25], least-squares/ robust (LM/RLM) and censored (TOBIT) [26], LOESS or splines with mclust [27], 5 new algorithms each with single, two features and batchwise modes (all also in subtraction/division versions) – KNN (caret package), decision tree (rpart), bagging tree (ipred), XGB (xgboost) and catboost (catboost) gradient boostings and QC-NORM for between batches correction [28]. QC-NORM was implemented in two versions: division-based [28] and subtraction-based [29]. Moreover, sequential combinations of methods without consideration of batch number and QC-NORM were also performed. De facto, random forest method is equal to bagging for the regression with only one variable, but due to some differences in trees growing and pruning, modeling results differ. Also, some other algorithms were tested: Cubist (Cubist package), conditional trees (partykit) and smoothing splines with auto-determination of knots (mgcv). Modeling results for these algorithms did not outperform existing approaches and were not included into the script.

5 new methods in single division mode were performed according to the equations (1-3):
![1-3](/vignette/1-3.png)
where, i – is the index of metabolite, QC – is the index of QC samples (all samples are denoted without this index), M – is a machine learning model, is fitted by function f, y – is an intensity vector (dependent variable), x – is a run order vector (independent variable), F – is a vector of predicted values by the model M (correction factor), I’ – is a vector of a corrected intensity values, I – is a vector of an original intensity values.

Equations (1-3) are one of the basic algorithms for all QC sample based algorithms for the signal drift correction and are most similar to the statTarget package [24]. The MetNormalizer package [25] also performs clustering of features via Spearman correlation and the same correction operation (division). Correction algorithm, which is described by the eq. 1-3, was marked as single division mode (“d”).
![4](/vignette/4.png)
where, i – is the index of metabolite, QC – is the index of QC samples (all samples are denoted without this index), y – is an intensity vector, F – is a vector of predicted values by the model M (correction factor), I’ – is a vector of a corrected intensity values, I – is a vector of an original intensity values, mean – the function, that generates the mean value from the input vector.

Equation (4) is the modification for equation (3) of previously described algorithm (eq. 1-3), which is called single subtraction mode (eq. 1-2,4; “s”). The same mode was implemented in the BatchCorrMetabolomics [26] and notame [23] packages. Other mode in notame utilizes a correction factor with some changes (including prediction for the 1st QC sample and raw data).

Some QC sample based algorithms consider batch number. Batchwise mode (“bw”) is identical to the two types of single mode (eq. 1-3 division and 1-2,4 subtraction), but each equation was repeated for each batch separately. This implementation is similar to pmp package [22] (the difference with division mode in the correction factor, which includes the median value of the feature) and is closely to batchCorr package [27] (the differences are in the correction factor and clustering features by the mclust algorithm). 

Other way is an inclusion batch number into the regression model equation (into eq. 1, as in package BatchCorrMetabolomics [26]). This implementation was called Two Features mode (“tf”) and was implemented both for subtraction and division modes. 

Thus, totally 30 algorithms were proposed: each of 5 regression algorithms were implemented in 6 modes (3 (single, “bw”, “tf”) multiplied by 2 (division, subtraction)). All new correction methods in all modes were written with apply family functions and progress bar.

15 methods for the quality checking and comparison of corrections methods [26,28,30,31] were performed, which include: guided principal component analysis (gPCA), principal variance component analysis (PVCA), fraction of features with p-value below the significance level in the ANOVA model, mean Pearson correlation coefficient for QC samples, mean Silhouette score in HCA on PCs scores for QC samples, fraction of features with p-value below the significance level in the One-Sample test (sequential implementation of Shapiro-Wilk normality test and corresponding type of test (t-test or Wilcoxon)) for QC samples, fraction of features under the criterion of 30% RSD reproducibility in QC samples, PCA, hierarchical cluster analysis (HCA), box-plots and scatter plot (total intensity versus run order) data projection and visualization, mean Bhattacharyya distance for QC samples between batches, mean dissimilarity score for QC samples, mean classification accuracy for QC samples after hierarchical clustering on PCs scores, Partial R-square (PC-PR2) method. The best correction method should demonstrate the maximum values of the fraction of features under the criterion of 30% RSD reproducibility in QC samples and the mean Pearson correlation coefficient for QC samples. The values of other metrics should be the lowest. PCA, box-plots, scatter plots and HCA dendrogram should show the absence of tendency for the samples to cluster according to batch number or run order, and simultaneously QC samples should be placed almost between study samples.

## 5.	Annotation
Features in peak table (the output of XCMS processing) were annotated via CAMERA [32], RAMClustR [33], xMSannotator [34] and mWISE [35]. The sr and st parameters in RAMClustR were optimized by functions from [5]. Other settings which include mass detection accuracy in ppm and m/z values and time deviation in seconds were derived from IPO optimized parameters and manual inspection respectively. Annotations should be performed for files in the same folder as for XCMS integration and alignment, xcms objects (“xcms obj … .RData”) should be obtained by user. A search of isotope peaks, adducts, neutral losses and in-source fragments were included in all algorithms. The xMSannotator also was considered retention time and database search (HMDB, KEGG, LipidMaps, T3DB, and ChemSpider).

Some functions in xMSannotator package are not available for R version ≥ 4.0.0. For this reason, you should load file "fix for R 4.0 get_peak_blocks_modulesvhclust.R" in your environment (based on https://github.com/yufree/xMSannotator, other instructions are listed in the script). mWISE package forked from b2slab/mWISE to plyush1993/mWISE and depends were manually changed to R (>= 4.0) (available in “Auxiliary files (RData)” folder).

Each algorithm was characterized by a fraction of annotated features. The xMSannotator and mWISE provided the maximum value of annotated ions.

Also, metID [36] can be used for simple databases search from peak table. Databases from (https://github.com/jaspershen/demoData/tree/master/inst/ms2_database) should be copied in "ms2_database" folder in metID library folder.

## 6.	Filtering
Several most common options are available [1,23,37]: RSD based filtering by QC samples, annotation based (remove non-annotated features), by descriptive statistics for biological samples (by mean, max, min, median values of the feature or sd), by the fraction of missing values in the feature for biological samples, by the cutoff between median values of feature in QS or biological samples and blanks injections and by the Pearson correlation coefficient for samples (mainly for repeated injections). The mean values for every feature for all repeats are also computed. The mean values for two repeats were obtained in our study.

## 7.	Normalization
Five normalization methods are performed [38,39]: mass spectrometry total useful signal (MSTUS), probabilistic quotient normalization (PQN), quantile, median fold change and sample-specific factor. Also, several scaling and transformation algorithms are introduced: log transformation, auto- , range- , pareto- , vast- , level- , power-, etc. scaling.

Biological and phenotyping metadata (classification class, age, sex, BMI, etc.) as well as experimental conditions (run order, batch number, etc.) can be adjusted by linear model (LM) [40] or linear mixed effect model (LMM) [41] fitting. Both algorithms were implemented in the study.

The LMM approach fits the model according to a user-defined formula. In this study, model was constructed for each metabolite intensity as an dependent variable; age, class and sex – as fixed effects and batch ID – as a random effect. Then, the original data was employed as input and predicted values from the models were the final output.
![5-6](/vignette/5-6.png)

where, i – is the index of metabolite, j – is the index of independent variable j = 1,2, … n (fixed effects), Φ – is a LMM model is fitted by function f, β0 – is the intercept constant, ε – is the random error, β – is the coefficient of the regression model, x – is an independent variable, I – is a dependent variable (a vector of an original intensity values), μ0 – is the random effect, I’ – is a vector of an adjusted intensity values, that was predicted by model Φ.

The LM factorwise approach works in a similar way. At first, linear model was constructed for each metabolite intensity as a dependent variable and other covariates – as independent variables (Age, Sex, Batch, Class, Order in our study). One of the covariates should be noted as a feature of interest (this biological factor was kept, Class in our study). Then, residuals of models were obtained and design matrix was constructed for the feature of interest and constant intercept. For each model, design matrix was multiplied by regression coefficients of intercept and feature of interest and then was summed with residuals. Obtained values were summed with the difference between their mean values and mean values of the original data. The resulted values were the final output.
![7-9](/vignette/7-9.png)

where, i – is the index of metabolite, j – is the index of independent variable j = 1,2, … n (fixed effects), Φ – is a LM model is fitted by function f, β0 – is the intercept constant, ε – is the random error, β – is the coefficient of the regression model, x – is an independent variable, I – is a dependent variable (a vector of an original intensity values), I’ – is a vector of a transformed intensity values, I’’ – is a vector of an adjusted intensity value, M – is a design matrix, that was constructed for the feature of interest and constant intercept, coef – is a vector of regression coefficients of intercept and feature of interest from the model Φ, res – is a vector of the residuals of model Φ, mean – the function, that generates the mean value from the input vector.

Other type of LM adjustment was implemented in similar way to LMM (eq. 5-6 without random effect).

GAM, GAMM [42], GBM, GBMM [43] adjustment were also realized (as in eq. 5-6).

Four algorithms were utilized for evaluation of normalization methods, including: MA- , Box- , relative log abundance (RLA)– plots and mean/range of relative abundance between all features [21,44,45]. Biological adjustment methods were tested by classification accuracy via four machine learning (ML) algorithms (RF, linear SVM, nearest shrunken centroids (PAM), PLS), the high value of accuracy indicates successful removing of biological variation. ML parameters are listed in Section 9. PCA modeling was used for visualization of final adjustment. LM and LMM models were fitted in a similar way as described above (eq. 7,5). The fraction of p-values for the target variable and/or other covariates in the model less than the threshold (0.05) can serve as other character of the adjustment quality (the fraction should be the lowest for covariates and the largest for the target variable).

## 8.	Grouping
Molecular features from data table after integration and alignment can be clustered (grouped) by different algorithms in order to determine an independent molecular feature and dimensionality reduction. The pmd package [46] performs hierarchical cluster analysis for this purpose (and also provides reactomics). The notame package [23] and CROP algorithm [47] utilize the Pearson correlation coefficient. Final signals represent the largest/mean/median values of features groups. 

All .RData files that are needed to implement the CROP algorithm can be loaded from the GitHub repository.

## 9.Statistics
Statistical analysis part is divided into several logical blocks. Basic information and principles behind this chapter are available in [48-52]. Calculations are performed with basic functions, if no package is specified.

Outlier detection is performed by calculation of Mahalanobis/Euclidean distance for PCA scores (packages: OutlierDetection, pcaMethods or ClassDiscovery).

Statistical filtration block includes numerous methods for selection of informative variables: 

1. A combination of four ML models with stable feature extraction (SFE) and recursive feature selection (RFS) [53] (packages: caret, tuple). Lists of the top N (50,100, etc) important features from each model were selected according to the model-specific metric in order to perform SFE. The unique variables are subsequently extracted from the set of important features from all models that are matched in at least n times (from 2 to 4 or equal to 50%-100% frequency). RFS with Naïve Bayes(NB) or RF algorithm is applied to subsets of variables after SFE for collection final variable set. The N and n are determined experimentally. Filtration can be performed by SFE only or SFE with RFS.
The following parameters were applied when building the models: initial dataset was divided into two subsets by the 10 folds cross-validation with 10 times repeats for internal validation. The hyper-parameters tune length was 10. The four base ML algorithms were implemented: RF, SVM with radial kernel function, PLS and PAM. Only the most important value was tuned in each model (SVM – cost of constraints violation, RF – number of variables at each split, PLS – number of components, PAM – shrinkage threshold). All models were optimized by maximum mean accuracy metric across all cross-validation resampling. RFS was performed by backwards selection between all variables in subset. ML parameters in RFS were equal to those described above.
2. Filtration by VIP value from PLS (Orthogonal-PLS) model (package: ropls). The filtration threshold is determined experimentally.
3. Filtration by permutation importance from RF model (packages: permimp, party). The filtration threshold is determined experimentally.
4. Filtration by penalized or stepwise logistic regression models (packages: glmnet, MASS). The filtration threshold is determined experimentally.
5. Area Under Receiver Operating Characteristic (AUROC) calculations are computed in caret package. The trapezoidal rule is used to compute the area under the ROC curve for each predictor in binary classification. This area is used as the measure of variable importance and for filtration. The problem is decomposed into all pairwise tasks for multiclass outcomes and the area under the curve is calculated for each class pair. The mean value of AUROC across all predictors is directly calculated for binary classification task, for multiclass – at first, sum of AUROC for all groups pairs is determined and is divided by number of class labels and then mean value of AUROC across all predictors is measured. The filtration threshold is determined experimentally.
6. Univariate filtering (UVF) is a subsequent implementation of Shapiro-Wilk normality test, Bartlett test of homoscedasticity and Wilcox, Welch, Student tests with Benjamini-Hochberg method for multiple comparison (significance level was set 5% by default). The feature is filtered if significant level is reached between two groups in binary classification and at least any two groups – in multilabel task.
7. Kruskal-Wallis and t-test can be performed separately. The feature is filtered if significant level is reached.
8. Filtration by fold change value. The value is set manually. Fold change calculation for multiclass dataset can be also performed.
9. Filtration by moderated t-test [54] (limma package). The significant level is set manually.
10. Filtration by LM and LMM models (eq. 7,5; packages: MetabolomicsBasics, lme4, lmerTest). The features are filtered by p-value (set manually) for the target variable.
11. Filtration by RUV2 algorithm ([55], package NormalizeMets). The p-value cutoff for the target variable is set manually.
12. Filtration by two-dimensional false discovery rate control [56].
13. Filtration by removing highly correlated features (package caret). The absolute values of pair-wise correlations are considered. If two variables have a high correlation, the function looks at the mean absolute correlation of each variable and removes the variable with the largest mean absolute correlation. The cutoff for pairwise absolute correlation is set manually.
14. Filtration by correlation coefficient. Each variable is compared with vector of target numeric variable. Correlation cutoff is set manually.

Statistical filtration output can be adapted for any combination of filtration algorithms by selection of intersected variables or all unique features.

Another part of “Statistics” script is a classification model task. At first, dataset is divided into training and validation subsets. Then, resampling procedure is determined (bootstrap or cross-validation) and classification metric (ROC or accuracy) is specified (internal validation). In the next part ML model is fitted to training set via caret wrapper function. Hyperparameter tune length is set to 10 and can be changed. The resulted ML model is used for prediction on validation set (external validation) and resulted model overview with resampling statistics are provided. Also, logistic, penalized and stepwise logistic regression (packages: glmnet, MASS) are constructed in similar way. Performance of the models are tested by construction of confusion matrix on training and predicted class labels and ROC analysis. Generalized Weighted Quantile Sum (gWQS) Regression [57] is also available as option for classification task. Features can be visualized via box-plot or scatter-plot by groups.

Regression model task is constructed in the same manner as classification task. The differences are as follows: metric is RMSE, MAE or R2; performance is tested by MAE, RMSE, R2 values; optimal subset of variables selection, penalized regression and stepwise regression are performed via leaps, glmnet and MASS packages; visualization is conducted by scatter plot of predicted and training data points. gWQS is also available for regression task.

Next block is testing of biomarker set. Features can be visualized via scatter plot by data points and box or violin plots by one or two groups. MANOVA and PERMANOVA (packages: vegan, pairwiseAdonis) procedures are available. Means comparison can be performed as: fold change calculation, t-/moderated t-/Wilcoxon/ Kruskal-Wallis tests, two-way ANOVA and one-way ANOVA (equal to UVF). In all cases, p-values are adjusted for multiple comparisons.

Metabolomic-Wide association study (MWAS) and analysis of covariance (ANCOVA) are carried out in the form of LM and LMM modeling (eq. 7,5) or similarly by GAM, GAMM and nonlinear dose-response curves (DRC).

Signal modeling can be implemented by LM, LMM, GAM, GAMM and some other nonlinear functions for Dose-Response curve (DRC) analysis.

N-factor analysis is implemented by LM/LMM/GAM/GAMM/DRC modeling equally to MWAS/ANCOVA and Signal modeling case. ANOVA-simultaneous component analysis (ASCA, package MetStaT, [58]) is also accessible as an option for multiway analysis, two-way ANOVA, two-dimensional false discovery rate control, PVCA and PC-PR2.

Analysis of repeated measures is performed by LM/LMM/GAM/GAMM/DRC modeling (eq. 7,5) and multilevel sPLS algorithm (package mixOmics, [59]).

Time series or longitudinal analysis is introduced by LM/LMM/GAM/GAMM/DRC modeling (eq. 7,5), ASCA and multivariate empirical Bayes statistics (package timecourse, [60]). Also, dose-response modeling is available (DRomics [61], TOXcms [62]), profile modeling by LMM Splines (timeOmics [63]) and PVCA, PC-PR2. 

Multivariate data visualization and projection is provided by unsupervised data analysis. PCA, HCA, k-means clustering (all in packages: factoextra, FactoMineR, dendextend), HCA on PCA scores, heatmap (package pheatmap), t-distributed stochastic neighbor embedding (package Rtsne) and validation with optimization of clustering (Dunn index, silhouette analysis, gap stats, Rand index, p-value for HCA, classification accuracy, etc.; packages: NbClust, clustertend, mclust, clValid, fpc, pvclust) are available as an options for unsupervised data projection.

Correlation analysis is represented by computing correlation matrix and correlograms (packages: Hmisc, corrplot, psych).

Distance analysis part allows to calculate distance matrix for observations or features by distance metrics (Euclidean, maximum, Manhattan, Canberra, binary or Minkowski) or correlation and plot heatmap.

In the next section, effect size (Cohen's d and Hedges'g) and power analysis with sample size calculation are available (packages: effectsize, pwr).

## Conclusion
The distribution of the number of strings in script file could be used for visualization of the relative contribution of each of nine processing steps into the overall workflow (Fig. 1). 
