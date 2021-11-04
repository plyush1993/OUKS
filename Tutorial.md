# Omics Untargeted Key Scripts

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

