<!-- badges: starts -->
[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)
![](https://img.shields.io/github/release/plyush1993/OUKS?color=teal&logoColor=teal)
[![](https://img.shields.io/badge/R≥4.1.2-5fb9ed.svg?style=flat&logo=r&logoColor=white?)](https://cran.r-project.org/index.html)
[![License](https://img.shields.io/badge/GPLv3-indianred.svg?style=flat&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![](https://img.shields.io/badge/📜-F8F8F8.svg)](https://doi.org/10.1021/acs.jproteome.1c00392)
[![](https://img.shields.io/badge/🌐-F8F8F8.svg)](https://plyush1993.github.io/OUKS/)
<!-- badges: end -->

# Omics Untargeted Key Script *(OUKS)* <img src="GH logo.gif" align="right" height="220" width="200"> 
### Brief Description :old_key:
R based open-source collection of scripts called :red_circle:*OUKS*:large_blue_circle: (*Omics Untargeted Key Script*) providing comprehensive nine step LC-MS untargeted metabolomic profiling data processing toolbox :toolbox:  

See [website](https://plyush1993.github.io/OUKS/) :computer:

---
### Table of contents :clipboard:

<details>
  <summary><b>Click to expand&nbsp;&nbsp;</b></summary>
  <br/>
  
- Instruction and introduction into the :red_circle:*OUKS* toolbox:large_blue_circle: are provided by [`Basic tutorial`](https://github.com/plyush1993/OUKS/blob/main/Basic%20tutorial.pdf) file. [`Session info`](https://github.com/plyush1993/OUKS/blob/main/Session%20Info.txt) and [`Used packages`](https://github.com/plyush1993/OUKS/blob/main/Used%20packages.pdf) are listed in corresponding files.
- [Scripts](https://github.com/plyush1993/OUKS/tree/main/Scripts%20(R)) with comments, notes and references are stored in Scripts folder at a previously defined order along with [code](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/Figures%20for%20OUKS.R) for plotting figures associated with [article](https://doi.org/10.1021/acs.jproteome.1c00392).
- [MS2 spectra](https://github.com/plyush1993/OUKS/tree/main/MS2%20spectra%20(mzXML)) for selected potential biomarkers of bladder cancer are stored in mzXML format at corresponding folder.
- [Datasets](https://github.com/plyush1993/OUKS/tree/main/Datasets%20(csv)) in .csv and [other files](https://github.com/plyush1993/OUKS/tree/main/Auxiliary%20files%20(RData)) (.RData, .R) are available for reproducibility from corresponding folders. Files descriptions are provived by [`Roadmap`](https://github.com/plyush1993/OUKS/blob/main/Roadmap.pdf) file. Raw data (.CDF format) are available from [Metabolomics Workbench Repository](https://www.metabolomicsworkbench.org/), study ID: [ST001682](http://doi.org/10.21228/M8ZT4C). [`Metadata`](https://github.com/plyush1993/OUKS/blob/main/Datasets%20(csv)/metadata.csv) table is also provided.
- [Report](https://github.com/plyush1993/OUKS/tree/main/Report%20(Rmd)) in [`.Rmd`](https://github.com/plyush1993/OUKS/blob/main/Report%20(Rmd)/Report%20example%20OUKS.Rmd), [`.pdf`](https://github.com/plyush1993/OUKS/blob/main/Report%20(Rmd)/Report-example-OUKS.pdf) and [`.docx`](https://github.com/plyush1993/OUKS/blob/main/Report%20(Rmd)/Report-example-OUKS.docx) formats were provided as an example to reproduce the *OUKS* code script.
- [`Changelog`](https://github.com/plyush1993/OUKS/blob/main/CHANGELOG.md) file is provided and is constantly updated. See also [releases](https://github.com/plyush1993/OUKS/releases) page.
- [Required packages](https://github.com/plyush1993/OUKS/tree/main/Required%20packages%20(archive)) section is for storing packages archives with strong version dependency.
- Discussions, suggestions and error reports are [welcome](https://github.com/plyush1993/OUKS/issues).

</details>

Script | Purpose
------------ | -------------
[`1. Randomization.R`](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/1.%20Randomization.R) | experimental design and sample randomization
[`2. Integration.R`](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/2.%20Integration.R) | peaks integration and time alignment
[`3. Imputation.R`](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/3.%20Imputation.R) | missing value imputation (MVI) and artifacts removal
[`4. Correction.R`](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/4.%20Correction.R) | signal drift correction and batch effect removal
[`5. Annotation.R`](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/5.%20Annotation.R) | feature annotation and tentative identification by database search
[`6. Filtering.R`](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/6.%20Filtering.R) | peaks filtering for quality checking and accounting of technical variation
[`7. Normalization.R`](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/7.%20Normalization.R) | data normalization and adjusting of biological variation
[`8. Grouping.R`](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/8.%20Grouping.R) | peaks grouping and molecular features clustering
[`9. Statistics.R`](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/9.%20Statistics.R) | statistical analysis and hypothesis testing

### Requirements :building_construction:
The only requirements are to be familiar with the basic syntax of the R language, PC with Internet connection and Windows OS (desirable), [RStudio](https://www.rstudio.com/products/rstudio/download/) and [R](https://cloud.r-project.org/) (≥ 4.1.2).

### Citation :link:
*OUKS* has been published in the [Journal of Proteome Research](https://pubs.acs.org/journal/jprobs). If you use this software to analyze your own data, please cite it as below, thanks:

> [Ivan V. Plyushchenko, Elizaveta S. Fedorova, Natalia V. Potoldykova, Konstantin A. Polyakovskiy, Alexander I. Glukhov, and Igor A. Rodin
> Journal of Proteome Research 2022 21 (3), 833-847. DOI: 10.1021/acs.jproteome.1c00392](https://doi.org/10.1021/acs.jproteome.1c00392)

*OUKS* builds on many open-source software tools and open data sources. Therefore, it is important to also cite their
work when using these algorithms via *OUKS*: [*1*](https://github.com/plyush1993/OUKS/blob/main/Session%20Info.txt), [*2*](https://github.com/plyush1993/OUKS/blob/gh-pages/Tutorial.md#references).

### Contact :memo:
Please send any comment, suggestion or question you may have to the author (👨‍🔬 Dr. Ivan Plyushchenko): 
<div> 
  <a href="mailto:plyushchenko.ivan-@gmail.com"><img src="https://img.shields.io/badge/-4a9edc?style=for-the-badge&logo=gmail" height="28" alt="Email" /></a>
  <a href="https://github.com/plyush1993"><img src="https://img.shields.io/static/v1?style=for-the-badge&message= &color=181717&logo=GitHub&logoColor=FFFFFF&label=" height="28" alt="GH" /></a>
  <a href="https://orcid.org/0000-0003-3883-4695"><img src="https://img.shields.io/badge/-A6CE39?style=for-the-badge&logo=ORCID&logoColor=white" height="28" alt="ORCID" /></a>
</div>
