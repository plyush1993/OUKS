<!-- badges: starts -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://img.shields.io/badge/R≥4.1.2-5fb9ed.svg?style=flat-square&logo=r&logoColor=white?)](https://cran.r-project.org/index.html)
[![License](https://img.shields.io/badge/license-GPLv3-2186f8.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![](https://img.shields.io/github/v/release/plyush1993/OUKS?include_prereleases&color=ff6500)](https://github.com/plyush1993/OUKS/releases)
[![Report](https://img.shields.io/badge/see-Report-078a8a.svg?maxAge=2678400&style=flat-square)](https://github.com/plyush1993/OUKS/tree/main/Report%20(Rmd))
[![Changelog](https://img.shields.io/badge/keep-Changelog-d30b0b.svg?maxAge=2678400&style=flat-square)](https://github.com/plyush1993/OUKS/blob/main/CHANGELOG.md)
[![](https://img.shields.io/badge/dependencies-R%20session-yellow.svg)](https://github.com/plyush1993/OUKS/blob/main/Session%20Info.txt)
![](https://img.shields.io/github/repo-size/plyush1993/OUKS)
![](https://img.shields.io/github/languages/code-size/plyush1993/OUKS)
[![](https://img.shields.io/badge/article-JPR.1c00392-blueviolet.svg)](https://doi.org/10.1021/acs.jproteome.1c00392)
[![](https://img.shields.io/badge/web-site-ff69b4.svg)](https://plyush1993.github.io/OUKS/)
![](https://img.shields.io/github/release-date/plyush1993/OUKS)
<!-- badges: end -->

# Omics Untargeted Key Script *(OUKS)* <img src="GH logo.gif" align="right" height="300" width="280"> 
## Brief Description :old_key:
R based open-source collection of scripts called :red_circle:*OUKS*:large_blue_circle: (*Omics Untargeted Key Script*) providing comprehensive nine step LC-MS untargeted metabolomic profiling data processing toolbox :toolbox:  

### **See [website](https://plyush1993.github.io/OUKS/)**
---

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

## Table of Contents :clipboard:
- Instruction and introduction into the :red_circle:*OUKS* toolbox:large_blue_circle: are provided by [`Basic tutorial`](https://github.com/plyush1993/OUKS/blob/main/Basic%20tutorial.pdf) file. [`Session info`](https://github.com/plyush1993/OUKS/blob/main/Session%20Info.txt) and [`Used packages`](https://github.com/plyush1993/OUKS/blob/main/Used%20packages.pdf) are listed in corresponding files.
- [Scripts](https://github.com/plyush1993/OUKS/tree/main/Scripts%20(R)) with comments, notes and references are stored in Scripts folder at a previously defined order along with [code](https://github.com/plyush1993/OUKS/blob/main/Scripts%20(R)/Figures%20for%20OUKS.R) for plotting figures associated with [article](https://doi.org/10.1021/acs.jproteome.1c00392).
- [MS2 spectra](https://github.com/plyush1993/OUKS/tree/main/MS2%20spectra%20(mzXML)) for selected potential biomarkers of bladder cancer are stored in mzXML format at corresponding folder.
- [Datasets](https://github.com/plyush1993/OUKS/tree/main/Datasets%20(csv)) in .csv and [other files](https://github.com/plyush1993/OUKS/tree/main/Auxiliary%20files%20(RData)) (.RData, .R) are available for reproducibility from corresponding folders. Files descriptions are provived by [`Roadmap`](https://github.com/plyush1993/OUKS/blob/main/Roadmap.pdf) file. Raw data (.CDF format) are available from [Metabolomics Workbench Repository](https://www.metabolomicsworkbench.org/), study ID: [ST001682](http://doi.org/10.21228/M8ZT4C). [`Metadata`](https://github.com/plyush1993/OUKS/blob/main/Datasets%20(csv)/metadata.csv) table is also provided.
- [Report](https://github.com/plyush1993/OUKS/tree/main/Report%20(Rmd)) in [`.Rmd`](https://github.com/plyush1993/OUKS/blob/main/Report%20(Rmd)/Report%20example%20OUKS.Rmd), [`.pdf`](https://github.com/plyush1993/OUKS/blob/main/Report%20(Rmd)/Report-example-OUKS.pdf) and [`.docx`](https://github.com/plyush1993/OUKS/blob/main/Report%20(Rmd)/Report-example-OUKS.docx) formats were provided as an example to reproduce the *OUKS* code script.
- [`Changelog`](https://github.com/plyush1993/OUKS/blob/main/CHANGELOG.md) file is provided and is constantly updated. See also [releases](https://github.com/plyush1993/OUKS/releases) page.
- [Required packages](https://github.com/plyush1993/OUKS/tree/main/Required%20packages%20(archive)) is for storing packages archives with strong version dependency.

## Requirements :building_construction:
The only requirements are to be familiar with the basic syntax of the R language, PC with Internet connection and Windows OS (desirable), [RStudio](https://www.rstudio.com/products/rstudio/download/) and [R](https://cloud.r-project.org/) (≥ 4.1.2).

## Citation :link:
*OUKS* has been published in the [Journal of Proteome Research](https://pubs.acs.org/journal/jprobs). If you use this software to analyze your own data, please cite it as below, thanks:

<img src="qrcode1.png" align="right" height="100" width="200"/> 

> [Ivan V. Plyushchenko, Elizaveta S. Fedorova, Natalia V. Potoldykova, Konstantin A. Polyakovskiy, Alexander I. Glukhov, and Igor A. Rodin
> Journal of Proteome Research 2022 21 (3), 833-847. DOI: 10.1021/acs.jproteome.1c00392](https://doi.org/10.1021/acs.jproteome.1c00392)

*OUKS* builds on many open-source software tools and open data sources. Therefore, it is important to also cite their
work when using these algorithms via *OUKS*: [*1*](https://github.com/plyush1993/OUKS/wiki/Tutorial#references), [*2*](https://github.com/plyush1993/OUKS/blob/main/Session%20Info.txt).

## Contact :memo:
Please send any comment, suggestion or question you may have to the author (:man_scientist: Dr. Ivan Plyushchenko): 

[:octocat:](https://github.com/plyush1993)
[<img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png">](https://orcid.org/0000-0003-3883-4695)
<a href="https://scholar.google.com/citations?user=Mz4nxtwAAAAJ&hl=ru&oi=ao"><img alt="Google Scholar" src="https://img.shields.io/badge/Scholar%20-%23F6F6F6.svg?&style=flat-square&logoColor=white&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCA1MTIgNTEyIj48cGF0aCBmaWxsPSIjNDI4NWY0IiBkPSJNMjU2IDQxMS4xMkwwIDIwMi42NjcgMjU2IDB6Ii8+PHBhdGggZmlsbD0iIzM1NmFjMyIgZD0iTTI1NiA0MTEuMTJsMjU2LTIwOC40NTNMMjU2IDB6Ii8+PGNpcmNsZSBmaWxsPSIjYTBjM2ZmIiBjeD0iMjU2IiBjeT0iMzYyLjY2NyIgcj0iMTQ5LjMzMyIvPjxwYXRoIGZpbGw9IiM3NmE3ZmEiIGQ9Ik0xMjEuMDM3IDI5OC42NjdjMjMuOTY4LTUwLjQ1MyA3NS4zOTItODUuMzM0IDEzNC45NjMtODUuMzM0czExMC45OTUgMzQuODgxIDEzNC45NjMgODUuMzM0SDEyMS4wMzd6Ii8+PC9zdmc+"></a>
[<img src = https://upload.wikimedia.org/wikipedia/commons/5/5e/ResearchGate_icon_SVG.svg height="15" width="15">](https://www.researchgate.net/profile/Ivan-Plyushchenko-2)
[<img src = https://upload.wikimedia.org/wikipedia/commons/thumb/2/26/Scopus_logo.svg/2560px-Scopus_logo.svg.png height="16.4" width="50">](https://www.scopus.com/authid/detail.uri?authorId=57202386632)
[<img src = https://res.cloudinary.com/crunchbase-production/image/upload/c_lpad,h_256,w_256,f_auto,q_auto:eco,dpr_1/v1496403635/dbpokhjtlxno6psfu5r0.png height="25" width="25">](https://publons.com/researcher/4095888/ivan-plyushchenko/)
