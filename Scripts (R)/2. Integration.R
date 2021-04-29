##############################################################################################################################################################
# Table of contents
##############################################################################################################################################################

# Installation
# IPO for XCMS params optimization
# XCMS with the best params
# Warpgroup for increase precision
# Description of peaks table
# References

##############################################################################################################################################################
# Installation
##############################################################################################################################################################

# setup environment
setwd("D:/...")

wd_1 <- c("D:/...") # folder with files for IPO optimization process
wd_2 <- c("D:/...") # folder with all study samples

files_QC <- list.files(wd_1, recursive = TRUE, full.names = TRUE, pattern = ".CDF") # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.
files_all <- list.files(wd_2, recursive = TRUE, full.names = TRUE, pattern = ".CDF") # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.

##############################################################################################################################################################
# IPO for XCMS params optimization
##############################################################################################################################################################

library(IPO)
library(xcms)
library(doParallel)

nCore=detectCores()-1

#PeakPickingParameters
peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')
peakpickingParameters$noise=c(100,1000)
peakpickingParameters$value_of_prefilter=c(3,800)
peakpickingParameters$min_peakwidth<- c(3,15)
peakpickingParameters$max_peakwidth<- c(30,40)
peakpickingParameters$ppm<- c(15,35)
param=SnowParam(workers = nCore)
nSlaves=1 # or nCore or more

resultPeakpicking <-
  optimizeXcmsSet(files = files_QC,
                  params = peakpickingParameters,
                  BPPARAM = param,
                  nSlaves = nSlaves,
                  subdir = NULL,
                  plot = TRUE)

optimizedXcmsSetObject <- resultPeakpicking$best_settings$xset
# save(resultPeakpicking, file = "IPO_optimiz_xcms_CDF_7QC.RData")

# Retention Time Alignment Optimization
retcorGroupParameters <- getDefaultRetGroupStartingParams()
retcorGroupParameters$profStep <- c(0.33,1)
retcorGroupParameters$gapExtend <- c(2.0,3.0)
retcorGroupParameters$minfrac=c(0.05,0.5)
retcorGroupParameters$response=c(9.0,18.0)
retcorGroupParameters$gapInit=c(0.10, 0.50)
retcorGroupParameters$mzwid=c(0.020, 0.050)
BiocParallel::register(BiocParallel::SerialParam())
nSlaves=1 # or nCore or more

resultRetcorGroup <-
  optimizeRetGroup(xset = optimizedXcmsSetObject,
                   params = retcorGroupParameters,
                   nSlaves = nSlaves,
                   subdir = NULL,
                   plot = TRUE)

resultRetcorGroup$best_settings
# save(resultRetcorGroup, file = "IPO_optimiz_ret_align_CDF_7QC.RData")

# Get final results
writeRScript(resultPeakpicking$best_settings$parameters, resultRetcorGroup$best_settings)
param <- c(resultPeakpicking$best_settings$parameters, resultRetcorGroup$best_settings)
#save(param,file = 'all params IPO CDF 7QC.RData')

funs_params <- capture.output(writeRScript(resultPeakpicking$best_settings$parameters, resultRetcorGroup$best_settings), type = "message")
#save(funs_params,file = 'funs params IPO CDF 7QC.RData')

##############################################################################################################################################################
# XCMS with the best params
##############################################################################################################################################################

library(xcms)
library(data.table)
library(dplyr)
library(stringr)

# load best parameters from IPO
# load("IPO_optimiz_xcms_CDF_7QC.RData")
# load("IPO_optimiz_ret_align_CDF_7QC.RData")

# Create a phenodata data.frame
pd <- data.frame(sample_name = sub(basename(files_all), pattern = ".CDF",
                                   replacement = "", fixed = TRUE), stringsAsFactors = FALSE) # download filenames

rname <- pd
all_id <- sapply(1:nrow(rname), function(y) unlist(str_split(rname[y,], " "))) # split info from rownames
all_id1 <- lapply(1:length(all_id), function(y) as.numeric(all_id[[y]][1])) # as numeric run order in [[x]][1]
all_id2 <- lapply(1:length(all_id), function(y) replace(all_id[[y]],1, all_id1[[y]][1]))

p1_id <- unlist(lapply(all_id, function(y) unlist(y[4]))) # obtain patient ID (every [4] element)
un_raw_l <- unique(gsub("[[:digit:]]", "", p1_id)) # obtain unique raw label from p1_id (every [4] element)
true_l <- c("QC", "TG", "CG") # write desired label in order as in unique raw label (should be analogical)
rbind(un_raw_l,true_l) # visual check agreement between true and raw labels
raw_l <- gsub("[[:digit:]]", "", p1_id) # obtain raw label from p1_id (every [4] element)
n_gr_t <- as.character(match(raw_l, un_raw_l)) # obtain index for raw label by unique raw label
for (i in 1:length(unique(n_gr_t))) {
  n_gr_t <- str_replace(n_gr_t, unique(n_gr_t)[i], true_l[i]) } # exchange raw label by true label via index
n_gr_t <- data.frame(n_gr_t)

vec_gr <- as.numeric(as.factor(n_gr_t$n_gr_t))
sample_gr <- unique(as.numeric(as.factor(n_gr_t$n_gr_t)))
n_gr <- sapply(1:length(sample_gr), function(y) length(vec_gr[vec_gr == y]))
min_frac_man <- min(round(n_gr/length(vec_gr), 1)) # calculate min_frac manually

# download files
raw_data <- readMSData(files = files_all, pdata = new("NAnnotatedDataFrame", n_gr_t), mode = "onDisk") # or use pd only as: pdata = new("NAnnotatedDataFrame", pd) or pdata = new("NAnnotatedDataFrame", n_gr_t)
raw_data <- filterRt(raw_data, c(0,2500)) # time range in sec for unified rt range

# parallel processing
cores = detectCores()-1
register(bpstart(SnowParam(cores)))
BiocParallel::register(BiocParallel::SerialParam())

# feature detection
cwp <- CentWaveParam(ppm = resultPeakpicking$best_settings$parameters$ppm, 
                     peakwidth = c(resultPeakpicking$best_settings$parameters$min_peakwidth, resultPeakpicking$best_settings$parameters$max_peakwidth),
                     snthresh = resultPeakpicking$best_settings$parameters$snthresh,
                     prefilter = c(resultPeakpicking$best_settings$parameters$prefilter, resultPeakpicking$best_settings$parameters$value_of_prefilter),
                     mzCenterFun = resultPeakpicking$best_settings$parameters$mzCenterFun,
                     integrate = resultPeakpicking$best_settings$parameters$integrate,
                     mzdiff = resultPeakpicking$best_settings$parameters$mzdiff,
                     fitgauss = resultPeakpicking$best_settings$parameters$fitgauss,
                     noise = resultPeakpicking$best_settings$parameters$noise)

feat_det <- findChromPeaks(raw_data, param = cwp)

# retention time correction
BiocParallel::register(BiocParallel::SerialParam())
ret_cor <- adjustRtime(feat_det, param = ObiwarpParam(
                          binSize = resultRetcorGroup$best_settings$profStep,
                          center = resultRetcorGroup$best_settings$center,
                          response = resultRetcorGroup$best_settings$response,
                          distFun = resultRetcorGroup$best_settings$distFunc,
                          gapInit = resultRetcorGroup$best_settings$gapInit,
                          gapExtend = resultRetcorGroup$best_settings$gapExtend,
                          factorDiag = resultRetcorGroup$best_settings$factorDiag,
                          factorGap = resultRetcorGroup$best_settings$factorGap,
                          localAlignment = ifelse(resultRetcorGroup$best_settings$localAlignment==0, F,T)))

# peak grouping
pgp <- PeakDensityParam(sampleGroups = as.numeric(as.factor(n_gr_t$n_gr_t)), # rep(1, length(fileNames(feat_det))) or as.numeric(as.factor(n_gr_t))
                        bw = resultRetcorGroup$best_settings$bw,
                        minFraction =  min_frac_man, # or resultRetcorGroup$best_settings$minfrac
                        minSamples = resultRetcorGroup$best_settings$minsamp, 
                        binSize = resultRetcorGroup$best_settings$mzwid,
                        maxFeatures = resultRetcorGroup$best_settings$max) 

pk_gr <- groupChromPeaks(ret_cor, param = pgp)

# peak filling
pk_fil <- fillChromPeaks(pk_gr)

# final feature table
ft_tbl <- featureValues(pk_fil, value = "into")

# final peak info table
ft_inf <- featureDefinitions(pk_fil)

# join peak table and save
ft_tbl_f <- data.frame(t(ft_tbl))
colnames(ft_tbl_f) <- paste(ft_inf$mzmed, ft_inf$rtmed, sep = " / ")
fwrite(ft_tbl_f, "xcms after IPO.csv", row.names = T)

# save all xcms objects
save(feat_det, file = "xcms obj feat_det.RData")
save(ret_cor, file = "xcms obj ret_cor.RData")
save(pk_gr, file = "xcms obj pk_gr.RData")
save(pk_fil, file = "xcms obj pk_fil.RData")

##############################################################################################################################################################
# Warpgroup for increase precision
##############################################################################################################################################################

library(warpgroup)
library(xcms)
library(data.table)
library(doParallel)

# load peak table
load("xcms obj pk_gr.RData") # filename and the same folder of raw data as in section 2
pk_gr <- as(pk_gr, "xcmsSet") # convert XCMSnExp to xcmsSet object

# load optimal params from IPO
load("IPO_optimiz_xcms_CDF_7QC.RData")
ppm_diff <- resultPeakpicking[["best_settings"]][["parameters"]][["ppm"]]

# Parallel Backend Setup
cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Check time of processing
detach("package:dplyr", unload = T)
xs2 <- pk_gr
xs2@groupidx <- sample(pk_gr@groupidx, 5)
xr.l = llply(pk_gr@filepaths, xcmsRaw, profstep=0) # adjust to your data
start = Sys.time()
wg <- group.warpgroup(xs2, xr.l = xr.l, smooth.n = 1, rt.max.drift = 15, ppm.max.drift = ppm_diff, rt.aligned.lim = 5) # adjust to your data
end = Sys.time()
paste("Full dataset will take", length(pk_gr@groupidx)/length(xs2@groupidx), "times as long to process.")
paste("Subset took", end - start)

#Perform Warpgroup
detach("package:dplyr", unload = T)
xr.l = llply(pk_gr@filepaths, xcmsRaw, profstep=0) # adjust to your data
wg = group.warpgroup(pk_gr, xr.l = xr.l, smooth.n = 1, rt.max.drift = 15, ppm.max.drift = ppm_diff, rt.aligned.lim = 5) # adjust to your data

# Save
save(wg, file = "warpgroup res.RData")
data <- xcms::groupval(wg, "medret", "into")
mz <- as.numeric(xcms::groups(wg)[, 1])
time <- as.numeric(xcms::groups(wg)[, 4])
rownames(data) <- paste(mz, time, sep = " / ")
load("xcms obj pk_fil.RData")
ft_tbl <- featureValues(pk_fil, value = "into")
colnames(data) <- colnames(ft_tbl)
data <- as.data.frame(t(data))
fwrite(data, "xcms after WG test.csv", row.names = T)

##############################################################################################################################################################
# Description of peaks table
##############################################################################################################################################################

# dataset
ds_pt <- ft_tbl_f

# Missing value %
tn <- nrow(ds_pt)*ncol(ds_pt)
mv_c <- sum(is.na(ds_pt)) # or length(which(is.a(...)))
pr_mv <- round(mv_c/tn*100,0)

# Number of peaks
ncol(ds_pt) # or nrow

##############################################################################################################################################################
# References
##############################################################################################################################################################

# 1. Tautenhahn, Ralf, Christoph Boettcher, and Steffen Neumann. "Highly sensitive feature detection for high resolution LC/MS." BMC bioinformatics 9.1 (2008): 504.
# 2. Libiseller, Gunnar, et al. "IPO: a tool for automated optimization of XCMS parameters." BMC bioinformatics 16.1 (2015): 118.
# 3. Fernández-Ochoa, Álvaro, et al. "A Case Report of Switching from Specific Vendor-Based to R-Based Pipelines for Untargeted LC-MS Metabolomics." Metabolites 10.1 (2020): 28.
# 4. Mahieu, Nathaniel G., Jonathan L. Spalding, and Gary J. Patti. "Warpgroup: increased precision of metabolomic data processing by consensus integration bound analysis." Bioinformatics 32.2 (2016): 268-275.