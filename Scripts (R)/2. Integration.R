##############################################################################################################################################################
# Table of contents
##############################################################################################################################################################

# Installation
# Autotuner for XCMS params optimization
# MetaboAnalystR for XCMS params optimization
# Paramounter for XCMS params optimization
# IPO for XCMS params optimization
# XCMS with the best params (from IPO) and other options of processing
# Old version of XCMS
# XCMS for improve integration and alignment
# Warpgroup for increase precision
# ncGTW for realignment
# cpc for filtering peaks
# Description of peaks table
# EIC peaks integration
# References

##############################################################################################################################################################
# Installation
##############################################################################################################################################################

library(data.table)
library(dplyr)
library(stringr)

# setup environment
setwd("D:/...")

wd_1 <- c("D:/...") # folder with files for IPO optimization process
wd_2 <- c("D:/...") # folder with all study samples

files_QC <- list.files(wd_1, recursive = TRUE, full.names = TRUE, pattern = ".CDF") # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.
files_all <- list.files(wd_2, recursive = TRUE, full.names = TRUE, pattern = ".CDF") # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.

# order files by run order
# load in format: 1st column: 1. qc1(s78_1) b1 QC1(KMS 53).X --- ro. nameN(nameN_NREPEAT) BatchIndex QCN(Label N)
# "150. qc30 b6 QC30.CDF" 
# "167. s64_2 b7 MS 45.CDF"
rname <- files_all # obtain all info from rownames
rname <- str_remove(rname, wd_2) # remove folder name
rname <- str_remove(rname, ".CDF") # remove some pattern from vendor-specific format
all_id <- lapply(1:length(rname), function(y) unlist(str_split(rname[y], " ")), simplify = F) # split info from rownames
ro_id <- as.numeric(unlist(lapply(all_id, function(y) unlist(y[1])))) # obtain run order ID (every [1] element)
files_all_df <- as.data.frame(cbind(files_all, ro = as.numeric(ro_id)))
files_all_df <- files_all_df[order(as.numeric(files_all_df$ro), decreasing = F),]
files_all <- files_all_df[,-2]

##############################################################################################################################################################
# Autotuner for XCMS params optimization
##############################################################################################################################################################

library(Autotuner)

# Raw data
rawPaths <- "D:/..." # folder for raw data files
files <- list.files(rawPaths, recursive = TRUE, full.names = TRUE, pattern = ".mzData") # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.

# Generate metadata
col_f <- paste("QC", 1:length(rawPaths), sep = " ") # generate "Factor.Value.genotype.", adjust to your data
metadata <- data.frame(cbind(files, col_f))
colnames(metadata) <- c("Raw.Spectral.Data.File", "Factor.Value.genotype.")

# Perform Parameters Optimization
Autotuner <- createAutotuner(files,
                             metadata,
                             file_col = "Raw.Spectral.Data.File",
                             factorCol = "Factor.Value.genotype.")

lag <- 100 # adjust to your data
threshold <- 5 # adjust to your data
influence <- 0.01 # adjust to your data
signals <- lapply(getAutoIntensity(Autotuner), 
                  ThresholdingAlgo, lag, threshold, influence)

plot_signals(Autotuner, 
             threshold, 
             ## index for which data files should be displayed
             sample_index = 1:3, # number of sample
             signals = signals)

Autotuner <- isolatePeaks(Autotuner = Autotuner, 
                          returned_peaks = 30, 
                          signals = signals)


for(i in 1:10) {
  plot_peaks(Autotuner = Autotuner, 
             boundary = 100, 
             peak = i)  }

# get parameters
eicParamEsts <- EICparams(Autotuner = Autotuner, 
                          massThresh = .005, 
                          verbose = T,
                          returnPpmPlots = FALSE,
                          useGap = T,
                          varExpThresh = 1)

returnParams(eicParamEsts, Autotuner)

best_param <- returnParams(eicParamEsts, Autotuner)
# save(best_param, file = "Autotuner_optimiz_xcms.RData")
# load("Autotuner_optimiz_xcms.RData")

##############################################################################################################################################################
# MetaboAnalystR for XCMS params optimization
##############################################################################################################################################################

# devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = F, build_manual =F)

library(MetaboAnalystR)

# Raw data
rawPaths <- "D:/..." # folder for raw data files

# Inspect the MS data via a 3D image. "res" are used to specify the resolution for the MS data.
PerformDataInspect(rawPaths,res = 1000, rt.range = c(1,800), mz.range = c(80,800), dimension = "3D") # if some error: specify one file and mzXML

# Data trimming
#raw_data <- PerformDataTrimming(rawPaths,rt.idx = 1) # The percentage of RT dimension retained is set as 20%
raw_data <- PerformDataTrimming(rawPaths, rt.idx = 0.2) # Time select, adjust to your data

# Initial platform specific parameters
param_initial <- SetPeakParam(platform = "general") # adjust to your data

# Perform Parameters Optimization
param_optimized <- PerformParamsOptimization(raw_data, param = param_initial, ncore = 1)

param_optimized$best_parameters # see parameters
# save(param_optimized, file = "MetaboAnalystR_optimiz_xcms.RData")
# load("MetaboAnalystR_optimiz_xcms.RData")

##############################################################################################################################################################
# Paramounter for XCMS params optimization
##############################################################################################################################################################

library(xcms)
library(MSnbase)
library(dplyr)
library(ggplot2)
library(gridExtra)

directory <- "D:/..." # folder with files for IPO optimization process

# set parameters
massSDrange <- 2 # adjust to your data
smooth <- 0 # adjust to your data
cutoff <- 0.95 # adjust to your data

setwd(directory)
filename <- list.files(pattern = ".CDF") # adjust to your data
start_time <- Sys.time()
mzDiff <- c()
ppm <- c()
mzDiff2D <- as.data.frame(matrix(ncol = 3, nrow = 1)) # adjust to your data
colnames(mzDiff2D) <- c("mz", "rt", "mzdiff") # adjust to your data
ppm2D <- as.data.frame(matrix(ncol = 3, nrow = 1)) # adjust to your data
colnames(ppm2D) <- c("mz", "rt", "ppm") # adjust to your data

# Perform Paramounter part1
source(paste0("D:/.../","\\Paramounter_part1 (V2).R"))

Software <- "XCMS" # adjust to your data
massSDrange <- 2 # adjust to your data
ppmCut <- 40 # select on Figure, adjust to your data by figure from part1
smooth <- 0 # adjust to your data

# Perform Paramounter part2
source(paste0("D:/.../","\\Paramounter_part2 (V2).R"))

# see parameters in Folders directory/"Parameters for XCMS" and/or directory/"Universal Parameters"

##############################################################################################################################################################
# IPO for XCMS params optimization
##############################################################################################################################################################

library(IPO)
library(xcms)
library(doParallel)
library(BiocParallel)

nCore=detectCores()-1

# PeakPickingParameters (adjust to your data)
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

# Retention Time Alignment Optimization (adjust to your data)
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
# XCMS with the best params (from IPO) and other options of processing
##############################################################################################################################################################

library(xcms)
library(data.table)
library(dplyr)
library(stringr)
library(BiocParallel)
library(doParallel)

# load best parameters from IPO:
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
raw_data <- MSnbase::filterRt(raw_data, c(0,2500)) # time range in sec for unified rt range, use also MSnbase::filterMz, MSnbase::filterPolarity, MSnbase::filterMsLevel
# use also MSnbase::filterMz, MSnbase::filterPolarity, MSnbase::filterMsLevel
               
# data smoothing   
# raw_data <- MSnbase::smooth(raw_data, method = "SavitzkyGolay", halfWindowSize = 4L) # see ?smooth for edit parameters            

# refinement of the centroided m/z values  
# raw_data <- MSnbase::pickPeaks(raw_data, refineMz = "descendPeak") # see ?pickPeaks for edit parameters               
               
# parallel processing
cores = detectCores()-1
register(bpstart(SnowParam(cores)))
BiocParallel::register(BiocParallel::SerialParam())

# feature detection
cwp <- xcms::CentWaveParam(ppm = resultPeakpicking$best_settings$parameters$ppm, # maximal tolerated m/z deviation in consecutive scans, in ppm
                     peakwidth = c(resultPeakpicking$best_settings$parameters$min_peakwidth, resultPeakpicking$best_settings$parameters$max_peakwidth), # min/max chromatographic peak width
                     snthresh = resultPeakpicking$best_settings$parameters$snthresh, # Signal/Noise threshold
                     prefilter = c(resultPeakpicking$best_settings$parameters$prefilter, resultPeakpicking$best_settings$parameters$value_of_prefilter), # Prefilter step for the first phase. Mass traces are only retained if they contain at least [prefilter peaks] peaks with intensity >= [prefilter intensity]
                     mzCenterFun = resultPeakpicking$best_settings$parameters$mzCenterFun,
                     integrate = resultPeakpicking$best_settings$parameters$integrate, # Integration method
                     mzdiff = resultPeakpicking$best_settings$parameters$mzdiff, # minimum difference in m/z for peaks with overlapping retention times, can be negative to allow overlap
                     fitgauss = resultPeakpicking$best_settings$parameters$fitgauss,
                     noise = resultPeakpicking$best_settings$parameters$noise) # optional argument which is useful for data that was centroided without any intensity threshold, centroids with intensity < noise are omitted from ROI detection

feat_det <- xcms::findChromPeaks(raw_data, param = cwp)

# retention time correction
BiocParallel::register(BiocParallel::SerialParam())
app <- xcms::ObiwarpParam(binSize = resultRetcorGroup$best_settings$profStep, # step size (in m/z) to use for profile generation from the raw data files
                    center = resultRetcorGroup$best_settings$center, # or use 1st sample , Sample as center sample
                    response = resultRetcorGroup$best_settings$response,
                    distFun = resultRetcorGroup$best_settings$distFunc,
                    gapInit = resultRetcorGroup$best_settings$gapInit,
                    gapExtend = resultRetcorGroup$best_settings$gapExtend,
                    factorDiag = resultRetcorGroup$best_settings$factorDiag,
                    factorGap = resultRetcorGroup$best_settings$factorGap,
                    localAlignment = ifelse(resultRetcorGroup$best_settings$localAlignment==0, F,T))

ret_cor <- xcms::adjustRtime(feat_det, param = app)

# peak grouping
pgp <- xcms::PeakDensityParam(sampleGroups = as.numeric(as.factor(n_gr_t$n_gr_t)), # rep(1, length(fileNames(feat_det))) or as.numeric(as.factor(n_gr_t))
                        bw = resultRetcorGroup$best_settings$bw, # Allowable retention time deviations, in seconds. In more detail: bandwidth (standard deviation or half width at half maximum) of gaussian smoothing kernel to apply to the peak density chromatogram
                        minFraction =  min_frac_man, # or resultRetcorGroup$best_settings$minfrac , minimum fraction of samples necessary in at least one of the sample groups for it to be a valid group
                        minSamples = resultRetcorGroup$best_settings$minsamp, # minimum number of samples necessary in at least one of the sample groups for it to be a valid group
                        binSize = resultRetcorGroup$best_settings$mzwid, # width of overlapping m/z slices to use for creating peak density chromatograms and grouping peaks across samples
                        maxFeatures = resultRetcorGroup$best_settings$max) # maximum number of groups to identify in a single m/z slice

pk_gr <- xcms::groupChromPeaks(ret_cor, param = pgp)

# peak filling
pk_fil <- xcms::fillChromPeaks(pk_gr) # see also params: "expandMz", "expandRt", "ppm"

# final feature table
ft_tbl <- featureValues(pk_fil, value = "into")  # "index" for height "into" for area

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
# Old version of XCMS
##############################################################################################################################################################

library(xcms)
library(BiocParallel)

register(BiocParallel::SerialParam())
p <- SerialParam()

# feature detection
fd <- xcmsSet(files_all,
              polarity = "positive",
              BPPARAM = p,
              method = "centWave",
              ppm = resultPeakpicking$best_settings$parameters$ppm, 
              peakwidth = c(resultPeakpicking$best_settings$parameters$min_peakwidth, resultPeakpicking$best_settings$parameters$max_peakwidth),
              snthresh = resultPeakpicking$best_settings$parameters$snthresh,
              prefilter = c(resultPeakpicking$best_settings$parameters$prefilter, resultPeakpicking$best_settings$parameters$value_of_prefilter),
              mzCenterFun = resultPeakpicking$best_settings$parameters$mzCenterFun,
              integrate = resultPeakpicking$best_settings$parameters$integrate,
              mzdiff = resultPeakpicking$best_settings$parameters$mzdiff,
              fitgauss = resultPeakpicking$best_settings$parameters$fitgauss,
              noise = resultPeakpicking$best_settings$parameters$noise,
              verbose.columns = resultPeakpicking$best_settings$parameters$verbose.columns)

# retention time correction
rc <- retcor(fd,
             method = resultRetcorGroup$best_settings$retcorMethod,
             plottype = resultRetcorGroup$best_settings$plottype,
             profStep = resultRetcorGroup$best_settings$profStep,
             center = resultRetcorGroup$best_settings$center,
             response = resultRetcorGroup$best_settings$response,
             distFunc = resultRetcorGroup$best_settings$distFunc,
             gapInit = resultRetcorGroup$best_settings$gapInit,
             gapExtend = resultRetcorGroup$best_settings$gapExtend,
             factorDiag = resultRetcorGroup$best_settings$factorDiag,
             factorGap = resultRetcorGroup$best_settings$factorGap,
             localAlignment = resultRetcorGroup$best_settings$localAlignment)

# reak grouping
pg <- group(rc,
            method = "density",
            minfrac = resultRetcorGroup$best_settings$minfrac, 
            minsamp = resultRetcorGroup$best_settings$minsamp, 
            bw = resultRetcorGroup$best_settings$bw,
            mzwid = resultRetcorGroup$best_settings$mzwid,
            max = resultRetcorGroup$best_settings$max)

# peak filling
fp <- fillPeaks(pg)

# final feature table
ft <- groupval(fp, "medret", "into") # "index" for height "into" for area
mz <- as.numeric(xcms::groups(fp)[, 1])
time <- as.numeric(xcms::groups(fp)[, 4])
rownames(ft) <- paste(mz, time, sep = " / ")
fwrite(as.data.frame(t(ft)), "xcms old funs.csv", row.names = T)

##############################################################################################################################################################
# XCMS for improve integration and alignment
##############################################################################################################################################################

library(xcms)

# load feature detection xcmsSet object
load("xcms obj feat_det.RData") # filename and the same folder of raw data as in section "XCMS with the best params (from IPO)"
               
# parallel processing
cores = detectCores()-1
register(bpstart(SnowParam(cores)))
BiocParallel::register(BiocParallel::SerialParam())

######################################################### Improve data after peak detection

# Remove all peaks with a width larger "rt_max" seconds
rt_max <- 60 # define the max width of peak, adjust to your data
data_rt_clean <- refineChromPeaks(feat_det, param = CleanPeaksParam(maxPeakwidth = rt_max))

# Remove all peaks with a maximal intensity below "int_min" (for xcms version 3.12.0)
int_min <- 50000 # define the min intensity of peak, adjust to your data
data_int_clean <- refineChromPeaks(feat_det, param = FilterIntensityParam(threshold = int_min))

# Merge peaks by rt, mz, proportion, see ?MergeNeighboringPeaksParam
data_merge <- refineChromPeaks(feat_det, MergeNeighboringPeaksParam(minProp = 0.05, expandRt = 4, expandMz = 0.005, ppm = 50)) # adjust to your data

######################################################### Grouping features
library(MsFeatures)

# By rt parameter
xdata <- groupFeatures(pk_gr, param = SimilarRtimeParam(10)) # adjust to your data
table(featureGroups(xdata))
plotFeatureGroups(xdata)
grid()

# By abundance correlation parameter
xdata <- groupFeatures(pk_gr, AbundanceSimilarityParam(threshold = 0.7, transform = log2), filled = TRUE) # adjust to your data
table(featureGroups(xdata))
plotFeatureGroups(xdata)
grid()

# By similarity of their EICs
xdata <- groupFeatures(pk_gr, EicSimilarityParam(threshold = 0.7, n = 2)) # adjust to your data
table(featureGroups(xdata))
plotFeatureGroups(xdata)
grid()

######################################################### Subset-based alignment
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

# Remove previous alignment results
feat_det <- dropAdjustedRtime(feat_det)
               
# Define the experimental layout
feat_det$sample_type <- "study"
feat_det$sample_type <- n_gr_t[,1]

# Initial peak grouping. Use sample_type as grouping variable
pdp_subs <- PeakDensityParam(sampleGroups = feat_det$sample_type,
                             minFraction = min_frac_man) # adjust to your data
xdata <- groupChromPeaks(feat_det, param = pdp_subs)
               
# Identify overlapping features
overlappingFeatures(xdata)

# Define subset-alignment options and perform the alignment
align_subs <- ObiwarpParam(
  binSize = resultRetcorGroup$best_settings$profStep,
  center = resultRetcorGroup$best_settings$center,
  response = resultRetcorGroup$best_settings$response,
  distFun = resultRetcorGroup$best_settings$distFunc,
  gapInit = resultRetcorGroup$best_settings$gapInit,
  gapExtend = resultRetcorGroup$best_settings$gapExtend,
  factorDiag = resultRetcorGroup$best_settings$factorDiag,
  factorGap = resultRetcorGroup$best_settings$factorGap,
  localAlignment = ifelse(resultRetcorGroup$best_settings$localAlignment==0, F,T),
  subset = which(feat_det$sample_type == "QC"), # adjust to your data, select base group for alignment
  subsetAdjust = "average")

xdata <- adjustRtime(xdata, param = align_subs) # final file with specific sample group alignment

# Define subset-alignment options and perform the alignment (other type)
pgp_subs <- PeakGroupsParam(minFraction = min_frac_man, # or resultRetcorGroup$best_settings$minfrac
                            subset = which(feat_det$sample_type == "QC"), # adjust to your data, select base group for alignment
                            subsetAdjust = "average", span = 0.4) # adjust to your data

xdata <- adjustRtime(xdata, param = pgp_subs) # final file with specific sample group alignment

##############################################################################################################################################################
# Warpgroup for increase precision
##############################################################################################################################################################

library(xcms)
library(warpgroup)
library(data.table)
library(doParallel)

# increase memory
memory.limit(999999)
memory.size(max = TRUE)
invisible(utils::memory.limit(999999))

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

# Perform Warpgroup
detach("package:dplyr", unload = T)
xr.l = llply(pk_gr@filepaths, xcmsRaw, profstep=0) # adjust to your data
wg = group.warpgroup(pk_gr, xr.l = xr.l, smooth.n = 1, rt.max.drift = 15, ppm.max.drift = ppm_diff, rt.aligned.lim = 5) # adjust to your data

# Save
save(wg, file = "warpgroup res.RData")
data <- xcms::groupval(wg, "medret", "into")
mz <- as.numeric(xcms::groups(wg)[, 1])
time <- as.numeric(xcms::groups(wg)[, 4])
rownames(data) <- paste(mz, time, sep = " / ")
load("xcms obj pk_fil.RData") # load xcms object
ft_tbl <- featureValues(pk_fil, value = "into")
colnames(data) <- colnames(ft_tbl)
data <- as.data.frame(t(data))
fwrite(data, "xcms after Warpgroup.csv", row.names = T)

##############################################################################################################################################################
# ncGTW for realignment
##############################################################################################################################################################

library(xcms)
library(ncGTW)
library(data.table)
library(dplyr)
library(stringr)
library(doParallel)

# increase memory
memory.limit(999999)
memory.size(max = TRUE)
invisible(utils::memory.limit(999999)) 

# parallel processing
cores = detectCores()-1
register(bpstart(SnowParam(cores)))
BiocParallel::register(BiocParallel::SerialParam())

# folder with files
wd_2 <- c("D:/...") # folder with all study samples
files_all <- list.files(wd_2, recursive = TRUE, full.names = TRUE, pattern = ".CDF") # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.

# order files by run order
rname <- files_all # obtain all info from rownames
rname <- str_remove(rname, wd_2) # remove folder name
rname <- str_remove(rname, ".CDF") # remove some pattern from vendor-specific format
all_id <- lapply(1:length(rname), function(y) unlist(str_split(rname[y], " ")), simplify = F) # split info from rownames
ro_id <- as.numeric(unlist(lapply(all_id, function(y) unlist(y[1])))) # obtain run order ID (every [1] element)
files_all_df <- as.data.frame(cbind(files_all, ro = as.numeric(ro_id)))
files_all_df <- files_all_df[order(as.numeric(files_all_df$ro), decreasing = F),]
files_all <- files_all_df[,-2]

# load optimal params from IPO
load("IPO_optimiz_xcms_CDF_7QC.RData")
load("IPO_optimiz_ret_align_CDF_7QC.RData")
ppm_ipo <- resultPeakpicking[["best_settings"]][["parameters"]][["ppm"]]
mzdiff_ipo <- resultPeakpicking[["best_settings"]][["parameters"]][["mzdiff"]]
bw_ipo <- resultRetcorGroup[["best_settings"]][["bw"]]
minSamples_ipo <- resultRetcorGroup$best_settings$minsamp
min_frac_man <- 0.2 # adjust to your data

# feature detection
feat_det <- xcmsSet(files_all, method="centWave", peakwidth=c(resultPeakpicking$best_settings$parameters$min_peakwidth, resultPeakpicking$best_settings$parameters$max_peakwidth),
                    ppm=resultPeakpicking$best_settings$parameters$ppm, noise= resultPeakpicking$best_settings$parameters$noise,
                    snthresh = resultPeakpicking$best_settings$parameters$snthresh, integrate = resultPeakpicking$best_settings$parameters$integrate,
                    prefilter = c(resultPeakpicking$best_settings$parameters$prefilter, resultPeakpicking$best_settings$parameters$value_of_prefilter))

# or:
# load(file = "xcms obj feat_det.RData")                              
# feat_det <- as(feat_det, "xcmsSet")

# peak grouping
pk_gr <- group(feat_det, mzwid=mzdiff_ipo, bw=bw_ipo)

# retention time correction
missing <- minSamples_ipo # number of missing samples to allow in retention time correction groups, adjust to your data
ret_cor <- retcor(pk_gr, missing=missing)

#  set min and max bw values
bw_min <- bw_ipo/5
bw_max <- bw_ipo*10

# perform ncGTW
xcms_max_bw <- group(ret_cor, mzwid=mzdiff_ipo, bw=bw_max) # group by max bw
xcms_min_bw <- group(ret_cor, mzwid=mzdiff_ipo, bw=bw_min, minfrac=min_frac_man) # group by min bw, set minFrac to your data
excluGroups <- misalignDetect(xcms_max_bw, xcms_min_bw, ppm_ipo)
show(excluGroups)
rtAdd <- 15 # rt range, adjust to your data
workers <- 4 # for parallel processing, adjust to your PC
ncGTWinputs <- loadProfile(files_all, excluGroups, mzAdd = mzdiff_ipo, rtAdd = rtAdd, BPPARAM = BiocParallel::SnowParam(workers = workers))
for (n in seq_along(ncGTWinputs))
  plotGroup(ncGTWinputs[[n]], slot(xcms_max_bw, 'rt')$corrected, ind=n)
ncGTWoutputs <- vector('list', length(ncGTWinputs))
ncGTWparam <- new("ncGTWparam")
parSamp <- 5 # for parallel processing, adjust to your PC
for (n in seq_along(ncGTWinputs))
  ncGTWoutputs[[n]] <- ncGTWalign(ncGTWinputs[[n]], xcms_max_bw, parSamp=parSamp,
                                  bpParam=SnowParam(workers=workers), ncGTWparam=ncGTWparam)

ncGTWres <- xcms_max_bw # Perform ncGTW alignment
ncGTWRt <- vector('list', length(ncGTWinputs))
for (n in seq_along(ncGTWinputs)){
  adjustRes <- adjustRT(ncGTWres, ncGTWinputs[[n]], ncGTWoutputs[[n]], ppm_ipo)
  peaks(ncGTWres) <- ncGTWpeaks(adjustRes)
  ncGTWRt[[n]] <- rtncGTW(adjustRes)
}

for (n in seq_along(ncGTWinputs))
  plotGroup(ncGTWinputs[[n]], ncGTWRt[[n]], ind = n)

groups(ncGTWres) <- excluGroups[ , 2:9]
groupidx(ncGTWres) <- groupidx(xcms_max_bw)[excluGroups[ , 1]]
rtCor <- vector('list', length(files_all)) # Only consider the misaligned features
for (n in seq_along(files_all)){
  rtCor[[n]] <- vector('list', dim(excluGroups)[1])
  for (m in seq_len(dim(xcms::groups(ncGTWres))[1]))
    rtCor[[n]][[m]] <- ncGTWRt[[m]][[n]]
}
slot(ncGTWres, 'rt')$corrected <- rtCor
XCMSres <- xcms_max_bw
groups(XCMSres) <- excluGroups[ , 2:9]
groupidx(XCMSres) <- groupidx(xcms_max_bw)[excluGroups[ , 1]]

assignInNamespace("fillPeaksChromPar", ncGTW:::fillPeaksChromPar, ns="xcms",
                  envir=as.environment("package:xcms"))

ncGTWresFilled <- fillPeaks(ncGTWres)
XCMSresFilled <- fillPeaks(XCMSres)

# join peak table and save
data <- xcms::groupval(ncGTWresFilled, "medret", "into")
mz <- as.numeric(xcms::groups(ncGTWres)[, 1])
time <- as.numeric(xcms::groups(ncGTWres)[, 4])
rownames(data) <- paste(mz, time, sep = " / ")
load("xcms obj pk_fil.RData") # load xcms object
ft_tbl <- featureValues(pk_fil, value = "into")
colnames(data) <- colnames(ft_tbl)
data <- as.data.frame(t(data))
fwrite(data, "xcms after ncGTW.csv", row.names = T)

# compute CV
compCV(XCMSresFilled)
compCV(ncGTWresFilled)

##############################################################################################################################################################
# cpc for filtering peaks
##############################################################################################################################################################

library(xcms)
library(cpc)

# load xcms object after peak detection
load("xcms obj feat_det.RData") # filename and the same folder of raw data as in section "XCMS with the best params (from IPO)" 

# perform
cpc <- cpc::filter_xcms_peaklist(xd = feat_det, return_type = "cpc", 
                                 param = cpc::cpcProcParam(min_sn = 10, # adjust to your data
                                                           min_pts = 10, # adjust to your data
                                                           min_intensity = 2000)) # adjust to your data

cpcPeaktable <- cpc::cpt(cpc)

# get the filtered xcms object using getFilteredXCMS()
xdFilt <- cpc::getFilteredXCMS(cpc)

# get the filtered peak table
xcmsPeaktableFilt <- data.frame(xcms::chromPeaks(xdFilt))

# get the filtered cpc peak table
cpcPeaktableFilt <- cpcPeaktable[row.names(xcmsPeaktableFilt), ]

##############################################################################################################################################################
# Description of peaks table
##############################################################################################################################################################

# load data set
peak_table <- as.data.frame(fread("xcms after IPO.csv"))
rownames(peak_table) <- peak_table[,1]
peak_table <- peak_table[,-1]

# dataset
ds_pt <- peak_table

# Missing value %
tn <- nrow(ds_pt)*ncol(ds_pt)
mv_c <- sum(is.na(ds_pt)) # or length(which(is.a(...)))
pr_mv <- round(mv_c/tn*100,0)
pr_mv

# Number of peaks
ncol(ds_pt) # or nrow

##############################################################################################################################################################
# EIC peaks integration
##############################################################################################################################################################

library(xcms)
library(RColorBrewer)

# Load data
pd <- data.frame(sample_name = sub(basename(files_QC), pattern = ".CDF", replacement = "", fixed = TRUE), stringsAsFactors = FALSE) # just example, adjust to your study
files <- readMSData(files_QC, mode = "onDisk", pdata = new("NAnnotatedDataFrame", pd)) # just example, adjust to your study

# Define m/z value
mz <- 312.2150 + c(-0.5, 0.5) # adjust to your study

# Define rt value
rt <- 936 + c(-100,100)

# Extract EIC
chr <- chromatogram(files, rt = rt, mz = mz, aggregationFun = "sum")

# Plot EICs
group_colors <- paste0(brewer.pal(length(files_QC), "Set1"), "60")
names(group_colors) <- files_QC
plot(chr, col = group_colors[files_QC])

# Plot EICs with mass
files %>%
  filterRt(rt = rt) %>%
  filterMz(mz = mz) %>%
  plot(type = "XIC")

# Peak detection
cwp <- CentWaveParam(peakwidth = c(5, 30), noise = 200, ppm = 80, prefilter = c(3, 400), snthresh = 10) # adjust to your data

xchr <- findChromPeaks(chr, param = cwp)
chromPeaks(xchr) # peak table

plot(xchr, col = group_colors, peakBg = group_colors[chromPeaks(xchr)[, "column"]]) # plot results

##############################################################################################################################################################
# References
##############################################################################################################################################################

# 1. Tautenhahn, Ralf, Christoph Boettcher, and Steffen Neumann. "Highly sensitive feature detection for high resolution LC/MS." BMC bioinformatics 9.1 (2008): 504.
# 2. McLean, Craig, and Elizabeth B. Kujawinski. "AutoTuner: high fidelity and robust parameter selection for metabolomics data processing." Analytical chemistry 92.8 (2020): 5724-5732.
# 3. Pang, Zhiqiang, et al. "MetaboAnalystR 3.0: Toward an optimized workflow for global metabolomics." Metabolites 10.5 (2020): 186.
# 4. Guo, Jian, Sam Shen, and Tao Huan. "Paramounter: Direct Measurement of Universal Parameters To Process Metabolomics Data in a "White Box"." Analytical Chemistry (2022).
# 5. Libiseller, Gunnar, et al. "IPO: a tool for automated optimization of XCMS parameters." BMC bioinformatics 16.1 (2015): 118.
# 6. FernÃ¡ndez-Ochoa, Ãlvaro, et al. "A Case Report of Switching from Specific Vendor-Based to R-Based Pipelines for Untargeted LC-MS Metabolomics." Metabolites 10.1 (2020): 28.
# 7. Mahieu, Nathaniel G., Jonathan L. Spalding, and Gary J. Patti. "Warpgroup: increased precision of metabolomic data processing by consensus integration bound analysis." Bioinformatics 32.2 (2016): 268-275.
# 8. Wu, Chiung-Ting, et al. "Targeted realignment of LC-MS profiles by neighbor-wise compound-specific graphical time warping with misalignment detection." Bioinformatics 36.9 (2020): 2862-2871.
# 9. Pirttilä, Kristian, et al. "Comprehensive Peak Characterization (CPC) in Untargeted LC-MS Analysis." Metabolites 12.2 (2022): 137.
