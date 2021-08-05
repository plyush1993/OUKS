##############################################################################################################################################################
# Table of contents
##############################################################################################################################################################

# Installation
# Metadata generating & repeated measurements calculation
# QC RSD filtering
# Combine and filter with annotation data
# Descriptive Statistics filtering
# Missing Value filtering
# Blank filtering
# Pearson Correlation filtering of repeats
# References

##############################################################################################################################################################
# Installation
##############################################################################################################################################################

# setup environment
library(data.table)
library(dplyr)
library(stringr)
setwd("D:/...")

# Typical Operation Order
# 1. RSD filtering 
# 2. (Optional) Descriptive Statistics/Missing Value/Blank/Pearson Correlation filtering 
# 3. Repeats calculation
# 4. Annotation filtering
# Note: see "Supplementary code" subtitle. Colnames format after corrections can differ from raw data colnames.

##############################################################################################################################################################
# Metadata generating & repeated measurements calculation
##############################################################################################################################################################

# load peak table in csv
# load in format: 1st column: 1. qc1(s78_1) b1 QC1(KMS 53).X --- ro. nameN(nameN_NREPEAT) BatchIndex QCN(Label N)
# "150. qc30 b6 QC30.CDF" 
# "167. s64_2 b7 MS 45.CDF" 
dsr <-as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter.csv", header=T)) # first column with all metadata
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]

########################################## METADATA GENERATING

rname <- rownames(dsr) # obtain all info from rownames
rname <- str_remove(rname, ".CDF") # remove some pattern from vendor-specific format
# qc_id <- grep(pattern = "QC", x = rname) # find by pattern in info from rownames
# dsr_qc <- dsr[qc_id,] # select only observations by pattern from peak table

# na_var <- apply(dsr,2,function(x) any(is.na(x))) # obtain features with NA
# dsr_no_na <- dsr[,na_var==F] # select only features without NA

all_id <- sapply(1:length(rname), function(y) unlist(str_split(rname[y], " "))) # split info from rownames
ro_id <- as.numeric(unlist(lapply(all_id, function(y) unlist(y[1])))) # obtain run order ID (every [1] element)
b_id <- unlist(lapply(all_id, function(y) unlist(y[3]))) # obtain batch ID (every [3] element)
s_id <- unlist(lapply(all_id, function(y) unlist(y[2]))) # obtain sample ID (every [2] element)
p1_id <- unlist(lapply(all_id, function(y) unlist(y[4]))) # obtain patient ID (every [4] element)
p2_id <- unlist(lapply(all_id, function(y) unlist(y[5]))) # obtain patient ID (every [5] element)
p_all_id <- sapply(1:length(p1_id), function(y) ifelse(is.na(p2_id[y]), p1_id[y] ,paste(p1_id[y], p2_id[y], sep = " "))) # join patient ID (every [4-5] element)
un_raw_l <- unique(gsub("[[:digit:]]", "", p1_id)) # obtain unique raw label from p1_id (every [4] element)
true_l <- c("QC", "TG", "CG") # write desired label in order as in unique raw label (should be analogical)
rbind(un_raw_l,true_l) # visual check agreement between true and raw labels
raw_l <- gsub("[[:digit:]]", "", p1_id) # obtain raw label from p1_id (every [4] element)
n_gr_t <- as.character(match(raw_l, un_raw_l)) # obtain index for raw label by unique raw label
for (i in 1:length(unique(n_gr_t))) {
  n_gr_t <- str_replace(n_gr_t, unique(n_gr_t)[i], true_l[i]) } # exchange raw label by true label via index
s_id_p <- sapply(1:length(s_id), function(y) unlist(str_split(s_id[y], "_"))) # obtain sample ID (every [2] element) with splitting by "_"
s_id_p2 <- unlist(lapply(s_id_p, function(y) unlist(y[1]))) # obtain sample ID (every [2] element) with only one part of sample ID
all_meta_data<- as.data.frame(cbind(n_gr_t, ro_id, b_id, s_id, p_all_id, raw_l, s_id_p2)) # obtain all metadata for all repeats

# join features and metadata and order
ds <- data.frame(cbind(all_meta_data, dsr))
ds$ro_id <- as.numeric(ds$ro_id)
ds <- ds[order(ds$ro_id, decreasing = F),] 

########################################## REPEATED MEASURMENTS CALCULATIONS

counting <- sapply(1:length(s_id_p2), function(i) subset(s_id_p2, s_id_p2==s_id_p2[i])) # counting all repeats for sample ID (every [2] element) 
n_rep <- sapply(1:length(counting), function(y) length(counting[[y]])) # count number of repeats for sample ID (every [2] element)
rep_id <- t(data.frame(unique(counting[which(n_rep>1)]))) # remain only repeated (>1) unique sample ID (every [2] element)
all_meta_data_rep <- lapply(c(1:nrow(rep_id)), function(y) filter(all_meta_data,  s_id_p2 == rep_id[y,1])[1,]) # obtain all metadata for repeats
all_meta_data_rep <- do.call("rbind",all_meta_data_rep) # obtain all metadata for repeats
dsr2 <- data.frame(cbind(all_meta_data$s_id_p2, dsr)) # join sample ID (every [2] element) and peak table
dsr3 <- data.frame(t(sapply(1:nrow(rep_id), function(t) apply(filter(dsr2[,-1], dsr2$all_meta_data.s_id_p2 == rep_id[t,1]), 2, function(x) mean(x, na.rm = T))))) # obtain mean value for all variables for all repeats
dsr3 <- cbind(all_meta_data_rep, dsr3) # join metadata for repeats and mean value for all variables for all repeats
dsr3 <- dsr3[order(as.numeric(dsr3$ro_id), decreasing = F),] # order by run order 
colnames(dsr3)[-c(1:7)] <- colnames(dsr) # names as in loaded files

# save output 
fwrite(dsr3, "xcms after IPO MVI QC-XGB filter repeats.csv", row.names = T)
# fwrite(all_meta_data, "all_meta_data.csv", row.names = T)

##############################################################################################################################################################
# QC RSD filtering
##############################################################################################################################################################

# Settings
setwd("D:/...") # load peak table in csv
# load in format: 1st column: 1. qc1(s78_1) b1 QC1(KMS 53).X --- ro. nameN(nameN_NREPEAT) BatchIndex QCN(Label N)
# "150. qc30 b6 QC30.CDF" 
# "167. s64_2 b7 MS 45.CDF" 
dsr <-as.data.frame(fread(input = "xcms after IPO MVI QC-XGB.csv", header=T))
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]

####### Supplementary code
# colnames as from xcms
raw_data <-as.data.frame(fread(input = "xcms after IPO MVI.csv", header=T)) # load data
rownames(raw_data) <- raw_data[,1] # load data
raw_data <- raw_data[,-1] # load data
colnames(dsr) <- colnames(raw_data) # rename colnames as in raw data from xcms

# metadata generating
rname_all <- rownames(dsr) # obtain all info from rownames
qc_id <- grep(pattern = "QC", x = rname_all) # find by pattern in info from rownames
dsr_qc <- dsr[qc_id,] # select only observations by pattern from peak table
rname <- rownames(dsr_qc) # obtain all info from rownames
rname <- str_remove(rname, ".CDF") # remove some pattern from vendor-specific format
all_id <- as.data.frame(sapply(1:length(rname), function(y) unlist(str_split(rname[y], " ")))) # split info from rownames
b_id <- data.frame(t(all_id[3,])) # obtain batch ID (3 row)

ds <- cbind(b_id, dsr_qc)
colnames(ds)[1] <- "batch"

# rename batch index to numeric value
ds[,1] <- str_remove(ds[,1], "b")
ds[,1] <- as.numeric(ds[,1])
ds <- ds[order(ds$batch, decreasing = F),]

# Select batch by pattern
tn <- ds$batch # batch variables
tnu <- unique(tn)
vb <- list()
for (i in (1:length(tnu))){
  vb[[i]] <- which(tn == tnu[i])}

# RSD calculation by batches
RSD_by_batch_results <- list()
RSD_by_batch <- function (x) { 
  for (i in (1:length(vb))){
    RSD_by_batch_results[[i]] <- apply(x[vb[[i]],], 2, function(y) sd(y, na.rm = T)/mean(y, na.rm = T)*100)}
  RSD_by_batch_results[[(length(vb)+1)]] <- apply(x, 2, function(y) sd(y, na.rm = T)/mean(y, na.rm = T)*100) # all batches
  for (i in (1:length(RSD_by_batch_results))){
    RSD_by_batch_results[[i]] <- data.frame(colnames(x),RSD_by_batch_results[[i]])}
  for (i in (1:length(RSD_by_batch_results))){
    colnames(RSD_by_batch_results[[i]]) <- c("name", "rsd")}
  n <- c(paste(c(1:length(vb)), "batch"), "all batches")
  names(RSD_by_batch_results) <- n
  RSD_by_batch_results
}

n <- c(paste(c(1:length(vb)), "batch"), "all batches") # names
ds1 <- ds[,-1] # sample in row only metabolites variables
rsd <- RSD_by_batch(ds1)
rsd_all_b <- rsd[[length(rsd)]]

# cutoff by RSD in QC
cutoff <- 30
cutoff_all <- lapply(1:length(rsd), function(y) nrow(filter(rsd[[y]], rsd < cutoff)))
cutoff_f <- filter(rsd_all_b, rsd <= cutoff)
nrow(cutoff_f) # number of retained features
round(nrow(cutoff_f)/ncol(dsr)*100, 0) # percent of retained features
nan_l <- length(which(is.nan(rsd_all_b$rsd))) # number of NaN (when 0 values) 
dsr_QC_RSD <- dsr[,cutoff_f$name]

# save
fwrite(dsr_QC_RSD, "xcms after IPO MVI QC-XGB filter.csv", row.names = T)

##############################################################################################################################################################
# Combine and filter with annotation data
##############################################################################################################################################################

library(data.table)
library(dplyr)

# load data after all steps and preprocess
dsr <-as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats.csv", header=T)) 
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]
dsr <- dsr[,-c(4,6,7)]
colnames(dsr)[1:4] <- c("Label", "Injection", "Batch", "PatientID")
dsr_t <- data.frame(t(dsr))

########################################## ANNOTATION
# load annotation data
# CAMERA
camera <-as.data.frame(fread(input = "xcms CAMERA.csv", header=T))[,-1]
camera <- filter(camera, camera$annot_camera.name %in% rownames(dsr_t)[-c(1:4)])
identical(rownames(dsr_t)[-c(1:4)], camera$annot_camera.name) # check identical names
# xMSannotator
xma <- as.data.frame(fread(input = "xcms xMSannotator.csv", header=T))[-1]
xma <- filter(xma, xma$`raw data name` %in% rownames(dsr_t)[-c(1:4)])
identical(rownames(dsr_t)[-c(1:4)], xma$`raw data name`) # check identical names
# RAMClustR
rcr <- as.data.frame(fread(input = "xcms RAMClust.csv", header=T))[-1]
rcr <- filter(rcr, rcr$`raw data name` %in% rownames(dsr_t)[-c(1:4)])
identical(rownames(dsr_t)[-c(1:4)], rcr$`raw data name`) # check identical names
# mWISE
mwise <- as.data.frame(fread(input = "xcms mWISE.csv", header=T))[-1]
mwise <- filter(mwise, mwise$`raw name` %in% rownames(dsr_t)[-c(1:4)])
identical(rownames(dsr_t)[-c(1:4)], mwise$`raw name`) # check identical names

# combine with annotation data
dsr_t_a <- cbind(camera[,3:4], rcr_adduct = rcr[,5], xma[,c(4,6:9)], mwise[,c(4,7:9)], dsr_t[-c(1:4),])
title <- cbind(do.call("rbind", replicate(4, rep(NA, ncol(dsr_t_a)-ncol(dsr_t)), simplify = F)), dsr_t[c(1:4),])
colnames(title) <- colnames(dsr_t_a)
dsr_t_a <- rbind(title, dsr_t_a)
dsr_a <- data.frame(t(dsr_t_a))

# transform to na
empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) # since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}

dsr_a <- mutate_each(dsr_a, funs(empty_as_na)) 

# save
colnames(dsr_a) <- colnames(dsr)
fwrite(dsr_a, "xcms after IPO MVI QC-XGB filter repeats annot.csv", row.names = T)

########################################## FILTRATION
# load data after all steps and preprocess
dsr_a <-as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats annot.csv", header=T))
rownames(dsr_a) <- dsr_a[,1]
dsr_a <- dsr_a[,-1]

# transform to na
empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) # since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}

dsr_a <- mutate_each(dsr_a, funs(empty_as_na)) 

# remove features without annotation
cna <- apply(dsr_a[,-c(1:4)], 2, function(x) length(which(is.na(x))))
ind <- which(cna==12)+4 # if integer(0) --> all metabolites are retained
dsr_a_f <- dsr_a[,-c(ind)]

# save
fwrite(dsr_a_f, "xcms after IPO MVI QC-XGB filter repeats annot+filtr.csv", row.names = T)

##############################################################################################################################################################
# Descriptive Statistics filtering
##############################################################################################################################################################

# Settings
library(data.table)
# load peak table in csv
# load in format: 1st column: 1. qc1(s78_1) b1 QC1(KMS 53).X --- ro. nameN(nameN_NREPEAT) BatchIndex QCN(Label N)
# "150. qc30 b6 QC30.CDF" 
# "167. s64_2 b7 MS 45.CDF" 
dsr <-as.data.frame(fread(input = "xcms after IPO MVI.csv", header=T)) # first column with all metadata
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]

####### Select only QC samples
# colnames as from xcms
raw_data <-as.data.frame(fread(input = "xcms after IPO MVI.csv", header=T)) # load data
rownames(raw_data) <- raw_data[,1] # load data
raw_data <- raw_data[,-1] # load data
colnames(dsr) <- colnames(raw_data) # rename colnames as in raw data from xcms

# metadata generating
rname_all <- rownames(dsr) # obtain all info from rownames
qc_id <- grep(pattern = "QC", x = rname_all) # find by pattern in info from rownames
dsr_qc <- dsr[qc_id,] # select only observations by pattern from peak table
rname <- rownames(dsr_qc) # obtain all info from rownames
rname <- str_remove(rname, ".CDF") # remove some pattern from vendor-specific format
all_id <- as.data.frame(sapply(1:length(rname), function(y) unlist(str_split(rname[y], " ")))) # split info from rownames
b_id <- data.frame(t(all_id[3,])) # obtain batch ID (3 row)
dsr <- dsr_qc

# Descriptive Statistics calculation
dsc <- apply(dsr, 2, function(y) mean(y, na.rm = T)) # try also "min", "max", "mean", "median", "sd" funs
dsc <- data.frame(names(dsc), as.numeric(dsc))
colnames(dsc) <- c("name", "stats")

# cutoff by Descriptive Statistics
cutoff <- 30000 # set
cutoff_f <- filter(dsc, dsc$stats >= cutoff)
dsr_DS <- dsr[,cutoff_f$name]

# save
fwrite(dsr_DS, "xcms after IPO MVI filter DS.csv", row.names = T)

##############################################################################################################################################################
# Missing Value filtering
##############################################################################################################################################################

# Settings
library(data.table)
# load peak table in csv
# load in format: 1st column: 1. qc1(s78_1) b1 QC1(KMS 53).X --- ro. nameN(nameN_NREPEAT) BatchIndex QCN(Label N)
# "150. qc30 b6 QC30.CDF" 
# "167. s64_2 b7 MS 45.CDF" 
dsr <-as.data.frame(fread(input = "xcms after IPO.csv", header=T)) # first column with all metadata
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]

####### Select only QC samples
# colnames as from xcms
raw_data <-as.data.frame(fread(input = "xcms after IPO MVI.csv", header=T)) # load data
rownames(raw_data) <- raw_data[,1] # load data
raw_data <- raw_data[,-1] # load data
colnames(dsr) <- colnames(raw_data) # rename colnames as in raw data from xcms

# metadata generating
rname_all <- rownames(dsr) # obtain all info from rownames
qc_id <- grep(pattern = "QC", x = rname_all) # find by pattern in info from rownames
dsr_qc <- dsr[qc_id,] # select only observations by pattern from peak table
rname <- rownames(dsr_qc) # obtain all info from rownames
rname <- str_remove(rname, ".CDF") # remove some pattern from vendor-specific format
all_id <- as.data.frame(sapply(1:length(rname), function(y) unlist(str_split(rname[y], " ")))) # split info from rownames
b_id <- data.frame(t(all_id[3,])) # obtain batch ID (3 row)
dsr <- dsr_qc

# replace 0 to NA
dsr[dsr == 0] <- NA

# remove column with NA
colna <- apply(dsr,2,function(x) any(is.na(x)))
dsr_no_NA <- dsr[,colna==F]

# MV total amount calculation
ds_mv <- dsr # remove label and other categorical data
tn <- nrow(ds_mv)*ncol(ds_mv)
mv_c <- sum(is.na(ds_mv)) # or lenght(which(is.a(...)))
pr_mv <- round(mv_c/tn*100,0)

# MV calculation and filtering
mv_f <- data.frame(apply(dsr, 2, function(y) round(sum(is.na(y))/length(y)*100, 0)))
colnames(mv_f) <- "mv"
cutoff <- 50 # set
mv_f_n <- rownames(filter(mv_f, mv < cutoff))
dsr_MV <- dsr[,mv_f_n]

# save
fwrite(dsr_MV, "xcms after IPO MVI filter MV.csv", row.names = T)

##############################################################################################################################################################
# Blank filtering
##############################################################################################################################################################

library(data.table)
library(dplyr)
library(stringr)

# Settings
# load peak table in csv
# load in format: 1st column: 1. qc1(s78_1) b1 QC1(KMS 53).X --- ro. nameN(nameN_NREPEAT) BatchIndex QCN(Label N)
# "150. qc30 b6 QC30.CDF" 
# "167. s64_2 b7 MS 45.CDF" 
dsr <-as.data.frame(fread(input = "xcms after IPO MVI.csv", header=T)) # first column with all metadata
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]

# Group Info
rname <- rownames(dsr) # obtain all info from rownames
rname <- str_remove(rname, ".CDF") # remove some pattern from vendor-specific format
all_id <- sapply(1:length(rname), function(y) unlist(str_split(rname[y], " "))) # split info from rownames
p1_id <- unlist(lapply(all_id, function(y) unlist(y[4]))) # obtain patient ID (every [4] element)
p2_id <- unlist(lapply(all_id, function(y) unlist(y[5]))) # obtain patient ID (every [5] element)
p_all_id <- sapply(1:length(p1_id), function(y) ifelse(is.na(p2_id[y]), p1_id[y] ,paste(p1_id[y], p2_id[y], sep = " "))) # join patient ID (every [4-5] element)
un_raw_l <- unique(gsub("[[:digit:]]", "", p1_id)) # obtain unique raw label from p1_id (every [4] element)
true_l <- c("QC", "TG", "CG") # write desired label in order as in unique raw label (should be analogical)
rbind(un_raw_l,true_l) # visual check agreement between true and raw labels
raw_l <- gsub("[[:digit:]]", "", p1_id) # obtain raw label from p1_id (every [4] element)
n_gr_t <- as.character(match(raw_l, un_raw_l)) # obtain index for raw label by unique raw label
for (i in 1:length(unique(n_gr_t))) {
  n_gr_t <- str_replace(n_gr_t, unique(n_gr_t)[i], true_l[i]) } # exchange raw label by true label via index

# Dataset
dsr_bl <- dsr

# Perform
identical(length(n_gr_t), nrow(dsr_bl)) # check dimensions
bl_id <- which(n_gr_t == "blank") # set blank label "blank" or which(n_gr_t != "QC") if no blank samples
cond_ind <- which(n_gr_t == "QC") # set label for specific group (for example -> QC) or cond_ind <- n_gr_t[-bl_id]
bl_med <- apply(dsr_bl[bl_id, ], 2, function(y) median(y, na.rm = T)) # median values for blanks
cond_med <- apply(dsr_bl[cond_ind, ], 2, function(y) median(y, na.rm = T)) # median values for samples
thresh <- 0.1 # set threshold
ratio <- bl_med/cond_med # calculate ratio blank/samples
ind_bl_tr <- which(ratio < thresh) # index of features with ratio < threshold

# Subset features by blank filtering
dsr_bl_tr <- dsr_bl[,ind_bl_tr]

# save
fwrite(dsr_bl_tr, "xcms after IPO MVI filter blanks.csv", row.names = T)

##############################################################################################################################################################
# Pearson Correlation filtering of repeats
##############################################################################################################################################################

library(data.table)

# after processing from all first section "Metadata generating & repeated measurements calculation"

# Pearson correlation calculation and filtering
dsr2 <- data.frame(cbind(s_id_p2, dsr)) # join sample ID (every [2] element) and peak table

# fun for Pearson cor coef if n repeats samples
filtered.cor <- function(x){ 
  x <- t(x)
  num_var <- apply(x, 2,  function(x) is.numeric(x))    
  cor_mat <- cor(x[,num_var], method = 'pearson')    
  diag(cor_mat) <- 100    
  return(cor_mat[which.min(abs(cor_mat))])}

# if n repeats
dsr3 <- sapply(1:nrow(rep_id), function(t) filtered.cor(filter(dsr2[,-1], dsr2$all_meta_data.s_id_p2 == rep_id[t,1])))

# filter by correlation
cor_v <- data.frame(dsr3)
rownames(cor_v) <- all_meta_data_rep$s_id_p2
colnames(cor_v) <- "cor"
cutoff <- 0.90 # set
pc_f <- rownames(subset(cor_v, cor > cutoff))
dsr_PR <- filter(dsr2, s_id_p2 %in% pc_f)
dsr_PR <- dsr_PR[,-1]

# save
fwrite(dsr_PR, "xcms after IPO MVI filter corr.csv", row.names = T)

##############################################################################################################################################################
# References
##############################################################################################################################################################

# 1. Mathé, Ewy A., et al. "Noninvasive urinary metabolomic profiling identifies diagnostic and prognostic markers in lung cancer." Cancer research 74.12 (2014): 3259-3270.
# 2. Martín-Blázquez, Ariadna, et al. "Untargeted LC-HRMS-based metabolomics to identify novel biomarkers of metastatic colorectal cancer." Scientific reports 9.1 (2019): 1-9.
# 3. Pezzatti, Julian, et al. "Implementation of liquid chromatography-high resolution mass spectrometry methods for untargeted metabolomic analyses of biological samples: A tutorial." Analytica Chimica Acta 1105 (2020): 28-44.
# 4. Chong, Jasmine, David S. Wishart, and Jianguo Xia. "Using MetaboAnalyst 4.0 for comprehensive and integrative metabolomics data analysis." Current protocols in bioinformatics 68.1 (2019): e86.
# 5. Klåvus, Anton, et al. ""Notame": Workflow for Non-Targeted LC-MS Metabolic Profiling." Metabolites 10.4 (2020): 135.
# 6. Helmus, Rick, et al. "patRoon: open source software platform for environmental mass spectrometry based non-target screening." Journal of Cheminformatics 13.1 (2021): 1-25.
