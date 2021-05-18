##############################################################################################################################################################
# Table of contents
##############################################################################################################################################################

# Installation
# Preprocessing
# Artifacts Removal
# Univariate MVI 
# Multivariate MVI
# Performance by NRMSE
# Perform best MVI
# References

##############################################################################################################################################################
# Installation
##############################################################################################################################################################

# setup environment
library(data.table)
library(stringr)
library(dplyr)
setwd("D:/...")

# load dataset
# load peak table in csv
# load in format: 1st column: 1. qc1(s78_1) b1 QC1(KMS 53).X --- ro. nameN(nameN_NREPEAT) BatchIndex QCN(Label N)
# "150. qc30 b6 QC30.CDF" 
# "167. s64_2 b7 MS 45.CDF"
dsr <-as.data.frame(fread(input = "xcms after IPO.csv", header=T))
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]

# obtain run order and sort
rname <- rownames(dsr) # obtain all info from rownames
rname <- str_remove(rname, ".CDF") # remove some pattern from vendor-specific format
all_id <- sapply(1:length(rname), function(y) unlist(str_split(rname[y], " "))) # split info from rownames
ro_id <- as.numeric(unlist(lapply(all_id, function(y) unlist(y[1])))) # obtain run order ID (every [1] element)
dsr <- cbind(ro_id, dsr)
dsr$ro_id <- as.numeric(dsr$ro_id)
dsr <- dsr[order(dsr$ro_id, decreasing = F),] # order by ro col
dsr <- dsr[,-1]

# prepare for mvi and artifacts removal
ds_mvi <- dsr # remove all column with non-numeric values
colnames(ds_mvi) <- paste0("V", c(1:ncol(ds_mvi))) # simplify colnames

##############################################################################################################################################################
# Preprocessing
##############################################################################################################################################################

########################################## amount of NA in %
tn <- nrow(ds_mvi)*ncol(ds_mvi)
mv_c <- sum(is.na(ds_mvi)) # or length(which(is.a(...)))
pr_mv <- round(mv_c/tn*100,0)

########################################## remove all column with NA
colna <- apply(ds_mvi,2,function(x) any(is.na(x)))
dsr_no_NA <- ds_mvi[,colna==F]

########################################## convert 0 into NA
ds_mvi[ds_mvi == 0] <- NA

########################################## One value MVI for all dataset
ds_mvi[is.na(ds_mvi)] <- 0
ds_mvi[is.na(ds_mvi)] <- min(ds_mvi, na.rm = T)

########################################## Random values (from noise) MVI
load("IPO_optimiz_xcms_CDF_7QC.RData") # load IPO parameters
noise <- resultPeakpicking$best_settings$parameters$noise # set noise level for random value generation. This can be obtained from IPO optimization: resultPeakpicking$best_settings$parameters$noise or manually, for example: 100
# noise <- as.numeric(quantile(1:noise)[2]) # reduce noise by quantile: quantile(1:noise)[1] = 0% ... quantile(1:noise)[4] = 100%
sd <- as.numeric(quantile(1:noise)[2]) # set sd value for random value generation as.numeric(quantile(1:noise)[2]) or noise*0.3
set.seed(1234)
ds_mvi[is.na(ds_mvi)] <- 0 # convert NA into 0
NAidx <- ds_mvi == 0 # NA index
imp.rand <- abs(rnorm(sum(NAidx), mean = noise, sd = sd)) # generate random values
ds_mvi[NAidx] <- imp.rand

########################################## check NA status
length(which(is.na(ds_mvi)))

##############################################################################################################################################################
# Artifacts Removal
##############################################################################################################################################################

library(MetProc)
library(stringr)
library(dplyr)

# prepare data
ds <- as.data.frame(t(ds_mvi))
or_vn <- rownames(ds) # original variable names
qc_ind <- which(str_detect(colnames(ds), "QC")==T)
colnames(ds)[qc_ind] <- paste0("QC", c(1:length(qc_ind)))
colnames(ds)[-qc_ind] <- c(1:(length(colnames(ds))-length(qc_ind)))
rownames(ds) <- paste0("M", c(1:nrow(ds)))
t_vn <- as.data.frame(cbind(or_vn, MetProc = rownames(ds))) # original variable names and for MetProc
ds <- cbind("metab" = rownames(ds), ds)
ds <- rbind(colnames(ds), ds)
fwrite(ds, "xcms after IPO ordered for MetProc.csv", row.names = T, col.names = T)

# run MetProc
metdata <- read.met("xcms after IPO ordered for MetProc.csv",
                    headrow=2, metidcol=2, fvalue=3, sep=",", ppkey="QC", ippkey="BPP")

results <- met_proc(metdata, 
                    numsplit=1, cor_rates=c(.95), runlengths=c(10), mincut=0.01, maxcut=0.95, scut=0.75,
                    plot=F, ppkey = "QC")

# filter artifacts
art <- rownames(results$remove)
art_name <- filter(t_vn, MetProc %in% art)
retained <- colnames(ds_mvi)[! colnames(ds_mvi) %in% art_name$or_vn]
ds_mp <- ds_mvi[,retained]
fwrite(ds_mp, "xcms after IPO remove artifacts.csv", row.names = T)

##############################################################################################################################################################
# Univariate MVI 
##############################################################################################################################################################

########################################## Simple MVI by column
########################################## mean
for (i in 1:ncol(ds_mvi)) {
  ds_mvi[,i][is.na(ds_mvi[,i])] <- mean(ds_mvi[,i], na.rm = T) }

########################################## median
for (i in 1:ncol(ds_mvi)) {
  ds_mvi[,i][is.na(ds_mvi[,i])] <- median(ds_mvi[,i], na.rm = T) }

########################################## minimum
for (i in 1:ncol(ds_mvi)) {
  ds_mvi[,i][is.na(ds_mvi[,i])] <- min(ds_mvi[,i], na.rm = T) }

########################################## half minimum
for (i in 1:ncol(ds_mvi)) {
  ds_mvi[,i][is.na(ds_mvi[,i])] <- (min(ds_mvi[,i], na.rm = T))/2 }

########################################## apply version example for half minimum
HM_apply <- function(data) {
  result <- as.data.frame(sapply(1:ncol(data), function(y) replace(data[,y],which(is.na(data[,y])), min(data[,y], na.rm = T)/2)))
  result
}

ds_mvi_hf <- HM_apply(ds_mvi)
colnames(ds_mvi_hf) <- colnames(ds_mvi)
rownames(ds_mvi_hf) <- rownames(ds_mvi)

##############################################################################################################################################################
# Multivariate MVI
##############################################################################################################################################################

########################################## KNN MVI
library(impute)  
ds_mvi_knn <- impute.knn(t(ds_mvi), k = 3, rowmax = 90, colmax = 90)[[1]]

########################################## CART MVI
library(mice)
mvi_cart <- mice(ds_mvi, method="cart")
ds_mvi_cart <- complete(mvi_cart)

########################################## RF MVI
# example mice package
library(mice)
mvi_rf <- mice(ds_mvi, method="rf")
ds_mvi_rf1 <- complete(mvi_rf)

# example missForest package
library(missForest)
library(doParallel)
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
ds_mvi_rf2 <- missForest(ds_mvi, maxiter = 10, ntree = 100, mtry = floor(sqrt(ncol(ds_mvi))), parallelize = "forests")[[1]]
stopCluster(cl)

# example StatTools package
library(StatTools)
library(doParallel)
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
ds_mvi_rf3 <- data.frame(rfImpWrap(MAT = ds_mvi, nCore = nCore))

########################################## PLS MVI
library(StatTools)
library(doParallel)
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
ds_mvi_pls <- data.frame(mvImpWrap(ds_mvi, method = "PLS", nCore = nCore))
stopCluster(cl)

########################################## PCA MVI
library(pcaMethods)
pc <- pca(ds_mvi, nPcs=3, method="ppca") # svd, ppca, bpca, nipals, etc, see listPcaMethods()
ds_mvi_pca <- data.frame(completeObs(pc))

########################################## QRILC MVI
library(imputeLCMD)
ds_mvi_qrilc <- impute.QRILC(ds_mvi)[[1]]

##############################################################################################################################################################
# Performance by NRMSE
##############################################################################################################################################################

# amount of NA in %
tn <- nrow(ds_mvi)*ncol(ds_mvi)
mv_c <- sum(is.na(ds_mvi)) # or length(which(is.a(...)))
pr_mv <- round(mv_c/tn*100,0)

# remove all column with NA
colna <- apply(ds_mvi,2,function(x) any(is.na(x)))
dsr_no_NA <- ds_mvi[,colna==F]

########################################## ALL GROUPS ALL METHODS

# compute NRMSE for all groups (labels) simultaneously
library(imputeR)
library(dplyr)
ds_true <- dsr_no_NA
ds_miss <- SimIm(ds_true, p = pr_mv/100) # perform NA randomly in percent as amount NA in data
ds_mvi <- ds_miss
# now generate values at previous sections Univariate MVI / Multivariate MVI

NRMSE_hm <- Rmse(imp = ds_mvi_hf, mis = ds_miss, true = ds_true, norm = T)
NRMSE_knn <- Rmse(imp = ds_mvi_knn, mis = ds_miss, true = ds_true, norm = T)
NRMSE_cart <- Rmse(imp = ds_mvi_cart, mis = ds_miss, true = ds_true, norm = T)
NRMSE_rf1 <- Rmse(imp = ds_mvi_rf1, mis = ds_miss, true = ds_true, norm = T)
NRMSE_rf2 <- Rmse(imp = ds_mvi_rf2, mis = ds_miss, true = ds_true, norm = T)
NRMSE_rf3 <- Rmse(imp = ds_mvi_rf3, mis = ds_miss, true = ds_true, norm = T)
NRMSE_pls <- Rmse(imp = ds_mvi_pls, mis = ds_miss, true = ds_true, norm = T)
NRMSE_pca <- Rmse(imp = ds_mvi_pca, mis = ds_miss, true = ds_true, norm = T)
NRMSE_qrilc <- Rmse(imp = ds_mvi_qrilc, mis = ds_miss, true = ds_true, norm = T)

NRMSE_results <- data.frame(round(c(NRMSE_hm,NRMSE_knn, NRMSE_cart, NRMSE_rf1, NRMSE_rf2, NRMSE_rf3, NRMSE_pls, NRMSE_pca, NRMSE_qrilc), 2))
rownames(NRMSE_results) <- c("HM","KNN", "CART", "RF1", "RF2", "RF3", "PLS", "PCA", "QRILC")
colnames(NRMSE_results) <- "NRMSE"
best_method <- filter(NRMSE_results, NRMSE_results$NRMSE == min(NRMSE_results$NRMSE))
rownames(best_method)

########################################## EACH GROUPS ONE METHODS

# create Label (group) column
rname <- rownames(dsr)
rname <- str_remove(rname, ".CDF")
rname <- data.frame(rname)
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
un_l <- unique(n_gr_t)[[1]]

# create list with separated data frames with NA for each label (group)
dsr_Label <- cbind(n_gr_t, dsr)
colnames(dsr_Label)[1] <- "Label"
l_dsr_Label <- lapply(c(1:length(un_l)), function(y) subset(dsr_Label, Label == un_l[y]))
l_dsr <- lapply(c(1:length(un_l)), function(y) l_dsr_Label[[y]][,-1])

# calculate NA in % for each group
calc_NA <- function(x) {
 tn <- nrow(x)*ncol(x)
 mv_c <- sum(is.na(x)) 
 pr_mv <- round(mv_c/tn*100,0) 
 pr_mv }

groups_NA <- sapply(1:length(l_dsr), function(y) calc_NA(l_dsr[[y]]))

# create list with separated data frames without NA for each label (group)
ds_true_Label <- cbind(n_gr_t, ds_true)
colnames(ds_true_Label)[1] <- "Label"
l_ds_true_Label <- lapply(c(1:length(un_l)), function(y) subset(ds_true_Label, Label == un_l[y]))
l_ds_true <- lapply(c(1:length(un_l)), function(y) l_ds_true_Label[[y]][,-1])

# perform NA randomly in percent as amount NA in data
library(imputeR)
l_ds_miss <- lapply(1:length(l_ds_true), function(y) SimIm(l_ds_true[[y]], p = groups_NA[y]/100))
 
# copy and perform best methods funs
library(missForest)
library(doParallel)
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
l_ds_mvi_rf2 <- lapply(1:length(l_ds_miss), function(y) missForest(l_ds_miss[[y]], maxiter = 10, ntree = 100, mtry = floor(sqrt(ncol(l_ds_miss[[y]]))), parallelize = "forests")[[1]])
stopCluster(cl)

# proportion of each group
n_r <- lapply(1:length(l_ds_miss), function(y) nrow(l_ds_miss[[y]]))
s_n_r <- sum(unlist(n_r))
pr_gr <- round(unlist(n_r)/s_n_r,2)

# compute weighted mean NRMSE
library(imputeR)
NRMSE_gr <- lapply(1:length(l_ds_miss), function(y) Rmse(imp = l_ds_mvi_rf2[[y]], mis = l_ds_miss[[y]], true = l_ds_true[[y]], norm = T))
NRMSE_gr_w_mean <- round(weighted.mean(unlist(NRMSE_gr), pr_gr), 2)

##############################################################################################################################################################
# Perform best MVI
##############################################################################################################################################################

# copy and perform best methods funs
library(missForest)
library(doParallel)
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
ds_mvi_best <- missForest(ds_mvi, maxiter = 10, ntree = 100, mtry = floor(sqrt(ncol(ds_mvi))), parallelize = "forests")[[1]] 
stopCluster(cl)

colnames(ds_mvi_best) <- colnames(dsr) # assign initial colnames
identical(colnames(dsr), colnames(ds_mvi_best)) # check colnames

fwrite(ds_mvi_best, "xcms after IPO MVI.csv", row.names = T)

##############################################################################################################################################################
# References
##############################################################################################################################################################

# 1. Chaffin, Mark D., et al. "MetProc: separating measurement artifacts from true metabolites in an untargeted metabolomics experiment." Journal of proteome research 18.3 (2018): 1446-1450.
# 2. Marr, Sue, et al. "LC-MS based plant metabolic profiles of thirteen grassland species grown in diverse neighbourhoods." Scientific data 8.1 (2021): 1-12.
# 3. Wei, Runmin, et al. "Missing value imputation approach for mass spectrometry-based metabolomics data." Scientific reports 8.1 (2018): 1-10.
# 4. Di Guida, Riccardo, et al. "Non-targeted UHPLC-MS metabolomic data processing methods: a comparative investigation of normalisation, missing value imputation, transformation and scaling." Metabolomics 12.5 (2016): 93.
