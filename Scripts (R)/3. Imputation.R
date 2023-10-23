##############################################################################################################################################################
# Table of contents
##############################################################################################################################################################

# Installation
# Preprocessing
# Artifacts Removal
# Univariate MVI 
# Multivariate MVI
# Performance by different metrics
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
all_id <- lapply(1:length(rname), function(y) unlist(str_split(rname[y], " "))) # split info from rownames
ro_id <- as.numeric(unlist(lapply(all_id, function(y) unlist(y[1])))) # obtain run order ID (every [1] element)
dsr <- cbind(ro_id, dsr)
dsr$ro_id <- as.numeric(dsr$ro_id)
dsr <- dsr[order(dsr$ro_id, decreasing = F),] # order by ro col
dsr <- dsr[,-1]

# prepare for mvi and artifacts removal
ds_mvi <- dsr # remove all column with non-numeric values, use for mvi and artifacts removal
colnames(ds_mvi) <- paste0("V", c(1:ncol(ds_mvi))) # simplify colnames

##############################################################################################################################################################
# Preprocessing
##############################################################################################################################################################

########################################## amount of NA in %
tn <- nrow(ds_mvi)*ncol(ds_mvi)
mv_c <- sum(is.na(ds_mvi)) # or length(which(is.a(...)))
pr_mv <- round(mv_c/tn*100,0)
pr_mv

########################################## remove all column with NA
colna <- apply(ds_mvi,2,function(x) any(is.na(x)))
dsr_no_NA <- ds_mvi[,colna==F]

########################################## convert 0 into NA
ds_mvi[ds_mvi == 0] <- NA

########################################## check NA status
length(which(is.na(ds_mvi)))

########################################## filter by MVI % by column (feature)
sel_id <- grep(pattern = "QC", x = rownames(ds_mvi)) # select only biological or QC samples by pattern in row names
dsr_subset <- ds_mvi[sel_id,] # select only observations by pattern from peak table
dsr_subset <- ds_mvi # or select all observations
mv_f <- data.frame(apply(dsr_subset, 2, function(y) round(sum(is.na(y))/length(y)*100, 0)))
colnames(mv_f) <- "mv"
cutoff <- 50 # set cut-off for filtering
mv_f_n <- rownames(dplyr::filter(mv_f, mv < cutoff))
dsr_MVF <- ds_mvi[,mv_f_n] # filter
               
##############################################################################################################################################################
# Artifacts Removal
##############################################################################################################################################################

library(MetProc)
library(stringr)
library(dplyr)

# prepare data
ds <- as.data.frame(t(dsr))
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
art_name <- dplyr::filter(t_vn, MetProc %in% art)
retained <- colnames(dsr)[! colnames(dsr) %in% art_name$or_vn]
ds_mp <- dsr[,retained]
fwrite(ds_mp, "xcms after IPO remove artifacts.csv", row.names = T)

##############################################################################################################################################################
# Automatic determination type of NA and sequential imputation 
##############################################################################################################################################################

library(MAI)
library(parallel)
library(doParallel)

fc <- as.numeric(detectCores(logical = F))
cl <- makePSOCKcluster(fc+1)
registerDoParallel(cl)

set.seed(1234)

Results = MAI(data_miss = ds_mvi, # The data with missing values
              MCAR_algorithm = "random_forest", # The MCAR algorithm to use: "BPCA", "Multi_nsKNN", "random_forest"
              MNAR_algorithm = "Single", # The MNAR algorithm to use: "nsKNN", "Single"
              n_cores = -1, # To use all cores specify n_cores = -1, adjust to your data
              verbose = T) # allows console message output

stopCluster(cl)
stopImplicitCluster()
registerDoSEQ()

ds_mai <- as.data.frame(Results$Imputed_data) # as data frame
colnames(ds_mai) <- colnames(dsr) # assign initial colnames
identical(colnames(dsr), colnames(ds_mai)) # check colnames
rownames(ds_mai) <- rownames(dsr) # assign initial rownames
identical(rownames(dsr), rownames(ds_mai)) # check rownames

fwrite(ds_mai, "xcms after IPO MAI.csv", row.names = T)

# Alternatively you can use "model.Selector" & "impute.MAR.MNAR" functions from imputeLCMD package               
               
##############################################################################################################################################################
# Univariate MVI 
##############################################################################################################################################################

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
imp.rand <- abs(rnorm(sum(NAidx), mean = noise, sd = sd)) # generate random values. other option: runif(sum(NAidx), noise-sd, noise+sd)
ds_mvi[NAidx] <- imp.rand

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

########################################## left-censored data MVI
# alternatively you can use "impute.MinDet" & "impute.MinProb" functions from imputeLCMD package for for left-censored values or QRILC MVI (see below)
                         
########################################## apply version example for half minimum
HM_apply <- function(data) {
  result <- as.data.frame(sapply(1:ncol(data), function(y) replace(data[,y],which(is.na(data[,y])), min(data[,y], na.rm = T)/2)))
  result
}

ds_mvi_hm <- HM_apply(ds_mvi)
colnames(ds_mvi_hm) <- colnames(ds_mvi)
rownames(ds_mvi_hm) <- rownames(ds_mvi)

##############################################################################################################################################################
# Multivariate MVI
##############################################################################################################################################################

########################################## KNN MVI
library(impute)  
ds_mvi_knn <- as.data.frame(t(impute.knn(t(ds_mvi), k = 3, rowmax = 90, colmax = 90)[[1]]))

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
ds_mvi_rf2 <- as.data.frame(missForest(ds_mvi, maxiter = 10, ntree = 100, mtry = floor(sqrt(ncol(ds_mvi))), parallelize = "forests")[[1]])
stopCluster(cl)

# example StatTools package
library(StatTools)
library(doParallel)
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
ds_mvi_rf3 <- as.data.frame(rfImpWrap(MAT = ds_mvi, nCore = nCore))
stopCluster(cl)

########################################## PLS MVI
library(StatTools)
library(doParallel)
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
ds_mvi_pls <- as.data.frame(mvImpWrap(ds_mvi, method = "PLS", nCore = nCore))
stopCluster(cl)

########################################## PCA MVI
library(pcaMethods)
pc <- pca(ds_mvi, nPcs=3, method="ppca") # svd, ppca, bpca, nipals, etc, see listPcaMethods(), set "nPcs"
ds_mvi_pca <- as.data.frame(completeObs(pc))

########################################## QRILC MVI
library(imputeLCMD)
ds_mvi_qrilc <- as.data.frame(t(impute.QRILC(t(ds_mvi))[[1]]))

########################################## tWLSA MVI
library(tWLSA)
ds_mvi_twlsa <- as.data.frame(wlsMisImp(ds_mvi,lamda=0.1,conRate=99,group=1, outlier=T)) # adjust "group" and other to your data

##############################################################################################################################################################
# Performance by different metrics
##############################################################################################################################################################

##########################################
########################################## Preparation
##########################################

# amount of NA in %
tn <- nrow(ds_mvi)*ncol(ds_mvi)
mv_c <- sum(is.na(ds_mvi)) # or length(which(is.a(...)))
pr_mv <- round(mv_c/tn*100,0)

# remove all column with NA
colna <- apply(ds_mvi,2,function(x) any(is.na(x)))
dsr_no_NA <- ds_mvi[,colna==F]

# generate datasets for testing MVI
library(imputeR)
library(dplyr)
ds_true <- dsr_no_NA
ds_miss <- SimIm(ds_true, p = pr_mv/100) # perform NA randomly in percent as amount NA in data. Also use missMethods package (delete_... functions) or imputeLCMD (insertMVs function)
ds_mvi <- ds_miss # use it for testing mvi

##########################################
########################################## NRMSE
##########################################

library(imputeR)

# generate ds_mvi_... data sets at previous sections Univariate MVI / Multivariate MVI
NRMSE_hm <- Rmse(imp = ds_mvi_hm, mis = ds_miss, true = ds_true, norm = T)
NRMSE_knn <- Rmse(imp = ds_mvi_knn, mis = ds_miss, true = ds_true, norm = T)
NRMSE_cart <- Rmse(imp = ds_mvi_cart, mis = ds_miss, true = ds_true, norm = T)
NRMSE_rf1 <- Rmse(imp = ds_mvi_rf1, mis = ds_miss, true = ds_true, norm = T)
NRMSE_rf2 <- Rmse(imp = ds_mvi_rf2, mis = ds_miss, true = ds_true, norm = T)
NRMSE_rf3 <- Rmse(imp = ds_mvi_rf3, mis = ds_miss, true = ds_true, norm = T)
NRMSE_pls <- Rmse(imp = ds_mvi_pls, mis = ds_miss, true = ds_true, norm = T)
NRMSE_pca <- Rmse(imp = ds_mvi_pca, mis = ds_miss, true = ds_true, norm = T)
NRMSE_qrilc <- Rmse(imp = ds_mvi_qrilc, mis = ds_miss, true = ds_true, norm = T)
NRMSE_twlsa <- Rmse(imp = ds_mvi_twlsa, mis = ds_miss, true = ds_true, norm = T)

NRMSE_results <- data.frame(round(c(NRMSE_hm,NRMSE_knn, NRMSE_cart, NRMSE_rf1, NRMSE_rf2, NRMSE_rf3, NRMSE_pls, NRMSE_pca, NRMSE_qrilc, NRMSE_twlsa), 2))
rownames(NRMSE_results) <- c("HM","KNN", "CART", "RF1", "RF2", "RF3", "PLS", "PCA", "QRILC", "tWLSA")
colnames(NRMSE_results) <- "NRMSE"
NRMSE_results

##########################################
########################################## Other metrics
##########################################

library(missMethods)

# generate ds_mvi_... data sets at previous sections Univariate MVI / Multivariate MVI (see argument "criterion")
VALUE_hm <- evaluate_imputed_values(ds_imp = as.data.frame(ds_mvi_hm), ds_orig = ds_true, criterion = "cor")
VALUE_knn <- evaluate_imputed_values(ds_imp = as.data.frame(ds_mvi_knn), ds_orig = ds_true, criterion = "cor")
VALUE_cart <- evaluate_imputed_values(ds_imp = as.data.frame(ds_mvi_cart), ds_orig = ds_true, criterion = "cor")
VALUE_rf1 <- evaluate_imputed_values(ds_imp = as.data.frame(ds_mvi_rf1), ds_orig = ds_true, criterion = "cor")
VALUE_rf2 <- evaluate_imputed_values(ds_imp = as.data.frame(ds_mvi_rf2), ds_orig = ds_true, criterion = "cor")
VALUE_rf3 <- evaluate_imputed_values(ds_imp = as.data.frame(ds_mvi_rf3), ds_orig = ds_true, criterion = "cor")
VALUE_pls <- evaluate_imputed_values(ds_imp = as.data.frame(ds_mvi_pls), ds_orig = ds_true, criterion = "cor")
VALUE_pca <- evaluate_imputed_values(ds_imp = as.data.frame(ds_mvi_pca), ds_orig = ds_true, criterion = "cor")
VALUE_qrilc <- evaluate_imputed_values(ds_imp = as.data.frame(ds_mvi_qrilc), ds_orig = ds_true, criterion = "cor")
VALUE_twlsa <- evaluate_imputed_values(ds_imp = as.data.frame(ds_mvi_twlsa), ds_orig = ds_true, criterion = "cor")

VALUE_results <- data.frame(round(c(VALUE_hm,VALUE_knn, VALUE_cart, VALUE_rf1, VALUE_rf2, VALUE_rf3, VALUE_pls, VALUE_pca, VALUE_qrilc, VALUE_twlsa), 2))
rownames(VALUE_results) <- c("HM","KNN", "CART", "RF1", "RF2", "RF3", "PLS", "PCA", "QRILC", "tWLSA")
colnames(VALUE_results) <- "VALUE"
VALUE_results

##########################################
########################################## Procrustes analysis on PCs scores
##########################################

library(vegan)

nPCs <- 2 # number of PCs

# generate ds_mvi_... data sets at previous sections Univariate MVI / Multivariate MVI (see argument "criterion")
PA_hm <- procrustes(prcomp(ds_true, scale. = T, center = T)$x[, 1:nPCs], prcomp(ds_mvi_hm, scale. = T, center = T)$x[, 1:nPCs], symmetric = T)$ss
PA_knn <- procrustes(prcomp(ds_true, scale. = T, center = T)$x[, 1:nPCs], prcomp(ds_mvi_knn, scale. = T, center = T)$x[, 1:nPCs], symmetric = T)$ss
PA_cart <- procrustes(prcomp(ds_true, scale. = T, center = T)$x[, 1:nPCs], prcomp(ds_mvi_cart, scale. = T, center = T)$x[, 1:nPCs], symmetric = T)$ss
PA_rf1 <- procrustes(prcomp(ds_true, scale. = T, center = T)$x[, 1:nPCs], prcomp(ds_mvi_rf1, scale. = T, center = T)$x[, 1:nPCs], symmetric = T)$ss
PA_rf2 <- procrustes(prcomp(ds_true, scale. = T, center = T)$x[, 1:nPCs], prcomp(ds_mvi_rf2, scale. = T, center = T)$x[, 1:nPCs], symmetric = T)$ss
PA_rf3 <- procrustes(prcomp(ds_true, scale. = T, center = T)$x[, 1:nPCs], prcomp(ds_mvi_rf3, scale. = T, center = T)$x[, 1:nPCs], symmetric = T)$ss
PA_pls <- procrustes(prcomp(ds_true, scale. = T, center = T)$x[, 1:nPCs], prcomp(ds_mvi_pls, scale. = T, center = T)$x[, 1:nPCs], symmetric = T)$ss
PA_pca <- procrustes(prcomp(ds_true, scale. = T, center = T)$x[, 1:nPCs], prcomp(ds_mvi_pca, scale. = T, center = T)$x[, 1:nPCs], symmetric = T)$ss
PA_qrilc <- procrustes(prcomp(ds_true, scale. = T, center = T)$x[, 1:nPCs], prcomp(ds_mvi_qrilc, scale. = T, center = T)$x[, 1:nPCs], symmetric = T)$ss
PA_twlsa <- procrustes(prcomp(ds_true, scale. = T, center = T)$x[, 1:nPCs], prcomp(ds_mvi_twlsa, scale. = T, center = T)$x[, 1:nPCs], symmetric = T)$ss

PA_results <- data.frame(round(c(PA_hm,PA_knn, PA_cart, PA_rf1, PA_rf2, PA_rf3, PA_pls, PA_pca, PA_qrilc, PA_twlsa), 2))
rownames(PA_results) <- c("HM","KNN", "CART", "RF1", "RF2", "RF3", "PLS", "PCA", "QRILC", "tWLSA")
colnames(PA_results) <- "PA_ss"
PA_results

##############################################################################################################################################################
# Perform best MVI
##############################################################################################################################################################

# copy and perform best methods funs with ds_mvi from 1st section "Installation"
library(missForest)
library(doParallel)
nCore=detectCores()-1
cl=makeCluster(nCore)
registerDoParallel(cl)
ds_mvi_best <- missForest(ds_mvi, maxiter = 10, ntree = 100, mtry = floor(sqrt(ncol(ds_mvi))), parallelize = "forests")[[1]] 
stopCluster(cl)

ds_mvi_best <- as.data.frame(ds_mvi_best) # as data frame
colnames(ds_mvi_best) <- colnames(dsr) # assign initial colnames
identical(colnames(dsr), colnames(ds_mvi_best)) # check colnames
rownames(ds_mvi_best) <- rownames(dsr) # assign initial rownames
identical(rownames(dsr), rownames(ds_mvi_best)) # check rownames

fwrite(ds_mvi_best, "xcms after IPO MVI.csv", row.names = T)

##############################################################################################################################################################
# References
##############################################################################################################################################################

# 1. Chaffin, Mark D., et al. "MetProc: separating measurement artifacts from true metabolites in an untargeted metabolomics experiment." Journal of proteome research 18.3 (2018): 1446-1450.
# 2. Dekermanjian, Jonathan P., et al. "Mechanism-aware imputation: a two-step approach in handling missing values in metabolomics." BMC Bioinformatics (2022): 23, 179.
# 3. Marr, Sue, et al. "LC-MS based plant metabolic profiles of thirteen grassland species grown in diverse neighbourhoods." Scientific data 8.1 (2021): 1-12.
# 4. Wei, Runmin, et al. "Missing value imputation approach for mass spectrometry-based metabolomics data." Scientific reports 8.1 (2018): 1-10.
# 5. Di Guida, Riccardo, et al. "Non-targeted UHPLC-MS metabolomic data processing methods: a comparative investigation of normalisation, missing value imputation, transformation and scaling." Metabolomics 12.5 (2016): 93.
# 6. Kumar, Nishith, Md Hoque, and Masahiro Sugimoto. "Kernel weighted least square approach for imputing missing values of metabolomics data." Scientific reports 11.1 (2021): 1-12.
