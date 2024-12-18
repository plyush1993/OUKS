##############################################################################################################################################################
# Table of contents
##############################################################################################################################################################

# Installation
# Outlier detection
# Multiple statistical filtration
# Classification Machine Learning task
# Regression Machine Learning task
# Testing set of biomarkers
# MWAS/ANCOVA
# Signal Modeling
# N-Factor Analysis
# Repeated Measures
# Time series
# Unsupervised Data Projection
# Correlation Analysis
# Distance Analysis
# Sample Size and Power Calculation
# References

##############################################################################################################################################################
# Installation
##############################################################################################################################################################

# setup environment
library(data.table)
setwd("D:/...")

# PEAK TABLE
# dataset with intensities, annotation and label column
ds <- as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats annot+filtr LMM adj KEGG.csv", header=T))
ds <- ds[-c(1:12),] # adjust to your data
rownames(ds) <- ds[,5] # adjust to your data
ds <- ds[,-c(1,3:5)] # adjust to your data
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# METADATA
meta <- as.data.frame(fread(input = "metadata.csv", header=T))
rownames(meta) <- meta[,5] # adjust to your data
order <- meta[,2] # adjust to your data
meta <- meta[,c(3,8:10)] # adjust to your data
colnames(meta) <- c("Batch", "Creatinine", "Age", "Sex")
# check identical rownames between peak table and metadata
identical(rownames(ds), rownames(meta))

# NEW WD FOR STATISTICAL ANALYSIS
setwd("D:/...")

##############################################################################################################################################################
# Outlier detection
##############################################################################################################################################################

###############################################
############################################### OutlierDetection package
############################################### 

# remotes::install_github("rubak/spatstat.revdep", subdir = "OutlierDetection")
library(OutlierDetection) # require from "rubak/spatstat.revdep"
library(pcaMethods)

# perform
pcod <- pca(ds[,-1]) # compute pca
co_od <- 0.95 # set value
maha <- maha(pcod@scores,cutoff=co_od) # Mahalanobis distance outlier detection
maha
nn <- nn(pcod@scores,k=(1-co_od)*nrow(ds)) # Euclidean distance outlier detection
nn
nnk <- nnk(pcod@scores,k=(1-co_od)*nrow(ds)) # Euclidean distance outlier detection
nnk

###############################################
############################################### ClassDiscovery package
############################################### 

library(ClassDiscovery)

thr_out <- 0.05 # set value
spca <- SamplePCA(t(ds[,-1])) # compute pca
mc <- mahalanobisQC(spca, 2)
mc[mc$p.value < thr_out,] # Mahalanobis distance outlier detection

plot(spca@scores)
points(spca@scores[which(mc$p.value < thr_out),], col = "red")

###############################################
############################################### pcaMethods package
############################################### 

library(pcaMethods)

pca_res <- pca(ds[,-1])
dmodx <- DModX(pca_res)
plot(dmodx) # overview
thr_out <- 1.2 # set value
plot(dmodx)
points(dmodx[dmodx > thr_out], col = "red")
dmodx[dmodx > thr_out]
with(pca_res, plot(DModX(pca_res)~ds$Label))

###############################################
############################################### HotellingEllipse package
############################################### 

library(HotellingEllipse)
library(FactoMineR)
library(ggplot2)
library(ggforce)
library(dplyr)
library(tidyverse)

n = 2 # number of PCs, adjust to your data
pca_res <- FactoMineR::PCA(ds[,-1], scale.unit = T, graph = F, ncp = n)
res_2PCs <- ellipseParam(data = as.data.frame(pca_res$ind$coord), k = 2, pcx = 1, pcy = 2) # adjust to your data

a1 <- pluck(res_2PCs, "Ellipse", "a.99pct")
b1 <- pluck(res_2PCs, "Ellipse", "b.99pct")
a2 <- pluck(res_2PCs, "Ellipse", "a.95pct")
b2 <- pluck(res_2PCs, "Ellipse", "b.95pct") 
T2 <- pluck(res_2PCs, "Tsquare", "value")

as.data.frame(pca_res$ind$coord) %>%
  ggplot(aes(x = Dim.1, y = Dim.2)) +
  geom_ellipse(aes(x0 = 0, y0 = 0, a = a1, b = b1, angle = 0), size = .8, linetype = "dotted", fill = "white", color = "darkred") +
  geom_ellipse(aes(x0 = 0, y0 = 0, a = a2, b = b2, angle = 0), size = .8, linetype = "dashed", fill = "white", color = "darkblue") +
  geom_point(aes(fill = T2), shape = 21, size = 3, color = "black") +
  scale_fill_viridis_c(option = "viridis") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = .2) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = .2) +
  labs(title = "Scatterplot of PCA scores", subtitle = "PC1 vs. PC2", x = "PC1", y = "PC2", fill = "T2", caption = "Figure: Hotelling's T2 ellipse obtained\n using the ellipseParam function") +
  theme_grey()

pc2 <- as.data.frame(pca_res$ind$coord)
pc2 <- pc2[order(pc2$Dim.1, decreasing = T),] # set Dim.1 or Dim.2
as.data.frame(res_2PCs$Ellipse)
pc2

n = 5 # number of PCs, adjust to your data
pca_res <- FactoMineR::PCA(ds[,-1], scale.unit = T, graph = F, ncp = n)
res_nPCs <- ellipseParam(data = as.data.frame(pca_res$ind$coord), k = n)
tibble(
  T2 = pluck(res_nPCs, "Tsquare", "value"), 
  obs = 1:nrow(as.data.frame(pca_res$ind$coord))) %>%
  ggplot() +
  geom_point(aes(x = obs, y = T2, fill = T2), shape = 21, size = 3, color = "black") +
  geom_segment(aes(x = obs, y = T2, xend = obs, yend = 0), size = .5) +
  scale_fill_gradient(low = "black", high = "red", guide = "none") +
  geom_hline(yintercept = pluck(res_nPCs, "cutoff.99pct"), linetype = "dashed", color = "darkred", size = .5) +
  geom_hline(yintercept = pluck(res_nPCs, "cutoff.95pct"), linetype = "dashed", color = "darkblue", size = .5) +
  annotate("text", x = 80, y = 13, label = "99% limit", color = "darkred") +
  annotate("text", x = 80, y = 9, label = "95% limit", color = "darkblue") +
  labs(x = "Observations", y = "Hotelling's T-squared (n PCs)", fill = "T2 stats", caption = "Figure: Hotelling's T-squared vs. Observations") +
  theme_bw()

thr_out <-  res_nPCs$cutoff.95pct # set value res_nPCs$cutoff.95pct or res_nPCs$cutoff.99pct
hot_res <- as.data.frame(res_nPCs$Tsquare)
hot_res <- hot_res$value
names(hot_res) <- rownames(ds)
hot_res[hot_res > thr_out]

##############################################################################################################################################################
# Multiple statistical filtration
##############################################################################################################################################################

# Content:
# Machine Learning + Variable Importance + Recursive Feature Selection
# Consensus Nested Cross-Validation for Feature Selection
# VIP from PLS models
# RF permutation
# Penalized/Stepwise Regression
# AUROC analysis
# Univariate Filtering
# t-Test / Kruskal-Wallis test
# Fold Change Calculation
# Moderated t-test
# Linear Modeling 
# Linear Mixed-Effects Modeling
# RUV2 Method
# Boruta Method
# TDFDR Method
# Correlation Analysis
# Correlation with smth
# Combine results
# Combine with annotation

###############################
############ NOTE! ############
###############################

# It is critically important to test feature selection (filtering) on independent dataset!
# For this purpose use "Data sampling" from "Classification Machine Learning task" section or "Consensus Nested Cross-Validation for Feature Selection" from current section.

###############################################
############################################### MULTIPLE ML WITH SFE AND RFS
############################################### 

library(parallel)
library(doParallel)
library(caret)
library(tuple)

# stop parallel
stopCluster(cl)
stopImplicitCluster()

# start parallel processing
fc <- as.numeric(detectCores(logical = T))
cl <- makePSOCKcluster(fc-1)
registerDoParallel(cl)

# repeated cross validation
trainControl <- trainControl(method="repeatedcv", number=10, repeats=10) # or bootstrap: trainControl(method="boot", number=100)
metric <- "Accuracy"

# PAM
set.seed(1234)
# library(pamr)
fit.pam <- train(Label~ ., data=ds, method="pam", metric=metric, trControl=trainControl, tuneLength = 10) # adjust tuneLength to your data

# SVM
set.seed(1234)
# library(kernlab)
fit.svm <- train(Label~ ., data=ds, method="svmRadial", metric=metric, trControl=trainControl, tuneLength = 10) # adjust tuneLength to your data

# PLS
set.seed(1234)
# library(pls)
fit.pls <- train(Label~ ., data=ds, method="pls", metric=metric, trControl=trainControl, tuneLength = 10) # adjust tuneLength to your data

# RF
set.seed(1234)
# library(randomForest)
fit.rf <- train(Label~ ., data=ds, method="rf", metric=metric, trControl=trainControl, tuneLength = 10) # adjust tuneLength to your data

# only accuracy for all models
results <- resamples(list(pam=fit.pam, svm=fit.svm, rf=fit.rf, pls=fit.pls), trControl = trainControl, metric=metric)
results_df <- as.data.frame(results)
results_df_fin <- apply(results_df[,-5], 2, mean)
results_df_fin

# by model specific value
set.seed(1234)
Imp.rf <- varImp(fit.rf, scale = FALSE)
Imp.rf <- Imp.rf$importance
rownames(Imp.rf) <- gsub("`", '', rownames(Imp.rf))

set.seed(1234)
Imp.pls <- varImp(fit.pls, scale = FALSE)
Imp.pls <- Imp.pls$importance
rownames(Imp.pls) <- gsub("`", '', rownames(Imp.pls))

set.seed(1234)
Imp.svm <- varImp(fit.svm, scale = FALSE)
Imp.svm <- Imp.svm$importance
rownames(Imp.svm) <- gsub("`", '', rownames(Imp.svm))

set.seed(1234)
Imp.pam <- varImp(fit.pam, scale = FALSE)
Imp.pam <- Imp.pam$importance
rownames(Imp.pam) <- gsub("`", '', rownames(Imp.pam))

# creating of list with top = n important features by model
n_model = 100
set.seed(1234)

Imp.rf[,c("sum")] <- apply(X = Imp.rf, 1, FUN = sum)
Imp.rf <- Imp.rf[order(Imp.rf$sum, decreasing = T),]
Imp.rf <- rownames(Imp.rf)[c(1:n_model)]

Imp.pls[,c("sum")] <- apply(X = Imp.pls, 1, FUN = sum)
Imp.pls <- Imp.pls[order(Imp.pls$sum, decreasing = T),]
Imp.pls <- rownames(Imp.pls)[c(1:n_model)]

Imp.svm[,c("sum")] <- apply(X = Imp.svm, 1, FUN = sum)
Imp.svm <- Imp.svm[order(Imp.svm$sum, decreasing = T),]
Imp.svm <- rownames(Imp.svm)[c(1:n_model)]

Imp.pam[,c("sum")] <- apply(X = Imp.pam, 1, FUN = sum)
Imp.pam <- Imp.pam[order(Imp.pam$sum, decreasing = T),]
Imp.pam <- rownames(Imp.pam)[c(1:n_model)]

# minimum n time of duplicated 
all <- c(Imp.rf, Imp.svm, Imp.pam, Imp.pls)
#all <- unique(c(Imp.rf, Imp.svm, Imp.pam, Imp.pls)) 
#all <- unique(c(Imp.svm,Imp.rf)) 
n_tuple <- 4
all1 <- all[which(tuplicated(all, n_tuple), T)]
all1 <- unique(all1)
ds_d <- cbind(ds[,1], ds[,all1])
#ds_d <- cbind(ds[,1], ds[,all]) 

# Recursive Feature Elimination
set.seed(1234)
subsets <- c(1:(ncol(ds_d)-1))
#ctrl_rfe <- rfeControl(functions = rfFuncs,method = "repeatedcv",number = 10, repeats = 10,verbose = FALSE)
ctrl_rfe <- rfeControl(functions = nbFuncs,method = "repeatedcv",number = 10, repeats = 10,verbose = FALSE)
# library(klaR) #(Version 0.6-14)
rfe <- rfe(ds_d[,-1], ds_d[,1], sizes = subsets, metric = "Accuracy", rfeControl = ctrl_rfe) # for RFE klaR package shuld be exactly 0.6-14 version
rfe_vi <- rfe[["optVariables"]]
ds_rfe <- data.frame(cbind(as.character(ds[,1]), ds[, rfe_vi]))

# features
sfs <- colnames(ds_d)[-1] # all1 or colnames(ds_d)[-1]
rfe <- rfe_vi # rfe_vi or colnames(ds_rfe)[-1]

# save outputs
# save(fit.pls, file = "pls.RData")
# save(fit.pam, file = "pam.RData")
# save(fit.rf, file = "rf.RData")
# save(fit.svm, file = "svm.RData")

# fwrite(ds_d, "ds_d tuple 2 top 5.csv", row.names = T)
# fwrite(ds_rfe, "ds_rfe nb.csv", row.names = T)
# fwrite(ds_rfe, "ds_rfe rf.csv", row.names = T)

###############################################
############################################### CONSENSUS NESTED CROSS-VALIDATION FOR FEATURE SELECTION
############################################### 

library(limma)
library(dplyr)
library(ropls)
library(caret)
library(parallel)
library(doParallel)

# stop parallel
stopCluster(cl)
stopImplicitCluster()

# start parallel processing
fc <- as.numeric(detectCores(logical = T))
cl <- makePSOCKcluster(fc-1)
registerDoParallel(cl)

# wrapper function with feature selection methods (in this example: intersection of 2 algorithms (p-value in moderated t-test, Fold-Change, VIP in PLS))
perform_tests <- function(inner_indecies,dataframe, cutoff_tt = 0.05, cutoff_fc = 1, cutoff_vip = 1){
  ds <- dataframe[inner_indecies,]
  
  # Moderated t-test
  mdl_mtrx <- model.matrix(~Label, ds) # adjust to your data
  lmf <- lmFit(t(log2(ds[,-1])), method = "ls", design = mdl_mtrx, maxit = 1000) # "robust" or "ls"
  efit <- eBayes(lmf)
  tableTop <- topTable(efit, coef = 2, adjust = "BH", number = ncol(ds), sort.by = "none") # select method
  lim_pval <- rownames(dplyr::filter(tableTop, as.numeric(adj.P.Val) <= cutoff_tt))

  # Fold Change
  ds_log <- as.data.frame(log2(ds[,-1]))
  ds_log <- cbind(ds[,1], ds_log)
  colnames(ds_log)[1] <- "Label"
  mdl_mtrx_log <- model.matrix(~Label, ds_log) # adjust to your data
  lmf_log <- lmFit(t(ds_log[,-1]), method = "ls", design = mdl_mtrx_log, maxit = 1000) # "robust" or "ls"
  efit_log <- eBayes(lmf_log)
  tableTop_log <- topTable(efit_log, coef = 2, adjust = "BH", number = ncol(ds), sort.by = "none") # select method
  fc <- rownames(dplyr::filter(tableTop_log, abs(as.numeric(logFC)) >= cutoff_fc))

  # VIP value from PLS
  pls <- opls(ds[,-1], ds[,1], orthoI = 0, predI = NA, crossvalI = 10, permI = 100, fig.pdfC = "none", info.txtC = "none") # orthoI -> 0/NA PLS/OPLS, predI -> No of components, crossvalI -> cross-validation, permI -> No of permutations
  vip <- as.data.frame(getVipVn(pls))
  vip <- cbind(name = rownames(vip), VIP = vip)
  colnames(vip)[2] <- "VIP"
  vip <- vip[order(vip$VIP, decreasing = T),]
  vip_th <- subset(vip, vip$VIP > cutoff_vip)
  vip_pls<- rownames(vip_th)
  
  combine <- Reduce(intersect, list(vip_pls, fc, lim_pval))
  return(combine)
}

# Set parameters
cutoff_tt <- 0.05 # set p-value for moderated t-test
cutoff_fc <- 1 # set cutoff for Fold Change
cutoff_vip <- 1 # set cutoff for VIP in PLS
outer <- createFolds(ds$Label, 10) # Type of resampling in outer loop, adjust to your data
tc_out <- trainControl(method = "cv", number = 10, savePredictions = T, classProbs = T) # Type of resampling in outer loop for machine learning, adjust to your data
l_r <- list()
results <- data.frame(fold = numeric(), par = numeric(), acc_test = numeric()) 

# Perform
for (i in 1:length(outer)) { # !!! if some errors try different cutoffs and/or tuneLength 
  inner <- createFolds(ds[-outer[[i]],]$Label, 10) # Type of resampling in inner loop, adjust to your data 
  is <- lapply(1:length(inner), function(y) perform_tests(as.numeric(unlist(inner[-y])), ds[-outer[[i]],], cutoff_tt = cutoff_tt, cutoff_fc = cutoff_fc, cutoff_vip = cutoff_vip)) # adjust cutoffs to your data
  ifs <- Reduce(intersect, is) # table(unlist(is)) or Reduce(intersect, is)
  l_r[[i]] <- ifs

  # PLS (as example) see http://topepo.github.io/caret/train-models-by-tag.html for other regression algorithm
  m_outer <- train(x = ds[-outer[[i]], l_r[[i]]], y = ds[-outer[[i]],]$Label, method = "pls", trControl = tc_out, tuneLength = 5) # adjust tuneLength to your data
  acc <- confusionMatrix(predict(m_outer, ds[outer[[i]], l_r[[i]]]), ds[outer[[i]],]$Label) 
  results <- dplyr::add_row(results, fold = i, par = as.numeric(m_outer$bestTune), acc_test = acc$overall["Accuracy"]) # adjust ti your data
}

nfs <- Reduce(intersect, l_r) # Selected features (intersect from all inner and outer loops)
results # par -> best hyperparameter
results$acc_test 
mean(results$acc_test)
plot(results$acc_test)

###############################################
############################################### VIP FROM PLS MODELS
############################################### 

library(ropls)

pls <- opls(ds[,-1], ds[,1], orthoI = 0, predI = NA, crossvalI = 10, permI = 100) # orthoI -> 0/NA PLS/OPLS, predI -> No of components, crossvalI -> cross-validation, permI -> No of permutations
# plot(pls)
vip <- as.data.frame(getVipVn(pls))
vip <- cbind(name = rownames(vip), VIP = vip) 
colnames(vip)[2] <- "VIP"
vip <- vip[order(vip$VIP, decreasing = T),]
th_vip <- 1.0 # set value for filtration
vip_th <- subset(vip, vip$VIP > th_vip)
vip_pls <- rownames(vip_th) # features

###############################################
############################################### RF PERMUTATION
############################################### 

library(permimp)
library(party)

set.seed(1234)
mtry <- round(sqrt(ncol(ds[,-1])),0) # as.numeric(fit.rf[["bestTune"]]) or round(sqrt(ncol(ds[,-1])),0)
ntree <- 500 # fit.rf[["finalModel"]][["ntree"]] or 500
rf_model <- cforest(Label ~ ., data = ds, control = party::cforest_unbiased(mtry = mtry, ntree = ntree)) # fit RF with hyperparameters: mtry and ntree
rf_perm <-  permimp(rf_model, threshold = 0.95, nperm = 100, scaled = F) # perform permutation importance: nperm -> No of permutations, threshold -> threshold value
# plot(rf_perm, horizontal = T)
rf_perm_df <- as.data.frame(rf_perm$values)
rf_perm_df <- as.data.frame(cbind(rownames(rf_perm_df), RFPERM = rf_perm_df$`rf_perm$values`))
rownames(rf_perm_df) <- rf_perm_df$V1
rf_perm_df$RFPERM <- abs(as.numeric(rf_perm_df$RFPERM))
rf_perm_df <- rf_perm_df[order(rf_perm_df$RFPERM, decreasing = T),]
th_rf <- 0 # set value for filtration
rf_th <- subset(rf_perm_df, rf_perm_df$RFPERM > th_rf)
rf_perm_f<- rownames(rf_th) # features

############################################### 
############################################### PENALIZED/STEPWISE REGRESSION
###############################################

################################################### Penalized Regression

library(glmnet)

# predictor variables
x <- model.matrix(Label~., ds)[,-1]
# outcome (Label) to a numerical variable
y <- ifelse(ds$Label == "TG", 1, 0) # set class for example: target -> "TG"

################################# perform

# alpha: the elasticnet mixing parameter. Allowed values include:
# "1": for lasso regression
# "0": for ridge regression
# a value between 0 and 1 for elastic net regression
set.seed(1234) 
alpha <- 0.5 # set type of penalized regression
cv.penal <- cv.glmnet(x, y, alpha = alpha, family = "binomial", nfolds = 10, type.measure="auc") # adjust type.measure for your data
# plot(cv.penal)

################################# Final model with optimal lambda

# lambda = cv.penal$lambda.1se or lambda = cv.penal$lambda.min
# lambda = lambda.1se produces a simpler model compared to lambda.min, but the model might be a little bit less accurate than the one obtained with lambda.min
penal.model <- glmnet(x, y, alpha = alpha, family = "binomial", lambda = cv.penal$lambda.1se)
coef_pr <- as.matrix(coef(penal.model))
rownames(coef_pr) <- gsub("`", '', rownames(coef_pr))
coef_pr_df <- as.data.frame(cbind(rownames(coef_pr), coef_pr))
rownames(coef_pr_df)[-c(1:2)] <- rownames(coef_pr)[-c(1:2)]
penal_f <- subset(coef_pr_df, as.numeric(coef_pr_df$s0) != 0)
penal_f <- rownames(penal_f[!grepl("Intercept", penal_f$V1),]) # features

################################################### Logistic Stepwise

library(MASS)
library(dplyr)

log_step_mod <- glm(Label ~ ., data = ds, family = binomial) %>% stepAIC(trace = F, direction = "both") # set direction
coef_st_log <- as.data.frame(coef(log_step_mod))
rownames(coef_st_log) <- gsub("`", '', rownames(coef_st_log))
coef_st_log <- as.data.frame(cbind(rownames(coef_st_log), coef_st_log))
colnames(coef_st_log) <- c("V1", "V2")
step_f <- rownames(coef_st_log[!grepl("Intercept", coef_st_log$V1),]) # features

############################################### 
############################################### AUROC ANALYSIS
###############################################

library(caret)

set.seed(1234)
Imp.ROC <- filterVarImp(x = ds[,-1], y = ds[,1])
Imp.ROC[,c("sum")] <- apply(X = Imp.ROC, 1, FUN = sum)
Imp.ROC <- Imp.ROC[order(Imp.ROC$sum, decreasing = T),]
Imp.ROC$sum <- Imp.ROC$sum/length(unique(ds[,1])) # mean AUROC value
th_roc <- 0.80 # set value for filtration
Imp.ROC_sel <- subset(Imp.ROC, sum > th_roc)

# features
roc <- rownames(Imp.ROC_sel)

###############################################
############################################### UNIVARIATE FILTERING
############################################### 

#----------------------------------
# designed for 2 groups comparison
#----------------------------------

# if some error in Shapiro normality test:
# use shapiro.wilk.test function from cwhmisc instead shapiro.test from stats
# library(cwhmisc)
# norm.test <- apply(xx, 2, function(t) cwhmisc::shapiro.wilk.test(t)$p)

uvf <- function(x, p.val.sig = 0.05, p.adjust = "BH"){
  
  norm_homog_tests <- function(x) {
    xx <- x[,-1]
    # normality test
    norm.test <- as.numeric(apply(xx, 2, function(t) shapiro.test(t)$p.value))
    # if some error in Shapiro normality test:
    # use shapiro.wilk.test function from cwhmisc instead shapiro.test from stats
    # library(cwhmisc)
    # norm.test <- apply(xx, 2, function(t) cwhmisc::shapiro.wilk.test(t)$p)
    
    # homogeneity test
    homog.test <- as.numeric(apply(xx, 2, function(t) bartlett.test(t,g = x[,1])$p.value))
    return(as.data.frame(cbind(norm.test, homog.test)))}
  
  res_tests <- norm_homog_tests(x)
  
   wilcox_test <- function(x,y) {
    xx <- x[,-1]
    wx.t <- as.vector(which(y[,1] < 0.05))
    wilcox_test <- list()
    ifelse(identical(wx.t, integer(0)), return (wilcox_test <- 1), wx.t)
    wilcox_test <- lapply(as.data.frame(xx[,wx.t]), function(t) as.numeric(pairwise.wilcox.test(x = t, g =  x[,1], paired=F)$p.value))
    wilcox_test <- as.list(p.adjust(unlist(wilcox_test), method = p.adjust))    
    names(wilcox_test) <- (colnames(x)[-1])[wx.t]
    return(as.list(wilcox_test))}
  
  wx.t.res <- wilcox_test(x, res_tests)
  
  welch_test <- function(x,y) {
    xx <- x[,-1]
    wl.t <- as.vector(which(y[,1] > 0.05 & y[,2] < 0.05))
    welch_test <- list()
    ifelse(identical(wl.t, integer(0)), return (welch_test <- 1), wl.t)
    welch_test <- lapply(as.data.frame(xx[,wl.t]), function(t) as.numeric(pairwise.t.test(x = t, g = x[,1], pool.sd = F)$p.value))
    welch_test <- as.list(p.adjust(unlist(welch_test), method = p.adjust))
    names(welch_test) <- (colnames(x)[-1])[wl.t]
    return(as.list(welch_test))}
  
  wl.t.res <- welch_test(x, res_tests)
  
  student_test <- function(x,y) {
    xx <- x[,-1]
    st.t <- as.vector(which(y[,1] > 0.05 & y[,2] > 0.05))
    student_test <- list()
    ifelse(identical(st.t, integer(0)), return (student_test <- 1), st.t)
    student_test <- lapply(as.data.frame(xx[,st.t]), function(t) as.numeric(pairwise.t.test(x = t, g = x[,1], pool.sd = T)$p.value))
    student_test <- as.list(p.adjust(unlist(student_test), method = p.adjust))
    names(student_test) <- (colnames(x)[-1])[st.t]
    return(as.list(student_test))}
  
  st.t.res <- student_test(x, res_tests)
  
  filt_p_val <- function(x, y, z, w){
    
    #x = ds
    #y = wx.t.res
    #z = wl.t.res
    #w = st.t.res
    
    wx.t.n <- names(y)
    wx.t.res2 <-as.data.frame(y)
    colnames(wx.t.res2) <- wx.t.n
    wx.t.res2[is.na(wx.t.res2)] <- max(wx.t.res2, na.rm = T)
    
    wl.t.n <- names(z)
    wl.t.res2 <- as.data.frame(z)
    colnames(wl.t.res2) <- wl.t.n
    wl.t.res2[is.na(wl.t.res2)] <- max(wl.t.res2, na.rm = T)
    
    st.t.n <- names(w)
    st.t.res2 <- as.data.frame(w)
    colnames(st.t.res2) <- st.t.n
    st.t.res2[is.na(st.t.res2)] <- max(st.t.res2, na.rm = T)
    
    wxx <- apply(wx.t.res2, 2, min)
    wll <- apply(wl.t.res2, 2, min)
    stt <- apply(st.t.res2, 2, min)
    wxx2 <- as.data.frame(wxx)
    wll2 <- as.data.frame(wll)
    stt2 <- as.data.frame(stt)
    
    wxx3 <- rownames(wxx2)[which(wxx2 <= p.val.sig)]
    wll3 <- rownames(wll2)[which(wll2 <= p.val.sig)]
    stt3 <- rownames(stt2)[which(stt2 <= p.val.sig)]
    aff <- c(wxx3, wll3, stt3)
    
    ds_fil <- cbind(x[,1], x[, aff])
    return(ds_fil)
  }
  return(filt_p_val(x, wx.t.res, wl.t.res, st.t.res))
}

# 1 st argument -> dataset with 1st "Label" column, 2nd -> p-value, 3rd -> the method of adjustment for multiple comparisons.
ds_uvf <- uvf(ds, p.val.sig = 0.05, p.adjust = "BH") 

# features
uvf <- colnames(ds_uvf)[-1]

############################################### 
############################################### # t-TEST / Kruskal-Wallis TEST
###############################################

#----------------------------------
# designed for 2 groups comparison
#----------------------------------

################################################### t-test

# prepare data
data_l <- lapply(1:length(unique(ds$Label)), function(y) subset(ds, Label==unique(ds$Label)[y])[,-1])
data_l <- lapply(1:length(data_l), function(y) sapply(data_l[[y]], as.numeric))
# perform
res.t.test.pval <- as.numeric(sapply(1:ncol(data_l[[1]]), function(y) t.test(x = data_l[[1]][,y], y = data_l[[2]][,y])$p.value)) # t-test
p_adj <- as.numeric(p.adjust(res.t.test.pval, method = "BH")) # Adjust P-values for Multiple Comparisons
pval.ttest <- as.data.frame(cbind(name = colnames(ds[,-1]), pval = as.numeric(p_adj))) 
pval.ttest$pval <- as.numeric(pval.ttest$pval)
rownames(pval.ttest) <- pval.ttest$name
th_ttest <- 0.05 # set value for filtration
pval_ttest <- rownames(subset(pval.ttest, pval.ttest$pval <= th_ttest))

################################################### Kruskal-Wallis test

res.kw.test <- lapply(2:ncol(ds), function(y) kruskal.test(ds[,y] ~ ds[,1])) # Kruskal test # Try also wilcox.test or pairwise.wilcox.test 
res.kw.test.pval <- sapply(2:ncol(ds), function(y) kruskal.test(ds[,y] ~ ds[,1])$p.value)
p_adj <-p.adjust(res.kw.test.pval, method = "BH") # Adjust P-values for Multiple Comparisons
pval.kw <- as.data.frame(cbind(name = colnames(ds[,-1]), pval = as.numeric(p_adj))) 
pval.kw$pval <- as.numeric(pval.kw$pval)
rownames(pval.kw) <- pval.kw$name
th_kw <- 0.05 # set value for filtration
pval_kw <- rownames(subset(pval.kw, pval.kw$pval <= th_kw))

############################################### 
############################################### FOLD CHANGE CALCULATION
###############################################

library(dplyr)

FOLD.CHANGE <- function(data) {
  ds_subsets <- lapply(1:length(unique(data[,1])), function(y) dplyr::filter(data[,-1], data$Label == unique(data[,1])[y])) # list of subsets by label
  mean_r_l <- lapply(1:length(ds_subsets), function(y) apply(ds_subsets[[y]], 2, mean, na.rm = T)) # calculate mean for feature
  foldchange <- log2((mean_r_l[[1]] / mean_r_l[[2]]))
  fc_res <- as.data.frame(foldchange)
  return(fc_res)
}

fc_t <- 1.0 # set threshold
fc_res <- FOLD.CHANGE(ds)
fc_res$foldchange <- as.numeric(fc_res$foldchange)
fc <- rownames(subset(fc_res, foldchange > fc_t | foldchange < -fc_t))

################################################### Multigroup Fold Change

FOLD.CHANGE.MG <- function(x, f, aggr_FUN = colMeans, combi_FUN = {function(x,y) "-"(x,y)}){
  f <- as.factor(f)
  i <- split(1:nrow(x), f)
  x <- sapply(i, function(i){ aggr_FUN(x[i,])})
  x <- t(x)
  x <- log2(x)
  j <- combn(levels(f), 2)
  ret <- combi_FUN(x[j[1,],], x[j[2,],])
  rownames(ret) <- paste(j[1,], j[2,], sep = '/')
  t(ret)
}

fdr <- FOLD.CHANGE.MG(ds[,-1], ds[,1])
fdt_tr <- 1.0 # set threshold
f <- as.data.frame(apply(fdr, 1, function(x) x> fdt_tr | x< -fdt_tr))
f1 <- as.data.frame(which(f == T, arr.ind = T))
colnames(f) <- colnames(ds[,-1])
fc <- colnames(f)[unique(f1$col)]

# by mean
fdr_mean <- apply(abs(fdr),1, mean, na.rm=T)
fdt_tr <- 1.0 # set threshold
fc <- names(which(fdr_mean > fdt_tr)) 

############################################### 
############################################### MODERATED T-TEST
###############################################

library(limma)
library(dplyr)

mdl_mtrx <- model.matrix(~Label, ds) # adjust to your data
lmf <- lmFit(t(log2(ds[,-1])), method = "robust", design = mdl_mtrx, maxit = 1000) # "robust" or "ls"
efit <- eBayes(lmf)
tableTop <- topTable(efit, coef = 2, adjust = "BH", number = ncol(ds), sort.by = "none") # select method
cutoff <- 0.05 # set cutoff
lim_pval <- rownames(dplyr::filter(tableTop, as.numeric(adj.P.Val) <= cutoff)) # features

############################################### 
############################################### LM MODELING OF BIOLOGICAL FACTORS
###############################################

library(MetabolomicsBasics)

# data
dat <- ds[,-1]
s_d <- cbind(Class = ds$Label, Order = order, meta)

# perform by package
model <- MetaboliteANOVA(dat=dat, sam=s_d, model="Batch+Class+Age+Sex", method = "BH") # "Batch+Class+Age+Sex" or "Batch+Class+Age+Sex+Order" # select method # adjust to your data

# features
# lm_c <- names(which(model[,"Batch"]>0.05 & model[,"Age"]>0.05 & model[,"Sex"]>0.05 & model[,"Class"]<0.05))
lm <- rownames(model[which(model[,"Class"]<0.05),]) # or lm_c

# perform manually
model="Batch+Class+Age+Sex" # "Batch+Class+Age+Sex" or "Batch+Class+Age+Sex+Order" # adjust to your data # select method 
mod_f <- as.formula(paste("y", model, sep = " ~ "))
lm_fit <- lapply(1:ncol(dat), function(y) lm(mod_f, data = cbind(s_d, y = dat[,y])))
lm_pval <- as.data.frame(t(sapply(1:ncol(dat), function(y) anova(lm_fit[[y]])[, "Pr(>F)"])))
lm_pval <- as.data.frame(lm_pval[,-ncol(lm_pval)])
colnames(lm_pval) <- rownames(anova(lm_fit[[1]]))[1:ncol(lm_pval)]
rownames(lm_pval) <- colnames(dat)

# features
# lm_c <- names(which(lm_pval[,"Batch"]>0.05 & lm_pval[,"Age"]>0.05 & lm_pval[,"Sex"]>0.05 & lm_pval[,"Class"]<0.05))
lm <- rownames(lm_pval[which(lm_pval[,"Class"]<0.05),]) # or lm_c

############################################### 
############################################### LMM MODELING OF BIOLOGICAL FACTORS
###############################################

library(lme4)
library(lmerTest)

# data
dat <- cbind(meta, ds)
n_meta <- 5 # adjust to your data
colnames(dat)[-c(1:n_meta)] <- paste0("X", c(1:ncol(dat))) # adjust to your data
dat$Batch <- as.factor(dat$Batch)
dat$Age <- as.integer(dat$Age)
dat$Sex <- as.integer(dat$Sex)
dat$Label <- as.factor(dat$Label)

# perform
n_start <- 6 # adjust to your data
lmm_fit <- lapply(n_start:ncol(dat), function(x) lmer(dat[,x] ~ Label + Sex + Age + (1|Batch), dat)) # adjust to your data # lmer from lmerTest package
lmm_fit_coef <- lapply(1:length(lmm_fit), function(x) summary(lmm_fit[[x]])$coefficients)
lmm_fit_pval <- sapply(1:length(lmm_fit_coef), function(x) lmm_fit_coef[[x]][,5][2]) # adjust to your data
lmm_fit_pval_all_df <- as.data.frame(t(sapply(1:length(lmm_fit_coef), function(x) lmm_fit_coef[[x]][,5]))) # select method
lmm_fit_pval_all_df <- as.data.frame(sapply(1:ncol(lmm_fit_pval_all_df), function(x) p.adjust(lmm_fit_pval_all_df[,x], method = "BH"))) # adjust to your data # select method
dat2 <- cbind(meta, ds)
rownames(lmm_fit_pval_all_df) <- colnames(dat2)[-c(1:n_meta)]
colnames(lmm_fit_pval_all_df) <- rownames(lmm_fit_coef[[1]]) # adjust to your data

# features
class_un <- sum(lmm_fit_pval_all_df[,"Age"]>0.05 & lmm_fit_pval_all_df[,"Sex"]>0.05 & lmm_fit_pval_all_df[,"Class"]<0.05)
class_un
class_er <- sum(lmm_fit_pval_all_df[,"Age"]>0.05 & lmm_fit_pval_all_df[,"Sex"]>0.05)
class_er
class_er2 <- sum(lmm_fit_pval_all_df[,"Age"]<0.05 & lmm_fit_pval_all_df[,"Sex"]<0.05)
class_er2
class <- length(rownames(lmm_fit_pval_all_df[which(lmm_fit_pval_all_df[,"Class"]<0.05),]))
class
# lmm_ind_c <- which(lmm_fit_pval_all_df[,"Age"]>0.05 & lmm_fit_pval_all_df[,"Sex"]>0.05 & lmm_fit_pval_all_df[,"Class"]<0.05)+n_meta
lmm_ind <- which(lmm_fit_pval <= 0.05)+n_meta
lmm <- colnames(dat2)[lmm_ind] # or lmm_ind_c

############################################### 
############################################### RUV2 METHOD
###############################################

library(NormalizeMets)
library(stringr)

# generate featuredata/sampledata for RUVs
featuredata <- as.data.frame(ds[,-1])
colnames(featuredata) <- colnames(ds[,-1])
sampledata <- cbind(Class = ds$Label, Order = order, meta)

# generate data of IS
# search m/z of IS
cn <- colnames(featuredata)
is_ind <- which(str_detect(string = cn, pattern = "340.15")==T) # write in pattern m/z value of IS. If many ISs -> repeat this and save in 1 vector

smart.corr.test <- function(x, is_idx){
  b <- vector(mode="numeric")
  res <- apply(x, 2, shapiro.test)
  for (i in 1:ncol(x)) b[i] = res[[i]]$p.value
  names(b) <- colnames(x)
  for (i in 1:ncol(x))
    if (b[i] < 0.05) {
      b[i] <- (cor.test(x = x[, i], y = x[, is_idx], method = "spearman")$estimate)
    } else{
      b[i] <- (cor.test(x = x[, i], y = x[, is_idx], method = "pearson")$estimate)
    }
  return(b) }

corr <- smart.corr.test(featuredata, is_idx = is_ind) # detect other signals with correlation > corr_tr by IS (QC signals)
corr_tr <- 0.95 
corr <- abs(corr)
iss <- names(corr[corr > corr_tr]) # names(corr[corr > corr_tr]) or cn[is_ind] if multiple ISs

qcmets <- sapply(1:length(iss), function(y) which(colnames(featuredata) == iss[y])) # which(colnames(featuredata) == iss) or sapply(1:length(iss), function(y) which(colnames(featuredata) == iss[y]))

# Perform
featuredata <- log10(featuredata)
factormat <-model.matrix(~ Sex + Age + Class + Order, sampledata) # adjust to your data
ruv2.analysis <-NormQcmets(featuredata, factormat=factormat,method = "ruv2", # adjust to your data
                      k=4, qcmets = qcmets)
# ruv2.analysis
pval_ruv2 <- as.data.frame(ruv2.analysis$adj.p.value)
cutoff <- 0.05 # set cutoff
ruv2_pval <- rownames(subset(pval_ruv2, pval_ruv2$ClassTG <= cutoff)) # adjust to your data

############################################### 
############################################### BORUTA METHOD
###############################################

library(Boruta)

fs_b <- Boruta(Label~., data = ds, maxRuns = 1000) # adjust to your data
plot(fs_b)
bor <- getSelectedAttributes(fs_b, withTentative = F) # adjust to your data "withTentative" parameter
bor <- gsub("`", '', bor)

############################################### 
############################################### TDFDR METHOD
###############################################
 
library(tdfdr)

# data
dat <- ds[,-1]
s_d <- cbind(Class = ds$Label, Order = order, meta) # adjust to your data
cutoff <- 0.05 # adjust to your data

# perform
tdfdr_f <- tdfdr(dat, s_d$Class, s_d[,c(3,5,6)], alpha = cutoff) # adjust to your data
tdfdr_f_sel <- tdfdr.select(tdfdr_f, fdr.level = cutoff) # adjust to your data
tdfdr_f_sel <- tdfdr_f_sel$pos
td_fdr <- names(tdfdr_f_sel[which(tdfdr_f_sel==T)])

############################################### 
############################################### REMOVE HIGHLY CORRELATED
###############################################

library(caret)

set.seed(1234)
cutoff <- 0.95 # set cut-off
dsc <- ds[,-1]
correlations <- cor(dsc) # method = c("pearson", "kendall", "spearman")
highlyCorrelated <- findCorrelation(correlations, cutoff=cutoff, names = F)
ds_corr <- cbind(Label = ds[,1], dsc[,-highlyCorrelated])
corr <- colnames(ds_corr)[-1]

############################################### 
############################################### CORRELATION WITH SMTH
###############################################

target <- meta$Creatinine # set target feature/component
corr.an <- apply(ds[,-1], 2, function(y) cor.test(x = y, y = target, method = "pearson")[["estimate"]][["cor"]]) # or "pearson", "kendall", "spearman"
corr.an <- data.frame(cbind(colnames(ds[,-1]), corr = abs(as.numeric(corr.an))))
rownames(corr.an) <- colnames(ds[,-1])
th_cor <- 0.9 # set value for filtration
corr_th <- subset(corr.an, corr.an$corr > th_cor)
corr_f <- rownames(corr_th) # features

################################################### Auto detect type of correlation for correlation with specific feature
smart.corr.test <- function(x, n){
  b <- vector(mode="numeric")
  res <- apply(x, 2, shapiro.test)
  for (i in 1:ncol(x)) b[i] = res[[i]]$p.value 
  names(b) <- colnames(x)
  for (i in 1:ncol(x))
    if (b[i] < 0.05) {
      b[i] <- (cor.test(x = x[, i], y = n, method = "spearman")$estimate)
    } else{
      b[i] <- (cor.test(x = x[, i], y = n, method = "pearson")$estimate)
    }
  return(b) }

corr_smart <- as.data.frame(smart.corr.test(ds[,-1], n = target)) # n - vector with target feature
colnames(corr_smart) <- "corr"
corr_th <- subset(corr_smart, abs(corr_smart$corr) > th_cor)
corr_f <- rownames(corr_th) # features

###############################################
############################################### COMBINE ALL RESULTS BY INTERSECTION OR REMAINING ALL UNIQUE VALUES
############################################### 

# intersect
combine <- Reduce(intersect, list(roc, uvf, lim_pval, fc, vip_pls)) # or use tuple package all => (rfe, sfs, nfs, roc, fc, uvf, lm, lmm, corr, lim_pval, pval_ttest, pval_kw, ruv2_pval, bor, vip_pls, rf_perm_f, td_fdr, corr_f, penal_f, step_f)
# all unique
combine <- unique(c(sfs, combine)) # or use tuple package all => (rfe, sfs, nfs, roc, fc, uvf, lm, lmm, corr, lim_pval, pval_ttest, pval_kw, ruv2_pval, bor, vip_pls, rf_perm_f, corr_f, penal_f, step_f, td_fdr)
# combination
combine_df <- cbind(Label = ds[,1], ds[,combine])
fwrite(combine_df, "8 peaks.csv", row.names = T)

###############################################
############################################### COMBINE WITH ANNOTATION
############################################### 

setwd("D:/...")

dsa <- as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats annot+filtr LMM adj KEGG.csv", header=T))[,-1]
dsa <- cbind(dsa[,1:4], dsa[,combine]) # adjust to your data
fwrite(dsa, "8 peaks.csv", row.names = T)

##############################################################################################################################################################
# Classification Machine Learning task
##############################################################################################################################################################

# Content:
# Load data and installation
# Data sampling
# Resampling
# Fit caret classification and make predictions
# Fit logistic/penalized classification and make predictions
# Performance
# Nested Cross-Validation 
# gWQS regression
# Plot results

###############################################
############################################### LOAD DATA AND INSTALLATION
############################################### 

# setup environment
library(data.table)
setwd("D:/...")

# PEAK TABLE
# dataset with intensities and label column only
ds <- as.data.frame(fread(input = "8 peaks.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <-"Label"
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# PEAK TABLE
# dataset with intensities, annotation and label column
ds <- as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats annot+filtr LMM adj KEGG.csv", header=T))
ds <- ds[-c(1:12),] # adjust to your data
rownames(ds) <- ds[,5] # adjust to your data
ds <- ds[,-c(1,3:5)] # adjust to your data
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# METADATA
setwd("D:/...")
meta <- as.data.frame(fread(input = "metadata.csv", header=T))
rownames(meta) <- meta[,5] # adjust to your data
order <- meta[,2] # adjust to your data
meta <- meta[,c(3,8:10)] # adjust to your data
colnames(meta) <- c("Batch", "Creatinine", "Age", "Sex")
# check identical rownames between peak table and metadata
identical(rownames(ds), rownames(meta))

# NEW WD FOR STATISTICAL ANALYSIS
setwd("D:/...")

###############################################
############################################### DATA SPLITTING
############################################### 

library(caret)

# split data into train and validation with balanced class 
set.seed(1234) 
trainIndex <- createDataPartition(ds$Label, p = 0.8, list = F, times = 1) # or simple unbalanced splitting: sample(2, length(ds$Label), replace = T, prob=c(0.8, 0.2))
dsTrain <- ds[ trainIndex,]
dsValid <- ds[-trainIndex,]

###############################################
############################################### RESAMPLING
############################################### 

library(caret)

# repeated cross validation, ROC metric (two class study)
trainControl <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs = T, summaryFunction = twoClassSummary) # or bootstrap: trainControl(method="boot", number=100); adjust to your data
metric <- "ROC" 

# repeated cross validation, Accuracy metric
trainControl <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs = T) # or bootstrap: trainControl(method="boot", number=100); adjust to your data
metric <- "Accuracy" 

###############################################
############################################### FIT CARET CLASSIFICATION AND MAKE PREDICTIONS
############################################### 

library(caret)
library(parallel)
library(doParallel)

# stop parallel
stopCluster(cl)
stopImplicitCluster()

# start parallel processing
fc <- as.numeric(detectCores(logical = T))
cl <- makePSOCKcluster(fc-1)
registerDoParallel(cl)

# PLS (as example) see http://topepo.github.io/caret/train-models-by-tag.html for other classification algorithm
set.seed(1234)
fit.cl <- train(Label~., data=dsTrain, method="pls", metric=metric, trControl=trainControl, tuneLength = 10) # adjust tuneLength to your data
# model and resampling summary
results <- resamples(list(m=fit.cl, m1=fit.cl), trControl = trainControl, metric=metric)$values # see "results" for resampling summary
# to compare 2 or more models (fit.cl and fit.cl1, for example) 
summary(resamples(list(m=fit.cl, m1=fit.cl1), trControl = trainControl, metric=metric)) # summary(resamples(list(m=fit.cl, m1=fit.cl1), trControl = trainControl, metric=metric)) or summary(diff(resamples(list(m=fit.cl, m1=fit.cl1), trControl = trainControl, metric=metric)))
# model info
print(fit.cl)
# make predictions on the validation dataset
predicted.classes <- predict(fit.cl, newdata=dsValid)
probabilities <- predict(fit.cl, newdata=dsValid, type = "prob")[,1]

###############################################
############################################### FIT LOGISTIC/PENALIZED CLASSIFICATION AND MAKE PREDICTIONS
############################################### 

#################################### Logistic

library(dplyr)

log_mod <- glm(Label ~ ., data = dsTrain, family = binomial) # for multiclass -> nnet::multinom
summary(log_mod)
print(log_mod)
# make predictions on the validation dataset
probabilities <- log_mod %>% predict(dsValid, type = "response")
contrasts(dsValid$Label)
predicted.classes <- ifelse(probabilities > 0.5, "TG", "CG") # adjust to your data

#################################### Logistic Stepwise

library(dplyr)
library(MASS)

log_step_mod <- glm(Label ~ ., data = dsTrain, family = binomial) %>% stepAIC(trace = F, direction = "both")
summary(log_step_mod)
print(log_step_mod)
# make predictions on the validation dataset
probabilities <- log_step_mod %>% predict(dsValid, type = "response")
contrasts(dsValid$Label)
predicted.classes <- ifelse(probabilities > 0.5, "TG", "CG") # adjust to your data

#################################### Penalized

library(glmnet)
library(dplyr)

# predictor variables
x <- model.matrix(Label~., dsTrain)[,-1]
# outcome (Label) to a numerical variable
y <- ifelse(dsTrain$Label == "TG", 1, 0) # set class for target -> "TG"

# perform

# alpha: the elasticnet mixing parameter. Allowed values include:
# "1": for lasso regression
# "0": for ridge regression
# a value between 0 and 1 for elastic net regression
set.seed(1234) 
alpha <- 0.5 # set type of penalized regression
cv.penal <- cv.glmnet(x, y, alpha = alpha, family = "binomial", nfolds = 10, type.measure="auc") # adjust type.measure for your data
# plot(cv.penal)

# Final model with optimal lambda

# lambda = cv.penal$lambda.1se or lambda = cv.penal$lambda.min
# lambda = lambda.1se produces a simpler model compared to lambda.min, but the model might be a little bit less accurate than the one obtained with lambda.min
penal.model <- glmnet(x, y, alpha = alpha, family = "binomial", lambda = cv.penal$lambda.min)
summary(penal.model)
print(penal.model)
coef_glm <- as.matrix(coef(penal.model)) # coefficients of regression
fs_glm <- subset(coef_glm[-1,], coef_glm[-1,] != 0) # selected variables
names(fs_glm) # variable names
# make predictions on the validation dataset
x.test <- model.matrix(Label ~., dsValid)[,-1]
probabilities <- penal.model %>% predict(newx = x.test, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "TG", "CG") # adjust to your data

###############################################
############################################### PERFORMANCE
###############################################

library(caret)
library(pROC)

# Confusion Matrix
confusionMatrix(predicted.classes, dsValid$Label) # confusionMatrix(as.factor(predicted.classes), as.factor(dsValid$Label)) or confusionMatrix(predicted.classes, dsValid$Label)

# ROC curve
res.roc <- roc(dsValid$Label, probabilities, levels = levels(dsValid$Label))
plot.roc(res.roc, print.auc = TRUE) # for compute only AUC: auc(res.roc)

###############################################
############################################### NESTED CROSS-VALIDATION
###############################################

library(tibble)
library(dplyr)
library(rsample)
library(caret)
library(parallel)
library(doParallel)

# stop parallel
stopCluster(cl)
stopImplicitCluster()

# start parallel processing
fc <- as.numeric(detectCores(logical = T))
cl <- makePSOCKcluster(fc-1)
registerDoParallel(cl)

# Set parameters
set.seed(1234)
metric <- "Accuracy" # set metric
outer <- vfold_cv(ds, v = 10, repeats = 1, strata = "Label") # Type of resampling in outer loop, adjust to your data
tc_out <- trainControl(method = "cv", number = 10, savePredictions = T, classProbs = T) # Type of resampling in inner loop, adjust to your data
results <- tibble(fold = numeric(), par = numeric(), acc_test = numeric()) 

# Perform
for (i in seq_along(outer$splits)) { # !!! if some errors try different tuneLength 
  trn_x <- outer$splits[[i]] %>% 
    analysis() %>% 
    dplyr::select(-Label) %>% 
    as.data.frame() 
  
  trn_y <- outer$splits[[i]] %>% 
    analysis() %>%
    pull(Label) 
  
  test <- outer$splits[[i]] %>% 
    assessment() 
  
  # PLS (as example) see http://topepo.github.io/caret/train-models-by-tag.html for other classification algorithm
  m_outer <- train(x = trn_x, y = trn_y, method = "pls", metric=metric, trControl = tc_out, tuneLength = 5) # adjust tuneLength to your data
  
  # predictions in held-out fold from outer loop 
  acc <- confusionMatrix(predict(m_outer, test), test$Label) 
  results <- add_row(results, fold = i, par = as.numeric(m_outer$bestTune), acc_test = acc$overall[metric]) # adjust to your data
}

# Results
as.data.frame(results) # par -> best hyperparameter
results$acc_test 
mean(results$acc_test)
plot(results$acc_test)

#################################### Only caret version
library(caret)
library(parallel)
library(doParallel)

# stop parallel
stopCluster(cl)
stopImplicitCluster()

# start parallel processing
fc <- as.numeric(detectCores(logical = T))
cl <- makePSOCKcluster(fc-1)
registerDoParallel(cl)

# Set parameters
set.seed(1234)
metric <- "Accuracy" # set metric
outer <- createFolds(ds$Label, k = 10) # Type of resampling in outer loop, adjust to your data
tc_out <- trainControl(method = "cv", number = 10, savePredictions = T, classProbs = T) # Type of resampling in inner loop, adjust to your data
results <- data.frame(fold = numeric(), par = numeric(), acc_test = numeric()) 

# Perform
for (i in 1:length(outer)) { 
  
  trn_x <- ds[-outer[[i]], -1]
  trn_y <- ds[-outer[[i]], 1]
  test <- ds[outer[[i]], ]
  
  # PLS (as example) see http://topepo.github.io/caret/train-models-by-tag.html for other classification algorithm
  m_outer <- train(x = trn_x, y = trn_y, method = "pls", metric=metric, trControl = tc_out, tuneLength = 5) # adjust tuneLength to your data
  
  # predictions in held-out fold from outer loop 
  acc <- confusionMatrix(predict(m_outer, test), test$Label) 
  results <- dplyr::add_row(results, fold = i, par = as.numeric(m_outer$bestTune), acc_test = acc$overall[metric]) # adjust to your data
}

# Results
as.data.frame(results) # par -> best hyperparameter
results$acc_test 
mean(results$acc_test)
plot(results$acc_test)

###############################################
############################################### gWQS REGRESSION
###############################################

library(gWQS)
library(caret)
library(pROC)

# generate data
meta <- cbind(Order = order, meta)
ds1 <- cbind(Label = ds$Label, meta, ds[,-1])
ds1$Label <- as.numeric(ds1$Label) # factors to numeric

# split data into train and validation with balanced class 
set.seed(1234) 
trainIndex <- createDataPartition(ds1$Label, p = 0.8, list = F, times = 1) # or simple unbalanced splitting: sample(2, length(ds$Label), replace = T, prob=c(0.8, 0.2))
dsTrain <- ds1[ trainIndex,]
dsValid <- ds1[-trainIndex,]

# perform
formula <- names(ds1)[-c(1:6)] # adjust to your data
results <- gwqs(Label ~ wqs, mix_name = formula, data = dsTrain, 
                q = 10, validation = 0.8, b = 50, b1_pos = F,  # adjust to your data
                b1_constr = F, family = gaussian, seed = 1234)

# summary
summary(results)

# performance
predictions <- predict(results, newdata = dsValid)
predictions <- round(predictions$df_pred$ypred,0)
confusionMatrix(as.factor(predictions), as.factor(dsValid$Label))

# visualization
# bar plot
gwqs_barplot(results)
# scatter plot y vs wqs
gwqs_scatterplot(results)
# scatter plot residuals vs fitted values
gwqs_fitted_vs_resid(results)
# ROC plot
probabilities <- predict(results, newdata = dsValid, type = "response")$df_pred[,2]
res.roc <- roc(as.factor(dsValid$Label), probabilities)
plot.roc(res.roc, print.auc = TRUE)

###############################################
############################################### PLOT RESULTS
###############################################

library(reshape2)
library(ggplot2)

df.m <- melt(ds, id.var = "Label") # reshape data frame

p <- ggplot(data = df.m, aes(x=variable, y=value)) + xlab("") + ylab("") +
  geom_boxplot(aes(fill=Label)) + theme(legend.position="bottom") + theme_classic() + scale_x_discrete(labels=c("")) 

pp <- p + facet_wrap( ~ variable, scales="free") + theme_classic() + theme(legend.position="bottom") 

pp

##############################################################################################################################################################
# Regression Machine Learning task
##############################################################################################################################################################

# Content:
# Load data and installation
# Data sampling
# Resampling
# Fit caret regression and make predictions
# Stepwise regression
# Penalized regression
# Performance
# Nested Cross-Validation 
# gWQS regression
# Plot results

###############################################
############################################### LOAD DATA AND INSTALLATION
############################################### 

# setup environment
library(data.table)
setwd("D:/...")

# PEAK TABLE
# dataset with intensities and label column only
ds <- as.data.frame(fread(input = "8 peaks.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <-"Label"
ds[,-1] <- sapply(ds[,-1], as.numeric)
# ds$Label should be numeric: 
ds$Label <- c(1:nrow(ds)) # toy example

# PEAK TABLE
# dataset with intensities, annotation and label column
ds <- as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats annot+filtr LMM adj KEGG.csv", header=T))
ds <- ds[-c(1:12),] # adjust to your data
rownames(ds) <- ds[,5] # adjust to your data
ds <- ds[,-c(1,3:5)] # adjust to your data
ds[,-1] <- sapply(ds[,-1], as.numeric)
# ds$Label should be numeric: 
ds$Label <- c(1:nrow(ds)) # toy example

# METADATA
setwd("D:/...")
meta <- as.data.frame(fread(input = "metadata.csv", header=T))
rownames(meta) <- meta[,5] # adjust to your data
order <- meta[,2] # adjust to your data
meta <- meta[,c(3,8:10)] # adjust to your data
colnames(meta) <- c("Batch", "Creatinine", "Age", "Sex")
# check identical rownames between peak table and metadata
identical(rownames(ds), rownames(meta))

###############################################
############################################### DATA SPLITTING
############################################### 

library(caret)

# split data into train and validation with balanced class 
set.seed(1234) 
trainIndex <- createDataPartition(ds$Label, p = 0.8, list = F, times = 1) # or simple unbalanced splitting: sample(2, length(ds$Label), replace = T, prob=c(0.8, 0.2))
dsTrain <- ds[ trainIndex,]
dsValid <- ds[-trainIndex,]

###############################################
############################################### RESAMPLING
############################################### 

library(caret)

# repeated cross validation
trainControl <- trainControl(method="repeatedcv", number=10, repeats=10) # or bootstrap: trainControl(method="boot", number=100); adjust to your data
metric <- "RMSE" # RMSE or MAE or Rsquared

###############################################
############################################### FIT CARET REGRESSION AND MAKE PREDICTIONS
############################################### 

library(caret)
library(parallel)
library(doParallel)

# stop parallel
stopCluster(cl)
stopImplicitCluster()

# start parallel processing
fc <- as.numeric(detectCores(logical = T))
cl <- makePSOCKcluster(fc-1)
registerDoParallel(cl)

# LM (as example) see http://topepo.github.io/caret/train-models-by-tag.html for other regression algorithm
set.seed(1234)
fit.regr <- train(Label~ ., data=dsTrain, method="lm", metric=metric, trControl=trainControl, tuneLength = 10) # adjust tuneLength to your data
# model summary
summary(fit.regr)
print(fit.regr)
# make predictions on the validation dataset
predictions <- predict(fit.regr, newdata=dsValid)

###############################################
############################################### STEPWISE REGRESSION
###############################################

library(leaps)
library(MASS)
library(dplyr)
library(caret)
library(parallel)
library(doParallel)

# Label should be numeric: 
ds$Label <- as.numeric(ds$Label)
dsTrain$Label <- as.numeric(dsTrain$Label)
dsValid$Label <- as.numeric(dsValid$Label)

# Full model
full_m <- lm(Label~ ., data=dsTrain)

#################################### Stepwise regression by MASS
step.model <- stepAIC(full_m, direction = "both", trace = F) # "both", "backward", "forward"
summary(step.model)
# print coefficients of selected model
coef(step.model)

# predictions
predictions <- step.model %>% predict(dsValid)

#################################### Stepwise regression by leaps
models_ss <- regsubsets(Label~ ., data=dsTrain, nvmax = ncol(dsTrain[,-1]),  method = "seqrep") # adjust nvmax, method=c("exhaustive","backward", "forward", "seqrep")
summary(models_ss)
res_sum <- summary(models_ss)
# print Adj.R2, CP and BIC stats for every model ID
data.frame(Adj.R2 = which.max(res_sum$adjr2), 
           CP = which.max(res_sum$cp),
           BIC = which.max(res_sum$bic))
# print coefficients of selected model
coefs <- coef(models_ss, 2) # enter selected model ID
coefs

# predict
valid.mat = model.matrix(Label ~ ., data = dsValid)
predictions <- valid.mat[, names(coefs)]%*%coefs

#################################### Stepwise regression by caret
# stop parallel
stopCluster(cl)
stopImplicitCluster()

# start parallel processing
fc <- as.numeric(detectCores(logical = T))
cl <- makePSOCKcluster(fc-1)
registerDoParallel(cl)

# repeated cross validation
trainControl <- trainControl(method="repeatedcv", number=10, repeats=10) # or bootstrap: trainControl(method="boot", number=100)
metric <- "RMSE" # RMSE or MAE or Rsquared

# see http://topepo.github.io/caret/train-models-by-tag.html for other regression algorithm
set.seed(1234)
fit.regr <- train(Label~ ., data=dsTrain, method="lmStepAIC", metric=metric, trControl=trainControl, tuneLength = 5) # leapBackward or leapSeq or leapForward or lmStepAIC, adjust tuneLength to your data

# results 
fit.regr$results
fit.regr$bestTune
summary(fit.regr)
print(fit.regr)
# print variables names for results of best tune for leaps methods only!
results_leaps <- as.data.frame(summary(fit.regr)$which)
colnames(results_leaps)[which(results_leaps[as.numeric(fit.regr$bestTune),] == T)]
# print coefficients and variable names for results of best tune for MASS method only!
summary(fit.regr)$coefficients

# predictions
predictions <- fit.regr %>% predict(dsValid)

###############################################
############################################### PENALIZED REGRESSION
############################################### 

library(glmnet)
library(dplyr)

# predictor variables
x <- model.matrix(Label~., dsTrain)[,-1]
# outcome (Label) to a numerical variable
y <- as.numeric(dsTrain$Label)

# perform

# alpha: the elasticnet mixing parameter. Allowed values include:
# "1": for lasso regression
# "0": for ridge regression
# a value between 0 and 1 for elastic net regression
set.seed(1234) 
alpha <- 0.5 # set type of penalized regression
cv.penal <- cv.glmnet(x, y, alpha = alpha, nfolds = 10, type.measure="mse") # adjust type.measure for your data
# plot(cv.penal)

# Final model with optimal lambda

# lambda = cv.penal$lambda.1se or lambda = cv.penal$lambda.min
# lambda = lambda.1se produces a simpler model compared to lambda.min, but the model might be a little bit less accurate than the one obtained with lambda.min
penal.model <- glmnet(x, y, alpha = alpha, lambda = cv.penal$lambda.min)
summary(penal.model)
print(penal.model)
coef_glm <- as.matrix(coef(penal.model)) # coefficients of regression
fs_glm <- subset(coef_glm[-1,], coef_glm[-1,] != 0) # selected variables
names(fs_glm) # variable names
# make predictions on the validation dataset
x.test <- model.matrix(Label ~., dsValid)[,-1]
predictions <- penal.model %>% predict(newx = x.test)

###############################################
############################################### PERFORMANCE
###############################################

library(caret)

# Quality metrics
rmse <- caret::RMSE(predictions, dsValid$Label)
r2 <- caret::R2(predictions, dsValid$Label)
mae <- caret::MAE(predictions, dsValid$Label)
tbl_regr <- c("RMSE" = rmse, "R2" = r2, "MAE" = mae)
tbl_regr

###############################################
############################################### NESTED CROSS-VALIDATION
###############################################

library(tibble)
library(dplyr)
library(rsample)
library(caret)
library(parallel)
library(doParallel)

# stop parallel
stopCluster(cl)
stopImplicitCluster()

# start parallel processing
fc <- as.numeric(detectCores(logical = T))
cl <- makePSOCKcluster(fc-1)
registerDoParallel(cl)

# Set parameters
set.seed(1234)
metric <- "RMSE" # set metric
outer <- vfold_cv(ds, v = 10, repeats = 1, strata = "Label") # Type of resampling in outer loop, adjust to your data
tc_out <- trainControl(method = "cv", number = 10, savePredictions = T, classProbs = F) # Type of resampling in inner loop, adjust to your data
results <- tibble(fold = numeric(), par = numeric(), perf_test = numeric()) 

# Perform
for (i in seq_along(outer$splits)) { # !!! if some errors try different tuneLength 
  trn_x <- outer$splits[[i]] %>% 
    analysis() %>% 
    dplyr::select(-Label) %>% 
    as.data.frame() 
  
  trn_y <- outer$splits[[i]] %>% 
    analysis() %>%
    pull(Label) 
  
  test <- outer$splits[[i]] %>% 
    assessment() 
  
  # LM (as example) see http://topepo.github.io/caret/train-models-by-tag.html for other regression algorithm
  m_outer <- train(x = trn_x, y = trn_y, method = "lm", metric=metric, trControl = tc_out, tuneLength = 5) # adjust tuneLength to your data
  
  # predictions in held-out fold from outer loop 
  rmse <- caret::RMSE(predict(m_outer, test), test$Label) # adjust to your data
  results <- add_row(results, fold = i, par = as.numeric(m_outer$bestTune), perf_test = rmse) # adjust to your data
}

# Results
as.data.frame(results) # par -> best hyperparameter
results$perf_test 
mean(results$perf_test)
plot(results$perf_test)

#################################### Only caret version
library(caret)
library(parallel)
library(doParallel)

# stop parallel
stopCluster(cl)
stopImplicitCluster()

# start parallel processing
fc <- as.numeric(detectCores(logical = T))
cl <- makePSOCKcluster(fc-1)
registerDoParallel(cl)

# Set parameters
set.seed(1234)
metric <- "RMSE" # set metric
outer <- createFolds(ds$Label, k = 10) # Type of resampling in outer loop, adjust to your data
tc_out <- trainControl(method = "cv", number = 10, savePredictions = T, classProbs = F) # Type of resampling in inner loop, adjust to your data
results <- data.frame(fold = numeric(), par = numeric(), perf_test = numeric()) 

# Perform
for (i in 1:length(outer)) { # !!! if some errors try different tuneLength 
  
  trn_x <- ds[-outer[[i]], -1]
  trn_y <- ds[-outer[[i]], 1]
  test <- ds[outer[[i]], ]
  
  # LM (as example) see http://topepo.github.io/caret/train-models-by-tag.html for other regression algorithm
  m_outer <- train(x = trn_x, y = trn_y, method = "lm", metric=metric, trControl = tc_out, tuneLength = 5) # adjust tuneLength to your data
  
  # predictions in held-out fold from outer loop 
  rmse <- caret::RMSE(predict(m_outer, test), test$Label) # adjust to your data
  results <- dplyr::add_row(results, fold = i, par = as.numeric(m_outer$bestTune), perf_test = rmse) # adjust to your data
}

# Results
as.data.frame(results) # par -> best hyperparameter
results$perf_test 
mean(results$perf_test)
plot(results$perf_test)

###############################################
############################################### gWQS REGRESSION
###############################################

library(gWQS)
library(caret)

# generate data
meta <- cbind(Order = order, meta)
ds1 <- cbind(Label = ds$Label, meta, ds[,-1])

# split data into train and validation with balanced class 
set.seed(1234) 
trainIndex <- createDataPartition(ds1$Label, p = 0.8, list = F, times = 1) # or simple unbalanced splitting: sample(2, length(ds$Label), replace = T, prob=c(0.8, 0.2))
dsTrain1 <- ds1[ trainIndex,]
dsValid1 <- ds1[-trainIndex,]

# perform
formula <- names(ds1)[-c(1:6)] # adjust to your data
results <- gwqs(Label ~ wqs, mix_name = formula, data = dsTrain1, 
                q = 10, validation = 0.8, b = 10, b1_pos = T,  # adjust to your data
                b1_constr = F, family = gaussian, seed = 1234)

# summary
summary(results)

# performance
predictions <- predict(results, newdata = dsValid1)
predictions <- predictions$df_pred$ypred
rmse <- caret::RMSE(predictions, dsValid1$Label)
r2 <- caret::R2(predictions, dsValid1$Label)
mae <- caret::MAE(predictions, dsValid1$Label)
tbl_regr <- c("RMSE" = rmse, "R2" = r2, "MAE" = mae)
tbl_regr

# visualization
# bar plot
gwqs_barplot(results)
# scatter plot y vs wqs
gwqs_scatterplot(results)
# scatter plot residuals vs fitted values
gwqs_fitted_vs_resid(results)

###############################################
############################################### PLOT RESULTS
###############################################

# plot original values versus predicted
x <- 1:length(dsValid$Label)
plot(x, dsValid$Label, col = "red", type = "l", lwd=2,
     main = "", xlab = "Observations", ylab = "Values")
lines(x, predictions, col = "blue", lwd=2)
legend("top",  legend = c("original", "predicted"), 
       fill = c("red", "blue"), col = 2:3,  adj = c(0, 0.6))

##############################################################################################################################################################
# Testing set of biomarkers
##############################################################################################################################################################

# Content:
# Load data and installation
# Plot results
# MANOVA
# PERMANOVA
# Moderated t-test
# Fold Change Calculation
# Comparing Means

###############################################
############################################### LOAD DATA AND INSTALLATION
############################################### 

# setup environment
library(data.table)
setwd("D:/...")

# PEAK TABLE
# dataset with intensities and label column only
ds <- as.data.frame(fread(input = "8 peaks.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <-"Label"
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# METADATA
setwd("D:/...")
meta <- as.data.frame(fread(input = "metadata.csv", header=T))
rownames(meta) <- meta[,5] # adjust to your data
order <- meta[,2] # adjust to your data
meta <- meta[,c(3,8:10)] # adjust to your data
colnames(meta) <- c("Batch", "Creatinine", "Age", "Sex")
# check identical rownames between peak table and metadata
identical(rownames(ds), rownames(meta))

###############################################
############################################### PLOT RESULTS
###############################################

library(reshape2)
library(ggplot2)

#################################### box plots / violin plots
df.m <- melt(ds, id.var = "Label") # reshape data frame

p <- ggplot(data = df.m, aes(x=variable, y=value)) + xlab("") + ylab("") +
  geom_boxplot(aes(fill=Label)) + theme(legend.position="bottom") + theme_classic() + scale_x_discrete(labels=c("")) 

# try also: geom_boxplot or geom_violin

pp <- p + facet_wrap( ~ variable, scales="free") + theme_classic() + theme(legend.position="bottom") 

pp

#################################### scatter plots

df.m <- melt(ds, id.var = "Label")
df.m <- cbind(1:nrow(df.m), df.m)
colnames(df.m)[1] <- "Patient"
p <- ggplot(data = df.m, aes(x=Patient, y=value)) + xlab("") + ylab("") +
  geom_point(aes(colour=Label)) + theme(legend.position="bottom") + theme_classic() + scale_x_discrete(labels=c("")) + geom_smooth(method = "lm") # or method = "loess" or no geom_smooth

pp <- p + facet_wrap( ~ variable, scales="free") + theme_classic() + theme(legend.position="bottom") 

pp
# or pairs(ds[,-1]) or psych::pairs.panels(ds[,-1], method = "pearson", hist.col = "#00AFBB",density = T, ellipses = T)

#################################### box plots / violin plots for 2 group variables

Sex <- as.factor(meta$Sex) # name and select factor from metadata
ds2 <- cbind(Sex, ds) 
df.m2 <- melt(ds2, id.var = c("Label", "Sex"))

p <- ggplot(data = df.m2, aes(x=variable, y=value)) + xlab("") + ylab("") +
  geom_boxplot(aes(fill = Label, linetype = Sex)) + theme(legend.position="bottom") + theme_classic() + scale_x_discrete(labels=c("")) 

# try also: geom_boxplot or geom_violin

pp <- p + facet_wrap( ~ variable, scales="free") + theme_classic() + theme(legend.position="bottom") 

pp

###############################################
############################################### MANOVA
###############################################

mnv <- sapply(2:ncol(ds), function(y) cbind(ds[,y]))
res.man <- manova(mnv ~ Label, data = ds) # adjust to data
summary(res.man)
summary.aov(res.man)
sum.manova <- summary.aov(res.man)
pval_manova <- lapply(1:length(sum.manova), function(y) sum.manova[[y]][,5]) # adjust to data 
pval_manova <- as.data.frame(do.call(rbind, pval_manova))
pval_manova <- as.data.frame(apply(pval_manova, 2, function(x) p.adjust(x, method = "BH")))
colnames(pval_manova) <- rownames(sum.manova[[1]])
pval_manova

###############################################
############################################### PERMANOVA
###############################################

library(vegan)
library(pairwiseAdonis)

pmnv <- ds[,-1]
res.perman <- adonis(pmnv ~ Label, data = ds, method = "euclidean", permutations = 1000) # perform PERMANOVA
res.perman
res.perman.ph <- pairwise.adonis2(ds[,-1]~Label, data = ds, p.adjust.m = "BH", perm = 1000) # perform PERMANOVA with multilevel comparison
res.perman.ph

###############################################
############################################### MODERATED T-TEST
###############################################

library(limma)

mdl_mtrx <- model.matrix(~Label, ds) # adjust to your data
lmf <- lmFit(t(log2(ds[,-1])), method = "robust", design = mdl_mtrx, maxit = 1000) # "robust" or "ls"
efit <- eBayes(lmf)
tableTop <- topTable(efit, coef = 2, adjust = "BH", p.value = 0.05, number = ncol(ds)) # select method
tableTop

############################################### 
############################################### FOLD CHANGE CALCULATION
###############################################

library(dplyr)

FOLD.CHANGE <- function(data) {
  ds_subsets <- lapply(1:length(unique(data[,1])), function(y) dplyr::filter(data[,-1], data$Label == unique(data[,1])[y])) # list of subsets by label
  mean_r_l <- lapply(1:length(ds_subsets), function(y) apply(ds_subsets[[y]], 2, mean, na.rm = T)) # calculate mean for feature
  foldchange <- log2((mean_r_l[[1]] / mean_r_l[[2]]))
  fc_res <- as.data.frame(foldchange)
  return(fc_res)
}

fc_res <- FOLD.CHANGE(ds)
fc_res

################################################### Multigroup Fold Change

FOLD.CHANGE.MG <- function(x, f, aggr_FUN = colMeans, combi_FUN = {function(x,y) "-"(x,y)}){
  f <- as.factor(f)
  i <- split(1:nrow(x), f)
  x <- sapply(i, function(i){ aggr_FUN(x[i,])})
  x <- t(x)
  x <- log2(x)
  j <- combn(levels(f), 2)
  ret <- combi_FUN(x[j[1,],], x[j[2,],])
  rownames(ret) <- paste(j[1,], j[2,], sep = '/')
  t(ret)
}

fdr <- FOLD.CHANGE.MG(ds[,-1], ds[,1])

# by mean
fdr_mean <- apply(abs(fdr),1, mean, na.rm=T)

############################################### 
############################################### COMPARING MEANS 
###############################################

#----------------------------------
# designed for 2 groups comparison
#----------------------------------

#################################### t-test

# prepare data
data_l <- lapply(1:length(unique(ds$Label)), function(y) subset(ds, Label==unique(ds$Label)[y])[,-1])
data_l <- lapply(1:length(data_l), function(y) sapply(data_l[[y]], as.numeric))
# perform
res.t.test <- lapply(1:ncol(data_l[[1]]), function(y) t.test(x = data_l[[1]][,y], y = data_l[[2]][,y])) # t-test # Try also pairwise.t.test
res.t.test.pval <- sapply(1:ncol(data_l[[1]]), function(y) t.test(x = data_l[[1]][,y], y = data_l[[2]][,y])$p.value)
p_adj <-p.adjust(res.t.test.pval, method = "BH") # Adjust P-values for Multiple Comparisons
p_adj

#################################### Wilcoxon test (aka Mann-Whitney test)

# prepare data
data_l <- lapply(1:length(unique(ds$Label)), function(y) subset(ds, Label==unique(ds$Label)[y])[,-1])
data_l <- lapply(1:length(data_l), function(y) sapply(data_l[[y]], as.numeric))
# perform
res.w.test <- lapply(1:ncol(data_l[[1]]), function(y) wilcox.test(x = data_l[[1]][,y], y = data_l[[2]][,y])) # wilcox test # Try also pairwise.wilcox.test
res.w.test.pval <- sapply(1:ncol(data_l[[1]]), function(y) wilcox.test(x = data_l[[1]][,y], y = data_l[[2]][,y])$p.value)
p_adj <-p.adjust(res.w.test.pval, method = "BH") # Adjust P-values for Multiple Comparisons
p_adj

#################################### Kruskal-Wallis test

# perform
res.kw.test <- lapply(2:ncol(ds), function(y) kruskal.test(ds[,y] ~ ds[,1])) # kruskal test
res.kw.test.pval <- sapply(2:ncol(ds), function(y) kruskal.test(ds[,y] ~ ds[,1])$p.value)
p_adj <-p.adjust(res.kw.test.pval, method = "BH") # Adjust P-values for Multiple Comparisons
p_adj

#################################### One-way/Two-way ANOVA

# prepare data
dat <- cbind(sex = as.factor(meta$Sex), ds) # adjust to your data
# perform
res.anova <- lapply(3:ncol(dat), function(y) aov(dat[,y] ~ sex*Label, data = dat)) # adjust to your data
res.anova
res.anova.sum <- lapply(3:ncol(dat), function(y) summary(aov(dat[,y] ~ sex*Label, data = dat)))
res.anova.sum
p_adj <- lapply(1:length(res.anova.sum), function(y) res.anova.sum[[y]][[1]][,5]) # adjust to your data 
p_adj <- do.call(rbind, p_adj)
p_adj <- as.data.frame(apply(p_adj, 2, function(x) p.adjust(x, "BH"))) # Adjust P-values for Multiple Comparisons
colnames(p_adj) <- rownames(res.anova.sum[[1]][[1]])
p_adj

library(multcomp)
glh <- lapply(1:length(res.anova), function(y) summary(glht(res.anova[[y]]))) # adjust to your data # linfct = mcp(Label = "Tukey") or adjust "linfct = mcp()" to your data
glh
tukey <- lapply(1:length(res.anova), function(y) TukeyHSD(res.anova[[y]])) # adjust to your data # which = "Label" or adjust "which" to your data
tukey

#################################### Sequential implementation of appropriate statistical test with automatic detection for normality and homogeneity, returns only adjusted p-value 

#----------------------------------
# designed for 2 groups comparison
#----------------------------------

# if some error in Shapiro normality test:
# use shapiro.wilk.test function from cwhmisc instead shapiro.test from stats
# library(cwhmisc)
# norm.test <- apply(xx, 2, function(t) cwhmisc::shapiro.wilk.test(t)$p)

uva <- function(x, p.adjust = "BH"){
  
  norm_homog_tests <- function(x) {
    xx <- x[,-1]
    # normality test
    norm.test <- apply(xx, 2, function(t) shapiro.test(t)$p.value)
    # if some error in Shapiro normality test:
    # use shapiro.wilk.test function from cwhmisc instead shapiro.test from stats
    # library(cwhmisc)
    # norm.test <- apply(xx, 2, function(t) cwhmisc::shapiro.wilk.test(t)$p)
    
    # homogeneity test
    homog.test <- apply(xx, 2, function(t) bartlett.test(t,g = x[,1])$p.value)
    return(as.data.frame(cbind(norm.test, homog.test)))}
  
  res_tests <- norm_homog_tests(x)
  
  wilcox_test <- function(x,y) {
    xx <- x[,-1]
    wx.t <- as.vector(which(y[,1] < 0.05))
    wilcox_test <- list()
    ifelse(identical(wx.t, integer(0)), return (wilcox_test <- 1), wx.t)
    wilcox_test <- sapply(as.data.frame(xx[,wx.t]),  function(t) pairwise.wilcox.test(x = t, g =  x[,1], paired=F)$p.value) # or as.numeric(as,vector(...)))$p.value))) or ...))
    wilcox_test <- p.adjust(wilcox_test, method = p.adjust)
    names(wilcox_test) <- (colnames(x)[-1])[wx.t]
    return(as.list(wilcox_test))}
  
  wx.t.res <- wilcox_test(x, res_tests)
    
  welch_test <- function(x,y) {
    xx <- x[,-1]
    wl.t <- as.vector(which(y[,1] > 0.05 & y[,2] < 0.05))
    welch_test <- list()
    ifelse(identical(wl.t, integer(0)), return (welch_test <- 1), wl.t)
    welch_test <- sapply(as.data.frame(xx[,wl.t]), function(t) pairwise.t.test(x = t, g = x[,1], pool.sd = F)$p.value) # or as.numeric(as,vector(...)))$p.value))) or ...))
    welch_test <- p.adjust(welch_test, method = p.adjust)
    names(welch_test) <- (colnames(x)[-1])[wl.t]
    return(as.list(welch_test))}
  
  wl.t.res <- welch_test(x, res_tests)
  
  student_test <- function(x,y) {
    xx <- x[,-1]
    st.t <- as.vector(which(y[,1] > 0.05 & y[,2] > 0.05))
    student_test <- list()
    ifelse(identical(st.t, integer(0)), return (student_test <- 1), st.t)
    student_test <- sapply(as.data.frame(xx[,st.t]), function(t) pairwise.t.test(x = t, g = x[,1], pool.sd = T)$p.value) # or as.numeric(as,vector(...)))$p.value))) or ...))
    student_test <- p.adjust(student_test, method = p.adjust)
    names(student_test) <- (colnames(x)[-1])[st.t]
    return(as.list(student_test))}
  
  st.t.res <- student_test(x, res_tests)
  
  return(list(wilcox_test = wx.t.res, welch_test = wl.t.res, student_test = st.t.res))
}

# 1 st argument -> dataset with 1st "Label" column, 2nd -> the method of adjustment for multiple comparisons.
ds_uva <- uva(ds, p.adjust = "BH")
ds_uva

##############################################################################################################################################################
# MWAS/ANCOVA
##############################################################################################################################################################

# setup environment
library(data.table)
setwd("D:/...")

# dataset with intensities and label column only
ds <- as.data.frame(fread(input = "8 peaks.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <-"Label"
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# METADATA
setwd("D:/...")
meta <- as.data.frame(fread(input = "metadata.csv", header=T))
rownames(meta) <- meta[,5] # adjust to your data
order <- meta[,2] # adjust to your data
meta <- meta[,c(3,8:10)] # adjust to your data
colnames(meta) <- c("Batch", "Creatinine", "Age", "Sex")
# check identical rownames between peak table and metadata
identical(rownames(ds), rownames(meta))

############################################### 
############################################### ANCOVA
###############################################

################################################### LM modeling
library(MetabolomicsBasics)

# data
dat <- ds[,-1]
s_d <- cbind(Class = ds$Label, Order = order, meta)

# perform
model <- MetaboliteANOVA(dat=dat, sam=s_d, model="Batch+Class+Age+Sex", method = "BH") # "Batch+Class+Age+Sex" or "Batch+Class+Age+Sex+Order" # select method and adjust formula
model

################################################### other type of LM modeling

library(multcomp)

# data
dat <- ds[,-1]
s_d <- cbind(Class = ds$Label, Order = order, meta)
ds_ancova <- as.data.frame(cbind(s_d, dat))

# perform
n_start <- 7 # adjust to your data
ancova <- lapply(n_start:ncol(ds_ancova), function(y) aov(ds_ancova[,y]~ Batch+Class+Age+Sex, ds_ancova)) # adjust formula to your data
sum_ancova <- lapply(1:length(ancova), function(y) summary(ancova[[y]]))
sum_ancova
res_ancova <- lapply(1:length(ancova), function(y) summary(glht(ancova[[y]]))) # linfct = mcp(Class = "Tukey") or adjust linfct = mcp() to your data
res_ancova

# adjusted p-value
res_ancova <- as.data.frame(t(sapply(1:length(ancova), function(y) summary(ancova[[y]])[[1]][,5]))) # adjust to your data 
res_ancova <- as.data.frame(apply(res_ancova, 2, function(x) p.adjust(x, "BH"))) # choose p.adjust method
colnames(res_ancova) <- rownames(summary(ancova[[1]])[[1]])
res_ancova

# perform other type
n_start <- 7 # adjust to your data
ancova <- lapply(n_start:ncol(ds_ancova), function(y) anova(lm(ds_ancova[,y]~ Batch+Class+Age+Sex, ds_ancova))) # adjust formula to your data
ancova

# adjusted p-value
res_ancova <- as.data.frame(t(sapply(1:length(ancova), function(y) ancova[[y]][,5]))) # adjust to your data 
res_ancova <- as.data.frame(apply(res_ancova, 2, function(x) p.adjust(x, "BH"))) # choose p.adjust method
colnames(res_ancova) <- rownames(summary(ancova[[1]])[[1]])
res_ancova

################################################### LMM modeling

library(lme4)
library(lmerTest)

# data
dat <- cbind(meta, ds)
n_meta <- 5 # adjust to your data
colnames(dat)[-c(1:n_meta)] <- paste0("X", c(1:ncol(dat))) # adjust to your data
dat$Batch <- as.factor(dat$Batch) # as.factor(dat$Batch) or as.numeric(dat$Batch) or as.numeric(stringr::str_remove(dat$Batch, "b"))
dat$Age <- as.integer(dat$Age)
dat$Sex <- as.integer(dat$Sex)
dat$Label <- as.factor(dat$Label)
dat[,-c(1:n_meta)] <- sapply(dat[,-c(1:n_meta)], as.numeric)

# perform
n_start <- 6 # adjust to your data
lmm_fit <- lapply(n_start:ncol(dat), function(x) lmer(dat[,x] ~ Label + Sex + Age + (1|Batch), dat)) # adjust to your data # lmer from lmerTest package
lmm_fit_coef <- lapply(1:length(lmm_fit), function(x) summary(lmm_fit[[x]])$coefficients)
lmm_fit_pval_all_df <- as.data.frame(t(sapply(1:length(lmm_fit_coef), function(x) lmm_fit_coef[[x]][,5])))
lmm_fit_pval_all_df <- as.data.frame(sapply(1:ncol(lmm_fit_pval_all_df), function(x) p.adjust(lmm_fit_pval_all_df[,x], method = "BH"))) # select method
dat2 <- cbind(meta, ds) 
rownames(lmm_fit_pval_all_df) <- colnames(dat2)[-c(1:n_meta)]
colnames(lmm_fit_pval_all_df) <- rownames(lmm_fit_coef[[1]]) # adjust to your data
lmm_fit_pval_all_df

############################################### 
############################################### MWAS
###############################################
                                              
# Prepare data for MWAS
library(stringr)
dat <- cbind(meta, ds)
n_meta <- 5 # adjust to your data
dat$Label <- as.numeric(dat$Label) # numeric values, adjust to your data
dat$Label <- ifelse(dat$Label == 2, 1, 0) # binary output, adjust to your data
dat$Batch <- as.numeric(str_remove(dat$Batch, "b")) # adjust to your data
dat$Creatinine <- as.numeric(dat$Creatinine) # numeric values, adjust to your data
dat$Age <- as.numeric(dat$Age) # numeric values, adjust to your data
dat$Sex <- as.numeric(dat$Sex) # numeric values, adjust to your data

library(MWASTools)
metabolic_data = as.matrix(dat[,-c(1:n_meta)])
clinical_data = as.matrix(dat[,c(1:n_meta)])
sample_type = ifelse(ds$Label == "QC", 1, 0) # adjust to your data, experimental sample = 0, QC sample = 1
data_SE = MWAS_SummarizedExperiment(metabolic_data, clinical_data, sample_type) 

################################################### Correlation

corr_mwas <- MWAS_stats(data_SE, disease_id = "Label",   # adjust to your data
                        confounder_ids = c("Age", "Sex"), # c("...", "...", ...) or NULL
                        assoc_method = "spearman", mt_method = "BH", output = "pvalues") 
corr_mwas

################################################### LM model

lm_model <- MWAS_stats(data_SE, disease_id = "Label", # adjust to your data
                       confounder_ids = c("Age", "Sex"), # c("...", "...", ...) or NULL
                       assoc_method = "linear", mt_method = "BH", output = "pvalues")
lm_model

################################################### LM model other type

n_start <- 6 # adjust to your data
lm_fit <- lapply(n_start:ncol(dat), function(x) lm(Label ~ dat[,x] + Sex + Age, dat)) # adjust to your data 
lm_fit_coef <- lapply(1:length(lm_fit), function(x) summary(lm_fit[[x]])$coefficients)
lm_fit_pval_all_df <- as.data.frame(t(sapply(1:length(lm_fit_coef), function(x) lm_fit_coef[[x]][,4]))) 
lm_fit_pval_all_df <- as.data.frame(sapply(1:ncol(lm_fit_pval_all_df), function(x) p.adjust(lm_fit_pval_all_df[,x], method = "BH"))) # select method
dat2 <- cbind(meta, ds)
rownames(lm_fit_pval_all_df) <- colnames(dat2)[-c(1:n_meta)]
colnames(lm_fit_pval_all_df) <- rownames(summary(lm_fit[[1]])$coefficients)  # adjust to your data
lm_fit_pval_all_df

################################################### GLM model

glm_model <- MWAS_stats(data_SE, disease_id = "Label", # adjust to your data
                        confounder_ids = c("Age", "Sex"), # c("...", "...", ...) or NULL
                        assoc_method = "logistic", mt_method = "BH", output = "pvalues")
glm_model

################################################### GLM model other type

n_start <- 6 # adjust to your data
glm_fit <- lapply(n_start:ncol(dat), function(x) glm(Label ~ dat[,x] + Sex + Age, dat, family=binomial())) # adjust to your data 
glm_fit_coef <- lapply(1:length(glm_fit), function(x) summary(glm_fit[[x]])$coefficients)
glm_fit_pval_all_df <- as.data.frame(t(sapply(1:length(glm_fit_coef), function(x) glm_fit_coef[[x]][,4]))) 
glm_fit_pval_all_df <- as.data.frame(sapply(1:ncol(glm_fit_pval_all_df), function(x) p.adjust(glm_fit_pval_all_df[,x], method = "BH"))) # select method
dat2 <- cbind(meta, ds)
rownames(glm_fit_pval_all_df) <- colnames(dat2)[-c(1:n_meta)]
colnames(glm_fit_pval_all_df) <- rownames(summary(glm_fit[[1]])$coefficients) # adjust to your data
glm_fit_pval_all_df

################################################### LMM model

library(lme4)
library(lmerTest)

n_start <- 6 # adjust to your data
lmm_fit <- lapply(n_start:ncol(dat), function(x) lmer(Label ~ dat[,x] + Sex + Age + (1|Batch), dat)) # adjust to your data # lmer from lmerTest package
lmm_fit_coef <- lapply(1:length(lmm_fit), function(x) summary(lmm_fit[[x]])$coefficients)
lmm_fit_pval_all_df <- as.data.frame(t(sapply(1:length(lmm_fit_coef), function(x) lmm_fit_coef[[x]][,5]))) 
lmm_fit_pval_all_df <- as.data.frame(sapply(1:ncol(lmm_fit_pval_all_df), function(x) p.adjust(lmm_fit_pval_all_df[,x], method = "BH"))) # select method
dat2 <- cbind(meta, ds)
rownames(lmm_fit_pval_all_df) <- colnames(dat2)[-c(1:n_meta)]
colnames(lmm_fit_pval_all_df) <- rownames(summary(lmm_fit[[1]])$coefficients) # adjust to your data
lmm_fit_pval_all_df

################################################### GLMM model

library(glmmsr)

n_start <- 6 # adjust to your data
dss <- lapply(n_start:ncol(dat), function(y) as.data.frame(dat[, c(y, 4, 5, 1, 3)])) # adjust to your data
dss <- lapply(1:length(dss), function(y) {colnames(dss[[y]])[1] <-"X" 
return(dss[[y]])}) # adjust to your data
glmm_fit <- lapply(1:length(dss), function(x) glmm(Label ~ X + Sex + Age + (1|Batch), data=dss[[x]], family=binomial, method = "Laplace")) # adjust to your data
glmm_fit_coef <- lapply(1:length(glmm_fit), function(x) summary(glmm_fit[[x]])$p_value)
glmm_fit_pval_all_df <- as.data.frame(t(sapply(1:length(glmm_fit_coef), function(x) glmm_fit_coef[[x]]))) # select method
glmm_fit_pval_all_df <- as.data.frame(sapply(1:ncol(glmm_fit_pval_all_df), function(x) p.adjust(glmm_fit_pval_all_df[,x], method = "BH"))) # select method
dat2 <- cbind(meta, ds)
rownames(glmm_fit_pval_all_df) <- colnames(dat2)[-c(1:n_meta)]
colnames(glmm_fit_pval_all_df) <- colnames(summary(glmm_fit[[1]])[["fit"]][["modfr"]][["fr"]]) # adjust to your data
glmm_fit_pval_all_df

############################################### 
############################################### Other Modeling (GAM, GAMM, DRC)
###############################################

# see "Signal Modeling" for modeling 

##############################################################################################################################################################
# Signal Modeling
##############################################################################################################################################################

# setup environment
library(data.table)
setwd("D:/...")

# dataset with intensities and label column only
ds <- as.data.frame(fread(input = "8 peaks.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <-"Label"
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# METADATA
setwd("D:/...")
meta <- as.data.frame(fread(input = "metadata.csv", header=T))
rownames(meta) <- meta[,5] # adjust to your data
order <- meta[,2] # adjust to your data
meta <- meta[,c(1,3,8:10)] # adjust to your data
colnames(meta) <- c("Class", "Batch", "Creatinine", "Age", "Sex")
# check identical rownames between peak table and metadata
identical(rownames(ds), rownames(meta))

############################################### 
############################################### LM MODELING
###############################################

library(MetabolomicsBasics)

# data
meta$Batch <- as.numeric(stringr::str_remove(meta$Batch, "b")) # as.factor(meta$Batch) or as.numeric(meta$Batch) or as.numeric(stringr::str_remove(meta$Batch, "b"))
meta$Age <- as.integer(meta$Age)
meta$Sex <- as.integer(meta$Sex)
meta$Class <- as.factor(meta$Class)
s_d <- cbind(Order = order, meta)
dat <- ds[,-1]

# perform
model <- MetaboliteANOVA(dat=dat, sam=s_d, model="Batch+Class+Age+Sex", method = "BH") # "Batch+Class+Age+Sex" or "Batch+Class+Age+Sex+Order" # select method and adjust formula
model

################################################### other type of LM modeling

library(multcomp)

# data
meta$Batch <- as.numeric(stringr::str_remove(meta$Batch, "b")) # as.factor(meta$Batch) or as.numeric(meta$Batch) or as.numeric(stringr::str_remove(meta$Batch, "b"))
meta$Age <- as.integer(meta$Age)
meta$Sex <- as.integer(meta$Sex)
meta$Class <- as.factor(meta$Class)
s_d <- cbind(Order = order, meta)
dat <- ds[,-1]
ds_ancova <- as.data.frame(cbind(s_d, dat))

# perform
n_start <- 7 # adjust to your data
ancova <- lapply(n_start:ncol(ds_ancova), function(y) aov(ds_ancova[,y]~ Batch+Class+Age+Sex, ds_ancova)) # adjust formula to your data
sum_ancova <- lapply(1:length(ancova), function(y) summary(ancova[[y]]))
sum_ancova
res_ancova <- lapply(1:length(ancova), function(y) summary(glht(ancova[[y]]))) # linfct = mcp(Class = "Tukey") or adjust linfct = mcp() to your data
res_ancova

# adjusted p-value
res_ancova <- t(sapply(1:length(ancova), function(y) summary(ancova[[y]])[[1]][,5])) # adjust to your data 
res_ancova <- apply(res_ancova, 2, function(x) p.adjust(x, "BH")) # select method
colnames(res_ancova) <- rownames(summary(ancova[[1]])[[1]])
rownames(res_ancova) <- colnames(ds_ancova[, c(n_start:ncol(ds_ancova))])
res_ancova

# perform other type
n_start <- 7 # adjust to your data
ancova <- lapply(n_start:ncol(ds_ancova), function(y) anova(lm(ds_ancova[,y]~ Batch+Class+Age+Sex, ds_ancova))) # adjust formula to your data
ancova

# adjusted p-value
res_ancova <- t(sapply(1:length(ancova), function(y) ancova[[y]][,5])) # adjust to your data 
res_ancova <- apply(res_ancova, 2, function(x) p.adjust(x, "BH")) # select method
colnames(res_ancova) <- rownames(ancova[[1]])
rownames(res_ancova) <- colnames(ds_ancova[, c(n_start:ncol(ds_ancova))])
res_ancova

############################################### 
############################################### LMM MODELING
###############################################

library(lme4)
library(lmerTest)

# data
dat <- cbind(meta, ds[,-1])
n_meta <- 5 # adjust to your data
colnames(dat)[-c(1:n_meta)] <- paste0("X", c(1:ncol(dat))) # adjust to your data
dat$Batch <- as.factor(dat$Batch) # as.factor(dat$Batch) or as.numeric(dat$Batch) or as.numeric(stringr::str_remove(dat$Batch, "b"))
dat$Age <- as.integer(dat$Age)
dat$Sex <- as.integer(dat$Sex)
dat$Class <- as.factor(dat$Class)
dat[,-c(1:n_meta)] <- sapply(dat[,-c(1:n_meta)], as.numeric)

# perform
n_start <- 6 # adjust to your data
lmm_fit <- lapply(n_start:ncol(dat), function(x) lmer(dat[,x] ~ Class + Sex + Age + (1|Batch), dat)) # adjust to your data # lmer from lmerTest package
lmm_fit_coef <- lapply(1:length(lmm_fit), function(x) summary(lmm_fit[[x]])$coefficients)
lmm_fit_pval_all_df <- as.data.frame(t(sapply(1:length(lmm_fit_coef), function(x) lmm_fit_coef[[x]][,5]))) 
lmm_fit_pval_all_df <- apply(lmm_fit_pval_all_df, 2, function(x) p.adjust(x, "BH")) # select method
dat2 <- cbind(meta, ds[,-1])
rownames(lmm_fit_pval_all_df) <- colnames(dat2)[-c(1:n_meta)]
colnames(lmm_fit_pval_all_df) <- rownames(lmm_fit_coef[[1]]) # adjust to your data
lmm_fit_pval_all_df

############################################### 
############################################### GAM MODELING
###############################################

library(pbapply)
library(mgcv)

# data
dat <- cbind(meta, ds[,-1])
n_meta <- 5 # adjust to your data
colnames(dat)[-c(1:n_meta)] <- paste0("X", c(1:ncol(dat))) # adjust to your data
dat$Batch <- as.numeric(stringr::str_remove(dat$Batch, "b")) # as.factor(dat$Batch) or as.numeric(dat$Batch) or as.numeric(stringr::str_remove(dat$Batch, "b"))
dat$Age <- as.integer(dat$Age)
dat$Sex <- as.integer(dat$Sex)
dat$Class <- as.factor(dat$Class)
dat[,-c(1:n_meta)] <- sapply(dat[,-c(1:n_meta)], as.numeric)

# perform
n_start <- 6 # adjust to your data
dss <- lapply(n_start:ncol(dat), function(y) as.data.frame(dat[, c(y, 4, 5, 1, 2)])) # adjust to your data
dss <- lapply(1:length(dss), function(y) {colnames(dss[[y]])[1] <-"Y" 
                             return(dss[[y]])}) # adjust to your data
gam_fit <- pblapply(1:length(dss), function(x) mgcv::gam(Y~s(Age)+Class+Sex+Batch, data = dss[[x]])) # adjust to your data (if no s() or lo() etc. -> as lm)
gam_res <- pblapply(1:length(gam_fit), function(x) summary(gam_fit[[x]]))
gam_pval <- sapply(1:length(gam_res), function(x) gam_res[[x]][["s.pv"]]) # adjust to your data
gam_pval <- as.data.frame(p.adjust(gam_pval, "BH")) # select method
rownames(gam_pval) <- colnames(ds[,-1])
colnames(gam_pval) <- names(gam_res[[1]][["chi.sq"]])
gam_pval2 <- as.data.frame(t(sapply(1:length(gam_res), function(x) gam_res[[x]][["p.pv"]]))) # adjust to your data
gam_pval2 <- apply(gam_pval2, 2, function(x) p.adjust(x, "BH")) # select method
rownames(gam_pval2) <- colnames(ds[,-1])
gam_pval
gam_pval2

############################################### 
############################################### GAMM MODELING
###############################################

library(pbapply)
library(gamm4)

# data
dat <- cbind(meta, ds[,-1])
n_meta <- 5 # adjust to your data
colnames(dat)[-c(1:n_meta)] <- paste0("X", c(1:ncol(dat))) # adjust to your data
dat$Batch <- as.numeric(stringr::str_remove(dat$Batch, "b")) # as.factor(dat$Batch) or as.numeric(dat$Batch) or as.numeric(stringr::str_remove(dat$Batch, "b"))
dat$Age <- as.integer(dat$Age)
dat$Sex <- as.integer(dat$Sex)
dat$Class <- as.factor(dat$Class)
dat[,-c(1:n_meta)] <- sapply(dat[,-c(1:n_meta)], as.numeric)

# perform
n_start <- 6 # adjust to your data
dss <- lapply(n_start:ncol(dat), function(y) as.data.frame(dat[, c(y, 4, 5, 1, 2)])) # adjust to your data
dss <- lapply(1:length(dss), function(y) {colnames(dss[[y]])[1] <-"Y" 
                             return(dss[[y]])}) # adjust to your data
gamm_fit <- pblapply(1:length(dss), function(x) gamm4(Y~s(Age)+Class+Sex, random=~(1|Batch), data = dss[[x]])) # adjust to your data (if no s() or lo() etc. -> as lm)
gamm_res <- pblapply(1:length(gamm_fit), function(x) summary(gamm_fit[[x]]$gam))
gamm_pval <- sapply(1:length(gamm_res), function(x) gamm_res[[x]][["s.pv"]]) # adjust to your data
gamm_pval <- as.data.frame(p.adjust((gamm_pval), method = "BH")) # select method
rownames(gamm_pval) <- colnames(ds[,-1])
colnames(gamm_pval) <- names(gamm_res[[1]][["chi.sq"]])
gamm_pval2 <- as.data.frame(t(sapply(1:length(gamm_res), function(x) gamm_res[[x]][["p.pv"]]))) # adjust to your data
gamm_pval2 <- apply(gamm_pval2, 2, function(x) p.adjust(x, "BH")) # select method
rownames(gamm_pval2) <- colnames(ds[,-1])
gamm_pval
gamm_pval2

############################################### 
############################################### DRC MODELING
###############################################

library(pbapply)
library(drc)

# data
dat <- cbind(meta, ds[,-1])
n_meta <- 5 # adjust to your data
colnames(dat)[-c(1:n_meta)] <- paste0("X", c(1:ncol(dat))) # adjust to your data
dat$Batch <- as.numeric(stringr::str_remove(dat$Batch, "b")) # as.factor(dat$Batch) or as.numeric(dat$Batch) or as.numeric(stringr::str_remove(dat$Batch, "b"))
dat$Age <- as.integer(dat$Age)
dat$Sex <- as.integer(dat$Sex)
dat$Class <- as.factor(dat$Class)
dat[,-c(1:n_meta)] <- sapply(dat[,-c(1:n_meta)], as.numeric)

# perform
n_start <- 6 # adjust to your data
dss <- lapply(n_start:ncol(dat), function(y) as.data.frame(dat[, c(y, 4, 5, 1, 2)])) # adjust to your data
dss <- lapply(1:length(dss), function(y) {colnames(dss[[y]])[1] <-"Y" 
                             return(dss[[y]])}) # adjust to your data
drm_fit <- pblapply(1:length(dss), function(x) drm(Y~Age+Class+Sex, data = dss[[x]], fct = LL.4())) # adjust to your data 
drm_res <- pblapply(1:length(drm_fit), function(x) summary(drm_fit[[x]]))
drm_pval <- as.data.frame(t(sapply(1:length(drm_res), function(x) drm_res[[x]][["coefficients"]][,4]))) # adjust to your data
drm_pval <- apply(drm_pval, 2, function(x) p.adjust(x, method = "BH")) # select method
rownames(drm_pval) <- colnames(ds[,-1])
drm_pval

############################################### 
############################################### Other Modeling (Signal as variable)
###############################################

# see "MWAS/ANCOVA" for modeling  

##############################################################################################################################################################
# N-Factor Analysis
##############################################################################################################################################################

# setup environment
library(data.table)
library(dplyr)
library(stringr)
setwd("D:/...")

# dataset with intensities and label column only
ds <- as.data.frame(fread(input = "8 peaks.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <-"Label"
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# METADATA
setwd("D:/...")
meta <- as.data.frame(fread(input = "metadata.csv", header=T))
rownames(meta) <- meta[,5] # adjust to your data
order <- meta[,2] # adjust to your data
meta <- meta[,c(3,8:10)] # adjust to your data
colnames(meta) <- c("Batch", "Creatinine", "Age", "Sex")
# check identical rownames between peak table and metadata
identical(rownames(ds), rownames(meta))

############################################### 
############################################### Other Modeling (LM, LMM, GAM, GAMM, DRC)
###############################################

# see "Signal Modeling" for modeling  

############################################### 
############################################### Two-way ANOVA
###############################################

# see"Two-way ANOVA" for N-way ANOVA in "Testing set of biomarkers" 

############################################### 
############################################### ASCA
###############################################

library(MetStaT)

# data
dat <- as.matrix(ds[,-1])
s_d <- cbind(Class = ds$Label, Order = order, meta)
factors <- cbind(as.numeric(s_d$Class),as.numeric(s_d$Sex)) # adjust to your data

# perform
ASCA <- ASCA.Calculate(dat, factors, scaling = F)
ASCA.Plot(ASCA)
ASCA.DoPermutationTest(ASCA, perm = 500)

############################################### 
############################################### PLS/sPLS
###############################################

library(mixOmics)

# data
dat <- ds[,-1]
s_d <- cbind(Class = ds$Label, Order = order, meta) # adjust to your data
Y <- as.matrix(cbind(Class = as.factor(s_d$Class), Sex = as.factor(s_d$Sex)))

################################################### PLS

pls.multilabel <- pls(X=dat, Y = Y, 
                            ncomp = 3) # adjust to your data

pls.multilabel
plotIndiv(pls.multilabel, group = as.factor(s_d$Class))

Q2.pls1 <- perf(pls.multilabel, validation = 'Mfold', folds = 10, nrepeat = 5)
plot(Q2.pls1, criterion = 'Q2')

################################################### sPLS

tune.spls <- tune.spls(X=dat, Y=Y, 
                           ncomp=3, # adjust to your data
                           validation = 'Mfold', folds = 10, nrepeat = 5, # adjust to your data
                           test.keepX = c(1:ncol(dat)),
                           progressBar = T)

tune.spls
tune.spls$choice.ncomp # selected ncomp
plot(tune.spls) 

spls.multilabel <- spls(X=dat,
                        Y = Y, 
                        ncomp = tune.spls$choice.ncomp) # adjust to your data tune.spls$choice.ncomp or N

spls.multilabel
plotIndiv(spls.multilabel, group = as.factor(s_d$Class))

Q2.spls1 <- perf(spls.multilabel, validation = 'Mfold', folds = 10, nrepeat = 5)
plot(Q2.spls1, criterion = 'Q2')
                                 
############################################### 
############################################### PVCA
###############################################

library(proBatch)

# generate batch data
batch <- meta$Batch
batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# generate group data
class <- as.character(ds$Label)
qc_id <- which(class == "QC")
class[-qc_id] <- "Subject"

# Sample data
s_d <- data.frame(FullRunName = rownames(ds), batch, class)  # adjust to your data
rownames(s_d) <- NULL

# Feature data
f_d <- as.matrix(t(ds[,-1]))
colnames(f_d) <- rownames(ds)

# Perform
pvca_df <- calculate_PVCA(f_d, s_d, 
                          factors_for_PVCA = c('batch', 'class'),  # adjust to your data
                          pca_threshold = .6, variance_threshold = .01, fill_the_missing = 0)

pvca_df

############################################### 
############################################### PC-PR2
###############################################

library(pcpr2)

# generate batch data
batch <- meta$Batch
batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# generate group data
class <- as.character(ds$Label)
qc_id <- which(class == "QC")
class[-qc_id] <- "Subject"

# Sample data
s_d <- data.frame(batch, class)

# Feature data
f_d <- as.matrix(ds[,-1])

# Perform
pct.threshold <- 0.8 # adjust to your data
PCPR2 <- runPCPR2(X = f_d, Z = s_d, pct.threshold = pct.threshold) 

PCPR2

############################################### 
############################################### TDFDR
###############################################

library(tdfdr)

# data
dat <- ds[,-1]
s_d <- cbind(Class = ds$Label, Order = order, meta) # adjust to your data
cutoff <- 0.05 # adjust to your data

# perform
tdfdr_f <- tdfdr(dat, s_d$Class, s_d[,c(3,5,6)], alpha = cutoff) # adjust to your data
tdfdr_f_sel <- tdfdr.select(tdfdr_f, fdr.level = cutoff) # adjust to your data
tdfdr_f_sel
pl <- tdfdr.plot(tdfdr_f)
pl

##############################################################################################################################################################
# Repeated Measures
##############################################################################################################################################################

# setup environment
library(data.table)
setwd("D:/...")

# dataset with intensities and label column only
ds <- as.data.frame(fread(input = "8 peaks.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <-"Label"
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# METADATA
setwd("D:/...")
meta <- as.data.frame(fread(input = "metadata.csv", header=T))
rownames(meta) <- meta[,5] # adjust to your data
order <- meta[,2] # adjust to your data
meta <- meta[,c(3,8:10)] # adjust to your data
colnames(meta) <- c("Batch", "Creatinine", "Age", "Sex")
# check identical rownames between peak table and metadata
identical(rownames(ds), rownames(meta))

############################################### 
############################################### LM MODELING
###############################################

# data
dat <- ds[,-1]
s_d <- cbind(Class = ds$Label, Order = order, meta)
ds_rep <- as.data.frame(cbind(s_d, dat))
id <- as.numeric(stringr::str_remove(ds_rep$Batch, "b")) # define repeats -> sample id, this is just toy example
ds_rep <- as.data.frame(cbind(id = id, ds_rep))

# perform
n_start <- 9 # adjust to your data
rep_aov <- lapply(n_start:ncol(ds_rep), function(y) aov(ds_rep[,y]~ Class+Age+Sex+Error(id), ds_rep)) # adjust formula to your data, Error(id) -> for repeated measures
sum_rep_aov <- lapply(1:length(rep_aov), function(y) summary(rep_aov[[y]]))
sum_rep_aov

# adjusted p-value
pval_rep_aov <- t(sapply(1:length(rep_aov), function(y) summary(rep_aov[[y]])[["Error: Within"]][[1]][,5])) # adjust  to your data
pval_rep_aov <- apply(pval_rep_aov, 2, function(x) p.adjust(x, "BH")) # select method
rownames(pval_rep_aov) <- colnames(ds_rep[,c(n_start:ncol(ds_rep))])
colnames(pval_rep_aov) <- rownames(summary(rep_aov[[1]])[["Error: Within"]][[1]])
pval_rep_aov

############################################### 
############################################### Other Modeling (LM, LMM, GAM, GAMM, DRC)
###############################################

# see "Signal Modeling" for modeling 

############################################### 
############################################### multilevel sPLS for repeated measure
###############################################

library(mixOmics)

# data
dat <- as.matrix(ds[,-1])
s_d <- cbind(Class = ds$Label, Order = order, meta)
Y <- data.frame(Class = as.factor(s_d$Class), Sex = as.factor(s_d$Sex)) # adjust to your data
id <- as.numeric(stringr::str_remove(s_d$Batch, "b")) # define repeats -> sample id. This just example

################################################### Two factor analysis

tune.splsda <- tune.splsda(X=dat, Y=Y, 
                              ncomp=3, # adjust to your data
                              multilevel = id, # use for repeated measure 
                              dist = 'mahalanobis.dist',
                              validation = 'Mfold', folds = 10, nrepeat = 5, # adjust to your data
                              test.keepX = c(1:ncol(dat)),
                              progressBar = T)

tune.splsda # select ncomp

splsda.multilevel <- splsda(X=dat,
                            Y = Y, 
                            multilevel = id, # use for repeated measure
                            ncomp = 3) # adjust to your data, select ncomp

splsda.multilevel

plotIndiv(splsda.multilevel, ind.names = id, group = as.factor(s_d$Class))

################################################### One factor analysis

tune.splsda <- tune.splsda(X=dat, Y=as.factor(s_d$Class),
                           ncomp=3, # adjust to your data
                           multilevel = id, # use for repeated measure 
                           dist = 'mahalanobis.dist',
                           validation = 'Mfold', folds = 10, nrepeat = 5, # adjust to your data
                           test.keepX = c(1:ncol(dat)),
                           progressBar = T)

tune.splsda # select ncomp

tune.splsda$choice.ncomp$ncomp
tune.splsda$choice.keepX
plot(tune.splsda)

splsda.multilevel <- splsda(X=dat,
                            Y = as.factor(s_d$Class), 
                            multilevel = id, # use for repeated measure
                            ncomp = 3) # adjust to your data tune.splsda$choice.ncomp$ncomp or N

splsda.multilevel

plotIndiv(splsda.multilevel, ind.names = id, group = as.factor(s_d$Class))

##############################################################################################################################################################
# Time series
##############################################################################################################################################################

# setup environment
library(data.table)
setwd("D:/...")

# dataset with intensities and label column only
ds <- as.data.frame(fread(input = "8 peaks.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <-"Label"
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# METADATA
setwd("D:/...")
meta <- as.data.frame(fread(input = "metadata.csv", header=T))
rownames(meta) <- meta[,5] # adjust to your data
order <- meta[,2] # adjust to your data
meta <- meta[,c(3,8:10)] # adjust to your data
colnames(meta) <- c("Batch", "Creatinine", "Age", "Sex")
# check identical rownames between peak table and metadata
identical(rownames(ds), rownames(meta))

############################################### 
############################################### Other Modeling (LM, LMM, GAM, GAMM, DRC)
###############################################

# see "Signal Modeling" for modeling  

############################################### 
############################################### ASCA, PVCA, PC-PR2
###############################################

# see "N-Factor Analysis" for ASCA, PVCA, PC-PR2

############################################### 
############################################### Multivariate Empirical Bayes Statistics
###############################################

library(timecourse)

# data
dat <- as.matrix(ds[,-1])
s_d <- cbind(Class = ds$Label, Order = order, meta)
ds_tc <- as.data.frame(cbind(s_d, dat))
ds_tc <- ds_tc[order(ds_tc$Class, decreasing = F),] # order by group variable
dat <- t(ds_tc[,-c(1:ncol(s_d))])
tp <- 4 # define time point, this is just toy example
assay <- rep(c(1:(ncol(dat)/tp)), each = tp) # define repeats -> sample id, this is just toy example. This 4 time points and 31 repeats
trt <- s_d$Class
l <- levels(as.factor(trt))
reps <- sapply(1:length(l), function(y) length(trt[trt==l[y]])/tp)

################################################# with condition group
size <- matrix(reps, nrow=nrow(dat), ncol=length(l), byrow=T)
long <- mb.long(dat, method="2", times=tp, reps=size, rep.grp=assay, condition.grp=trt)
ht2 <- long$HotellingT2 # Hotelling T2
ht2

plotProfile(long, type="b", ranking = 1) # plot 1 feature by rank
names <- as.character(1:nrow(dat))
plotProfile(long, gid="8", type="b", gnames=names)

################################################# no condition group
size <- matrix(length(unique(assay)), nrow=nrow(dat), ncol=1, byrow=T)
long <- mb.long(dat, method="1", times=tp, reps=size, rep.grp=assay)
ht2 <- as.data.frame(long$HotellingT2) # Hotelling T2

# plot
plotProfile(long, type="b", ranking = 1) # plot 1 feature by rank
names <- as.character(1:nrow(dat))
plotProfile(long, gid="8", type="b", gnames=names)

############################################### 
############################################### Dose-Response Modeling (DRomics)
###############################################

library(DRomics)
library(NormalizeMets)

# data
dat <- cbind(Time = meta$Batch, ds[,-1]) # adjust to your data # toy example
dat$Time <- as.numeric(stringr::str_remove(dat$Time, "b")) # as.factor(dat$Batch) or as.numeric(dat$Batch) or as.numeric(stringr::str_remove(dat$Batch, "b"))

######################### Perform for log data
dat1 <- dat
dat1[,-1] <- LogTransform(dat1[,-1])$featuredata
ds_t <- as.data.frame(t(dat1))
ds_t <- as.data.frame(sapply(ds_t, as.numeric))
ds_t <- as.data.frame(cbind(colnames(dat), ds_t))
o <- continuousomicdata(ds_t, backgrounddose = dat$Time[1]) # adjust "backgrounddose" argument to your data
s_mod <- itemselect(o, select.method = "linear", FDR = 0.05) # adjust to your data "quadratic", "linear", "ANOVA"
f <- drcfit(s_mod, progressbar = T)
res_dr <- f$fitres
res_dr
plot(f)

######################### Perform for raw data
ds_t1 <- as.data.frame(t(dat))
ds_t1 <- as.data.frame(sapply(ds_t1, as.numeric))
ds_t1 <- as.data.frame(cbind(colnames(dat), ds_t1))
o1 <- continuousomicdata(ds_t1, backgrounddose = dat$Time[1]) # adjust "backgrounddose" argument to your data
s_mod1 <- itemselect(o1, select.method = "linear", FDR = 0.05) # adjust to your data "quadratic", "linear", "ANOVA"
f1 <- drcfit(s_mod1, progressbar = T)
res_dr1 <- f1$fitres
res_dr1
plot(f1)

# Try also parallel version of "drcfit":
# library(parallel)
# library(doParallel)
# f <- drcfit(s_mod, progressbar = T, parallel = "snow", ncpus = 3)  # adjust ncores (ncpus) to your PC            
               
######################### Calculation of BMD
r <- bmdcalc(f, z = 1, x = 10)
r$res
plot(r, BMDtype = "zSD", plottype = "ecdf")
plot(r, BMDtype = "zSD", plottype = "density", by = "trend")
bmdplotwithgradient(r$res, BMDtype = "zSD", facetby = "trend", shapeby = "model", line.size = 1.2) 

######################### Bootstrap calculation
b <- bmdboot(r, niter = 50, progressbar = F)
b$res
plot(b, BMDtype = "zSD", by = "trend")

############################################### 
############################################### Dose-Response Modeling (TOXcms)
###############################################

library(toxcms)
library(stringr)
library(batchCorr)

# data
dat <- ds[,-1] # adjust to your data # toy example
time <- as.numeric(stringr::str_remove(meta$Batch, "b")) # as.factor(dat$Batch) or as.numeric(dat$Batch) or as.numeric(stringr::str_remove(dat$Batch, "b"))
dat <- as.data.frame(t(dat)) 
colnames(dat) <- paste(colnames(dat), sep = "_", time, "time") # add time index, adjust to your data
level_time <- paste0("_", unique(time), "_time") 
rn <- rownames(dat) # add mz and rt
rn <- gsub(pattern = " / ", replacement = "_", x = rn, fixed = T)
rn2 <- t(data.frame(rn))
colnames(rn2) <- rn
peakIn <- peakInfo(PT = rn2, sep = "_", start = 1)
dat <- as.data.frame(cbind(peakIn, dat))
name <- colnames(ds[,-1]) # add names and index
X <- c(1:nrow(dat))
dat <- as.data.frame(cbind(X, name, dat))
rownames(dat) <- colnames(ds[,-1])

# perform
etom_dosestat <- calcdosestat(Feature = dat, Dose_Levels = level_time, multicomp = "none", p.adjust.method = "none",projectName = "dataset") # adjust to your data
etom_dosestat$pvalue
# mono trend
etom_drreport_mono <- trendfilter(etom_dosestat, pval_cutoff = 0.05, pval_thres = 1, anova_cutoff = 0.05, trend = "mono", relChange_cutoff = 50, export = F) # adjust to your data trend = c("increase","decrease","mono","reverse","all"), pval, etc.
etom_drreport_mono$pvalue
etom_drreport_fit <- fitdrc(DoseResponse_report=etom_drreport_mono, Dose_values=unique(time), ED=0.5, export = T, mz_tag = "mz", rt_tag = "rt", plot=T) # adjust to your data
etom_drreport_clust <- clusttrend(etom_drreport_mono, reference_index = NULL, sort.method =c("clust","layer"), sort.thres = 20, dist.method = "euclidean", hclust.method = "average", mztag = "mz", rttag = "rt", heatmap.on = T, plot.all = T, filename = "testdataset_hclust.pdf") # adjust to your data
plotpca(DoseResponse_report = etom_drreport_fit, DoseStat = etom_dosestat, EDrange = c(0,max(time))) # adjust to your data
plottrend(etom_drreport_mono, Dose_conditions = level_time, y_transform = T, mz_tag = "mz", rt_tag = "rt") # adjust to your data # set file format .pdf in file name in the folder (working directory by default)
# reverse trend
etom_drreport_reverse <- trendfilter(etom_dosestat, pval_cutoff = 0.05, pval_thres = 1, anova_cutoff = 0.05, trend = "reverse", relChange_cutoff = 0.05, export = T) # adjust to your data trend = c("increase","decrease","mono","reverse","all"), pval, etc.
etom_drreport_reverse$pvalue
plottrend(etom_drreport_reverse, Dose_conditions = level_time, y_transform = T, mz_tag = "mz", rt_tag = "rt") # adjust to your data # set file format .pdf in file name in the folder (working directory by default)

############################################### 
############################################### LMM Splines Modeling (timeOmics)
###############################################

library(timeOmics)
library(lmms)
library(parallel)
library(tidyverse)
library(ggplot2)
library(dplyr)

# data
dat <- cbind(ds[,-1]) # adjust to your data # toy example
time <- as.numeric(stringr::str_remove(meta$Batch, "b")) # toy example # as.factor(dat$Batch) or as.numeric(dat$Batch) or as.numeric(stringr::str_remove(dat$Batch, "b"))

# perform
lmms.output <- lmms::lmmSpline(data = dat, time = time, # adjust to your
                               sampleID = rownames(dat), deri = FALSE, # adjust to your 
                               basis = "p-spline", numCores = (detectCores(logical = F)-1), timePredict = 1:max(time), # adjust to your
                               keepModels = TRUE) # adjust to your

modelled.data <- t(slot(lmms.output, 'predSpline')) # modeling results

# gather data
data.gathered <- modelled.data %>% as.data.frame() %>% 
  rownames_to_column("time") %>%
  mutate(time = as.numeric(time)) %>%
  pivot_longer(names_to="feature", values_to = 'value', -time)

# plot profiles
ggplot(data.gathered, aes(x = time, y = value, color = feature)) + geom_line() +
  theme_bw() + ggtitle("`lmms` profiles") + ylab("Feature expression") +
  xlab("Time")

# filter results
filter.res <- lmms.filter.lines(data = dat, lmms.obj = lmms.output, time = time, homoskedasticity.cutoff = 0.05) # adjust to your
profile.filtered <- filter.res$filtered
colnames(profile.filtered)

# PCA
pca.res <- mixOmics::pca(X = profile.filtered, ncomp = 2, scale = F, center=F) # adjust to your
pca.ncomp <- getNcomp(pca.res, max.ncomp = 5, X = profile.filtered, 
                      scale = FALSE, center=FALSE) # adjust to your data
pca.ncomp$choice.ncomp # optimal comp
plotIndiv(pca.res)
plotLong(pca.res, scale = F, center = F, title = "PCA longitudinal clustering")

# sPCA
tune.spca.res <- tuneCluster.spca(X = profile.filtered, ncomp = 2, test.keepX = c(2:ncol(profile.filtered))) # adjust to your data
spca.res <- spca(X = profile.filtered, ncomp = 2, keepX = c(1,2), scale = F) # keepX = tune.spca.res$choice.keepX or c(n, m), adjust to your data
plotLong(spca.res, title = "s-PLS longitudinal clustering")

# PLS
X <- profile.filtered
Y <- cbind(unique(time), unique(time)) # adjust to your data, just toy example
pls.res <- pls(X,Y, ncomp = 5, scale = FALSE)
pls.ncomp <- getNcomp(pls.res, max.ncomp = 5, X=X, Y=Y, scale = FALSE) # adjust to your data
pls.ncomp$choice.ncomp # optimal comp
pls.res <- pls(X,Y, ncomp = pls.ncomp$choice.ncomp, scale = FALSE)
head(getCluster(pls.res))
plotLong(pls.res, title = "PLS longitudinal clustering", legend = TRUE)

# sPLS
X <- profile.filtered
Y <- cbind(unique(time), unique(time)) # adjust to your data, just toy example
tune.spls <- tuneCluster.spls(X, Y, ncomp = 2, test.keepX = c(4:7), test.keepY <- c(1,2)) # adjust to your data
spls.res <- spls(X,Y, ncomp = 2, scale = FALSE, 
                 keepX = tune.spls$choice.keepX, keepY = tune.spls$choice.keepY)
spls.cluster <- getCluster(spls.res)
spls.cluster
plotLong(spls.res, title = "sPLS clustering")

# Multiblock PLS
X <- list("X" = profile.filtered, "Z" = meta[1:nrow(profile.filtered),2:4]) # adjust to your data, just toy example
rownames(X$Z) <- rownames(X$X) # identical names
Y <- as.matrix(cbind(unique(time), unique(time))) # adjust to your data, just toy example
block.pls.res <- block.pls(X=X, Y=Y, ncomp = 5, scale = FALSE, mode = "canonical")
block.ncomp <- getNcomp(block.pls.res,X=X, Y=Y, 
                        scale = FALSE, mode = "canonical", max.ncomp = 4)
block.ncomp$choice.ncomp
block.pls.res <- block.pls(X=X, Y=Y, ncomp = 1, scale = FALSE, mode = "canonical")
block.pls.cluster <- getCluster(block.pls.res)
block.pls.cluster
plotLong(block.pls.res)

# Multiblock sPLS
X <- list("X" = profile.filtered, "Z" = meta[1:nrow(profile.filtered),2:4]) # adjust to your data, just toy example
rownames(X$Z) <- rownames(X$X) # identical names
Y <- as.matrix(cbind(unique(time), unique(time))) # adjust to your data, just toy example
test.list.keepX <- list("X" = 1:2, "Z" = c(1,3)) # adjust to your data, just toy example
test.keepY <- c(1,2) # adjust to your data, just toy example
tune.block.res <- tuneCluster.block.spls(X= X, Y= Y, 
                                         test.list.keepX=test.list.keepX, 
                                         test.keepY= test.keepY, 
                                         scale=FALSE, 
                                         mode = "canonical", ncomp = 1)
block.pls.res <- block.spls(X=X, Y=Y, 
                            ncomp = 1, 
                            scale = FALSE, 
                            mode = "canonical", 
                            keepX = tune.block.res$choice.keepX, 
                            keepY = tune.block.res$choice.keepY)

head(getCluster(block.pls.res))
plotLong(block.pls.res)

res <- proportionality(block.pls.res)
pval.propr <- res$pvalue
pval.propr
plot(res)

############################################### 
############################################### Pharmacokinetics (polyPK)
###############################################

library(polyPK)
library(batchCorr)
library(stringr)
library(dplyr)
library(openxlsx)

# data
cn <- colnames(ds[,-1]) # adjust to your data # toy example
cn <- gsub(pattern = " / ", replacement = "_", x = cn, fixed = T)
cn2 <- t(data.frame(cn))
colnames(cn2) <- cn
peakIn <- peakInfo(PT = cn2, sep = "_", start = 1)
data_pk <- as.data.frame(t(ds)) # 0 values confused algorithm, replace by smth
group <- as.numeric(str_remove(meta$Batch, "b")) # adjust to your data, just toy example, should be from 0 to N in order (0,1,2,...N)
group <- group-1 # adjust to your data, just toy example, should be from 0 to N in order (0,1,2,...N)
rownames(data_pk)[1] <- "timepoints(h)" # adjust to your data, just toy example
data_pk[1,] <- group # adjust to your data, just toy example
gender <- meta$Sex # adjust to your data, just toy example
gender[which(gender==0)] <- 2 # adjust to your data, just toy example
data_pk <- rbind(gender=gender, group=group, data_pk)
ID <- c(NA,NA,NA, colnames(ds[,-1]))
m.z <- c(NA,NA,NA, as.numeric(peakIn[,1]))
RT.min. <- c(NA,NA,NA, as.numeric(peakIn[,2]))
data_pk <- cbind(ID=ID, m.z=m.z, RT.min.=RT.min., data_pk)
Name <- rownames(data_pk)
data_pk <- cbind(Name=Name, data_pk)

# perform
pks <- PKs(data_pk, d.point="mean", d.ebar="SE", filepath=getwd(), design=FALSE)
pk_res <- paste0(getwd(), "/PKresults/PKresults(all)/PK-parameters.xlsx")
pk_res_df <- read.xlsx(pk_res, sheet = 1, skipEmptyRows = F, colNames = T)
pk_res_df

##############################################################################################################################################################
# Unsupervised Data Projection
##############################################################################################################################################################

# Content:
# PCA
# HCA
# KM
# HCA on PCA
# Heatmap
# t-SNE
# Other types of clustering
# Validation Clustering

# LOAD DATA
library(data.table)
setwd("D:/...")

# PEAK TABLE
# dataset with intensities and label column only
ds <- as.data.frame(fread(input = "8 peaks.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <-"Label"
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# PEAK TABLE
# dataset with intensities, annotation and label column
ds <- as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats annot+filtr LMM adj KEGG.csv", header=T))
ds <- ds[-c(1:12),] # adjust to your data
rownames(ds) <- ds[,5] # adjust to your data
ds <- ds[,-c(1,3:5)] # adjust to your data
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# METADATA
meta <- as.data.frame(fread(input = "metadata.csv", header=T))
rownames(meta) <- meta[,5] # adjust to your data
order <- meta[,2] # adjust to your data
meta <- meta[,c(3,8:10)] # adjust to your data
colnames(meta) <- c("Batch", "Creatinine", "Age", "Sex")
# check identical rownames between peak table and metadata
identical(rownames(ds), rownames(meta))

# Settings for projection
library(factoextra)
library(FactoMineR)
library(pcaMethods)
library(dendextend)
library(rafalib)
library(RSEIS)
library(ggsci)
library(pheatmap)
library(Rtsne)
library(dbscan)
library(umap)
library(NMF)
library(Rdimtools)
library(dimRed)
library(cluster)
library(NbClust)
library(clustertend)
library(mclust)
library(clValid)
library(fpc)
library(pvclust)
library(parallel)
library(doParallel)

# dataset
base1 <- ds # dataset
mtrx1 <- ds[,-1] # numeric data
grp1 <- as.character(base1[,1]) # label of dataset
k <- length(unique(grp1)) # groups in clustering

###############################################
############################################### PRINCIPAL COMPONENT ANALYSIS
############################################### 

palette_pca <- pal_gsea("default", n = length(unique(grp1)), alpha = 0.6, reverse = T)(length(unique(grp1))) # color: "lancet" or "category20" or pal_gsea("default", n = length(unique(grp1)), alpha = 0.6, reverse = T)(length(unique(grp1)))
# palette_pca <- JGRAY(length(unique(grp1))) # grey

n = 5 # number of PCs, adjust to your data
pca.ds1 <- FactoMineR::PCA(mtrx1, scale.unit = T, graph = F, ncp = n)
pca <- fviz_pca_ind(pca.ds1,
                    title = "",
                    geom.ind = "point", # show points only 
                    col.ind = grp1, # color by groups
                    palette = palette_pca, # color "jco" , "lancet", gray JGRAY(length(unique(grp1)))
                    addEllipses = T, # Concentration ellipses
                    ellipse.level = 0.95, # size of the ellipse in Normal probability
                    legend.title = "Groups")
pca # add: "+scale_shape_manual(values=rep(0:length(unique(grp1))))" if shape palette error 

# Scree plot
fviz_eig(pca.ds1)

# Variables plot
fviz_pca_var(pca.ds1, col.var = "black")

# Loadings 
pc <- pcaMethods::pca(mtrx1, nPcs=n)
pc@loadings

###############################################
############################################### HIERARCHICAL CLUSTER ANALYSIS
############################################### 

# number of groups
k <- length(unique(grp1)) # groups in HC

# color
Cols = function(vec, ord){
  cols = pal_lancet(palette = c("lanonc"), alpha = 1)(length(unique(vec))) # or other palette from ggsci
  return(cols[as.fumeric(vec)[ord]])}

# grey
#Cols = function(vec, ord){
# cols = JGRAY(length(unique(vec)))
# return(cols[as.fumeric(vec)[ord]])}

mtrx1_1 <- mtrx1
#mtrx1_1 <- data.frame(scale(mtrx1, center = T, scale = T))
rownames(mtrx1_1) = make.names(grp1, unique=TRUE)
res.dist1 <- dist(mtrx1_1, method = "manhattan") #{euclidean}, {maximum}, {manhattan}, {canberra}, {binary}, {minkowski}
res.hc1 <- hclust(d = res.dist1, method = "ward.D2") #{ward (ward.D), (ward.D2)}, {single}, {complete}, {average}, {mcquitty},{median}, {centroid}
hca <- fviz_dend(res.hc1, k = k, # Cut in k groups
                 cex = 0.3, # label size 0.3/0.7
                 k_colors = unique(Cols(grp1,res.hc1$order)), # "lancet" color "jco" gray JGRAY(k_hc)
                 color_labels_by_k = F, # color labels by groups
                 label_cols = Cols(grp1,res.hc1$order),#Cols(ds[,1])[res.hc1$order], #as.fumeric(ds[,1])[res.hc1$order]
                 rect = T, # Add rectangle around groups
                 rect_fill = T,
                 rect_border = unique(Cols(grp1,res.hc1$order)), #"lancet"# color "jco" gray JGRAY(k_hc)
                 horiz = F,
                 lwd = 0.3, # lines size 0.3/0.7
                 show_labels = T,
                 main = "",
                 ylab = "")
hca

########################################### other HCA

# number of groups
k <- length(unique(grp1)) # groups in HC

# color
Cols = function(vec, ord){
  cols = pal_lancet(palette = c("lanonc"), alpha = 1)(length(unique(vec))) # or other palette from ggsci
  return(cols[as.fumeric(vec)[ord]])}

# grey
#Cols = function(vec, ord){
# cols = JGRAY(length(unique(vec)))
# return(cols[as.fumeric(vec)[ord]])}

mtrx <- mtrx1
#mtrx <- data.frame(scale(mtrx1, center = T, scale = T))
rownames(mtrx) = make.names(grp1, unique=TRUE)
res.dist1 <- dist(mtrx, method = "manhattan")
res.hc1 <- hclust(d = res.dist1, method = "ward.D2")
dend1 <- as.dendrogram(res.hc1, hang =70)
dend1 <- color_branches(dend1, k=k, groupLabels = F, col =unique(Cols(grp1,res.hc1$order))) 
labels_colors(dend1) <-Cols(grp1,res.hc1$order)
dend1 <- assign_values_to_leaves_nodePar(dend1, 0.6, "lab.cex")
dend1 <- dendextend::set(dend1, "branches_lwd", 2.5)
plot(dend1)
dend1 <- rect.dendrogram(dend1, k=k, border = 1, lty = 1, lwd = 1, col=rgb(0.1, 0.2, 0.4, 0.1))
legend("topright", legend = unique(grp1), fill = pal_lancet(palette = c("lanonc"), alpha = 1)(length(unique(grp1))))


########################################### other HCA

# number of groups
k <- length(unique(grp1)) # groups in HCA

grp1_1 <- as.numeric(as.factor(grp1))
dend1 <- as.dendrogram(res.hc1, hang =70)
k_1 <- cutree(dend1, k = k, order_clusters_as_data = F)
bars <- cbind(grp1_1, k_1)
bars <- as.data.frame(cbind(hcl.colors(k, palette = "ag_Sunset")[grp1_1], hcl.colors(k, palette = "ag_Sunset")[k_1])) # or use "hcl.colors" function with different "palette" argument from hcl.pals()
colnames(bars) <- c("Label", "HCA")

dend1 %>% set("labels", "") %>% plot
colored_bars(colors = bars, dend = dend1, sort_by_labels_order = F,
             y_shift = -0.5, y_scale = 2.5, text_shift = -0.1, cex = 1.2)
               
###############################################
############################################### k-MEANS CLUSTERING
############################################### 

k <- length(unique(grp1)) # groups in KM
km.res1 <- kmeans(mtrx1, centers = k, nstart = 25)
km.res1$cluster # cluster membership
fviz_cluster(list(data = mtrx1, cluster = km.res1$cluster), repel = T,
             ellipse.type = "euclid", geom = "point", stand = FALSE,
             palette = "jco", ggtheme = theme_classic()) # or other palette from ggsci

# Hierarchical K-Means Clustering
res.hk <- hkmeans(mtrx1, k)
fviz_dend(res.hk, cex = 0.6, palette = "jco", # or other palette from ggsci
          rect = TRUE, rect_border = "jco", rect_fill = TRUE)
hkmeans_tree(res.hk, cex = 0.6)

###############################################
############################################### HCA ON PCA
###############################################

k <- length(unique(grp1)) # groups in HC
res.pca <- FactoMineR::PCA(mtrx1, ncp = 3, graph = F)
res.hcpc <- HCPC(res.pca, nb.clust = k, graph = F)

fviz_dend(res.hcpc,
          cex = 0.7, # Label size
          palette = "jco", # Color palette see ?ggpubr::ggpar or other palette from ggsci
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco", # Rectangle color
          labels_track_height = 0.8 # Augment the room for labels
)

fviz_cluster(res.hcpc,
             repel = TRUE, # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco", # Color palette see ?ggpubr::ggpar or other palette from ggsci
             ggtheme = theme_minimal(),
             main = "PCA")

plot(res.hcpc, choice = "3D.map")

###############################################
############################################### HEATMAP
###############################################

# Correlation Distance by sample
rows.cor <- cor(t(mtrx1), use = "pairwise.complete.obs", method = "pearson") # use any: "pearson", "kendall", "spearman"

# Correlation Distance by feature
cols.cor <- cor(mtrx1, use = "pairwise.complete.obs", method = "pearson") # use any: "pearson", "kendall", "spearman"

# Plot the heatmap of correlation
pheatmap(mtrx1, 
         clustering_distance_cols = as.dist(1 - cols.cor),
         clustering_distance_rows = as.dist(1 - rows.cor)
)

# Plot the heatmap by sample
pheatmap(mtrx1, clustering_rows = T, cluster_cols = F,
         clustering_distance_rows = "manhattan", clustering_distance_cols	= "manhattan", # use any: {euclidean}, {maximum}, {manhattan}, {canberra}, {binary}, {minkowski}
         clustering_method = "ward.D2") # use any: {ward (ward.D), (ward.D2)}, {single}, {complete}, {average}, {mcquitty},{median}, {centroid}

# Plot the heatmap by feature
pheatmap(t(mtrx1), clustering_rows = T, cluster_cols = F,
         clustering_distance_rows = "manhattan", clustering_distance_cols	= "manhattan", # use any: {euclidean}, {maximum}, {manhattan}, {canberra}, {binary}, {minkowski}
         clustering_method = "ward.D2") # use any: {ward (ward.D), (ward.D2)}, {single}, {complete}, {average}, {mcquitty},{median}, {centroid}

###############################################
############################################### t-DISTRIBUTED STOCHASTIC NEIGHBOR EMBEDDING
############################################### 

# in Rtsne function matrix or dist(matrix or df) can be used and pca = T or F, perplexity number should be checked
set.seed(1234)
ds_ul_tsne <- as.matrix(unique(mtrx1))
tsne_out <- Rtsne(ds_ul_tsne, pca = T, perplexity = 10)
Cols = function(vec){
  cols = rainbow(length(unique(vec))) # or use "hcl.colors" function with different "palette" argument from hcl.pals()
  return(cols[as.numeric(as.factor(vec))])}
plot(tsne_out$Y, col = Cols(grp1),pch = 19, xlab = "", ylab = "", main = "")
legend(-20, 20, unique(grp1), col = rainbow(length(unique(grp1))), # or use "hcl.colors" function with different "palette" argument from hcl.pals()
       text.col = "black", pch = c(19), bty = "Y", merge = F, bg = "white")

########################################### other type

tsne <- do.tsne(as.matrix(mtrx1), ndim=2, perplexity=5) 
plot(tsne$Y, pch=19, col=hcl.colors(length(unique(grp1)), palette = "Blue-Red")[as.integer(as.factor(grp1))], main="t-SNE", xlab = "", ylab = "") # or use "hcl.colors" function with different "palette" argument from hcl.pals()

###############################################
############################################### OTHER TYPES OF CLUSTERING
###############################################

# dataset
base1 <- ds # dataset
mtrx1 <- ds[,-1] # numeric data
grp1 <- as.character(base1[,1]) # label of dataset
k <- length(unique(grp1)) # groups in clustering

########################################### DBSCAN, HDBSCAN

res <- dbscan::dbscan(mtrx1, eps = 0.7, minPts = 5)
fviz_cluster(list(data = mtrx1, cluster = res$cluster), repel = T,
             ellipse.type = "euclid", geom = "point", stand = FALSE,
             palette = "jco", ggtheme = theme_classic()) # or other palette from ggsci

res <- dbscan::hdbscan(mtrx1, minPts = 5)
fviz_cluster(list(data = mtrx1, cluster = res$cluster), repel = T,
             ellipse.type = "euclid", geom = "point", stand = FALSE,
             palette = "jco", ggtheme = theme_classic()) # or other palette from ggsci

########################################### Spectral Clustering

k <- length(unique(grp1)) # groups in clustering
sc <- speccCBI(mtrx1, k = k)
fviz_cluster(list(data = mtrx1, cluster = sc$partition), repel = T,
             ellipse.type = "euclid", geom = "point", stand = FALSE,
             palette = "jco", ggtheme = theme_classic()) # or other palette from ggsci

########################################### UMAP

plot.umap = function(x, labels, main="UMAP", colors=hcl.colors(length(unique(labels)), "Fall"),
                     pad=0.1, cex=0.95, pch=19, add= F, legend.suffix="", cex.main=1.5, cex.legend=1) {
  
  layout = x
  if (is(x, "umap")) {
    layout = x$layout
  } 
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    #par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", xlab = "", ylab = "")
    #rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)  
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  
  labels.u = unique(labels)
  legend.pos = "topleft"
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomleft"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  
  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

umap = umap(mtrx1)
plot.umap(umap, labels = as.factor(grp1))

########################################### MCLUST

BIC <- mclustBIC(mtrx1)
plot(BIC)
mod1 <- Mclust(mtrx1, x = BIC)
fviz_mclust(mod1, "BIC", palette = "jco")
fviz_mclust(mod1, "classification", geom = "point", pointsize = 1.5, palette = "jco")
fviz_mclust(mod1, "uncertainty", palette = "jco")
fviz_cluster(list(data = mtrx1, cluster = mod1$classification), repel = T,
             ellipse.type = "euclid", geom = "point", stand = FALSE,
             palette = "jco", ggtheme = theme_classic()) # or other palette from ggsci

# other type
mod2 <- Mclust(mtrx1)
fviz_mclust(mod2, "BIC", palette = "jco")
fviz_mclust(mod2, "classification", geom = "point", pointsize = 1.5, palette = "jco")
fviz_mclust(mod2, "uncertainty", palette = "jco")
fviz_cluster(list(data = mtrx1, cluster = mod2$classification), repel = T,
             ellipse.type = "euclid", geom = "point", stand = FALSE,
             palette = "jco", ggtheme = theme_classic()) # or other palette from ggsci

# other type
mod3 <- densityMclust(mtrx1, plot = F)
fviz_mclust(mod3, "BIC", palette = "jco")
fviz_mclust(mod3, "classification", geom = "point", pointsize = 1.5, palette = "jco")
fviz_mclust(mod3, "uncertainty", palette = "jco")
fviz_cluster(list(data = mtrx1, cluster = mod3$classification), repel = T,
             ellipse.type = "euclid", geom = "point", stand = FALSE,
             palette = "jco", ggtheme = theme_classic()) # or other palette from ggsci

# other type
mod1dr <- MclustDR(mod1) # try mod1,mod2,mod3
plot(mod1dr, what = "boundaries", ngrid = 200)
plot(mod1dr, what = "scatterplot")

########################################### MDS
mds <- do.mds(as.matrix(mtrx1), ndim=2)
plot(mds$Y, pch=19, col=hcl.colors(length(unique(grp1)), palette = "Blue-Red")[as.integer(as.factor(grp1))], main="MDS") # or use "hcl.colors" function with different "palette" argument from hcl.pals()

########################################### LLE
lle <- do.lle(as.matrix(mtrx1),ndim=2,type=c("proportion",0.20))
plot(lle$Y, pch=19, col=hcl.colors(length(unique(grp1)), palette = "Blue-Red")[as.integer(as.factor(grp1))], main="LLE") # or use "hcl.colors" function with different "palette" argument from hcl.pals()

########################################### IsoMap
isomap <- do.isomap(as.matrix(mtrx1),ndim=2,type=c("proportion",0.25),weight=FALSE)
plot(isomap$Y, pch=19, col=hcl.colors(length(unique(grp1)), palette = "Blue-Red")[as.integer(as.factor(grp1))], main="IsoMap") # or use "hcl.colors" function with different "palette" argument from hcl.pals()

########################################### Laplacian Score
ls <- do.lscore(as.matrix(mtrx1), t=0.1)
plot(ls$Y, pch=19, col=hcl.colors(length(unique(grp1)), palette = "Blue-Red")[as.integer(as.factor(grp1))], main="Laplacian Score") # or use "hcl.colors" function with different "palette" argument from hcl.pals()

########################################### Diffusion Maps
dm <- do.dm(as.matrix(mtrx1), bandwidth=10)
plot(dm$Y, pch=19, col=hcl.colors(length(unique(grp1)), palette = "Blue-Red")[as.integer(as.factor(grp1))], main="Diffusion Maps") # or use "hcl.colors" function with different "palette" argument from hcl.pals()

########################################### kernel PCA
kpca <- do.kpca(as.matrix(mtrx1),ndim = 2)
plot(kpca$Y, pch=19, col=hcl.colors(length(unique(grp1)), palette = "Blue-Red")[as.integer(as.factor(grp1))], main="kernel PCA") # or use "hcl.colors" function with different "palette" argument from hcl.pals()

########################################### sparse PCA
spca <- do.spca(as.matrix(mtrx1),ndim = 2)
plot(spca$Y, pch=19, col=hcl.colors(length(unique(grp1)), palette = "Blue-Red")[as.integer(as.factor(grp1))], main="sparse PCA") # or use "hcl.colors" function with different "palette" argument from hcl.pals()

########################################### ICA
ica <- do.ica(as.matrix(mtrx1),ndim=2,type="poly")
plot(ica$Y, pch=19, col=hcl.colors(length(unique(grp1)), palette = "Blue-Red")[as.integer(as.factor(grp1))], main="ICA") # or use "hcl.colors" function with different "palette" argument from hcl.pals()

########################################### FA
fa <- do.fa(as.matrix(mtrx1),ndim=2)
plot(fa$Y, pch=19, col=hcl.colors(length(unique(grp1)), palette = "Blue-Red")[as.integer(as.factor(grp1))], main="FA") # or use "hcl.colors" function with different "palette" argument from hcl.pals()

########################################### tSNE
tsne <- do.tsne(as.matrix(mtrx1), ndim=2, perplexity=40)
plot(tsne$Y, pch=19, col=hcl.colors(length(unique(grp1)), palette = "Blue-Red")[as.integer(as.factor(grp1))], main="t-SNE", xlab = "", ylab = "") # or use "hcl.colors" function with different "palette" argument from hcl.pals()
 
########################################### NMF

data_set <- as(mtrx1, "dimRedData")
data_set@meta <- as.data.frame(grp1, meta.prefix = "meta.", data.prefix = "")
labels_v <- as.factor(grp1)

nmf <- embed(data_set, "NNMF")
plot(nmf@data@data, col=hcl.colors(length(unique(grp1)), palette = "Blue-Red")[as.integer(as.factor(grp1))]) # or use "hcl.colors" function with different "palette" argument from hcl.pals()

########################################### PAM

k <- length(unique(grp1)) # groups in clustering
pam <- pam(mtrx1, k, metric = "euclidean", stand = FALSE)
fviz_cluster(list(data = mtrx1, cluster = pam$clustering), repel = T,
             ellipse.type = "euclid", geom = "point", stand = FALSE,
             palette = "jco", ggtheme = theme_classic()) # or other palette from ggsci

########################################### CLARA

k <- length(unique(grp1)) # groups in clustering
clara <- clara(mtrx1, k, metric = "euclidean", stand = FALSE, pamLike = FALSE)
fviz_cluster(list(data = mtrx1, cluster = clara$clustering), repel = T,
             ellipse.type = "euclid", geom = "point", stand = FALSE,
             palette = "jco", ggtheme = theme_classic()) # or other palette from ggsci

########################################### Fuzzy Clustering

k <- length(unique(grp1)) # groups in clustering
fc <- fanny(mtrx1, k, metric = "euclidean", stand = FALSE)
fviz_cluster(list(data = mtrx1, cluster = fc$clustering), repel = T,
             ellipse.type = "euclid", geom = "point", stand = FALSE,
             palette = "jco", ggtheme = theme_classic()) # or other palette from ggsci

###############################################
############################################### VALIDATION CLUSTERING
###############################################

# "fviz_nbclust" function
c <- fviz_nbclust(mtrx1, hcut, linecolor = "red", method = "gap_stat", nboot = 100)+ # method: "silhouette" or "wss" or "gap_stat" FUNcluster: "means", "cluster::pam", "cluster::clara", "cluster::fanny", "hcut"
  labs(subtitle = "")
c

# "NbClust" function
nb <- NbClust(mtrx1, distance = "manhattan", diss=NULL, min.nc = 2, max.nc = 20, method = "ward.D2", index ="all") # see index in documentation
d <- fviz_nbclust(nb)
d

# hopkins statistics
set.seed(1234)
hopkins(mtrx1, n = nrow(mtrx1)-1)

# mclust
mc <- Mclust(mtrx1) # mclust method
summary(mc)

mc$modelName # Optimal selected model
mc$G # Optimal number of cluster
head(mc$z, 30) # Probality to belong to a given cluster
head(mc$classification, 30) # Cluster assignement of each observation
plot(mc, what = "classification")

BIC <- mclustBIC(mtrx1) # mclust method with BIC
mc <- Mclust(mtrx1, x = BIC, G = 1:20)
summary(mc)

# clValid
clmethods <- c("hierarchical","kmeans")
stab <- clValid(mtrx1, nClust = 2:6, clMethods = clmethods, validation = "stability")
optimalScores(stab)

# fpc
k <- length(unique(grp1)) # groups for clusters
km.res2 <- eclust(mtrx1, "kmeans", k = k, nstart = 25, graph = FALSE)
hc.res2 <- eclust(mtrx1, "hclust", k = k, hc_metric = "euclidean", hc_method = "ward.D2", graph = FALSE)
# Silhouette information
silinfo <- km.res2$silinfo
names(silinfo)
# Silhouette widths of each observation
head(silinfo$widths[, 1:3], 10)
# Average silhouette width of each cluster
silinfo$clus.avg.widths
# The total average (mean of all individual silhouette widths)
silinfo$avg.width
# The size of each clusters
km.res2$size
# Silhouette width of observation
sil <- km.res2$silinfo$widths[, 1:3]
# Objects with negative silhouette
neg_sil_index <- which(sil[, "sil_width"] < 0)
sil[neg_sil_index, , drop = F]
# Dunn index
km_stats <- cluster.stats(dist(mtrx1), km.res2$cluster)
km_stats$dunn
# Compute cluster stats
species <- as.numeric(as.factor(grp1)) 
clust_stats <- cluster.stats(d = dist(mtrx1), species, km.res2$cluster)
# Corrected Rand index
clust_stats$corrected.rand
# VI
clust_stats$vi

# Classification Accuracy by clustering
misclass <- function(pred, obs) {
  tbl <- table(pred, obs)
  sum <- colSums(tbl)
  dia <- diag(tbl)
  msc <- (sum - dia)/sum * 100
  m.m <- mean(msc)
  cat("Classification table:", "\n")
  print(tbl)
  cat("Misclassification errors:", "\n")
  (round(msc, 1))}
# for n groups
k <- length(unique(grp1)) # groups in HC
res.hc <- eclust(mtrx1, "hclust", k = k, graph = FALSE, hc_metric = "manhattan", hc_method = "ward.D2")
pred_cl <- cutree(res.hc, k=k)
for (i in 1:length(unique(pred_cl))) {
  pred_cl[pred_cl==i] <- as.character(as.data.frame(unique(grp1))[i,1])} # use <- levels(ds_group) or <- as.character(as.data.frame(unique(ds_group))[i,1]) function or as.numeric(as.factor(ds_group)) or as.fumeric(ds_group) function from rafalib package
misclass(pred_cl, grp1)

# Computing P-value for HCA
set.seed(1234)
# start parallel processing
fc <- as.numeric(detectCores(logical = T))
cl <- makePSOCKcluster(fc-1)
registerDoParallel(cl)
# parallel version of pvclust
res.pv <- parPvclust(cl, mtrx1, nboot=100)
plot(res.pv, hang = -1, cex = 0.5)
pvrect(res.pv)
clusters <- pvpick(res.pv)
clusters
# stop parallel
stopCluster(cl)
stopImplicitCluster()

##############################################################################################################################################################
# Correlation analysis
##############################################################################################################################################################

# LOAD DATA
library(data.table)
setwd("D:/...")

# PEAK TABLE
# dataset with intensities and label column only
ds <- as.data.frame(fread(input = "8 peaks.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <-"Label"
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# Library
library(Hmisc)
library(corrplot)
library(psych)

# Correlation matrix
# for features
corr <- cor(as.matrix(ds[,-1]), method = "spearman") # method: "pearson", "kendall", "spearman"
corr
# for samples
corr <- cor(as.matrix(t(ds[,-1])), method = "spearman") # method: "pearson", "kendall", "spearman"
corr
# other type
res_corr <- rcorr(as.matrix(ds[,-1]), type = "spearman") # type: "pearson", "spearman". For correlation between samples use: as.matrix(t(ds[,-1]))
res_corr$r # matrix of correlations
res_corr$P # p values

# Correlograms
corrplot(corr, diag = F, order = "FPC", tl.pos = "td", tl.cex = 0.5, method = "color", type = "upper")

pairs.panels(ds[,-1], method = "pearson", # correlation method: "pearson","spearman","kendall"
       hist.col = "#00AFBB",
       density = T,  # show density plots
       ellipses = T )# show correlation ellipses

# Auto detect type of correlation for correlation with specific feature
smart.corr.test <- function(x, n){
  b <- vector(mode="numeric")
  res <- apply(x, 2, shapiro.test)
  for (i in 1:ncol(x)) b[i] = res[[i]]$p.value 
  names(b) <- colnames(x)
  for (i in 1:ncol(x))
    if (b[i] < 0.05) {
      b[i] <- (cor.test(x = x[, i], y = x[, n], method = "spearman")$estimate)
    } else{
      b[i] <- (cor.test(x = x[, i], y = x[, n], method = "pearson")$estimate)
    }
  return(b) }

corr.2.f <- as.data.frame(smart.corr.test(ds[,-1], n = 1)) # n - target feature
corr.2.f

##############################################################################################################################################################
# Distance analysis
##############################################################################################################################################################

# LOAD DATA
library(data.table)
setwd("D:/...")

# PEAK TABLE
# dataset with intensities and label column only
ds <- as.data.frame(fread(input = "8 peaks.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <-"Label"
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# Library
library(pheatmap)

# Distance by sample
raw.dist <- dist(ds[,-1], method = "euclidian") # use any: {euclidean}, {maximum}, {manhattan}, {canberra}, {binary}, {minkowski}

# Distance by feature
col.dist <- dist(t(ds[,-1]), method = "euclidian") # use any: {euclidean}, {maximum}, {manhattan}, {canberra}, {binary}, {minkowski}

# Correlation Distance by sample
rows.cor <- cor(t(ds[,-1]), use = "pairwise.complete.obs", method = "pearson") # use any: "pearson", "kendall", "spearman"

# Correlation Distance by feature
cols.cor <- cor(ds[,-1], use = "pairwise.complete.obs", method = "pearson") # use any: "pearson", "kendall", "spearman"

# Plot the heatmap of correlation
pheatmap(ds[,-1], 
  clustering_distance_cols = as.dist(1 - cols.cor),
  clustering_distance_rows = as.dist(1 - rows.cor))

# Plot the heatmap by sample
pheatmap(ds[,-1], clustering_rows = T, cluster_cols = F,
         clustering_distance_rows = "manhattan", clustering_distance_cols	= "manhattan", # use any: {euclidean}, {maximum}, {manhattan}, {canberra}, {binary}, {minkowski}
         clustering_method = "ward.D2") # use any: {ward (ward.D), (ward.D2)}, {single}, {complete}, {average}, {mcquitty},{median}, {centroid}

# Plot the heatmap by feature
pheatmap(t(ds[,-1]), clustering_rows = T, cluster_cols = F,
         clustering_distance_rows = "manhattan", clustering_distance_cols	= "manhattan", # use any: {euclidean}, {maximum}, {manhattan}, {canberra}, {binary}, {minkowski}
         clustering_method = "ward.D2") # use any: {ward (ward.D), (ward.D2)}, {single}, {complete}, {average}, {mcquitty},{median}, {centroid}

##############################################################################################################################################################
# Sample Size and Power Calculation
##############################################################################################################################################################

# LOAD DATA
library(data.table)
setwd("D:/...")

# PEAK TABLE
# dataset with intensities and label column only
ds <- as.data.frame(fread(input = "8 peaks.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <-"Label"
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# Library
library(effectsize)
library(pwr)
library(ggplot2)

# Effect size
es <- mean(sapply(2:ncol(ds), function(y) cohens_d(ds[,y], ds[,1])$Cohens_d)) # Cohen's d
es <- mean(sapply(2:ncol(ds), function(y) hedges_g(ds[,y], ds[,1])$Hedges_g)) # Hedges' g

# Other values
n <- nrow(ds) # Number of observations
n1 <- length(which(ds[,1] == unique(ds[,1])[1])) # Number of observations in the 1 sample
n2 <- length(which(ds[,1] == unique(ds[,1])[2])) # Number of observations in the 2 sample

# Chi squared power calculation
pwr <- pwr.chisq.test(w = es, df = 2, sig.level = 0.05, power = 0.8) # adjust to your study

# t test power calculation
pwr <- pwr.t.test(d = es, sig.level = 0.05, power = 0.8, type = "two.sample") # adjust to your study
pwr <- pwr.t2n.test(n1 = length(which(ds[,1] == unique(ds[,1])[1])), # adjust to your study
             n2 = length(which(ds[,1] == unique(ds[,1])[2])), sig.level = 0.05, power = 0.8)

# ANOVA power calculation
pwr <- pwr.anova.test(k = length(unique(ds[,1])), n = nrow(ds)/length(unique(ds[,1])), # adjust to your study
                sig.level = 0.05, power = 0.8)

# results
pwr

# plot
plot(pwr) + theme_minimal()

##############################################################################################################################################################
# References
##############################################################################################################################################################

# 1. Li, Shuzhao, ed. Computational Methods and Data Analysis for Metabolomics. Humana Press, 2020.
# 2. Kuhn, Max, and Kjell Johnson. Applied predictive modeling. Vol. 26. New York: Springer, 2013.
# 3. James, Gareth, et al. An introduction to statistical learning. Vol. 112. New York: springer, 2013.
# 4. Plyushchenko, Ivan, et al. "An approach for feature selection with data modelling in LC-MS metabolomics." Analytical Methods 12.28 (2020): 3582-3591.
# 5. Vabalas, Andrius, et al. "Machine learning algorithm validation with a limited sample size." PloS one 14.11 (2019): e0224365.
# 6. Parvandeh, Saeid, et al. "Consensus features nested cross-validation." Bioinformatics 36.10 (2020): 3093-3098.
# 7. Lewis, Nigel Da Costa. 100 Statistical Tests in R: What to Choose, how to Easily Calculate, with Over 300 Illustrations and Examples. Heather Hills Press, 2013.
# 8. Kabacoff, Robert. R in Action. Shelter Island, NY, USA: Manning publications, 2011.
# 9. Smyth, Gordon K. "Limma: linear models for microarray data." Bioinformatics and computational biology solutions using R and Bioconductor. Springer, New York, NY, 2005. 397-420.
# 10. De Livera, Alysha M., et al. "Normalizing and integrating metabolomics data." Analytical chemistry 84.24 (2012): 10768-10776.
# 11. Kursa, Miron B., and Witold R. Rudnicki. "Feature selection with the Boruta package." Journal of statistical software 36 (2010): 1-13.
# 12. Tanner, Eva M., Carl-Gustaf Bornehag, and Chris Gennings. "Repeated holdout validation for weighted quantile sum regression." MethodsX 6 (2019): 2855-2860.
# 13. Wood, Simon N. Generalized additive models: an introduction with R. CRC press, 2017.
# 14. Larras, Floriane, et al. "DRomics: A turnkey tool to support the use of the dose-response framework for Omics data in ecological risk assessment." Environmental science & technology 52.24 (2018): 14461-14468.
# 15. Yao, Cong-Hui, et al. "Dose-response metabolomics to understand biochemical mechanisms and off-target drug effects with the TOXcms software." Analytical chemistry 92.2 (2019): 1856-1864.
# 16. Bodein, Antoine, et al. "timeOmics: an R package for longitudinal multi-omics data integration." Bioinformatics (2021).
# 17. Li, Mengci, et al. "polyPK: an R package for pharmacokinetic analysis of multi-component drugs using a metabolomics approach." Bioinformatics 34.10 (2018): 1792-1794.
# 18. Smilde, Age K., et al. "ANOVA-simultaneous component analysis (ASCA): a new tool for analyzing designed metabolomics data." Bioinformatics 21.13 (2005): 3043-3048.
# 19. Sanchez-Illana, Angel, et al. "Evaluation of batch effect elimination using quality control replicates in LC-MS metabolite profiling." Analytica chimica acta 1019 (2018): 38-48.
# 20. Fages, Anne, et al. "Investigating sources of variability in metabolomic data in the EPIC study: the Principal Component Partial R-square (PC-PR2) method." Metabolomics 10.6 (2014): 1074-1083.
# 21. Yi, Sangyoon, et al. "2dFDR: a new approach to confounder adjustment substantially increases detection power in omics association studies." Genome biology 22.1 (2021): 1-18.
# 22. Rodriguez-Martinez, Andrea, et al. "MWASTools: an R/bioconductor package for metabolome-wide association studies." Bioinformatics 34.5 (2018): 890-892.
# 23. Liquet, Benoit, et al. "A novel approach for biomarker selection and the integration of repeated measures experiments from two assays." BMC bioinformatics 13.1 (2012): 1-14.
# 24. Tai, Yu Chuan, and Terence P. Speed. "A multivariate empirical Bayes statistic for replicated microarray time course data." The Annals of Statistics (2006): 2387-2412.
