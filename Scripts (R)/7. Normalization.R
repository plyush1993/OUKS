##############################################################################################################################################################
# Table of contents
##############################################################################################################################################################

# Installation
# Normalization
# Scaling
# Wrapper function
# Removing and adjustment of biological variation
# Evaluation
# Save with all info
# References

##############################################################################################################################################################
# Installation
##############################################################################################################################################################

# setup environment
library(data.table)
setwd("D:/...")

# load dataset
ds <- as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats annot+filtr.csv", header=T))
ds <- ds[-c(1:8),]
rownames(ds) <- ds[,5]
ds <- ds[,-c(1:5)]
dsr <- sapply(ds, as.numeric)
rownames(dsr) <- rownames(ds)
ds_norm <- dsr # remove all column with non-numeric values
ds_norm_t <- t(ds_norm)

# load metadata
meta <- as.data.frame(fread(input = "metadata.csv", header=T))
rownames(meta) <- meta[,5]
order <- meta[,2]
meta <- meta[,c(1,3,8:10)]
colnames(meta) <- c("Class", "Batch", "Creatinine", "Age", "Sex")
# check identical rownames between peak table and metadata
identical(rownames(ds), rownames(meta))

# save without normalization for evaluation
# ds_no_norm <- as.data.frame(ds_norm)
# fwrite(ds_no_norm, "xcms after IPO MVI filter repeats annot+filtr no norm.csv", row.names = T)

##############################################################################################################################################################
# Normalization
##############################################################################################################################################################

########################################## MSTUS

# function
MSTUS <- function(data) {
  data_sum <- matrix(colSums(data, na.rm = T), nrow = 1)
  uni <- matrix(rep(1, nrow(data)), ncol = 1)
  area.uni <- uni %*% data_sum
  MSTUS <- data/area.uni
  return(MSTUS)
}

# perform
ds_norm_mstus <- as.data.frame(t(MSTUS(ds_norm_t)))

# save
fwrite(ds_norm_mstus, "xcms after IPO MVI QC-XGB filter repeats annot+filtr MSTUS.csv", row.names = T)

########################################## PQN

# function
PQN <- function(data) {
  reference <- apply(data, 1, function(x) median(x, na.rm = T))
  quotient <- data / reference
  quotient.median <- apply(quotient, 2, function(x) median(x, na.rm = T))
  pqn.data <- t(t(data) / quotient.median)
  return(pqn.data)
}

# perform
ds_norm_pqn <- as.data.frame(t(PQN(ds_norm_t)))

# save
fwrite(ds_norm_pqn, "xcms after IPO MVI QC-XGB filter repeats annot+filtr PQN.csv", row.names = T)

########################################## Quantile

library(affy)

# function
QUANTILE <- function(data) {
  normalize.quantile <- get("normalize.quantiles", en = asNamespace("affy"))
  quantile.data <- normalize.quantile(data)
  rownames(quantile.data) <- rownames(data)
  colnames(quantile.data) <- colnames(data)
  return(quantile.data)
}

# perform
ds_norm_quant <- as.data.frame(t(QUANTILE(ds_norm_t)))

# save
fwrite(ds_norm_quant, "xcms after IPO MVI QC-XGB filter repeats annot+filtr QUANTILE.csv", row.names = T)

########################################## Median Fold Change

# function
medFC <- function(mat) {
  # Perform median fold change normalisation
  #          X - data set [Variables & Samples]
  medSam <- apply(mat, 1, median)
  medSam[which(medSam==0)] <- 0.0001
  mat <- apply(mat, 2, function(mat, medSam){
    medFDiSmpl <- mat/medSam
    vec<-mat/median(medFDiSmpl)
    return(vec)
  }, medSam)
  return (mat)
}

# perform
ds_norm_medfc <- as.data.frame(t(medFC(ds_norm_t)))

# save
fwrite(ds_norm_medfc, "xcms after IPO MVI QC-XGB filter repeats annot+filtr MedFC.csv", row.names = T)

########################################## Sample-specific Factor

# Sample factor
# ssf <- c(1:nrow(ds_norm)) # just simple example of sample-specific factor
ssf <- meta$Creatinine # adjust to your data: select variable from metadata
uni <- matrix(rep(1, ncol(ds_norm)), ncol = 1)
ssfm <- t(uni %*% ssf)

# Perform
ds_norm_ssf <- as.data.frame(ds_norm/ssfm)

# save
fwrite(ds_norm_ssf, "xcms after IPO MVI QC-XGB filter repeats annot+filtr CREATININE.csv", row.names = T)

##############################################################################################################################################################
# Scaling
##############################################################################################################################################################

########################################## log transformation

# perform
ds_norm[ds_norm == 0] <- NA
ds_norm_log <- log2(ds_norm)
ds_norm_log <- as.data.frame(ds_norm_log)

# save
fwrite(ds_norm_log, "xcms after IPO MVI QC-XGB filter repeats annot+filtr LOG.csv", row.names = T)

########################################## scaling

# perform
ds_norm_sc <- data.frame(scale(ds_norm, center = F, scale = apply(ds_norm, 2, sd, na.rm = T))) # for centering: center = T

# save
fwrite(ds_norm_sc, "xcms after IPO MVI QC-XGB filter repeats annot+filtr scal.csv", row.names = T)

##############################################################################################################################################################
# Wrapper function
##############################################################################################################################################################

##########################################

library(DiffCorr)

# perform
ds_norm_dc <- scalingMethods(ds_norm_t, "pareto") # or "auto", "range", "pareto", "vast", "level", "power"
ds_norm_dc <- as.data.frame(t(ds_norm_dc))

# save
fwrite(ds_norm_dc, "xcms after IPO MVI QC-XGB filter repeats annot+filtr dc.csv", row.names = T)

##########################################

library(clusterSim)

# perform
ds_norm_cs <- data.Normalization(ds_norm, type = "n6", normalization = "column") # or n0 - n13
ds_norm_cs <- as.data.frame(ds_norm_cs)

# save
fwrite(ds_norm_cs, "xcms after IPO MVI QC-XGB filter repeats annot+filtr cs.csv", row.names = T)

##############################################################################################################################################################
# Removing and adjustment of biological variation
##############################################################################################################################################################

############################################### LMM MODELLING OF BIOLOGICAL FACTORS

library(lmm2met)

# data
n_meta <- 5 # adjust to your data
dat <- cbind(meta, ds_norm)
colnames(dat)[-c(1:n_meta)] <- paste0("X", c(1:ncol(dat))) # adjust to your data
dat$Batch <- as.factor(dat$Batch)
dat$Age <- as.integer(dat$Age)
dat$Sex <- as.integer(dat$Sex)
dat$Class <- as.factor(dat$Class)
# dat[,-c(1:n_meta)] <- sapply(dat[,-c(1:n_meta)], as.numeric)

# perform
fitMet <- fitLmm(fix=c("Sex","Age","Class"), random="(1|Batch)", data=dat, start=6) # adjust to your data

# save adjusted dataset
ds_lmm_fit <- fitMet$fittedDat
dat2 <- cbind(meta, ds)
colnames(ds_lmm_fit) <- colnames(dat2)
ds_lmm_fit <- ds_lmm_fit[,-c(2:n_meta)]
fwrite(ds_lmm_fit, "QC-XGB after dual filt LMM adj.csv", row.names = T)

############################################### LM MODELLING OF BIOLOGICAL FACTORS (FACTORWISE)

library(MetabolomicsBasics)

# data
dat <- ds_norm
s_d <- cbind(Order = order, meta)

# perform
fmod <- "Batch+Class+Age+Sex+Order" # or "Batch+Class+Age+Sex+Order" or "Batch+Class+Age+Sex" adjust to your data
kmod <- "Class"
rawc <- RemoveFactorsByANOVA(y = dat, sam = s_d, fmod = fmod, kmod = kmod)

# save adjusted dataset
ds_lm_fit <- data.frame(cbind(Class = meta$Class, data.frame(rawc)))
colnames(ds_lm_fit)[-1] <- colnames(rawc)
fwrite(ds_lm_fit, "QC-XGB after dual filt LM adj+order.csv", row.names = T)

############################################### LM MODELLING OF BIOLOGICAL FACTORS

# data
dat <- ds_norm
s_d <- cbind(Order = order, meta)

# perform
model="Batch+Class+Age+Sex" # adjust to your data
mod_f <- as.formula(paste("y", model, sep = " ~ "))
lm_fit <- lapply(1:ncol(dat), function(y) lm(mod_f, data = cbind(s_d, y = dat[,y])))
lm_adj <- as.data.frame(sapply(1:ncol(dat), function(y) predict(lm_fit[[y]], cbind(s_d, y = dat[,y]))))
colnames(lm_adj) <- colnames(dat)

# save adjusted dataset
ds_lm_fit <- data.frame(cbind(Class = meta$Class, data.frame(lm_adj)))
fwrite(ds_lm_fit, "QC-XGB after dual filt LM simple adj+order.csv", row.names = T)

##############################################################################################################################################################
# Evaluation
##############################################################################################################################################################

############################################### Evaluation of biological variation adjustment

###############################################
###############################################

# load data
library(data.table)
setwd("D:/...")

# load dataset
ds <- as.data.frame(fread(input = "QC-XGB after dual filt LMM adj.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-c(1,2)] # only numeric variables
# ds <- sapply(ds, as.numeric)

# load metadata
meta <- as.data.frame(fread(input = "metadata.csv", header=T))
rownames(meta) <- meta[,5]
meta <- meta[,c(1,3,8:10)]
colnames(meta) <- c("Class", "Batch", "Creatinine", "Age", "Sex")
# check identical rownames between peak table and metadata
identical(rownames(ds), rownames(meta))

# create dataset
ds_abv <- as.data.frame(cbind(Label = meta$Class, ds))

############################################### ML accuracy

library(parallel)
library(doParallel)
library(caret)
library(tuple)

# stop parallel
stopCluster(cl)
stopImplicitCluster()
registerDoSEQ()

# start parallel processing
fc <- as.numeric(detectCores(logical = F))
cl <- makePSOCKcluster(fc+1)
registerDoParallel(cl)

# repeated cross validation
trainControl <- trainControl(method="repeatedcv", number=10, repeats=10) # or bootstrap: trainControl(method="boot", number=100)
metric <- "Accuracy"

# PAM
set.seed(1234)
# library(pamr)
fit.pam <- train(Label~ ., data=ds_abv, method="pam", metric=metric, trControl=trainControl, tuneLength = 10)

# SVM
set.seed(1234)
# library(kernlab)
fit.svm <- train(Label~ ., data=ds_abv, method="svmRadial", metric=metric, trControl=trainControl, tuneLength = 10)

# PLS
set.seed(1234)
# library(pls)
fit.pls <- train(Label~ ., data=ds_abv, method="pls", metric=metric, trControl=trainControl, tuneLength = 10)

# RF
set.seed(1234)
# library(randomForest)
fit.rf <- train(Label~ ., data=ds_abv, method="rf", metric=metric, trControl=trainControl, tuneLength = 10)

# only accuracy for all models
results <- resamples(list(pam=fit.pam, svm=fit.svm, rf=fit.rf, pls=fit.pls), trControl = trainControl, metric=metric)
results_df <- as.data.frame(results)
results_df_fin <- apply(results_df[,-5], 2, mean)
results_df_fin

############################################### PCA plot

library(factoextra)
library(FactoMineR)
library(RSEIS)
library(ggsci)

# dataset
base1 <- ds_abv # dataset
mtrx1 <- ds_abv[,-1] # numeric data
grp1 <- as.character(base1[,1]) # label of dataset

palette_pca <- "lancet" # color: "lancet" or "category20" or pal_gsea("default", n = length(unique(grp1)), alpha = 0.6, reverse = T)(length(unique(grp1)))
# palette_pca <- JGRAY(length(unique(grp1))) # grey

pca.ds1 <- PCA(mtrx1, scale.unit = T, graph = FALSE)
pca <- fviz_pca_ind(pca.ds1,
                    title = "",
                    geom.ind = "point", # show points only 
                    col.ind = grp1, # color by groups
                    palette = palette_pca, # color "jco" gray JGRAY(length(unique(grp1)))
                    addEllipses = T, # Concentration ellipses
                    legend.title = "Groups")

pca 

############################################### LM Modelling

library(MetabolomicsBasics)

# data
dat <- ds
s_d <- cbind(Order = order, meta)

# perform
model <- MetaboliteANOVA(dat=dat, sam=s_d, model="Batch+Class+Age+Sex") # "Batch+Class+Age+Sex" or "Batch+Class+Age+Sex+Order" # adjust to your data

# features
class_un <- sum(model[,"Batch"]>0.05 & model[,"Age"]>0.05 & model[,"Sex"]>0.05 & model[,"Class"]<0.05) # adjust to your data
class_un
class_er <- sum(model[,"Batch"]>0.05 & model[,"Age"]>0.05 & model[,"Sex"]>0.05) # adjust to your data
class_er
class_er2 <- sum(model[,"Batch"]<0.05 & model[,"Age"]<0.05 & model[,"Sex"]<0.05) # adjust to your data
class_er2
class <- length(rownames(model[which(model[,"Class"]<0.05),])) # adjust to your data
class

############################################### LMM Modelling

library(lme4)
library(lmerTest)

# data
dat <- cbind(meta, ds)
n_meta <- 5 # adjust to your data
colnames(dat)[-c(1:n_meta)] <- paste0("X", c(1:ncol(dat))) # adjust to your data
dat$Batch <- as.factor(dat$Batch)
dat$Age <- as.integer(dat$Age)
dat$Sex <- as.integer(dat$Sex)
dat$Class <- as.factor(dat$Class)
# dat[,-c(1:n_meta)] <- sapply(dat[,-c(1:n_meta)], as.numeric)

# perform
n_start <- 6 # adjust to your data
lmm_fit <- lapply(n_start:ncol(dat), function(x) lmer(dat[,x] ~ Class + Sex + Age + (1|Batch), dat)) # adjust to your data # lmer from lmerTest package
lmm_fit_coef <- lapply(1:length(lmm_fit), function(x) summary(lmm_fit[[x]])$coefficients)
lmm_fit_pval <- sapply(1:length(lmm_fit_coef), function(x) lmm_fit_coef[[x]][,5][2]) # adjust to your data
lmm_fit_pval_all_df <- as.data.frame(t(sapply(1:length(lmm_fit_coef), function(x) lmm_fit_coef[[x]][,5]))) # adjust to your data
dat2 <- cbind(meta, ds)
rownames(lmm_fit_pval_all_df) <- colnames(dat2)[-c(1:5)] # adjust to your data
colnames(lmm_fit_pval_all_df)[2:4] <- c("Class", "Sex", "Age") # adjust to your data

# features
class_un <- sum(lmm_fit_pval_all_df[,"Age"]>0.05 & lmm_fit_pval_all_df[,"Sex"]>0.05 & lmm_fit_pval_all_df[,"Class"]<0.05) # adjust to your data
class_un
class_er <- sum(lmm_fit_pval_all_df[,"Age"]>0.05 & lmm_fit_pval_all_df[,"Sex"]>0.05) # adjust to your data
class_er
class_er2 <- sum(lmm_fit_pval_all_df[,"Age"]<0.05 & lmm_fit_pval_all_df[,"Sex"]<0.05) # adjust to your data
class_er2
class <- length(rownames(lmm_fit_pval_all_df[which(lmm_fit_pval_all_df[,"Class"]<0.05),])) # adjust to your data
class

lmm_ind_c <- which(lmm_fit_pval_all_df[,"Age"]>0.05 & lmm_fit_pval_all_df[,"Sex"]>0.05 & lmm_fit_pval_all_df[,"Class"]<0.05)+n_meta # adjust to your data
lmm_ind <- which(lmm_fit_pval <= 0.05)+n_meta
lmm <- colnames(dat2)[lmm_ind]

############################################### Evaluation of Normalization and Scaling

###############################################
###############################################

# load data
library(data.table)
setwd("D:/...")

# load dataset
ds <- as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats annot+filtr LMM adj KEGG.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]

# load metadata
meta <- as.data.frame(fread(input = "metadata.csv", header=T))
rownames(meta) <- meta[,5]
meta <- meta[,c(1,3,8:10)]
colnames(meta) <- c("Class", "Batch", "Creatinine", "Age", "Sex")
# check identical rownames between peak table and metadata
identical(rownames(ds), rownames(meta))

# dataset for evaluation
ds_ra <- cbind(Class = meta$Class, ds) # use Class or Batch

######################################## PLOT RLA

library(NormalizeMets) # try also enviGCMS package (function: plotridges)

# Log transform
log_data <- LogTransform(featuredata = ds_ra[,-1])
# Single RLA
RlaPlots(featuredata = log_data$featuredata, groupdata = ds_ra[,1], ylim = c(2,-2), cex.axis = 0.6, saveinteractiveplot = F, interactiveplot = F)
# Multiple RLA
rla_data<-list(version1=log_data$featuredata, version2=log_data$featuredata)
CompareRlaPlots(rla_data ,groupdata=ds[,1], normmeth = c("version1","version2"), saveinteractiveplot = F)

######################################## PLOT BOXPLOT

library(RColorBrewer)
library(NormalizeMets)
library(ggplot2)

# Log transform
log_data <- LogTransform(featuredata = ds_ra[,-1])
lg_data <- log_data$featuredata

# boxplot
pal <-brewer.pal(n =9,name ="Set1")
cols <-pal[as.factor(ds_ra$Class)]
boxplot(t(lg_data), col =cols, main ="box plot")
legend("topright",levels(ds_ra$Class),fill =pal,bty ="n",xpd =T,  legend = c(unique(ds_ra$Class)))

# other type of boxplot
log_data <- LogTransform(featuredata = ds_ra[,-1])
dat <- stack(as.data.frame(t(log_data$featuredata)))
dat$Label <- ifelse(str_detect(dat$ind, "KMS"), "GG", "TG") # adjust to your data: create label for every value
bp <- ggplot(dat) + 
  geom_boxplot(aes(x = ind, y = values, fill = Label), outlier.size = 0.3) + theme_classic() + xlab("Sample") + ylab("Log Intensity")+theme(legend.position="none")+theme(axis.text.x = element_blank())

# bp

######################################## PLOT MA PLOT (two class study)

library(limma) # try plotMA function. Also use maplot from rafalib package or edgeR package
library(NormalizeMets)

# Log transform
log_data <- LogTransform(featuredata = ds_ra[,-1])
lg_data <- log_data$featuredata
lg_data <- as.data.frame(cbind(Class = ds_ra$Class, lg_data))
lg_data_l <- lapply(1:length(unique(ds_ra$Class)), function(y) subset(lg_data, Class==unique(ds_ra$Class)[y])[,-1])
lg_data_l <- lapply(1:length(lg_data_l), function(y) sapply(lg_data_l[[y]], as.numeric))

# MA plot 
maplot <- function(X1, X2,pch =21,main ="MA-plot",xlab ="Average log-intensities",ylab ="intensities log-ratio",lpars =list(col ="blue",lwd =2), ...){
 X <-(rowMeans(X2)+rowMeans(X1))/2
 Y <-rowMeans(X2) -rowMeans(X1)
  # plot
  scatter.smooth(x =X,y =Y,main =main,pch =pch,xlab =xlab,ylab =ylab,lpars =lpars, ...)
  abline(h =c(-1,0,1),lty =c(2,1,2))
}

maplot(t(lg_data_l[[1]]), t(lg_data_l[[2]]),main ="")

######################################## RA TOTAL

RLA <- function(data) { 
  # load dataset in form row-samples, columns-feature; 1st column-group variable
  data[data == 0] <- NA
  data_l <- data[,-1]
  fct <- apply(data_l, 2, function(y) median(y, na.rm = T)) # median value
  ra <- sapply(1:ncol(data_l), function(y) data_l[,y]-fct[y])
  radf <- as.data.frame(ra)
  return(radf)
}

rla <- RLA(ds_ra)
rsd_rla <- apply(rla, 1, function(y) sd(y, na.rm = T)/mean(y, na.rm = T))
round(median(rsd_rla, na.rm = T),0)
round(range(rsd_rla, na.rm = T)[2]-range(rsd_rla, na.rm = T)[1],0)

######################################## RA GROUPWISE

RLA_GROUPWISE <- function(data) { 
  # load dataset in form row-samples, columns-feature; 1st column-group variable
  data[data == 0] <- NA
  data_l <- data[,-1]
  gr <- split(data_l, as.factor(data$Class))
  fct <- lapply(1:length(gr), function(x) apply(gr[[x]], 2, function(y) median(y, na.rm = T))) # median value
  ra <- lapply(1:length(gr), function(x) sapply(1:ncol(gr[[x]]), function(y) gr[[x]][,y]-fct[[x]][y]))
  radf <- lapply(1:length(ra), function(x) as.data.frame(ra[[x]]))
  return(radf)
}

rla_gw <- RLA_GROUPWISE(ds_ra)
rsd_rla <- lapply(1:length(rla_gw), function(x) apply(rla_gw[[x]], 1, function(y) sd(y, na.rm = T)/mean(y, na.rm = T)))
med_rla <- lapply(1:length(rsd_rla), function(x) round(median(rsd_rla[[x]], na.rm = T),0))
ran_rla <- lapply(1:length(rsd_rla), function(x) round(range(rsd_rla[[x]], na.rm = T),0))

##############################################################################################################################################################
# Save with all info
##############################################################################################################################################################

# setup environment
library(data.table)
library(dplyr)
setwd("D:/...")

# load dataset used for normalization and/or adjustment
ds1 <- as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats annot+filtr LMM adj KEGG.csv", header=T))

ds <- ds1[-c(1:8),]
rownames(ds) <- ds[,5]
ds <- ds[,-c(1:5)]
dsr <- sapply(ds, as.numeric)
rownames(dsr) <- rownames(ds)

# load dataset after normalization and/or adjustment
dsf <- as.data.frame(fread(input = "QC-XGB after dual filt LMM adj.csv", header=T))
rownames(dsf) <- dsf[,1]
dsf <- dsf[,-c(1:2)]
identical(ds1[,5][-c(1:8)], rownames(dsf)) # check

# combine
dsf_c <- cbind(ds1[-c(1:8),2:5], dsf)
dsf_c <- rbind(ds1[1:8,-1], dsf_c)

empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) # since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}
dsf_c <- mutate_each(dsf_c, funs(empty_as_na)) 

# save
fwrite(dsf_c, "xcms after IPO MVI QC-XGB filter repeats annot+filtr LMM adj KEGG.csv", row.names = T)

##############################################################################################################################################################
# References
##############################################################################################################################################################

# 1. Wu, Yiman, and Liang Li. "Sample normalization methods in quantitative metabolomics." Journal of Chromatography A 1430 (2016): 80-95.
# 2. Li, Bo, et al. "Performance evaluation and online realization of data-driven normalization methods used in LC/MS based untargeted metabolomics analysis." Scientific reports 6 (2016): 38881.
# 3. António, Carla, ed. Plant metabolomics: Methods and protocols. Humana Press, 2018.
# 4. Wanichthanarak, Kwanjeera, et al. "Accounting for biological variation with linear mixed-effects modelling improves the quality of clinical metabolomics data." Computational and structural biotechnology journal 17 (2019): 611-618.
# 5. Gromski, Piotr S., et al. "The influence of scaling metabolomics data on model classification accuracy." Metabolomics 11.3 (2015): 684-695.
# 6. Cuevas-Delgado, Paula, et al. "Data-dependent normalization strategies for untargeted metabolomics-a case study." Analytical and Bioanalytical Chemistry (2020): 1-15.
# 7. Livera, Alysha M. De, et al. "Statistical methods for handling unwanted variation in metabolomics data." Analytical chemistry 87.7 (2015): 3606-3615.
