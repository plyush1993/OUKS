##############################################################################################################################################################
# Table of contents
##############################################################################################################################################################

# Installation
# Metadata generating
# Perform correction
# Evaluate correction
# References

##############################################################################################################################################################
# Installation
##############################################################################################################################################################

# setup environment
library(data.table)
library(dplyr)
library(stringr)
setwd("D:/...")

##############################################################################################################################################################
# Metadata generating
##############################################################################################################################################################

# load peak table in csv
# load in format: 1st column: 1. qc1(s78_1) b1 QC1(KMS 53).X --- ro. nameN(nameN_NREPEAT) BatchIndex QCN(Label N)
# "150. qc30 b6 QC30.CDF" 
# "167. s64_2 b7 MS 45.CDF"
dsr <-as.data.frame(fread(input = "xcms after IPO MVI.csv", header=T)) # first column with all metadata
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

# dataset with all metadata
ds <- data.frame(cbind(all_meta_data, dsr)) # 1-7 => metadata columns

# order by run order
ds$ro_id <- as.numeric(ds$ro_id)
ds <- ds[order(ds$ro_id, decreasing = F),] # order by ro col

##############################################################################################################################################################
# Perform correction
##############################################################################################################################################################

# Content:
# WaveICA (WaveICA, WaveICA 2.0)
# EigenMS
# Ber`s and Combat`s (Ber-Bagging, Ber, NP Combat, P Combat)
# RUVSeq`s (RUVs, RUVg, RUVr)
# QC-SPLINE (pmp)
# QC-SPLINE (notame)
# QC-LOESS
# QC-RF
# QC-SVR
# QC-LM/RLM/TOBIT
# QC-MCLUST
# QC-BT/DT/KNN
# QC-GBM (xgboost/catboost)
# QC-norm
# RUVs`s (ruvrand, ruvrandclust)
# ISs`s (SIS, NOMIS, CCMN, BMIS)

###############################################################################
########################################## WaveICA
###############################################################################

library(WaveICA)
library(WaveICA2.0)

# generate batch data
batch <- ds$b_id
batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# generate run order data
ro <- as.numeric(ds$ro_id)

######################### Perform WaveICA
# perform
ds_WaveICA<-WaveICA(data=ds[,-c(1:7)],wf="haar",batch=batch,group=NULL,K=20,t=0.05,t2=0.05,alpha=0) # adjust "batch", "group" arguments to your study
ds_WaveICA <- data.frame(ds_WaveICA$data_wave)

# save
fwrite(ds_WaveICA, "xcms after IPO MVI WaveICA.csv", row.names = T)

######################### Perform WaveICA 2.0
# perform
ds_WaveICA2<-WaveICA_2.0(data=ds[,-c(1:7)],wf="haar",Injection_Order = ro,alpha=0, Cutoff=0.1, K=10) # adjust arguments to your study
ds_WaveICA2 <- data.frame(ds_WaveICA2$data_wave)

# save
fwrite(ds_WaveICA2, "xcms after IPO MVI WaveICA 2.0.csv", row.names = T)

###############################################################################
########################################## EigenMS
###############################################################################

library(ProteoMM)

# batch data
batch <- as.factor(ds$b_id) 

# prepare data
data <- ds[,-c(1:7)]
m_logInts = as.data.frame(t(data))
m_logInts = convert_log2(m_logInts)
m_prot.info = cbind(colnames(data), colnames(data))

# perform
m_ints_eig1 = eig_norm1(m=m_logInts,treatment=batch,prot.info=m_prot.info) # for two-factor based correction see: ?eig_norm1 ("treatment" argument)
m_ints_norm1 = eig_norm2(rv=m_ints_eig1)
ds_EigenMS <- as.data.frame(t(m_ints_norm1$norm_m))

# save
fwrite(ds_EigenMS, "xcms after IPO MVI EigenMS.csv", row.names = T)

###############################################################################
########################################## Ber`s and Combat`s methods
###############################################################################

library(dbnorm)

# generate batch data
batch <- ds$b_id
batch <- str_remove(batch, "b")
batch <- as.factor(as.numeric(batch))

# data
dbn_d <- data.frame(batch,ds[,-c(1:7)])

######################### Perform Ber-Bagging correction
dbnormBagging(dbn_d)
bagg <- paste(tempdir(),'mydata_ber_baggingCorrected.csv', sep='/')
bagg <- as.data.frame(fread(input = bagg, header=T))
rownames(bagg) <- bagg[,1]
bagg <- bagg[,-c(1,2)]

# save
fwrite(bagg, "xcms after IPO MVI Ber-Bagging.csv", row.names = T)

######################### Perform Ber correction
dbnormBer(dbn_d)
ber <- paste(tempdir(),'mydata_berCorrected.csv', sep='/')
ber <- as.data.frame(fread(input = ber, header=T))
rownames(ber) <- ber[,1]
ber <- ber[,-c(1,2)]

# save
fwrite(ber, "xcms after IPO MVI Ber.csv", row.names = T)

######################### Perform NP Combat correction
dbnormNPcom(dbn_d)
npcomb <- paste(tempdir(),'mydata_nonparametricComBatCorrected.csv', sep='/')
npcomb <- as.data.frame(fread(input = npcomb, header=T))
rownames(npcomb) <- npcomb[,1]
npcomb <- npcomb[,-c(1,2)]

# save
fwrite(npcomb, "xcms after IPO MVI NP Combat.csv", row.names = T)

######################### Perform P Combat correction
dbnormPcom(dbn_d)
pcomb <- paste(tempdir(),'mydata_parametricComBatCorrected.csv', sep='/')
pcomb <- as.data.frame(fread(input = pcomb, header=T))
rownames(pcomb) <- pcomb[,1]
pcomb <- pcomb[,-c(1,2)]

# save
fwrite(pcomb, "xcms after IPO MVI P Combat.csv", row.names = T)

###############################################################################
########################################## RUVSeq`s based methods
###############################################################################

library(RUVSeq)
library(FactoMineR)
library(fpc)
# library(edgeR) #(Version 3.12.1)

# generate data for RUVSeq
dat <- ds[,-c(1:7)]
dat <- log10(dat)
idxQC <- which(ds$n_gr_t == "QC")
replicates.ind <- matrix(-1, nrow(dat) - length(idxQC) + 1, length(idxQC))
replicates.ind[1,] <- idxQC
replicates.ind[-1,1] <- (1:nrow(dat))[-idxQC]
class <- as.factor(ds$n_gr_t) # for RUVr

######################### Perform RUVs correction

NumOfComp <- 15 # Max number of components

selectNcomp <- sapply(1:NumOfComp, function(NumOfComp) {
  sprintf("Doing Batch correction using %d number of components\n",
          NumOfComp)
  doRUV <- t(RUVs(x = t(dat),
                  cIdx = 1:ncol(dat),
                  k = NumOfComp,
                  scIdx = replicates.ind,
                  round = FALSE,
                  isLog = TRUE)$normalizedCounts)
  
  pca.ds <- PCA(doRUV[idxQC,], scale.unit = T, graph = F)
  df.pca.x <- as.data.frame(pca.ds$ind$coord[,1:2]) # or pca.ds$X[,1:n], where n - number of components
  
  pca.class <- ds$b_id
  pca.class <- pca.class[idxQC]
  scores <- cbind(pca.class, df.pca.x)
  means <- lapply(1:length(unique(scores[,1])), function(y) colMeans(scores[scores$pca.class == unique(scores[,1])[y],-1] , na.rm = T))
  covs <- lapply(1:length(unique(scores[,1])), function(y) cov(scores[scores$pca.class == unique(scores[,1])[y],-1] , ))
  
  batch.dist_b <- matrix(0, length(unique(scores[,1])), length(unique(scores[,1])))
  for (i in 2:length(unique(scores[,1]))) {
    for (j in 1:(i-1)) {
      batch.dist_b[j, i] <- bhattacharyya.dist(means[[j]],
                                               means[[i]],
                                               covs[[j]],
                                               covs[[i]])
    }}
  
  m_b_d <- round(mean(batch.dist_b[col(batch.dist_b) > row(batch.dist_b)]), 2)
  
  return(m_b_d)
})

best_ncomp <- which.min(selectNcomp)

RUVs <- RUVs(x = t(dat),
             cIdx = 1:ncol(dat),
             k = best_ncomp,
             scIdx = replicates.ind,
             round = F,
             isLog = T)

ruvs <- as.data.frame(t(RUVs$normalizedCounts))

# save
fwrite(ruvs, "xcms after IPO MVI RUVs.csv", row.names = T)

######################### Perform RUVg correction

NumOfComp <- 15 # Max number of components

selectNcomp <- sapply(1:NumOfComp, function(NumOfComp) {
  sprintf("Doing Batch correction using %d number of components\n",
          NumOfComp)
  doRUV <- t(RUVg(x = t(dat),
                  cIdx = 1:ncol(dat),
                  k = NumOfComp,
                  round = FALSE,
                  isLog = TRUE)$normalizedCounts)
  
  pca.ds <- PCA(doRUV[idxQC,], scale.unit = T, graph = F)
  df.pca.x <- as.data.frame(pca.ds$ind$coord[,1:2]) # or pca.ds$X[,1:n], where n - number of components
  
  pca.class <- ds$b_id
  pca.class <- pca.class[idxQC]
  scores <- cbind(pca.class, df.pca.x)
  means <- lapply(1:length(unique(scores[,1])), function(y) colMeans(scores[scores$pca.class == unique(scores[,1])[y],-1] , na.rm = T))
  covs <- lapply(1:length(unique(scores[,1])), function(y) cov(scores[scores$pca.class == unique(scores[,1])[y],-1] , ))
  
  batch.dist_b <- matrix(0, length(unique(scores[,1])), length(unique(scores[,1])))
  for (i in 2:length(unique(scores[,1]))) {
    for (j in 1:(i-1)) {
      batch.dist_b[j, i] <- bhattacharyya.dist(means[[j]],
                                               means[[i]],
                                               covs[[j]],
                                               covs[[i]])
    }}
  
  m_b_d <- round(mean(batch.dist_b[col(batch.dist_b) > row(batch.dist_b)]), 2)
  
  return(m_b_d)
})

best_ncomp <- which.min(selectNcomp)


RUVg <- RUVg(x = t(dat),
             cIdx = 1:ncol(dat),
             k = best_ncomp,
             round = F,
             isLog = T)

ruvg <- as.data.frame(t(RUVg$normalizedCounts))

# save
fwrite(ruvg, "xcms after IPO MVI RUVg.csv", row.names = T)

######################### Perform RUVr correction

NumOfComp <- 15 # Max number of components

# Parameters for RUVr
design <- model.matrix(~class)
y <- DGEList(counts=t(dat), group=class)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

selectNcomp <- sapply(1:NumOfComp, function(NumOfComp) {
  sprintf("Doing Batch correction using %d number of components\n",
          NumOfComp)
  doRUV <- t(RUVr(x = t(dat),
                  cIdx = 1:ncol(dat),
                  k = NumOfComp,
                  round = FALSE,
                  residuals = res,
                  isLog = TRUE)$normalizedCounts)
  
  pca.ds <- PCA(doRUV[idxQC,], scale.unit = T, graph = F)
  df.pca.x <- as.data.frame(pca.ds$ind$coord[,1:2]) # or pca.ds$X[,1:n], where n - number of components
  
  pca.class <- ds$b_id
  pca.class <- pca.class[idxQC]
  scores <- cbind(pca.class, df.pca.x)
  means <- lapply(1:length(unique(scores[,1])), function(y) colMeans(scores[scores$pca.class == unique(scores[,1])[y],-1] , na.rm = T))
  covs <- lapply(1:length(unique(scores[,1])), function(y) cov(scores[scores$pca.class == unique(scores[,1])[y],-1] , ))
  
  batch.dist_b <- matrix(0, length(unique(scores[,1])), length(unique(scores[,1])))
  for (i in 2:length(unique(scores[,1]))) {
    for (j in 1:(i-1)) {
      batch.dist_b[j, i] <- bhattacharyya.dist(means[[j]],
                                               means[[i]],
                                               covs[[j]],
                                               covs[[i]])
    }}
  
  m_b_d <- round(mean(batch.dist_b[col(batch.dist_b) > row(batch.dist_b)]), 2)
  
  return(m_b_d)
})

best_ncomp <- which.min(selectNcomp)

RUVr <- RUVr(x = t(dat),
             cIdx = 1:ncol(dat),
             k = best_ncomp,
             residuals = res,
             round = F,
             isLog = T)

ruvr <- as.data.frame(t(RUVr$normalizedCounts))

# save
fwrite(ruvr, "xcms after IPO MVI RUVr.csv", row.names = T)

###############################################################################
########################################## QC-SPLINE (pmp)
###############################################################################

library(pmp)

# generate batch data
batch <- ds$b_id
batch <- str_remove(batch, "b")
batch <- as.character(batch)

# generate group data
group <- as.character(ds$n_gr_t)

# generate run order data
ro <- as.numeric(ds$ro_id)

# perform
ds_rsc <- QCRSC(df = ds[,-c(1:7)], order = ro, batch = batch,
             classes = group, spar = 0, minQC = 4, qc_label = "QC", log = F) # set minQC according to your data
ds_rsc <- data.frame(t(ds_rsc))

# save
fwrite(ds_rsc, "xcms after IPO MVI QC-SPLINE.csv", row.names = T)

###############################################################################
########################################## QC-SPLINE (notame)
###############################################################################

library(notame)
library(batchCorr)
library(openxlsx)

# extract m/z and rt from peak table
cn <- colnames(ds[,-c(1:7)])
cn <- str_remove(cn, "X")
cn0 <- colnames(ds[,-c(1:7)])
cn0 <- str_remove(cn0, "X")
cn <- gsub(pattern = "...", replacement = "_", x = cn, fixed = T)
cn2 <- t(data.frame(cn))
colnames(cn2) <- cn
peakIn <- peakInfo(PT = cn2, sep = "_", start = 1)
peakIn <- as.data.frame(cbind(peakIn, Name = cn0))
colnames(peakIn)[1:2] <- c("Mass", "RT")

# construct data table for notame
f_d <- ds[,-c(1:7)]
dsr2 <- as.data.frame(t(f_d))
Group <- as.character(ds$n_gr_t)
QC <- as.character(ds$n_gr_t)
Injection_order <- as.numeric(ds$ro_id)
qc_id <- which(QC == "QC")
QC[-qc_id] <- "Sample"
dsr_ntm <- rbind(Group, QC, Injection_order, Name = rname, dsr2)
dsr_ntm <- cbind(rownames(dsr_ntm), dsr_ntm)
nas <- as.data.frame(replicate(3, rep(NA, 4), simplify = F))
colnames(nas) <- colnames(peakIn)
add <- rbind(nas,peakIn)
dsr_ntm <- cbind(add, dsr_ntm)
dsr_ntm[1:4,4] <- c("Group", "QC", "Injection_order", "Name")
dsr_ntm[4,1:3] <- colnames(dsr_ntm)[1:3]
openxlsx::write.xlsx(dsr_ntm, file = "for notame.xlsx", col.names = F, row.names = F)

# perform preprocess
ntm <- read_from_excel(file = "for notame.xlsx", sheet = 1, name = "ntm", skip_checks =F)
ntm2 <- construct_metabosets(exprs = ntm$exprs, pheno_data = ntm$pheno_data, feature_data = ntm$feature_data, group_col = "Group")

######################### QC-SPLINE notame log=T
dc_T <- dc_cubic_spline(ntm2$ntm, log_transform = T, spar = NULL, spar_lower = 0.5, spar_upper = 1.5)
ntm_cs_T <- assayData(dc_T)
ntm_cs_T <- as.data.frame(t(exprs(ntm_cs_T$object)))
rownames(ntm_cs_T) <- rownames(f_d)
colnames(ntm_cs_T) <- colnames(f_d)
# save
fwrite(ntm_cs_T, "xcms after IPO MVI QC-SPLINE log=T.csv", row.names = T)

######################### QC-SPLINE notame log=F
dc_F <- dc_cubic_spline(ntm2$ntm, log_transform = F, spar = NULL, spar_lower = 0.5, spar_upper = 1.5)
ntm_cs_F <- assayData(dc_F)
ntm_cs_F <- as.data.frame(t(exprs(ntm_cs_F$object)))
rownames(ntm_cs_F) <- rownames(f_d)
colnames(ntm_cs_F) <- colnames(f_d)
# save
fwrite(ntm_cs_F, "xcms after IPO MVI QC-SPLINE log=F.csv", row.names = T)

###############################################################################
########################################## QC-LOESS
###############################################################################

library(statTarget)

# generate batch data
batch <- ds$b_id
batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# generate group data
group <- as.character(ds$n_gr_t)
qc_id <- which(group == "QC") # adjust to your QC label
group[qc_id] <- "NA"
group[-qc_id] <- as.numeric(as.factor(group[-qc_id]))
class <- group

# generate run order data
order <- as.numeric(ds$ro_id)

# Sample data
s_d <- data.frame(cbind(sample = paste(1:length(ds$n_gr_t), ds$n_gr_t), batch, class, order))
rownames(s_d) <- rownames(ds)

# Feature data
f_d <- data.frame(t(ds[,-c(1:7)]))
f_d <- cbind(name = colnames(ds[,-c(1:7)]), f_d)
colnames(f_d) <- c("name", paste(1:length(ds$n_gr_t), ds$n_gr_t))

# save data
fwrite(s_d, "s_d.csv", row.names = F)
fwrite(f_d, "f_d.csv", row.names = F)

# perform
datpath <- getwd()
samPeno <- paste(datpath,'s_d.csv', sep='/')
samFile <- paste(datpath,'f_d.csv', sep='/')
shiftCor(samPeno,samFile, MLmethod = 'QCRLSC', QCspan = 0, degree = 2, imputeM = 'min', coCV = 16000, Frule = 0.001) # set coCV as max as possible

# save
ds_rlc <- fread(input = "statTarget/shiftCor/After_shiftCor/shift_all_cor.csv", header = T)
ds_rlc <- ds_rlc[,-c(1:2)]
rownames(ds_rlc) <- rownames(ds)
fwrite(ds_rlc, "xcms after IPO MVI QC-LOESS.csv", row.names = T)

###############################################################################
########################################## QC-RF
###############################################################################

library(statTarget)

# generate batch data
batch <- ds$b_id
batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# generate group data
group <- as.character(ds$n_gr_t)
qc_id <- which(group == "QC") # adjust to your QC label
group[qc_id] <- "NA"
group[-qc_id] <- as.numeric(as.factor(group[-qc_id]))
class <- group

# generate run order data
order <- as.numeric(ds$ro_id)

# Sample data
s_d <- data.frame(cbind(sample = paste(1:length(ds$n_gr_t), ds$n_gr_t), batch, class, order))
rownames(s_d) <- rownames(ds)

# Feature data
f_d <- data.frame(t(ds[,-c(1:7)]))
f_d <- cbind(name = colnames(ds[,-c(1:7)]), f_d)
colnames(f_d) <- c("name", paste(1:length(ds$n_gr_t), ds$n_gr_t))

# save data
fwrite(s_d, "s_d.csv", row.names = F)
fwrite(f_d, "f_d.csv", row.names = F)

# perform
datpath <- getwd()
samPeno <- paste(datpath,'s_d.csv', sep='/')
samFile <- paste(datpath,'f_d.csv', sep='/')
shiftCor(samPeno,samFile, MLmethod = 'QCRFSC', ntree = 500, imputeM = 'min', coCV = 16000, Frule = 0.001) # set coCV as max as possible

# save
ds_rfc <- fread(input = "statTarget/shiftCor/After_shiftCor/shift_all_cor.csv", header = T)
ds_rfc <- ds_rfc[,-c(1:2)]
rownames(ds_rfc) <- rownames(ds)
fwrite(ds_rfc, "xcms after IPO MVI QC-RF.csv", row.names = T)

###############################################################################
########################################## QC-SVR
###############################################################################

library(MetNormalizer)

# generate group data
class <- as.character(ds$n_gr_t)
qc_id <- which(class == "QC")
class[-qc_id] <- "Subject"

# generate run order data
injection.order <- as.numeric(ds$ro_id)

# Sample data
s_d <- data.frame(cbind(sample.name = rownames(ds), injection.order, class))
rownames(s_d) <- rownames(ds)

# Feature data
f_d <- ds[,-c(1:7)]
f_d <- rbind(c(1:ncol(f_d)), colnames(f_d), colnames(f_d), f_d)
f_d <- cbind(c("name","mz", "rt", rownames(ds)), f_d)
f_d <- data.frame(t(f_d))

# save data
fwrite(s_d, "sample_data.csv", row.names = F, col.names = T)
fwrite(f_d, "feature_data.csv", row.names = F, col.names = F)

# perform
datpath <- getwd()
nd <- dir.create("MetNormalizer")
samPeno <- paste(datpath,'sample_data.csv', sep='/')
samFile <- paste(datpath,'feature_data.csv', sep='/')
file.copy(from = samPeno, to = "MetNormalizer/", overwrite = TRUE, recursive = TRUE)
file.copy(from = samFile, to = "MetNormalizer/", overwrite = TRUE, recursive = TRUE)
new.path <- paste(datpath,'MetNormalizer', sep='/')

######################### perform QC-SVR
metNor(
  ms1.data.name = "feature_data.csv",
  sample.info.name = "sample_data.csv",
  minfrac.qc = 0,
  minfrac.sample = 0,
  optimization = T,
  multiple = 1, # try 1 (based on inj ord) or 2-5 (top correlated features)
  threads = 3, # adjust to your system
  path = new.path
)

# save
ds_svr <- fread(input = "MetNormalizer/svr_normalization_result/data_svr_normalization.csv", header = F)
ds_svr <- data.frame(t(ds_svr[,-c(1:3,5,6)]))
rownames(ds_svr) <- str_remove(ds_svr[,1], "Sample")
ds_svr <- ds_svr[,-1]
colnames(ds_svr) <- ds_svr[1,]
ds_svr <- ds_svr[-1,]

rname1 <- rownames(ds_svr) # obtain all info from rownames
rname1 <- str_remove(rname1, ".CDF") # remove some pattern from vendor-specific format
all_id1 <- sapply(1:length(rname1), function(y) unlist(str_split(rname1[y], " "))) # split info from rownames
ro_id1 <- as.numeric(unlist(lapply(all_id1, function(y) unlist(y[1])))) # obtain run order ID (every [1] element)
ds_svr <- cbind(ro = ro_id1, ds_svr)
ds_svr <- ds_svr[order(ds_svr$ro, decreasing = F),] # order by run order
ds_svr <- ds_svr[,-1]
fwrite(ds_svr, "xcms after IPO MVI QC-SVR.csv", row.names = T)

######################### perform QC-SVR top n correlated
metNor(
  ms1.data.name = "feature_data.csv",
  sample.info.name = "sample_data.csv",
  minfrac.qc = 0,
  minfrac.sample = 0,
  optimization = T,
  multiple = 5, # try 1 (based on inj ord) or 2-5 (top correlated features)
  threads = 3, # adjust to your system
  path = new.path
)

# save
ds_svr <- fread(input = "MetNormalizer/svr_normalization_result/data_svr_normalization.csv", header = F)
ds_svr <- data.frame(t(ds_svr[,-c(1:3,5,6)]))
rownames(ds_svr) <- str_remove(ds_svr[,1], "Sample")
ds_svr <- ds_svr[,-1]
colnames(ds_svr) <- ds_svr[1,]
ds_svr <- ds_svr[-1,]

rname1 <- rownames(ds_svr) # obtain all info from rownames
rname1 <- str_remove(rname1, ".CDF") # remove some pattern from vendor-specific format
all_id1 <- sapply(1:length(rname1), function(y) unlist(str_split(rname1[y], " "))) # split info from rownames
ro_id1 <- as.numeric(unlist(lapply(all_id1, function(y) unlist(y[1])))) # obtain run order ID (every [1] element)
ds_svr <- cbind(ro = ro_id1, ds_svr)
ds_svr <- ds_svr[order(ds_svr$ro, decreasing = F),]
ds_svr <- ds_svr[,-1]
fwrite(ds_svr, "xcms after IPO MVI QC-SVR top 5.csv", row.names = T)

###############################################################################
########################################## QC-LM/RLM/TOBIT
###############################################################################

library(BatchCorrMetabolomics)

# generate group data
class <- as.character(ds$n_gr_t)
qc_id <- which(class == "QC")
class[qc_id] <- "ref"
class[-qc_id] <- "Subject"
class <- as.factor(class)

# generate batch data
batch <- as.character(ds$b_id)

# generate run order data
injection.order <- as.numeric(ds$ro_id)

# Feature data
f_d <- as.matrix(ds[,-c(1:7)])

######################### perform QC-RLM
ref.idx <- which(class == "ref")
ds_rlm <- doBC(f_d, ref.idx = ref.idx,
               batch.idx = batch, method = "rlm",
               seq.idx = injection.order, minBsamp = 3) # set minBsamp according to your data

# save
fwrite(data.frame(ds_rlm), "xcms after IPO MVI QC-RLM.csv", row.names = T)

######################### perform QC-LM
ref.idx <- which(class == "ref")
ds_lm <- doBC(f_d, ref.idx = ref.idx,
               batch.idx = batch, method = "lm",
               seq.idx = injection.order, minBsamp = 3) # set minBsamp according to your data

# save
fwrite(data.frame(ds_lm), "xcms after IPO MVI QC-LM.csv", row.names = T)

######################### perform QC-TOBIT
ref.idx <- which(class == "ref")
ds_tobit <- doBC(f_d, ref.idx = ref.idx,
               batch.idx = batch, method = "tobit",
               seq.idx = injection.order, imputeVal = 0, minBsamp = 3) # set minBsamp according to your data

# save
fwrite(data.frame(ds_tobit), "xcms after IPO MVI QC-TOBIT.csv", row.names = T)

###############################################################################
########################################## QC-MCLUST
###############################################################################

library(batchCorr)

# generate group data
class <- as.character(ds$n_gr_t)
qc_id <- which(class == "QC")
class[-qc_id] <- "Case"

# generate batch data
batch <- ds$b_id

# generate run order data
injection.order <- as.numeric(ds$ro_id)

# Sample data
s_d <- data.frame(cbind(sample_name = rownames(ds), sample_group = class, batch = batch, inj = injection.order))
rownames(s_d) <- rownames(ds)

# Feature data
f_d <- as.matrix(ds[,-c(1:7)])

# perform
# Within-batch intensity drift correction
batches_extract <- lapply(1:length(unique(s_d$batch)), function(y) getBatch(peakTable = f_d, meta = s_d, batch = s_d$batch, select = unique(s_d$batch)[y]))

batches_correct <- lapply(1:length(batches_extract), function(y) correctDrift(peakTable = batches_extract[[y]]$peakTable, injections = as.numeric(batches_extract[[y]]$meta$inj), sampleGroups = batches_extract[[y]]$meta$sample_group, 
                                                                              QCID = 'QC', G = seq(5,35,by=3), modelNames = c('VVE', 'VEE'), CVlimit = 30)) # set high value of CVlimit  
# Between-batch intensity drift correction
mergedData <- mergeBatches(batches_correct)
normData <- normalizeBatches(peakTableCorr = mergedData$peakTableCorr, peakTableOrg = mergedData$peakTableOrg,
                             batches = s_d$batch, 
                             sampleGroup = s_d$sample_group, 
                             population = 'all', CVlimit = 30) # set high value of CVlimit 

# save
ds_mclust <- data.frame(normData$peakTable)
rownames(ds_mclust) <- rownames(ds)
fwrite(ds_mclust, "xcms after IPO MVI QC-MCLUST.csv", row.names = T)

###############################################################################
########################################## QC-BT/DT/KNN
###############################################################################

# generate batch data
batch <- ds$b_id
batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# generate group data
class <- as.character(ds$n_gr_t)

# generate run order data
order <- as.numeric(ds$ro_id)

# Sample data
s_d <- data.frame(cbind(sample = rownames(ds), batch, class, order))
rownames(s_d) <- rownames(ds)

# Feature data
f_d <- data.frame(ds[,-c(1:7)])

################################################## Single mode
QC.SM <- function(int_data, order, class, qc_label, model = "bt", mode = "division", nbagg = 10, k = 2) {
  
  library(pbapply)
  
  qc_id <-  grep(qc_label, class)
  int_data_ro <- cbind(order, int_data)
  int_data_ro_qc <- int_data_ro[qc_id,]
  print("Start Correction")
  
  if (model == "knn") {
    library(caret)
    fit_sm <- pblapply(2:ncol(int_data_ro_qc), function(t) knnreg(int_data_ro_qc[,t] ~ order, data = int_data_ro_qc, k = k)) # KNN
  } else if (model == "dt") {
    library(rpart)
    fit_sm <- pblapply(2:ncol(int_data_ro_qc), function(t) rpart(int_data_ro_qc[,t] ~ order, data = int_data_ro_qc)) # Decision tree
  } else { (model == "bt")
    library(ipred)
    fit_sm <- pblapply(2:ncol(int_data_ro_qc), function(t) bagging(int_data_ro_qc[,t] ~ order, data = int_data_ro_qc, coob=T, nbagg = nbagg)) # Bagging tree
  }
  print("Done correction")
  predict_sm <- pblapply(1:length(fit_sm), function(t) predict(fit_sm[[t]], data.frame(order)))
  
  if (mode == "division") {
    val_sm <- pblapply(1:ncol(int_data), function(t) int_data[,t]/predict_sm[[t]])
    res_sm <- as.data.frame(t(do.call(rbind, val_sm)))
    rownames(res_sm) <- rownames(int_data)
    colnames(res_sm) <- colnames(int_data)
    res_sm <- res_sm*1000
  return(res_sm)
    } else { (mode == "subtraction")
    val_sm <- pblapply(1:ncol(int_data), function(t) int_data[,t]-predict_sm[[t]]+mean(int_data[qc_id,t], na.rm = T))
    res_sm <- as.data.frame(t(do.call(rbind, val_sm)))
    rownames(res_sm) <- rownames(int_data)
    colnames(res_sm) <- colnames(int_data)
    res_sm <- res_sm
  return(res_sm)
    }
  }

# Parameters of QC.SM function:  
# int_data -> numeric feature table, order -> numeric vector of the run order, class -> sample group variable, 
# qc_label -> label for QC samples in group, 
# model -> type of the model ("bt" for bagging tree, "dt" - decision tree, "knn" - KNN regression), k = number of neighbours in "knn" model, nbagg -> an integer giving the number of bootstrap replications in "bt" model
# mode -> type of correction factor ("division" or "subtraction")

######################### division mode
######################### perform bagging tree in single mode
qc_sm <- QC.SM(int_data = f_d, order = order, class = class, qc_label = "QC", model = "bt", mode = "division", nbagg = 10, k = 2)
# save
fwrite(qc_sm, "xcms after IPO MVI QC-BT d.csv", row.names = T)

######################### perform decision tree in single mode
qc_sm <- QC.SM(int_data = f_d, order = order, class = class, qc_label = "QC", model = "dt", mode = "division", nbagg = 10, k = 2)
# save
fwrite(qc_sm, "xcms after IPO MVI QC-DT d.csv", row.names = T)

######################### perform knn in single mode
qc_sm <- QC.SM(int_data = f_d, order = order, class = class, qc_label = "QC", model = "knn", mode = "division", nbagg = 10, k = 2)
# save
fwrite(qc_sm, "xcms after IPO MVI QC-KNN d.csv", row.names = T)

######################### subtraction mode
######################### perform bagging tree in single mode
qc_sm <- QC.SM(int_data = f_d, order = order, class = class, qc_label = "QC", model = "bt", mode = "subtraction", nbagg = 10, k = 2)
# save
fwrite(qc_sm, "xcms after IPO MVI QC-BT s.csv", row.names = T)

######################### perform decision tree in single mode
qc_sm <- QC.SM(int_data = f_d, order = order, class = class, qc_label = "QC", model = "dt", mode = "subtraction", nbagg = 10, k = 2)
# save
fwrite(qc_sm, "xcms after IPO MVI QC-DT s.csv", row.names = T)

######################### perform knn in single mode
qc_sm <- QC.SM(int_data = f_d, order = order, class = class, qc_label = "QC", model = "knn", mode = "subtraction", nbagg = 10, k = 2)
# save
fwrite(qc_sm, "xcms after IPO MVI QC-KNN s.csv", row.names = T)

################################################## Batchwise mode
QC.BW <- function(int_data, order, class, batch, qc_label, model = "bt", mode = "division", nbagg = 10, k = 2) {
  
  library(pbapply)
  
  print("Preprocessing")
  qc_id <-  grep(qc_label, class)
  int_data_ro <- cbind(order, int_data)
  int_data_ro_qc <- int_data_ro[qc_id,]
  int_data_ro_b <- cbind(batch, order, int_data)
  int_data_ro_b_qc <- int_data_ro_b[qc_id,]
  
  wr_batch <- pblapply(1:length(unique(int_data_ro_b_qc$batch)), function(t) filter(int_data_ro_b_qc[,-1], int_data_ro_b_qc$batch == unique(int_data_ro_b_qc$batch)[t]))
  wr_batch_a <- pblapply(1:length(unique(int_data_ro_b$batch)), function(t) filter(int_data_ro_b[,-1], int_data_ro_b$batch == unique(int_data_ro_b$batch)[t]))
  wr_batch_an <- pblapply(1:length(unique(int_data_ro_b$batch)), function(t) filter(int_data_ro_b[,-c(1,2)], int_data_ro_b$batch == unique(int_data_ro_b$batch)[t]))
  
  print("Start Correction")
  if (model == "knn") {
    library(caret)
    fit_wr_bw <- pblapply(1:length(wr_batch), function(x) pblapply(2:ncol(wr_batch[[x]]), function(t) knnreg(wr_batch[[x]][,t] ~ order, data = wr_batch[[x]], k = k)))  
  } else if (model == "dt") {
    library(rpart)
    fit_wr_bw <- pblapply(1:length(wr_batch), function(x) pblapply(2:ncol(wr_batch[[x]]), function(t) rpart(wr_batch[[x]][,t] ~ order, data = wr_batch[[x]])))
  } else { (model == "bt")
    library(ipred)
    fit_wr_bw <- pblapply(1:length(wr_batch), function(x) pblapply(2:ncol(wr_batch[[x]]), function(t) bagging(wr_batch[[x]][,t] ~ order, data = wr_batch[[x]], coob=T, nbagg = nbagg)))
  }
  print("Done Correction")
  predict_wr_bw <- pblapply(1:length(fit_wr_bw), function(x) pblapply(1:length(fit_wr_bw[[1]]), function(t) predict(fit_wr_bw[[x]][[t]], data.frame(order = wr_batch_a[[x]]$order))))
  
  if (mode == "division") {
    val_wr_bw <- pblapply(1:length(predict_wr_bw), function(x) pblapply(1:length(predict_wr_bw[[1]]), function(t) wr_batch_an[[x]][,t]/c(predict_wr_bw[[x]][[t]])))
    res_wr_bw <- pblapply(1:length(val_wr_bw), function(x) as.data.frame(t(do.call(rbind, val_wr_bw[[x]]))))
    res_wr_bw <- as.data.frame(do.call(rbind, res_wr_bw))
    res_wr_bw <- res_wr_bw*1000
    rownames(res_wr_bw) <- rownames(f_d)
    colnames(res_wr_bw) <- colnames(f_d)
  return(res_wr_bw)
    } else { (mode == "subtraction")
    val_wr_bw <- pblapply(1:length(predict_wr_bw), function(x) pblapply(1:length(predict_wr_bw[[1]]), function(t) wr_batch_an[[x]][,t]-c(predict_wr_bw[[x]][[t]])+mean(wr_batch_an[[x]][qc_id,t], na.rm = T)))
    res_wr_bw <- pblapply(1:length(val_wr_bw), function(x) as.data.frame(t(do.call(rbind, val_wr_bw[[x]]))))
    res_wr_bw <- as.data.frame(do.call(rbind, res_wr_bw))
    res_wr_bw <- res_wr_bw
    rownames(res_wr_bw) <- rownames(f_d)
    colnames(res_wr_bw) <- colnames(f_d)
  return(res_wr_bw) 
    }
}

# Parameters of QC.BW function:  
# int_data -> numeric feature table, order -> numeric vector of the run order, batch -> numeric vector of the batch,
# class -> sample group variable, qc_label -> label for QC samples in group, 
# model -> type of the model ("bt" for bagging tree, "dt" - decision tree, "knn" - KNN regression), k = number of neighbours in "knn" model, nbagg -> an integer giving the number of bootstrap replications in "bt" model
# mode -> type of correction factor ("division" or "subtraction")

######################### division mode
######################### perform bagging tree in batchwise mode
qc_bw <- QC.BW(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "bt", mode = "division", nbagg = 10, k = 2)
# save
fwrite(qc_bw, "xcms after IPO MVI QC-BT bw d.csv", row.names = T)

######################### perform decision tree in batchwise mode
qc_bw <- QC.BW(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "dt", mode = "division", nbagg = 10, k = 2)
# save
fwrite(qc_bw, "xcms after IPO MVI QC-DT bw d.csv", row.names = T)

######################### perform knn in batchwise mode
qc_bw <- QC.BW(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "knn", mode = "division", nbagg = 10, k = 2)
# save
fwrite(qc_bw, "xcms after IPO MVI QC-KNN bw d.csv", row.names = T)

######################### subtraction mode
######################### perform bagging tree in batchwise mode
qc_bw <- QC.BW(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "bt", mode = "subtraction", nbagg = 10, k = 2)
# save
fwrite(qc_bw, "xcms after IPO MVI QC-BT bw s.csv", row.names = T)

######################### perform decision tree in batchwise mode
qc_bw <- QC.BW(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "dt", mode = "subtraction", nbagg = 10, k = 2)
# save
fwrite(qc_bw, "xcms after IPO MVI QC-DT bw s.csv", row.names = T)

######################### perform knn in batchwise mode
qc_bw <- QC.BW(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "knn", mode = "subtraction", nbagg = 10, k = 2)
# save
fwrite(qc_bw, "xcms after IPO MVI QC-KNN bw s.csv", row.names = T)

################################################## Two Features mode
QC.TF <- function(int_data, order, batch, class, qc_label, model = "bt", mode = "division", nbagg = 10, k = 2) {
  
  library(pbapply)
  
  qc_id <-  grep(qc_label, class)
  int_data_ro_b <- cbind(order, batch, int_data)
  int_data_ro_b_qc <- int_data_ro_b[qc_id,]
  print("Start Correction")
  
  if (model == "knn") {
    library(caret)
    fit_tf <- pblapply(3:ncol(int_data_ro_b_qc), function(t) knnreg(int_data_ro_b_qc[,t] ~ order+batch, data = int_data_ro_b_qc, k = k)) # KNN
  } else if (model == "dt") {
    library(rpart)
    fit_tf <- pblapply(3:ncol(int_data_ro_b_qc), function(t) rpart(int_data_ro_b_qc[,t] ~ order+batch, data = int_data_ro_b_qc)) # Decision tree
  } else { (model == "bt")
    library(ipred)
    fit_tf <- pblapply(3:ncol(int_data_ro_b_qc), function(t) bagging(int_data_ro_b_qc[,t] ~ order+batch, data = int_data_ro_b_qc, coob=T, nbagg = nbagg)) # Bagging tree
  }
  print("Done correction")
  predict_tf <- pblapply(1:length(fit_tf), function(t) predict(fit_tf[[t]], data.frame(order))) 
  
  if (mode == "division") {
    val_tf <- pblapply(1:ncol(int_data), function(t) int_data[,t]/predict_tf[[t]])
    res_tf <- as.data.frame(t(do.call(rbind, val_tf)))
    rownames(res_tf) <- rownames(int_data)
    colnames(res_tf) <- colnames(int_data)
    res_tf <- res_tf*1000
  return(res_tf)
    } else { (mode == "subtraction")
    val_tf <- pblapply(1:ncol(int_data), function(t) int_data[,t]-predict_tf[[t]]+mean(int_data[qc_id,t], na.rm = T))
    res_tf <- as.data.frame(t(do.call(rbind, val_tf)))
    rownames(res_tf) <- rownames(int_data)
    colnames(res_tf) <- colnames(int_data)
    res_tf <- res_tf
  return(res_tf)
    }
}

# Parameters of QC.TF function:  
# int_data -> numeric feature table, order -> numeric vector of the run order, batch -> numeric vector of the batch, 
# class -> sample group variable, qc_label -> label for QC samples in group, 
# model -> type of the model ("bt" for bagging tree, "dt" - decision tree, "knn" - KNN regression), k = number of neighbours in "knn" model, nbagg -> an integer giving the number of bootstrap replications in "bt" model
# mode -> type of correction factor ("division" or "subtraction")

######################### division mode
######################### perform bagging tree in two features mode
qc_tf <- QC.TF(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "bt", mode = "division", nbagg = 10, k = 2)
# save
fwrite(qc_tf, "xcms after IPO MVI QC-BT tf d.csv", row.names = T)

######################### perform decision tree in two features mode
qc_tf <- QC.TF(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "dt", mode = "division", nbagg = 10, k = 2)
# save
fwrite(qc_tf, "xcms after IPO MVI QC-DT tf d.csv", row.names = T)

######################### perform knn in two features mode
qc_tf <- QC.TF(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "knn", mode = "division", nbagg = 10, k = 2)
# save
fwrite(qc_tf, "xcms after IPO MVI QC-KNN tf d.csv", row.names = T)

######################### subtraction mode
######################### perform bagging tree in two features mode
qc_tf <- QC.TF(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "bt", mode = "subtraction", nbagg = 10, k = 2)
# save
fwrite(qc_tf, "xcms after IPO MVI QC-BT tf s.csv", row.names = T)

######################### perform decision tree in two features mode
qc_tf <- QC.TF(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "dt", mode = "subtraction", nbagg = 10, k = 2)
# save
fwrite(qc_tf, "xcms after IPO MVI QC-DT tf s.csv", row.names = T)

######################### perform knn in two features mode
qc_tf <- QC.TF(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "knn", mode = "subtraction", nbagg = 10, k = 2)
# save
fwrite(qc_tf, "xcms after IPO MVI QC-KNN tf s.csv", row.names = T)

###############################################################################
########################################## QC-GBM (xgboost/catboost)
###############################################################################

# generate batch data
batch <- ds$b_id
batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# generate group data
class <- as.character(ds$n_gr_t)

# generate run order data
order <- as.numeric(ds$ro_id)

# Sample data
s_d <- data.frame(cbind(sample = rownames(ds), batch, class, order))
rownames(s_d) <- rownames(ds)

# Feature data
f_d <- data.frame(ds[,-c(1:7)])

# Parameters for catboost
fit_params <- list(
  iterations = 100,
  loss_function = 'RMSE',
  border_count = 32,
  depth = 2,
  learning_rate = 0.03,
  l2_leaf_reg = 3.5,
  train_dir = 'train_dir')

################################################## Single mode
QC.SM.GB <- function(int_data, order, class, qc_label, model = "xgboost", mode = "division", max.depth = 2, nrounds = 100, params) {
  
  library(pbapply)
  
  qc_id <-  grep(qc_label, class)
  int_data_ro <- cbind(order, int_data)
  int_data_ro_qc <- int_data_ro[qc_id,]
  ro_qc <- int_data_ro[qc_id,1]
  ro <- order
  print("Start Correction")
  
  if (model == "catboost") {
    library(catboost)
    catb_train <- pblapply(2:ncol(int_data_ro_qc), function(t) catboost.load_pool(data = as.matrix(as.numeric(ro_qc)), label = as.matrix(as.numeric(int_data_ro_qc[,t]))))
    catb_test <- pblapply(2:ncol(int_data_ro), function(t) catboost.load_pool(data = as.matrix(as.numeric(ro)), label = as.matrix(as.numeric(int_data_ro[,t]))))
    
    fit_catb <- pblapply(1:length(catb_train), function(t) catboost.train(catb_train[[t]], params = params))
    predict_gb <- pblapply(1:length(fit_catb), function(t) catboost.predict(fit_catb[[t]], catb_test[[t]], verbose = F)) 
  } else {
    library(xgboost)
    xgb_train <- pblapply(2:ncol(int_data_ro_qc), function(t) xgb.DMatrix(data = as.matrix(as.numeric(ro_qc)), label = as.matrix(as.numeric(int_data_ro_qc[,t]))))
    xgb_test <- pblapply(2:ncol(int_data_ro), function(t) xgb.DMatrix(data = as.matrix(as.numeric(ro)), label = as.matrix(as.numeric(int_data_ro[,t]))))
    
    fit_xgb <- pblapply(1:length(xgb_train), function(t) xgboost(data = xgb_train[[t]], max.depth = max.depth, nrounds = nrounds, verbose = 0))
    predict_gb <- pblapply(1:length(fit_xgb), function(t) predict(fit_xgb[[t]], xgb_test[[t]]))
  } 
  print("Done Correction")
  
  if (mode == "division") {
    val_sm <- pblapply(1:ncol(int_data), function(t) int_data[,t]/predict_gb[[t]])
    res_sm <- as.data.frame(t(do.call(rbind, val_sm)))
    rownames(res_sm) <- rownames(int_data)
    colnames(res_sm) <- colnames(int_data)
    res_sm <- res_sm*1000
  return(res_sm)
  } else { (mode == "subtraction")
    val_sm <- pblapply(1:ncol(int_data), function(t) int_data[,t]-predict_gb[[t]]+mean(int_data[qc_id,t], na.rm = T))
    res_sm <- as.data.frame(t(do.call(rbind, val_sm)))
    rownames(res_sm) <- rownames(int_data)
    colnames(res_sm) <- colnames(int_data)
    res_sm <- res_sm
  return(res_sm)
  }
}

# Parameters of QC.SM.GB function:  
# int_data -> numeric feature table, order -> numeric vector of the run order, class -> sample group variable, 
# qc_label -> label for QC samples in group, 
# model -> type of the model ("xgboost" for xgboost, "catboost" - catboost), 
# max_depth -> maximum depth of a tree in xgboost model, nrounds	-> max number of boosting iterations in xgboost, params -> parameters for the catboost model
# mode -> type of correction factor ("division" or "subtraction")

######################### division mode
######################### perform xgboost in single mode
qc_sm <- QC.SM.GB(int_data = f_d, order = order, class = class, qc_label = "QC", model = "xgboost", mode = "division", max.depth = 2, nrounds = 100, params = fit_params)
# save
fwrite(qc_sm, "xcms after IPO MVI QC-XGB d.csv", row.names = T)

######################### perform catboost in single mode
qc_sm <- QC.SM.GB(int_data = f_d, order = order, class = class, qc_label = "QC", model = "catboost", mode = "division", max.depth = 2, nrounds = 100, params = fit_params)
# save
fwrite(qc_sm, "xcms after IPO MVI QC-CTB d.csv", row.names = T)

######################### subtraction mode
######################### perform xgboost in single mode
qc_sm <- QC.SM.GB(int_data = f_d, order = order, class = class, qc_label = "QC", model = "xgboost", mode = "subtraction", max.depth = 2, nrounds = 100, params = fit_params)
# save
fwrite(qc_sm, "xcms after IPO MVI QC-XGB s.csv", row.names = T)

######################### perform catboost in single mode
qc_sm <- QC.SM.GB(int_data = f_d, order = order, class = class, qc_label = "QC", model = "catboost", mode = "subtraction", max.depth = 2, nrounds = 100, params = fit_params)
# save
fwrite(qc_sm, "xcms after IPO MVI QC-CTB s.csv", row.names = T)

################################################## Batchwise mode
QC.BW.GB <- function(int_data, order, class, batch, qc_label, model = "xgboost", mode = "division", max.depth = 2, nrounds = 100, params) {
  
  library(pbapply)
  
  print("Preprocessing")
  qc_id <-  grep(qc_label, class)
  int_data_ro <- cbind(order, int_data)
  int_data_ro_qc <- int_data_ro[qc_id,]
  int_data_ro_b <- cbind(batch, order, int_data)
  int_data_ro_b_qc <- int_data_ro_b[qc_id,]
  ro_qc <- int_data_ro[qc_id,1]
  ro <- order
  ro_qc_b <- pblapply(1:length(unique(int_data_ro_b_qc$batch)), function(t) filter(int_data_ro_b_qc[,-1], int_data_ro_b_qc$batch == unique(int_data_ro_b_qc$batch)[t])[,1])
  ro_b <- pblapply(1:length(unique(int_data_ro_b$batch)), function(t) filter(int_data_ro_b[,-1], int_data_ro_b$batch == unique(int_data_ro_b$batch)[t])[,1])
  wr_batch <- pblapply(1:length(unique(int_data_ro_b_qc$batch)), function(t) filter(int_data_ro_b_qc[,-1], int_data_ro_b_qc$batch == unique(int_data_ro_b_qc$batch)[t]))
  wr_batch_a <- pblapply(1:length(unique(int_data_ro_b$batch)), function(t) filter(int_data_ro_b[,-1], int_data_ro_b$batch == unique(int_data_ro_b$batch)[t]))
  wr_batch_an <- pblapply(1:length(unique(int_data_ro_b$batch)), function(t) filter(int_data_ro_b[,-c(1,2)], int_data_ro_b$batch == unique(int_data_ro_b$batch)[t]))
  
  print("Start Correction")
  if (model == "catboost") {
    library(catboost)
    ctb_train <- pblapply(1:length(ro_qc_b), function(x) pblapply(2:ncol(wr_batch[[1]]), function(t) catboost.load_pool(data = as.matrix(as.numeric(ro_qc_b[[x]])), label = as.matrix(as.numeric(wr_batch[[x]][,t])))))
    ctb_test <- pblapply(1:length(ro_b), function(x) pblapply(2:ncol(wr_batch_a[[1]]), function(t) catboost.load_pool(data = as.matrix(as.numeric(ro_b[[x]])), label = as.matrix(as.numeric(wr_batch_a[[x]][,t])))))
    
    fit_wr_bw <- pblapply(1:length(wr_batch), function(x) pblapply(1:length(ctb_train[[1]]), function(t) catboost.train(ctb_train[[x]][[t]], params = params)))
    predict_wr_bw <- pblapply(1:length(wr_batch), function(x) pblapply(1:length(ctb_test[[1]]), function(t) catboost.predict(fit_wr_bw[[x]][[t]], ctb_test[[x]][[t]], verbose = F)))
  } else {
    library(xgboost)
    xgb_train <- pblapply(1:length(ro_qc_b), function(x) pblapply(2:ncol(wr_batch[[1]]), function(t) xgb.DMatrix(data = as.matrix(as.numeric(ro_qc_b[[x]])), label = as.matrix(as.numeric(wr_batch[[x]][,t])))))
    xgb_test <- pblapply(1:length(ro_b), function(x) pblapply(2:ncol(wr_batch_a[[1]]), function(t) xgb.DMatrix(data = as.matrix(as.numeric(ro_b[[x]])), label = as.matrix(as.numeric(wr_batch_a[[x]][,t])))))
    
    fit_wr_bw <- pblapply(1:length(wr_batch), function(x) pblapply(1:length(xgb_train[[1]]), function(t) xgboost(data = xgb_train[[x]][[t]], max.depth = max.depth, nrounds = nrounds, verbose = 0)))
    predict_wr_bw <- pblapply(1:length(wr_batch), function(x) pblapply(1:length(xgb_test[[1]]), function(t) predict(fit_wr_bw[[x]][[t]], xgb_test[[x]][[t]])))
  } 
  print("Done Correction")
  
  if (mode == "division") {
    val_wr_bw <- pblapply(1:length(predict_wr_bw), function(x) pblapply(1:length(predict_wr_bw[[1]]), function(t) wr_batch_an[[x]][,t]/c(predict_wr_bw[[x]][[t]])))
    res_wr_bw <- pblapply(1:length(val_wr_bw), function(x) as.data.frame(t(do.call(rbind, val_wr_bw[[x]]))))
    res_wr_bw <- as.data.frame(do.call(rbind, res_wr_bw))
    res_wr_bw <- res_wr_bw*1000
    rownames(res_wr_bw) <- rownames(int_data)
    colnames(res_wr_bw) <- colnames(int_data)
  return(res_wr_bw)
    } else { (mode == "subtraction")
    val_wr_bw <- pblapply(1:length(predict_wr_bw), function(x) pblapply(1:length(predict_wr_bw[[1]]), function(t) wr_batch_an[[x]][,t]-c(predict_wr_bw[[x]][[t]])+mean(wr_batch_an[[x]][qc_id,t], na.rm = T)))
    res_wr_bw <- pblapply(1:length(val_wr_bw), function(x) as.data.frame(t(do.call(rbind, val_wr_bw[[x]]))))
    res_wr_bw <- as.data.frame(do.call(rbind, res_wr_bw))
    res_wr_bw <- res_wr_bw
    rownames(res_wr_bw) <- rownames(int_data)
    colnames(res_wr_bw) <- colnames(int_data)
  return(res_wr_bw)
    }
}

# Parameters of QC.BW.GB function:  
# int_data -> numeric feature table, order -> numeric vector of the run order, class -> sample group variable, 
# batch -> numeric vector of the batch, qc_label -> label for QC samples in group, 
# model -> type of the model ("xgboost" for xgboost, "catboost" - catboost), 
# max_depth -> maximum depth of a tree in xgboost model, nrounds	-> max number of boosting iterations in xgboost, params -> parameters for the catboost model
# mode -> type of correction factor ("division" or "subtraction")

######################### division mode
######################### perform xgboost in batchwise mode
qc_bw <- QC.BW.GB(int_data = f_d, order = order, class = class, batch = batch, qc_label = "QC", model = "xgboost", mode = "division", max.depth = 2, nrounds = 100, params = fit_params)
# save
fwrite(qc_bw, "xcms after IPO MVI QC-XGB bw d.csv", row.names = T)

######################### perform catboost in batchwise mode
qc_bw <- QC.BW.GB(int_data = f_d, order = order, class = class, batch = batch, qc_label = "QC", model = "catboost", mode = "division", max.depth = 2, nrounds = 100, params = fit_params)
# save
fwrite(qc_bw, "xcms after IPO MVI QC-CTB bw d.csv", row.names = T)

######################### subtraction mode
######################### perform xgboost in batchwise mode
qc_bw <- QC.BW.GB(int_data = f_d, order = order, class = class, batch = batch, qc_label = "QC", model = "xgboost", mode = "subtraction", max.depth = 2, nrounds = 100, params = fit_params)
# save
fwrite(qc_bw, "xcms after IPO MVI QC-XGB bw s.csv", row.names = T)

######################### perform catboost in batchwise mode
qc_bw <- QC.BW.GB(int_data = f_d, order = order, class = class, batch = batch, qc_label = "QC", model = "catboost", mode = "subtraction", max.depth = 2, nrounds = 100, params = fit_params)
# save
fwrite(qc_bw, "xcms after IPO MVI QC-CTB bw s.csv", row.names = T)

################################################## Two Features mode
QC.TF.GB <- function(int_data, order, batch, class, qc_label, model = "xgboost", mode = "division", max.depth = 2, nrounds = 100, params) {
  
  library(pbapply)
  
  qc_id <-  grep(qc_label, class)
  int_data_ro <- cbind(order, int_data)
  int_data_ro_qc <- int_data_ro[qc_id,]
  ro_qc <- int_data_ro[qc_id,1]
  ro <- order
  batch_qc <- batch[c(qc_id)]
  batch <- batch
  print("Start Correction")
  
  if (model == "catboost") {
    library(catboost)
    catb_train <- pblapply(2:ncol(int_data_ro_qc), function(t) catboost.load_pool(data = cbind(as.numeric(ro_qc), as.numeric(batch_qc)), label = as.matrix(as.numeric(int_data_ro_qc[,t]))))
    catb_test <- pblapply(2:ncol(int_data_ro), function(t) catboost.load_pool(data = cbind(as.numeric(ro), as.numeric(batch)), label = as.matrix(as.numeric(int_data_ro[,t]))))
    
    fit_catb <- pblapply(1:length(catb_train), function(t) catboost.train(catb_train[[t]], params = params))
    predict_gb <- pblapply(1:length(fit_catb), function(t) catboost.predict(fit_catb[[t]], catb_test[[t]], verbose = F)) 
  } else {
    library(xgboost)
    xgb_train <- pblapply(2:ncol(int_data_ro_qc), function(t) xgb.DMatrix(data = cbind(as.numeric(ro_qc), as.numeric(batch_qc)), label = as.matrix(as.numeric(int_data_ro_qc[,t]))))
    xgb_test <- pblapply(2:ncol(int_data_ro), function(t) xgb.DMatrix(data = cbind(as.numeric(ro), as.numeric(batch)), label = as.matrix(as.numeric(int_data_ro[,t]))))
    
    fit_xgb <- pblapply(1:length(xgb_train), function(t) xgboost(data = xgb_train[[t]], max.depth = max.depth, nrounds = nrounds, verbose = 0))
    predict_gb <- pblapply(1:length(fit_xgb), function(t) predict(fit_xgb[[t]], xgb_test[[t]]))
  } 
  print("Done Correction")
  
  if (mode == "division") {
    val_tf <- pblapply(1:ncol(int_data), function(t) int_data[,t]/predict_gb[[t]])
    res_tf <- as.data.frame(t(do.call(rbind, val_tf)))
    rownames(res_tf) <- rownames(int_data)
    colnames(res_tf) <- colnames(int_data)
    res_tf <- res_tf*1000
  return(res_tf)
    } else { (mode == "subtraction")
    val_tf <- pblapply(1:ncol(int_data), function(t) int_data[,t]-predict_gb[[t]]+mean(int_data[qc_id,t], na.rm = T))
    res_tf <- as.data.frame(t(do.call(rbind, val_tf)))
    rownames(res_tf) <- rownames(int_data)
    colnames(res_tf) <- colnames(int_data)
    res_tf <- res_tf
  return(res_tf)
    }
}

# Parameters of QC.TF.GB function:  
# int_data -> numeric feature table, order -> numeric vector of the run order, batch -> numeric vector of the batch,
# class -> sample group variable, qc_label -> label for QC samples in group, 
# model -> type of the model ("xgboost" for xgboost, "catboost" - catboost), 
# max_depth -> maximum depth of a tree in xgboost model, nrounds	-> max number of boosting iterations in xgboost, params -> parameters for the catboost model
# mode -> type of correction factor ("division" or "subtraction")

######################### division mode
######################### perform xgboost in two features mode
qc_tf <- QC.TF.GB(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "xgboost", mode = "division", max.depth = 2, nrounds = 100, params = fit_params)
# save
fwrite(qc_tf, "xcms after IPO MVI QC-XGB tf d.csv", row.names = T)

######################### perform catboost in two features mode
qc_tf <- QC.TF.GB(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "catboost", mode = "division", max.depth = 2, nrounds = 100, params = fit_params)
# save
fwrite(qc_tf, "xcms after IPO MVI QC-CTB tf d.csv", row.names = T)

######################### subtraction mode
######################### perform xgboost in two features mode
qc_tf <- QC.TF.GB(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "xgboost", mode = "subtraction", max.depth = 2, nrounds = 100, params = fit_params)
# save
fwrite(qc_tf, "xcms after IPO MVI QC-XGB tf s.csv", row.names = T)

######################### perform catboost in two features mode
qc_tf <- QC.TF.GB(int_data = f_d, order = order, batch = batch, class = class, qc_label = "QC", model = "catboost", mode = "subtraction", max.depth = 2, nrounds = 100, params = fit_params)
# save
fwrite(qc_tf, "xcms after IPO MVI QC-CTB tf s.csv", row.names = T)

###############################################################################
########################################## Median Between Batches (QC-norm)
###############################################################################

# generate batch data
batch <- ds$b_id
batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# generate group data
class <- as.character(ds$n_gr_t)
qc_id <- which(class == "QC")
class[-qc_id] <- "Subject"

# generate data
ds_bbc <- cbind(batch = ds$b_id, ds[,-c(1:7)])
ds_bbc$batch <- batch

######################### Perform division-based QC-norm 
QC.NORM <- function(data, qc_id){
  library(dplyr)
  ds_bbc <- data
  ds_bbc_qc <- ds_bbc[qc_id,]
  b_b_c_subsets <- lapply(1:length(unique(ds_bbc[,1])), function(y) filter(ds_bbc[,-1], ds_bbc$batch == unique(ds_bbc[,1])[y])) # list of subsets by batches for all samples
  b_b_c_subsets_qc <- lapply(1:length(unique(ds_bbc_qc[,1])), function(y) filter(ds_bbc_qc[,-1], ds_bbc_qc$batch == unique(ds_bbc_qc[,1])[y])) # list of subsets by batches for QC samples
  b_b_c_factor <- lapply(1:length(unique(ds_bbc_qc[,1])), function(y) sapply(1:ncol(b_b_c_subsets_qc[[1]]), function(z) median(b_b_c_subsets_qc[[y]][,z], na.rm = T)/median(ds_bbc_qc[,(z+1)], na.rm = T))) # calculate factor for every feature on QC samples
  b_b_c_result_l <- lapply(1:length(unique(ds_bbc[,1])), function(y) sapply(1:ncol(b_b_c_subsets[[1]]), function(z) b_b_c_subsets[[y]][,z]/b_b_c_factor[[y]][z])) # results by batch
  b_b_c_result <- data.frame(do.call("rbind",b_b_c_result_l)) # results in data frame
  rownames(b_b_c_result) <- rownames(ds_bbc)
  colnames(b_b_c_result) <- colnames(ds_bbc[,-c(1)])
  return(b_b_c_result)
}

ds_qc_norm <- QC.NORM(data = ds_bbc, qc_id = qc_id)

# save
fwrite(ds_qc_norm, "xcms after IPO MVI QC-RF + QC-NORM d.csv", row.names = T)

######################### Perform subtraction-based QC-norm 
QC.NORM <- function(data, qc_id){
  library(dplyr)
  ds_bbc <- data
  ds_bbc_qc <- ds_bbc[qc_id,]
  b_b_c_subsets <- lapply(1:length(unique(ds_bbc[,1])), function(y) filter(ds_bbc[,-1], ds_bbc$batch == unique(ds_bbc[,1])[y])) # list of subsets by batches for all samples
  b_b_c_subsets_qc <- lapply(1:length(unique(ds_bbc_qc[,1])), function(y) filter(ds_bbc_qc[,-1], ds_bbc_qc$batch == unique(ds_bbc_qc[,1])[y])) # list of subsets by batches for QC samples
  b_b_c_factor <- lapply(1:length(unique(ds_bbc_qc[,1])), function(y) sapply(1:ncol(b_b_c_subsets_qc[[1]]), function(z) median(b_b_c_subsets_qc[[y]][,z], na.rm = T)-median(ds_bbc_qc[,(z+1)], na.rm = T))) # calculate factor for every feature on QC samples
  b_b_c_result_l <- lapply(1:length(unique(ds_bbc[,1])), function(y) sapply(1:ncol(b_b_c_subsets[[1]]), function(z) b_b_c_subsets[[y]][,z]-b_b_c_factor[[y]][z])) # results by batch
  b_b_c_result <- data.frame(do.call("rbind",b_b_c_result_l)) # results in data frame
  rownames(b_b_c_result) <- rownames(ds_bbc)
  colnames(b_b_c_result) <- colnames(ds_bbc[,-c(1)])
  return(b_b_c_result)
}

ds_qc_norm <- QC.NORM(data = ds_bbc, qc_id = qc_id)

# save
fwrite(ds_qc_norm, "xcms after IPO MVI QC-RF + QC-NORM s.csv", row.names = T)

###############################################################################
########################################## RUVs metabolites based methods
###############################################################################

library(NormalizeMets)
library(stringr)

# generate batch data
batch <- ds$b_id
batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# generate featuredata/sampledata for RUVs
featuredata <- as.data.frame(ds[,-c(1:7)])
colnames(featuredata) <- colnames(dsr)
sampledata <- data.frame(type = as.factor(ds$n_gr_t), experiment = as.factor(ds$s_id_p2), batch = as.factor(batch))

# generate data of ISs/RUVs
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

######################### Perform ruvrand correction
logdata <- LogTransform(featuredata) # Log transformation
ruvrand <- NormQcmets(logdata$featuredata, method = "ruvrand", qcmets = qcmets, k = 1) # adjust k for your data manually (set high k and you can readjust)
ruvrand <- as.data.frame(ruvrand$featuredata)

# save
fwrite(ruvrand, "xcms after IPO MVI ruvrand.csv", row.names = T)

######################### Perform ruvrandclust correction
logdata <- LogTransform(featuredata) # Log transformation
ruvrandclust <- NormQcmets(logdata$featuredata, method = "ruvrandclust", qcmets = qcmets, k = 1, p = length(unique(sampledata$type))) # adjust k (you can readjust) for your data and also maxIter, nUpdate and p
ruvrandclust <- as.data.frame(ruvrandclust$featuredata)

# save
fwrite(ruvrandclust, "xcms after IPO MVI ruvrandclust.csv", row.names = T)

###############################################################################
########################################## ISs metabolites based methods
###############################################################################

library(NormalizeMets)
library(stringr)
library(tibble)

# generate batch data
batch <- ds$b_id
batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# generate featuredata/sampledata
featuredata <- as.data.frame(ds[,-c(1:7)])
colnames(featuredata) <- colnames(dsr)
sampledata <- data.frame(type = as.factor(ds$n_gr_t), experiment = as.factor(ds$s_id_p2), batch = as.factor(batch))

# generate data of ISs/RUVs
# search m/z of IS
cn <- colnames(featuredata)
is_ind <- which(str_detect(string = cn, pattern = "340.15")==T) # write in pattern m/z value of IS. If many ISs -> repeat this and save in 1 vector
isvec <- featuredata[,is_ind]

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

######################### Perform SIS correction
logdata <- LogTransform(featuredata) # Log transformation
sis <- NormQcmets(logdata$featuredata, method = "is", isvec=isvec) 
sis <- as.data.frame(sis$featuredata)

# save
fwrite(sis, "xcms after IPO MVI SIS.csv", row.names = T)

######################### Perform NOMIS correction
logdata <- LogTransform(featuredata) # Log transformation
nomis <- NormQcmets(logdata$featuredata, method = "nomis", qcmets = qcmets) 
nomis <- as.data.frame(nomis$featuredata)

# save
fwrite(nomis, "xcms after IPO MVI NOMIS.csv", row.names = T)

######################### Perform CCMN correction
logdata <- LogTransform(featuredata) # Log transformation
ccmn <- NormQcmets(logdata$featuredata, method = "ccmn", qcmets = qcmets, ncomp = 2, factors = sampledata$type) # adjust ncomp and factors = sampledata$type or sampledata$batch
ccmn <- as.data.frame(ccmn$featuredata)

# save
fwrite(ccmn, "xcms after IPO MVI CCMN.csv", row.names = T)

######################### Perform BMIS correction
#set columns that contain IS information. IS columns must be the first response value containing columns in the dataset. The value "0" has to be kept in the vector/dataframe as it represents non-normalization (use of raw data) for B-MIS.
#here,e.g.: 0 stands for non-normalization, values 1-5 indicate 5 IS in the first 5 rows.
ISvect<-as.vector(qcmets)
IS<-data.frame(t(ISvect))    
logdata <- LogTransform(featuredata) # Log transformation

#compute normalization with all IS
#add column with value 1 for all samples for non-normalization
data.raw <- cbind(Class = sampledata$type, as.data.frame(logdata$featuredata))
data.raw.BMIS<-add_column(data.raw,"noIS"=1:length(sampledata$type),.after="Class")
data.raw.BMIS$noIS<-1
data.raw.BMIS$Class <- as.numeric(data.raw.BMIS$Class) 

#calculate CVs for all potential ISs (via mean and standard deviation)
#class-specific mean after normalizatin with all IS
mean.BMIS<-data.frame(matrix(nrow=length(IS), ncol=ncol(data.raw.BMIS)-1))
for (i in 1:length(IS))
{mean.BMIS.norm <- paste("IS", i, sep = "")
assign(mean.BMIS.norm, data.frame(c(data.raw.BMIS[1],data.raw.BMIS[2:ncol(data.raw.BMIS)]/data.raw.BMIS[,IS[,i]+2]), check.names = FALSE))
Subset<-subset(get(mean.BMIS.norm), Class==2) #change value for Class here if other class than QC or for Class = QC (Class=1 or 2 or 3, etc.) should be processed
mean.BMIS[i,]<-(apply(Subset[,2:ncol(Subset)], 2, FUN=mean))}
#class-specific standard deviation after normalizatin with all IS
sd.BMIS<-data.frame(matrix(nrow=length(IS), ncol=ncol(data.raw.BMIS)-1))
for (i in 1:length(IS))
{sd.BMIS.norm <- paste("IS", i, sep = "")
assign(sd.BMIS.norm, data.frame(c(data.raw.BMIS[1],data.raw.BMIS[2:ncol(data.raw.BMIS)]/data.raw.BMIS[,IS[,i]+2]), check.names = FALSE))
#change value for Class here if other class than QC (Class=3) should be processed
Subset<-subset(get(sd.BMIS.norm), Class==3)
sd.BMIS[i,]<-(apply(Subset[,2:ncol(Subset)], 2, FUN=sd))}
CV.norm.BMIS<-sd.BMIS/mean.BMIS*100
colnames(CV.norm.BMIS)<-colnames(data.raw.BMIS[-1])

#Find IS that produces minimum class-specific CV
BMIS.assignment<-data.frame(matrix(nrow=2, ncol=ncol(data.raw.BMIS)-1))
BMIS.assignment[1,]<-apply(CV.norm.BMIS[,1:ncol(CV.norm.BMIS)], 2, FUN=min)
for (i in 1:ncol(BMIS.assignment))
{BMIS.assignment[2,i]<-rownames(CV.norm.BMIS)[which.min(apply(CV.norm.BMIS[i],MARGIN=1,min))]}
rownames(BMIS.assignment) <- c("minimum CV in %","IS ID")
colnames(BMIS.assignment) <- colnames(data.raw.BMIS[-1])

#Normalize each feature with assigned B-MIS
BMIS.assignment<-cbind(1,BMIS.assignment)
colnames(BMIS.assignment)<-colnames(data.raw.BMIS)
data.BMIS<-data.matrix(rbind(BMIS.assignment[2,],data.raw.BMIS))

for (i in length(IS):(length(data.BMIS[1,])-2))
{for (z in 1:length(IS))
{
  if(data.BMIS[1,i+2]==z)
  {data.BMIS[,i+2]<-data.BMIS[,i+2]/data.BMIS[,z+1]}}}
data.BMIS<-data.frame(data.BMIS[-1,-2])
colnames(data.BMIS)<-colnames(data.raw)

bmis <- as.data.frame(data.BMIS[,-1])

# save
fwrite(bmis, "xcms after IPO MVI BMIS.csv", row.names = T)

##############################################################################################################################################################
# Evaluate correction
##############################################################################################################################################################

# Content:
# Generate data for evaluation
# RSD calculation for QC
# Pearson Correlation of QCs
# mean Bhattacharyya distance
# mean Dissimilarity Score
# mean Silhouette Score
# Classif. Acc. by HCA on PCA
# PVCA
# gPCA
# PC-PR2
# ANOVA
# One-Sample Test
# all numeric metrics
# Data Projection (PCA, HCA, scatter plot, box-plot)

###############################################################################
########################################## Generate data for evaluation
###############################################################################

# setup environment
library(data.table) 
library(dplyr)
library(stringr)
setwd("D:/...")

###################### generate feature table with metadata (as in first section)
# load data
dsr <-as.data.frame(fread(input = "xcms after IPO MVI QC-XGB.csv", header=T)) # first column with all metadata
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]
# metadata generating
rname <- rownames(dsr) 
rname <- str_remove(rname, ".CDF")
all_id <- sapply(1:length(rname), function(y) unlist(str_split(rname[y], " "))) 
ro_id <- as.numeric(unlist(lapply(all_id, function(y) unlist(y[1])))) 
b_id <- unlist(lapply(all_id, function(y) unlist(y[3]))) 
s_id <- unlist(lapply(all_id, function(y) unlist(y[2]))) 
p1_id <- unlist(lapply(all_id, function(y) unlist(y[4]))) 
p2_id <- unlist(lapply(all_id, function(y) unlist(y[5]))) 
p_all_id <- sapply(1:length(p1_id), function(y) ifelse(is.na(p2_id[y]), p1_id[y] ,paste(p1_id[y], p2_id[y], sep = " "))) 
un_raw_l <- unique(gsub("[[:digit:]]", "", p1_id)) 
true_l <- c("QC", "TG", "CG") 
rbind(un_raw_l,true_l) 
raw_l <- gsub("[[:digit:]]", "", p1_id) 
n_gr_t <- as.character(match(raw_l, un_raw_l)) 
for (i in 1:length(unique(n_gr_t))) {
  n_gr_t <- str_replace(n_gr_t, unique(n_gr_t)[i], true_l[i]) } 
s_id_p <- sapply(1:length(s_id), function(y) unlist(str_split(s_id[y], "_"))) 
s_id_p2 <- unlist(lapply(s_id_p, function(y) unlist(y[1]))) 
all_meta_data<- as.data.frame(cbind(n_gr_t, ro_id, b_id, s_id, p_all_id, raw_l, s_id_p2)) 
# join features and metadata and order
ds <- data.frame(cbind(all_meta_data, dsr))
ds$ro_id <- as.numeric(ds$ro_id)
ds <- ds[order(ds$ro_id, decreasing = F),] 

###############################################################################
########################################## RSD calculation for QC
###############################################################################

# after processing from "Generate data for evaluation"
ds_rsd <- subset(ds, n_gr_t == "QC")
ds_rsd <- ds_rsd[,-c(1,2,4,5,6,7)]
colnames(ds_rsd)[1] <- "batch"

# rename batch index to numeric value
ds_rsd[,1] <- str_remove(ds_rsd[,1], "b")
ds_rsd[,1] <- as.numeric(ds_rsd[,1])
ds_rsd <- ds_rsd[order(ds_rsd$batch, decreasing = F),]

# Select batch by pattern
tn <- ds_rsd$batch # batch variables
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
ds_rsd1 <- ds_rsd[,-1] # sample in row only metabolites variables
rsd <- RSD_by_batch(ds_rsd1) 

# count in % by cut-off
# cutoff 1
cutoff1 <- 30 # adjust to your data
RSD_by_batch_cutoff1 <- c()
for (i in (1:length(rsd))){
  RSD_by_batch_cutoff1[i] <- round(length(which(rsd[[i]]$rsd <= cutoff1))/length(rsd[[1]]$rsd)*100,0)}
RSD_by_batch_cutoff1 <- data.frame("batch" = n,RSD_by_batch_cutoff1)
colnames(RSD_by_batch_cutoff1)[2] <- c(paste(cutoff1, "% RSD"))

# cutoff 2 
cutoff2 <- 50 # adjust to your data
RSD_by_batch_cutoff2 <- c()
for (i in (1:length(rsd))){
  RSD_by_batch_cutoff2[i] <- round(length(which(rsd[[i]]$rsd <= cutoff2))/length(rsd[[1]]$rsd)*100,0)}
RSD_by_batch_cutoff2 <- data.frame("batch" = n,RSD_by_batch_cutoff2)
colnames(RSD_by_batch_cutoff2)[2] <- c(paste(cutoff2, "% RSD"))

# result
result_RSD <- rbind(t(RSD_by_batch_cutoff1), t(RSD_by_batch_cutoff2))
result_RSD

###############################################################################
########################################## Pearson Correlation of QCs
###############################################################################

# after processing from "Generate data for evaluation"
ds_cor <- subset(ds, n_gr_t == "QC")
ds_cor <- data.frame(t(ds_cor[,-c(1:7)]))

cor_coef <- cor(ds_cor, method = "pearson") # adjust to your data
diag(cor_coef) <- 0
m_c_c <- round(mean(cor_coef, na.rm = T),3) 
m_c_c

###############################################################################
########################################## mean Bhattacharyya distance
###############################################################################

library(fpc)
library(FactoMineR)

# feature data
ds_bd <- subset(ds, n_gr_t == "QC")

# generate group data by batch
pca.class <- as.character(ds_bd$b_id)

# perform PCA
pca.ds <- PCA(ds_bd[,-c(1:7)], scale.unit = T, graph = F)
df.pca.x <- as.data.frame(pca.ds$ind$coord[,1:2]) # or pca.ds$X[,1:n], where n - number of components

# perform
scores <- cbind(pca.class, df.pca.x)
means <- lapply(1:length(unique(scores[,1])), function(y) colMeans(scores[scores$pca.class == unique(scores[,1])[y],-1] , na.rm = T))
covs <- lapply(1:length(unique(scores[,1])), function(y) cov(scores[scores$pca.class == unique(scores[,1])[y],-1] , ))

batch.dist_b <- matrix(0, length(unique(scores[,1])), length(unique(scores[,1])))
for (i in 2:length(unique(scores[,1]))) {
  for (j in 1:(i-1)) {
    batch.dist_b[j, i] <- bhattacharyya.dist(means[[j]],
                                             means[[i]],
                                             covs[[j]],
                                             covs[[i]])
  }}

m_b_d <- round(mean(batch.dist_b[col(batch.dist_b) > row(batch.dist_b)]), 2)
m_b_d

###############################################################################
########################################## mean Dissimilarity Score
###############################################################################

library(cluster)
library(FactoMineR)

# feature data
ds_bd <- subset(ds, n_gr_t == "QC")

# perform PCA
pca.ds <- PCA(ds_bd[,-c(1:7)], scale.unit = T, graph = F)
df.pca.x <- as.data.frame(pca.ds$ind$coord[,1:2]) # or pca.ds$X[,1:n], where n - number of components

# perform
dis <- daisy(df.pca.x, metric = "euclidean", stand = F)
m_d_s <- round(mean(dis),2)
m_d_s

###############################################################################
########################################## mean Silhouette Score
###############################################################################

library(cluster)
library(FactoMineR)

# feature data
ds_bd <- subset(ds, n_gr_t == "QC")

# perform PCA
pca.ds <- PCA(ds_bd[,-c(1:7)], scale.unit = T, graph = F)
df.pca.x <- as.data.frame(pca.ds$ind$coord[,1:2]) # or pca.ds$X[,1:n], where n - number of components

# generate group data by batch
batch <- ds_bd$b_id
batch <- str_remove(batch, "b")
pca.class <- as.numeric(batch)

# calculate Silhouette
sil <- silhouette(pca.class, daisy(df.pca.x, metric = "euclidean", stand = F))
sil_m <- round(mean(sil),2)
sil_m

###############################################################################
########################################## Classif. Acc. by HCA on PCA
###############################################################################

library(FactoMineR)

# feature data
ds_bd <- subset(ds, n_gr_t == "QC")

# perform PCA
pca.ds <- PCA(ds_bd[,-c(1:7)], scale.unit = T, graph = F)
df.pca.x <- as.data.frame(pca.ds$ind$coord[,1:2]) # or pca.ds$X[,1:n], where n - number of components

# generate group data by batch
batch <- ds_bd$b_id
batch <- str_remove(batch, "b")
pca.class <- as.numeric(batch)

# n groups
k <- length(unique(pca.class))

# perform HCA
res.dist <- dist(df.pca.x, method = "euclidean")
res.hc <- hclust(d = res.dist, method = "ward.D2")
pred_cl <- cutree(res.hc, k=k)

# perform
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

cl_ac_hca <- misclass(pred_cl, pca.class)
m_ca_h <- round(mean(cl_ac_hca),0)
m_ca_h

###############################################################################
########################################## PVCA
###############################################################################

library(proBatch)

# generate batch data
batch <- ds$b_id
batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# generate group data
class <- as.character(ds$n_gr_t)
qc_id <- which(class == "QC")
class[-qc_id] <- "Subject"

# Sample data
s_d <- data.frame(FullRunName = rownames(ds), batch, class)
rownames(s_d) <- NULL

# Feature data
f_d <- as.matrix(t(ds[,-c(1:7)]))
colnames(f_d) <- rownames(ds)

# Perform
pvca_df <- calculate_PVCA(f_d, s_d, 
                          factors_for_PVCA = c('batch', 'class'), # adjust to your data
                          pca_threshold = .6, variance_threshold = .01, fill_the_missing = 0)

pvca_w <- round(pvca_df$weights[2],2)
pvca_w

###############################################################################
########################################## gPCA
###############################################################################

library(gPCA)

# generate batch data
batch <- ds$b_id
batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# perform
gPCA_out<-gPCA.batchdetect(x=ds[,-c(1:7)],batch=batch,center=F,nperm=1000)
gPCA_d <- round(gPCA_out$delta, 2)
gPCA_d

###############################################################################
########################################## PC-PR2
###############################################################################

library(pcpr2)

# generate batch data
batch <- ds$b_id
batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# generate group data
class <- as.character(ds$n_gr_t)
qc_id <- which(class == "QC")
class[-qc_id] <- "Subject"

# Sample data
s_d <- data.frame(batch, class) # adjust to your data

# Feature data
f_d <- as.matrix(ds[,-c(1:7)])

# Perform
pct.threshold <- 0.8
PCPR2 <- runPCPR2(X = f_d, Z = s_d, pct.threshold = pct.threshold)

pcpr2 <- round(PCPR2$pR2[1],2)
pcpr2

###############################################################################
########################################## ANOVA
###############################################################################

library(MetabolomicsBasics)

# generate batch data
batch <- as.character(ds$b_id)

# generate run order data
order <- as.numeric(ds$ro_id)

# genetate class data
class <- as.character(ds$n_gr_t)

# Sample data
s_d <- data.frame(cbind(ID = rownames(ds), batch, class, order))

# feature data
dat <- as.matrix(ds[,-c(1:7)])

# perform
model <- MetaboliteANOVA(dat=dat, sam=s_d, model="batch") # "batch+class" or "batch+class+order" or "batch" # adjust to your data 
anova_per <- round(sum(model[,"batch"]<0.05)/nrow(model)*100, 0) # or round(sum(model[,"batch"]<0.05)/nrow(model)*100, 0) # round(sum(model[,"batch"]<0.05 & model[,"class"]>0.05)/nrow(model)*100, 0) or round(length(which(model[,"batch"]<0.05 & model[,"class"]>0.05))/nrow(model)*100, 0)
anova_per

###############################################################################
########################################## One-Sample Test
###############################################################################

# feature data
ds_bd <- subset(ds, n_gr_t == "QC")               

# perform
ost <- function(x, p.val.sig = 0.05, p.adjust = "BH"){
  
  norm_tests <- function(x) {
    xx <- x
    # normality test
    norm.test <- as.numeric(apply(xx, 2, function(t) shapiro.test(t)$p.value))
    # if some error in Shapiro normality test:
    # use shapiro.wilk.test function from cwhmisc instead shapiro.test from stats
    # library(cwhmisc)
    # norm.test <- apply(xx, 2, function(t) cwhmisc::shapiro.wilk.test(t)$p)
  }
  
  res_tests <- norm_tests(x)
  
  wilcox_test <- function(x,y) {
    xx <- x
    wx.t <- as.vector(which(y < 0.05))
    wilcox_test <- list()
    ifelse(identical(wx.t, integer(0)), return (wilcox_test <- 1), wx.t)
    wilcox_test <- apply(as.data.frame(xx[,wx.t]), 2, function(t) as.numeric(as.vector(p.adjust(wilcox.test(t, mu = mean(t, na.rm = T), alternative = "two.sided")$p.value, p.adjust))))
    names(wilcox_test) <- (colnames(x))[wx.t]
    return(wilcox_test)}
  
  wx.t.res <- wilcox_test(x, res_tests)
  
  student_test <- function(x,y) {
    xx <- x
    st.t <- as.vector(which(y > 0.05))
    student_test <- list()
    ifelse(identical(st.t, integer(0)), return (student_test <- 1), st.t)
    student_test <- apply(as.data.frame(xx[,st.t]), 2, function(t) as.numeric(as.vector(p.adjust(t.test(t, mu = mean(t, na.rm = T), alternative = "two.sided")$p.value, p.adjust))))
    names(student_test) <- (colnames(x))[st.t]
    return(student_test)}
  
  st.t.res <- student_test(x, res_tests)
  
  filt_p_val <- function(x, y, w){
    
    #x = ds
    #y = wx.t.res
    #w = st.t.res
    
    wx.t.n <- names(y)
    wx.t.res2 <-as.data.frame(y)
    rownames(wx.t.res2) <- wx.t.n
    
    st.t.n <- names(w)
    st.t.res2 <- as.data.frame(w)
    rownames(st.t.res2) <- st.t.n
    
    wxx <- rownames(wx.t.res2)[which(wx.t.res2 <= p.val.sig)]
    stt <- rownames(st.t.res2)[which(st.t.res2 <= p.val.sig)]
    aff <- c(wxx, stt)
    
    ds_fil <- x[, aff]
    return(ds_fil)
  }
  return(filt_p_val(x, wx.t.res, st.t.res))
}

# 1 st argument -> dataset with only numeric values, 2nd -> p-value, 3rd -> the method of adjustment for multiple comparisons.
ds_ost <- as.data.frame(ost(ds_bd[,-c(1:7)], p.val.sig = 0.05, p.adjust = "BH"))
os_test <- round(ncol(ds_ost)/ncol(ds[,-c(1:7)])*100, 0)
os_test

###############################################################################
########################################## all numeric metrics
###############################################################################

# RSD
result_RSD[2,ncol(result_RSD)]
# correlation
m_c_c
# Bhat Dist
m_b_d
# Diss Score
m_d_s
# Silhouette
sil_m
# HCA on PCA classif acc
m_ca_h
# PVCA
pvca_w
# gPCA
gPCA_d
# PC-PR2
pcpr2
# ANOVA
anova_per
# One-Sample Test
os_test

###############################################################################
########################################## Data Projection (PCA, HCA, scatter plot, box-plot)
###############################################################################

# after processing from "Generate data for evaluation"

# Settings for projection
library(factoextra)
library(FactoMineR)
library(rafalib)
library(RSEIS)
library(ggsci)
library(ggplot2)

############################################### PRINCIPAL COMPONENT ANALYSIS

# dataset
base1 <- ds
mtrx1 <- ds[,-c(1:7)] 
grp1 <- base1[,3] # 3 for batch; 1 for label

palette_pca <- "category20" # color: "lancet" or "category20" or pal_gsea("default", n = length(unique(grp1)), alpha = 0.6, reverse = T)(length(unique(grp1)))
# palette_pca <- JGRAY(length(unique(grp1))) # grey

pca.ds1 <- PCA(mtrx1, scale.unit = T, graph = F)
pca <- fviz_pca_ind(pca.ds1,
                  title = "",
                  geom.ind = "point", # show points only 
                  col.ind = grp1, # color by groups
                  palette = palette_pca, # color "jco" gray JGRAY(length(unique(grp1)))
                  addEllipses = F, # Concentration ellipses
                  legend.title = "Groups")

pca  +scale_shape_manual(values=rep(0:length(unique(grp1)))) # number of groups

################### other type of PCA

# generate batch data
batch <- ds$b_id
batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# generate group data
class <- as.numeric(as.integer(as.factor(ds$n_gr_t)))

# dataset
base1 <- ds
mtrx1 <- ds[,-c(1:7)] 
grp1 <- batch # batch or class

pca.ds <- PCA(mtrx1, scale.unit = T, graph = F)
pca <-fviz_pca_ind(pca.ds, col.ind=grp1, geom = "point", 
             gradient.cols = c("yellow", "green", "darkblue" ), legend.title = "Time", title = "")

pca+scale_color_gradient2(low="yellow", mid="green", high="darkblue", midpoint = c(as.numeric(quantile(grp1)))[3], breaks = c(as.numeric(quantile(grp1))))

############################################### HIERARCHICAL CLUSTER ANALYSIS

# number of groups
k <- length(unique(grp1)) # groups in HC

# color
Cols = function(vec, ord){
  cols = pal_lancet(palette = c("lanonc"), alpha = 1)(length(unique(vec)))
  return(cols[as.fumeric(vec)[ord]])}

# grey
#Cols = function(vec, ord){
# cols = JGRAY(length(unique(vec)))
# return(cols[as.fumeric(vec)[ord]])}

# dataset
base1 <- ds
mtrx1 <- ds[,-c(1:7)] 
grp1 <- base1[,3] # 3 for batch 1 for label

mtrx1_1 <- mtrx1
#mtrx1_1 <- data.frame(scale(mtrx1, center = T, scale = T))
rownames(mtrx1_1) = make.names(grp1, unique=TRUE)
res.dist1 <- dist(mtrx1_1, method = "manhattan") #{euclidean}, {maximum}, {manhattan}, {canberra}, {binary}, {minkowski}
res.hc1 <- hclust(d = res.dist1, method = "ward.D2") #{ward (ward.D)}, {single}, {complete}, {average}, {mcquitty},{median}, {centroid}
hca <- fviz_dend(res.hc1, k = k, # Cut in k groups
               cex = 0.3, # label size 0.3/0.7
               k_colors = unique(Cols(grp1,res.hc1$order)), # "lancet" color "jco" gray JGRAY(k_hc)
               color_labels_by_k = F, # color labels by groups
               label_cols = Cols(grp1,res.hc1$order),#Cols(ds_rfe[,1])[res.hc1$order], #as.fumeric(ds_rfe[,1])[res.hc1$order],
               rect = T, # Add rectangle around groups
               rect_fill = T,
               rect_border = unique(Cols(grp1,res.hc1$order)), #"lancet"# color "jco" gray JGRAY(k_hc)
               horiz = F,
               lwd = 0.3, # lines size 0.3/0.7
               show_labels = T,
               main = "",
               ylab = "") 

hca +scale_shape_manual(values=rep(0:length(unique(grp1)))) # number of groups

############################################### SCATTER PLOT BY BATCH NUMBER

dat <- ds[,-c(1:7)]
ds_s <- sapply(1:nrow(dat), function(x) sum(dat[x,], na.rm = T))
ds$b_id <- as.numeric(stringr::str_remove(ds$b_id, "b"))
ds$n_gr_t[-which(stringr::str_detect(ds$n_gr_t, "QC"))] <- "Sample"
ds_plot <- as.data.frame(cbind(ds[,1:7], ds_s))
# ds_plot <- subset(ds_plot, n_gr_t == "Sample") # subset type of samples
sp <- ggplot(ds_plot, aes(x=ro_id, y=ds_s, color=b_id, shape = n_gr_t)) + # or x=as.numeric(ro_id)
  geom_point(size = 2) + theme_bw() + ggtitle("") +
  xlab("Run order") + ylab("Total intensity") + labs(color = "Batch", shape = "Label") + theme(legend.position="top") +
  scale_color_gradient(low="blue", high="red", breaks = c(as.numeric(quantile(ds$b_id))))  +  geom_smooth(method=lm) # lm or loess

sp

######################################## PLOT BOXPLOT

# feature data
ds_bd <- subset(ds, n_gr_t == "QC")

# generate group data by batch
batch <- ds_bd$b_id
batch <- str_remove(batch, "b")
batch <- as.factor(as.numeric(batch))

# clean and reshape data
dat <- ds_bd[,-c(1:7)]
ds_s <- sapply(1:nrow(dat), function(x) sum(dat[x,], na.rm = T))
ds_bd <- as.data.frame(cbind(batch = batch, ds_s))
ds_bd$batch <- as.factor(ds_bd$batch)
colnames(ds_bd) <- c("batch", "value")

# boxplot by group
bp <- ggplot(data = ds_bd, aes(x=batch, y=value)) + xlab("") + ylab("") +
  geom_boxplot(aes(fill=batch)) + theme(legend.position="bottom") + theme_classic() 

bp

# boxplot by sample
dat2 <- stack(as.data.frame(t(dat)))
batch_long <- rep(batch, each = nrow(dat2)/nrow(dat)) 
dat2$Label <- batch_long
bp <- ggplot(dat2) + 
  geom_boxplot(aes(x = ind, y = values, fill = Label), outlier.size = 0.3) + theme_classic() + xlab("Sample") + ylab("Intensity")+theme(legend.position="top")+theme(axis.text.x = element_blank())

bp

##############################################################################################################################################################
# References
##############################################################################################################################################################

# 1. Deng, Kui, et al. "WaveICA: A novel algorithm to remove batch effects for large-scale untargeted metabolomics data based on wavelet analysis." Analytica chimica acta 1061 (2019): 60-69.
# 2. Deng, Kui, et al. "WaveICA 2.0: a novel batch effect removal method for untargeted metabolomics data without using batch information." Metabolomics 17.10 (2021): 1-8.
# 3. Karpievitch, Yuliya V., et al. "Metabolomics data normalization with EigenMS." PloS one 9.12 (2014): e116221.
# 4. Bararpour, Nasim, et al. "DBnorm as an R package for the comparison and selection of appropriate statistical methods for batch effect correction in metabolomic studies." Scientific reports 11.1 (2021): 1-13.
# 5. Risso, Davide, et al. "Normalization of RNA-seq data using factor analysis of control genes or samples." Nature biotechnology 32.9 (2014): 896-902.
# 6. Marr, Sue, et al. "LC-MS based plant metabolic profiles of thirteen grassland species grown in diverse neighbourhoods." Scientific data 8.1 (2021): 1-12.
# 7. Kirwan, J. A., et al. "Characterising and correcting batch variation in an automated direct infusion mass spectrometry (DIMS) metabolomics workflow." Analytical and bioanalytical chemistry 405.15 (2013): 5147-5157.
# 8. Klavus, Anton, et al. ""Notame": Workflow for Non-Targeted LC-MS Metabolic Profiling." Metabolites 10.4 (2020): 135.
# 9. Luan, Hemi, et al. "statTarget: A streamlined tool for signal drift correction and interpretations of quantitative mass spectrometry-based omics data." Analytica chimica acta 1036 (2018): 66-72.
# 10. Shen, Xiaotao, et al. "Normalization and integration of large-scale metabolomics data using support vector regression." Metabolomics 12.5 (2016): 89.
# 11. Wehrens, Ron, et al. "Improved batch correction in untargeted MS-based metabolomics." Metabolomics 12.5 (2016): 88.
# 12. Brunius, Carl, Lin Shi, and Rikard Landberg. "Large-scale untargeted LC-MS metabolomics data correction using between-batch feature alignment and cluster-based within-batch signal intensity drift correction." Metabolomics 12.11 (2016): 173.
# 13. Broadhurst, David, et al. "Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies." Metabolomics 14.6 (2018): 1-17.
# 14. Livera, Alysha M. De, et al. "Statistical methods for handling unwanted variation in metabolomics data." Analytical chemistry 87.7 (2015): 3606-3615.
# 15. Drotleff, Bernhard, and Michael Lammerhofer. "Guidelines for selection of internal standard-based normalization strategies in untargeted lipidomic profiling by LC-HR-MS/MS." Analytical chemistry 91.15 (2019): 9836-9843.
# 16. Sanchez-Illana, Angel, et al. "Evaluation of batch effect elimination using quality control replicates in LC-MS metabolite profiling." Analytica chimica acta 1019 (2018): 38-48.
# 17. Caesar, Lindsay K., Olav M. Kvalheim, and Nadja B. Cech. "Hierarchical cluster analysis of technical replicates to identify interferents in untargeted mass spectrometry metabolomics." Analytica chimica acta 1021 (2018): 69-77.
# 17. Fages, Anne, et al. "Investigating sources of variability in metabolomic data in the EPIC study: the Principal Component Partial R-square (PC-PR2) method." Metabolomics 10.6 (2014): 1074-1083.
