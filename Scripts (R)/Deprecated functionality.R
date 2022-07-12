################################################### Multigroup Fold Change (Stats)

library(structToolbox)

# transform data into log2 base
ds_log <- as.data.frame(log2(ds[,-1]))
ds_log <- cbind(Label = ds[,1], ds_log)

D1 <- DatasetExperiment(data = data.frame(ds_log[,-1]),
                        sample_meta = data.frame(Label = ds$Label),
                        variable_meta = data.frame(colnames(ds[,-1])))
M1 = fold_change(factor_name="Label", alpha = 0.05, threshold = 2)
M2 = model_apply(M1,D1)
fdr <- M2@fold_change@value
detach("package:structToolbox", unload = TRUE)

# by mean
fdr_mean <- apply(abs(fdr),1, mean, na.rm=T)

###############################################################################
########################################## Median Between Batches (QC-norm correction)
###############################################################################

# generate batch data
batch <- ds$b_id
batch <- str_remove(batch, "b")
batch <- as.numeric(batch)

# generate data
ds_bbc <- cbind(batch = ds$b_id, ds[,-c(1:7)])
ds_bbc$batch <- batch

# perform
QC.NORM <- function(data){
library(dplyr)
ds_bbc <- data
b_b_c_subsets <- lapply(1:length(unique(ds_bbc[,1])), function(y) filter(ds_bbc[,-1], ds_bbc$batch == unique(ds_bbc[,1])[y])) # list of subsets by batches
b_b_c_factor <- lapply(1:length(unique(ds_bbc[,1])), function(y) sapply(1:ncol(b_b_c_subsets[[1]]), function(z) median(b_b_c_subsets[[y]][,z], na.rm = T)/median(ds_bbc[,(z+1)], na.rm = T))) # calculate factor for every feature
b_b_c_result_l <- lapply(1:length(unique(ds_bbc[,1])), function(y) sapply(1:ncol(b_b_c_subsets[[1]]), function(z) b_b_c_subsets[[y]][,z]/b_b_c_factor[[y]][z])) # results by batch
b_b_c_result <- data.frame(do.call("rbind",b_b_c_result_l)) # results in data frame
rownames(b_b_c_result) <- rownames(ds_bbc)
colnames(b_b_c_result) <- colnames(ds_bbc[,-c(1)])
return(b_b_c_result)
}

ds_qc_norm <- QC.NORM(data = ds_bbc)

# save
fwrite(ds_qc_norm, "... QC-NORM.csv", row.names = T)

########################################## EACH GROUPS ONE METHODS (NRMSE FOR MVI)

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

############################################### GBMM MODELLING OF BIOLOGICAL FACTORS

# data
dat <- cbind(meta, ds_norm)
n_meta <- 5 # set non numeric columns, adjust to your data
colnames(dat)[-c(1:n_meta)] <- paste0("X", c(1:ncol(dat))) # adjust to your data
dat$Batch <- as.numeric(stringr::str_remove(dat$Batch, "b")) # as.factor(dat$Batch) or as.numeric(stringr::str_remove(dat$Batch, "b")) or as.numeric(dat$Batch)
dat$Age <- as.integer(dat$Age)
dat$Sex <- as.integer(dat$Sex)
dat$Class <- as.numeric(as.factor(dat$Class))
dat[,-c(1:n_meta)] <- sapply(dat[,-c(1:n_meta)], as.numeric)

library(pbapply)
library(gpboost)

# perform
n_start <- 6 # adjust to your data
dss <- lapply(n_start:ncol(dat), function(y) as.data.frame(dat[, c(y, 4, 5, 1, 2)])) # adjust to your data
dss <- lapply(1:length(dss), function(y) {colnames(dss[[y]])[1] <-"Y" 
return(dss[[y]])}) # adjust to your data
gbmm_mod <- pblapply(1:length(dss), function(x) GPModel(group_data = dss[[x]][,5])) # adjust to your data
gbmm_fit <- pblapply(1:length(dss), function(y) gpboost(data = as.matrix(dss[[y]][,c(2,3,4)]), label = as.matrix(dss[[y]][,1]), # adjust to your data
                                                        gp_model = gbmm_mod[[y]],
                                                        nrounds = 16,
                                                        learning_rate = 0.05,
                                                        max_depth = 6,
                                                        min_data_in_leaf = 5,
                                                        objective = "regression_l2",
                                                        verbose = -1)) # use 1 / -1

gbmm_pred <- as.data.frame(pbsapply(1:length(gbmm_fit), function(x) predict(gbmm_fit[[x]], data = as.matrix(dss[[x]][,c(2,3,4)]), group_data_pred = dss[[x]][,5], predict_var= TRUE)))
gbmm_adj <- as.data.frame(pbsapply(1:length(gbmm_fit), function(x) gbmm_pred[[x]]$fixed_effect))
colnames(gbmm_adj) <- colnames(dsr)
rownames(gbmm_adj) <- rownames(dsr)

# save adjusted dataset
ds_gbmm_fit <- data.frame(cbind(Class = meta$Class, data.frame(gbmm_adj)))
fwrite(ds_gbmm_fit, "QC-XGB after dual filt GBMM adj.csv", row.names = T)

############################################### GBM MODELLING OF BIOLOGICAL FACTORS

# data
dat <- cbind(meta, ds_norm)
n_meta <- 5 # set non numeric columns, adjust to your data
colnames(dat)[-c(1:n_meta)] <- paste0("X", c(1:ncol(dat))) # adjust to your data
dat$Batch <- as.numeric(stringr::str_remove(dat$Batch, "b")) # as.factor(dat$Batch) or as.numeric(stringr::str_remove(dat$Batch, "b")) or as.numeric(dat$Batch)
dat$Age <- as.integer(dat$Age)
dat$Sex <- as.integer(dat$Sex)
dat$Class <- as.numeric(as.factor(dat$Class))
dat[,-c(1:n_meta)] <- sapply(dat[,-c(1:n_meta)], as.numeric)

library(pbapply)
library(gpboost)

# perform
n_start <- 6 # adjust to your data
dss <- lapply(n_start:ncol(dat), function(y) as.data.frame(dat[, c(y, 4, 5, 1, 2)])) # adjust to your data
dss <- lapply(1:length(dss), function(y) {colnames(dss[[y]])[1] <-"Y" 
return(dss[[y]])}) # adjust to your data
gbm_fit <- pblapply(1:length(dss), function(y) gpboost(data = as.matrix(dss[[y]][,c(2,3,4)]), label = as.matrix(dss[[y]][,1]), # adjust to your data
                                                        nrounds = 16,
                                                        learning_rate = 0.05,
                                                        max_depth = 6,
                                                        min_data_in_leaf = 5,
                                                        objective = "regression_l2",
                                                        verbose = -1)) # use 1 / -1

gbm_adj <- as.data.frame(pbsapply(1:length(gbm_fit), function(x) predict(gbm_fit[[x]], data = as.matrix(dss[[x]][,c(2,3,4)]))))
colnames(gbm_adj) <- colnames(dsr)
rownames(gbm_adj) <- rownames(dsr)

# save adjusted dataset
ds_gbm_fit <- data.frame(cbind(Class = meta$Class, data.frame(gbm_adj)))
fwrite(ds_gbm_fit, "QC-XGB after dual filt GBM adj.csv", row.names = T)
