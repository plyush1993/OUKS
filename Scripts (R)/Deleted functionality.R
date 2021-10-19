################################################### Multigroup Fold Change

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
########################################## Median Between Batches (QC-norm)
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
fwrite(ds_qc_norm, "xcms after IPO MVI QC-RF + QC-NORM.csv", row.names = T)
