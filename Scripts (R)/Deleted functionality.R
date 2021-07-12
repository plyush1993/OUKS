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

