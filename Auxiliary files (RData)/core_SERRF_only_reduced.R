options(warn=-1)

path = paste0(workingdirectory,"\\",filename)
if(!"pacman" %in% rownames(installed.packages())){
  install.packages("pacman")
}
cat("Checking required packages (auto-installing if missing).\n")
pacman::p_load("randomForest", "affy", "e1071", "data.table", "parallel", "xlsx")
setwd(workingdirectory)
# load sources
source("normalizations.R")
source("utils.R")
source("evaluationMethods.R")
# read data.
cat("Reading Data.\n")
# metaData = read.csv("P20 dataset1 2016_11-30.csv")
# p = fread("p-pos.csv")
# p <- merge(p,metaData,by.x="Subject ID", by.y = "GBID", all.x = T,sort = F)
# f = fread("f-pos.csv")
# e = fread("e-pos.csv")
# e = as.matrix(e)
# p$`Acq. Date-Time`  = gsub("/13 ", "/2013 ", p$`Acq. Date-Time`)
# p$`Acq. Date-Time` = as.numeric(strptime(p$`Acq. Date-Time`, "%m/%d/%Y %H:%M"))
#
#
# path = "C:\\Users\\Sili Fan\\Desktop\\WORK\\WCMC\\projects\\mayo_depression_SSRIs_2013_normalization\\mayo_depression_SSRIs_2013_with_time_stamp.xlsx"

# for xlsx input
data = readData(path)
p = data$p
p$`Acq. Date-Time` = p$time
p$`Stat Level 1` = p$type
p$`Stat Level 1`[p$`Stat Level 1`=='validate'] = "NIST"
f = data$f
e = as.matrix(data$e)

if(sum(is.na(e)) > 0){
  cat(paste0("NOTE: ",sum(is.na(e)), " missing values detected in the data. They will be replaced by the half-minimum for each compound."))
  missing_compounds = which(is.na(e), arr.ind = T)[,1]
  for(i in missing_compounds){
    e[i, is.na(e[i,])] = 1/2 * min(e[i,!is.na(e[i,])])
  }
}


e = data.matrix(e)



# batch = rep(diff(sort(c(0,order(diff(p$time), decreasing = T)[c(1:3,5:6)], nrow(p)))), diff(sort(c(0,order(diff(p$time), decreasing = T)[c(1:3,5:6)], nrow(p)))))
# batch = revalue(as.character(batch), c("157"="Batch1", "140" = "Batch2", "167" = "Batch3", "179" = "Batch4", "204" = "Batch5", "178" = "Batch6"))

# png("batch defination.png", width = 800, height = 800)
# plot(p$time, col = factor(batch))
# dev.off()
# define batches according to time interval.
# cat("Visualize the batch vs injection order.\n")
batch = p$batch
batch = matrix(rep(batch,nrow(f)), nrow = nrow(f), byrow = T, ncol = nrow(p))
# dta = data.frame(injection.order = 1:nrow(p), Date = p$`Acq. Date-Time`, batch = batch[1,])
# ggplot(dta, aes(x=injection.order, y=Date, color=batch)) +
#   geom_point()+
#   labs(title="Batch Defination \n according to time intervals",
#        x="Injection Order", y = "Date") + theme.scatter

# set parallel computing if necessary.
cat("Initializing the parallel procedure, using all cores.\n")
cl = makeCluster(detectCores())

# results will be saved in the result_norm.
result_norm = list()

cat("\n<========== Normalizations Started! ==========>\n");

# SERRF
cat("\n<========== SERRF Normalization Started! ==========>\n")
cat("This may take some time. \n")
qc = rep(F, nrow(p))
qc[QC.index] = T
e. = e
result_norm[['SERRF']] = SERRF_norm(e., f, p, batch, QC.index, time = "Acq. Date-Time")

dta = data$original
dta[5:nrow(dta),3:ncol(dta)] = result_norm[['SERRF']]$e
dir.create("normalized-data-sets")
write.csv(dta, file="normalized-data-sets\\normalization-result-SERRF-normalization.csv",row.names=FALSE)
