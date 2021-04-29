##############################################################################################################################################################
# Table of contents
##############################################################################################################################################################

# Installation
# CROP
# Notame
# pmd
# References

##############################################################################################################################################################
# Installation
##############################################################################################################################################################

# setup environment
library(data.table)
library(stringr)
library(batchCorr)
setwd("D:/...")

# PEAK TABLE
# dataset with intensities and label column
ds <- as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats annot+filtr LMM adj KEGG.csv", header=T))
gr_nm <- ds[-c(1:8),c(2,5)]
ds <- ds[-c(1:8),]
rn <- ds[,5]
ds <- ds[,-c(1,3:5)]
ds <- as.data.frame(sapply(ds[,-1], as.numeric))
rownames(ds) <- rn

# extract m/z and rt from peak table
cn <- colnames(ds)
cn0 <- colnames(ds)
cn <- gsub(pattern = " / ", replacement = "_", x = cn, fixed = T)
cn2 <- t(data.frame(cn))
colnames(cn2) <- cn
peakIn <- peakInfo(PT = cn2, sep = "_", start = 1)
peakIn <- as.data.frame(cbind(peakIn, id = cn0))

##############################################################################################################################################################
# CROP
##############################################################################################################################################################

# Prepare
# files from https://github.com/rendju/CROP
load(file = "CROP.RData")
load(file = "mltpl_depict.RData")
load(file = "mltpl_detect.RData")
load(file = "mltpl_remove.RData")

library(gplots)
library(tidyverse)
library(igraph)
library(ggraph)
library(robCompositions)
library(compositions)
library(gridExtra)
library(stringi)
library(corrr)
library(data.table)

# perform
ds_crop <- as.data.frame(cbind(name = colnames(ds), rt = peakIn[,2], t(ds)))
fwrite(ds_crop, "for CROP.csv", row.names = F)
crop <- CROP(mscsv="for CROP.csv", name="CROP", ccth=0.85, rtw=5, maxrtw = 10, mcs=100, rtunit="s", funit="m/z") # adjust to your data
res_crop <- fread("CROP_CROPped_ccth0.85_rtw+-5_data_without_multiplicities.csv", sep = ";") # load results

##############################################################################################################################################################
# Notame
##############################################################################################################################################################

# Prepare
library(notame)
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

# perform
conn <- find_connections(data = ds, features = peakIn,
                         corr_thresh = 0.4, rt_window = 5, # adjust to your data
                         name_col = "id", mz_col = "mz", rt_col = "rt")

clusters <- find_clusters(connections = conn, d_thresh = 0.6)

features_clustered <- assign_cluster_id(ds, clusters, peakIn, name_col = "id")

pulled <- pull_clusters(ds, features_clustered, name_col = "id")
cluster_data <- pulled$cdata
cluster_features <- pulled$cfeatures

# save
fwrite(cluster_data, "after notame clust.csv", row.names = T)

##############################################################################################################################################################
# pmd
##############################################################################################################################################################

# Prepare
library(pmd)
library(igraph)

# Perform
dat_pmd <- list(data = t(ds), group = gr_nm, mz = as.numeric(peakIn$mz), rt = as.numeric(peakIn$rt))
pmd <- getpaired(dat_pmd, rtcutoff = 10, ng = 10) # adjust to your data
plotrtg(pmd)
plotpaired(pmd)

# show the unique pmd found by getpaired function
for(i in 1:length(unique(pmd$paired$diff2))){
  diff <- unique(pmd$paired$diff2)[i]
  index <- pmd$paired$diff2 == diff
  plotpaired(pmd,index)
}
std <- getstd(pmd) # use correlation cut-off -> corcutoff

# extract pseudospectra for std peak N
N <- 71 # specify peak
stdcluster <- getcluster(std)
idx <- unique(stdcluster$cluster$largei[stdcluster$cluster$i==N])
plot(stdcluster$cluster$mz[stdcluster$cluster$largei==idx],stdcluster$cluster$ins[stdcluster$cluster$largei==idx],type = 'h',xlab = 'm/z',ylab = 'intensity',main = 'pseudospectra for GlobalStd peak 71')

# export peaks with the highest intensities in each GlobalStd peaks groups
data <- stdcluster$data[stdcluster$stdmassindex2,] 
data <- as.data.frame(t(data))
fwrite(data, "after pmd clust.csv", row.names = T)

# reactomics
sda <- getsda(std)
plotstdsda(sda)
df <- sda$sda
net <- graph_from_data_frame(df,directed = F)
pal <- (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdYlBu")
)))(length(unique(E(net)$diff2)))
plot(net,vertex.label=round(as.numeric(V(net)$name)),vertex.size = 7,edge.width = 5,edge.color = pal[as.numeric(as.factor(E(net)$diff2))],main = 'PMD network')
legend("topright",bty = "n",
       legend=unique(E(net)$diff2),
       fill=unique(pal[as.numeric(as.factor(E(net)$diff2))]), border=NA,horiz = F)
deg <- degree(net, mode="all")
degree_distribution(net)
plot(net, vertex.size=deg/2,vertex.label=NA,vertex.size = 7, edge.width = 5)
ceb <- cluster_edge_betweenness(net,weights = abs(E(net)$cor), directed = F) 
plot(ceb, net,vertex.label=NA,) 

##############################################################################################################################################################
# References
##############################################################################################################################################################

# 1. Kouril, Stepán, et al. "CROP: Correlation-based reduction of feature multiplicities in untargeted metabolomic data." Bioinformatics 36.9 (2020): 2941-2942.
# 2. Klåvus, Anton, et al. ""Notame": Workflow for Non-Targeted LC-MS Metabolic Profiling." Metabolites 10.4 (2020): 135.
# 3. Yu, Miao, Mariola Olkowicz, and Janusz Pawliszyn. "Structure/reaction directed analysis for LC-MS based untargeted analysis." Analytica chimica acta 1050 (2019): 16-24.
# 4. Yu, Miao, and Lauren Petrick. "Untargeted high-resolution paired mass distance data mining for retrieving general chemical relationships." Communications Chemistry 3.1 (2020): 1-6.