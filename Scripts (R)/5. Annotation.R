##############################################################################################################################################################
# Table of contents
##############################################################################################################################################################

# Installation
# CAMERA
# RAMClustR
# xMSannotator
# References

##############################################################################################################################################################
# Installation
##############################################################################################################################################################

# setup environment
library(data.table)
library(stringr)
setwd("D:/...")

##############################################################################################################################################################
# CAMERA
##############################################################################################################################################################

library(CAMERA)
library(xcms)

# if data with MS2 scans:
# get only mslevel = 1 from your data
# data <- pk_fil@featureData@data
# cor_data <- data[which(data$msLevel == 1),]
# pk_fil@featureData@data <- cor_data
# pk_fil <- as(pk_fil, "xcmsSet")

# load optimal params from IPO
load("IPO_optimiz_xcms_CDF_7QC.RData")
ppm_diff <- resultPeakpicking[["best_settings"]][["parameters"]][["ppm"]]
mz_diff <- 0.02  # resultPeakpicking[["best_settings"]][["parameters"]][["mzdiff"]]
rt_diff <- 15 # adjust to your data

# load xcms object after peak filling
load("xcms obj pk_fil.RData") # filename and the same folder of raw data as in section 2
pk_fil <- as(pk_fil, "xcmsSet") # convert XCMSnExp to xcmsSet object

# perform
xsa <- xsAnnotate(pk_fil) # Create an xsAnnotate object 
xsaF <- groupFWHM(xsa, perfwhm=0.6) # Group after RT value of the xcms grouped peak
xsaC <- groupCorr(xsaF) # Verify grouping

# Annotate isotopes
xsaFI <- findIsotopes(xsaC, ppm = ppm_diff, mzabs=mz_diff, minfrac = 0.2) # set params according to your data

# Annotate adducts
xsaFA <- findAdducts(xsaFI, polarity= "positive", ppm=ppm_diff, mzabs=mz_diff) # set polarity and params according to your data

# Get annotation info
annot_camera <- getPeaklist(xsaFA)
vec_an <- paste(annot_camera$adduct, annot_camera$isotopes)
count_camera <- round(length(which(str_length(vec_an) != 1))/nrow(annot_camera)*100, 0) # percent of annotated data
count_camera

# Get annotation table
annot_camera$name <- paste(annot_camera$mz, annot_camera$rt, sep = " / ")
annot_camera$index <- c(1:nrow(annot_camera))
annot_camera <- data.frame(annot_camera$index, annot_camera$name, annot_camera$isotopes, annot_camera$adduct, annot_camera$mz, annot_camera$rt)

# save
fwrite(annot_camera, "xcms CAMERA.csv", row.names = T)

##############################################################################################################################################################
# RAMClustR
##############################################################################################################################################################

library(xcms)
library(RAMClustR)
library(BiocParallel)
library(tidyverse)
library(dplyr)

# load optimal params from IPO
load("IPO_optimiz_xcms_CDF_7QC.RData")
ppm_diff <- resultPeakpicking[["best_settings"]][["parameters"]][["ppm"]]
mz_diff <- 0.02  # resultPeakpicking[["best_settings"]][["parameters"]][["mzdiff"]]
rt_diff <- 15 # adjust to your data

# load xcms object after peak filling
load("xcms obj pk_fil.RData") # filename and the same folder of raw data as in section 2

# if data with MS2 scans:
# get only mslevel = 1 from your data
# data <- pk_fil@featureData@data
# cor_data <- data[which(data$msLevel == 1),]
# pk_fil@featureData@data <- cor_data
# pk_fil <- as(pk_fil, "xcmsSet")

# Function for optimization st parameter
getRamSt <- function(XObj) {
  featInfo <- featureDefinitions(XObj)
  hist((featInfo$rtmax-featInfo$rtmin)/2)
  st <- round(median(featInfo$rtmax-featInfo$rtmin)/2, digits = 2)
  abline(v=st)
  return(st)}

# Define optimal value of the st parameter
st_val <- getRamSt(pk_fil)

# Function for optimization sr parameter
plotClust=function(ram,clustnr,xcmsData,samps,dtime=5,dmz=.05) { # adjust to your data
  if(missing(samps)) {
    nSamp=nrow(ram$SpecAbund)
    samps=1:nSamp
  } else nSamp=length(samps)
  whichFeats=which(ram$featclus==clustnr)
  peakMeta=cbind(ram$fmz,ram$frt)
  pkMetaGrp=peakMeta[whichFeats,]
  rtr=ram$clrt[clustnr]+c(-dtime,dtime)
  rtr[rtr<0]=0
  mzr=cbind(ram$fmz[whichFeats]-dmz,ram$fmz[whichFeats]+dmz)
  chr <- chromatogram(xcmsData, mz = mzr, rt = rtr)
  plot(0:1,0:1,type='n',axes=F,xlab='Retention time (s)', ylab='Intensity (AU)',
       main=paste0('RAM cluster ',clustnr,'; RT ',signif(ram$clrt[clustnr],5),'s'))
  box(bty='l')
  for (pk in 1:length(whichFeats)) {
    rts=ints=list()
    for (samp in 1:nSamp) {
      rts[[samp]]=chr[pk,samps[samp]]@rtime
      ints[[samp]]=chr[pk,samps[samp]]@intensity
    }
    nrts=min(sapply(rts,length))
    rts=sapply(rts,function(x) x[1:nrts])
    rts=rowMeans(rts)
    ints=sapply(ints,function(x) x[1:nrts])
    ints=rowMeans(ints,na.rm=T)
    par(new=T)
    plot(rts,ints,type='l',col=pk+1,ylim=c(0,max(ints,na.rm=T)),axes=F,xlab='',ylab='')
  }
  axis(1)
  legend('topright',legend = paste0('F',whichFeats,'@mz',
                                    signif(pkMetaGrp[,1],5)), lty=1,col=(1:length(whichFeats))+1,bty='n')
}

# Define optimal value of the sr parameter
expDes=defineExperiment(force.skip = T)
sr=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9)
st=st_val
maxt=c(10) # adjust retention time similarity for your instrument
par=expand.grid(st=st,sr=sr,maxt=maxt)

nClust=nSing=sizeMax=sizeMed=sizeMean=numeric(nrow(par))
nFeat=list()
samps=sample(1:nrow(pk_fil@phenoData),20)
register(bpstart(SnowParam(4))) # Choose as many cores as you can/want
for (i in 1:nrow(par)) {
  RRP=ramclustR(ms='xcms after IPO MVI.csv', st = par$st[i], sr=par$sr[i], 
                maxt = par$maxt, timepos = 2, sampNameCol = 1, featdelim = ' / ', ExpDes = expDes, mspout = F) 
  nClust[i]=length(RRP$cmpd)
  nSing[i]=RRP$nsing
  sizeMax[i]=max(RRP$nfeat)
  sizeMed[i]=median(RRP$nfeat)
  sizeMean[i]=mean(RRP$nfeat)
  nFeat[[i]]=RRP$nfeat
  # pdf(file=paste0('clusts_par',i,'.pdf'),width=15,height=8)
  par(mfrow=c(4,5),mar=c(4,4,2,0)+.5)
  #clusts=round(c(2:6,seq(7,max(RRP$featclus),length.out = 15)))
 # for (c in clusts) {
 #   plotClust(ram = RRP,clustnr = c,xcmsData = pk_fil,samps = samps)
 # }
 # dev.off()
}

res_sr <- cbind(par,nClust,nSing,sizeMax,sizeMean,sizeMed) # nClust+nSing value should be max
res_sr$nclust_nSing <- res_sr$nClust + res_sr$nSing
res_sr
sr_val <- subset(res_sr, nClust == max(res_sr$nClust))$sr
save(res_sr, file = "optimiz sr RumClustR.RData")

# perform
expDes=defineExperiment(force.skip = T)
RC1 <- ramclustR(ms='xcms after IPO MVI.csv',
                 featdelim = " / ",
                 st = st_val,
                 sr = sr_val,
                 ExpDes=expDes,
                 sampNameCol = 1, mspout = F)

RC1 <- do.findmain(RC1, mode = "positive", mzabs.error = mz_diff, ppm.error = ppm_diff) # adjust for your instrument

# annotation table generating
raw_data <-as.data.frame(fread(input = "xcms after IPO MVI.csv", header=T)) # load data
rownames(raw_data) <- raw_data[,1] # load data
raw_data <- raw_data[,-1] # load data
cn <- colnames(raw_data)
raw_data_mz_rt <- data.frame(str_split_fixed(cn, pattern = " / ", 2))
raw_data_mz_rt$X1 <- round(as.numeric(raw_data_mz_rt$X1), 3) # round mz to 3
raw_data_mz_rt$X2 <- round(as.numeric(raw_data_mz_rt$X2), 1) # round rt to 1
raw_data_mz_rt$X3 <- paste(round(raw_data_mz_rt$X1, 3), round(raw_data_mz_rt$X2, 1), sep = " / ")
colnames(raw_data_mz_rt) <- c("mz", "time", "mz_time_round")

# generate all annotated data info
clustered_mz<- lapply(RC1$M.ann, function (x) round(x$mz,3)) # round mz to 3
clustered_rt<- lapply(1:length(RC1$M.ann), function (x) rep(round(RC1$clrt[x],1), length(clustered_mz[[x]]))) # round rt to 1
clustered_label <- lapply(1:length(RC1$M.ann), function (x) RC1$M.ann[[x]]$adduct)
clustered_name <- lapply(1:length(RC1$M.ann), function(y) paste(clustered_mz[[y]], clustered_rt[[y]], sep = " / "))
annot_rcr <- data.frame(cbind(unlist(clustered_name), unlist(clustered_label)))
annot_rcr_match <- filter(annot_rcr, annot_rcr$X1 %in% raw_data_mz_rt$mz_time)

# Get annotation info
count_rcr <- round(nrow(annot_rcr_match)/length(cn)*100, 0) # percent of annotated data
count_rcr

# Get annotation table
an_t_rcr <- data.frame(cn,cn,cn,cn,cn)
an_t_rcr[,3:5] <- NA
an_t_rcr$cn <- c(1:nrow(an_t_rcr))
cnr <- data.frame(raw_data_mz_rt$mz_time)
an_t_rcr[,3] <- cnr
colnames(an_t_rcr) <- c("index", "raw data name", "raw data name round", "RAMClustR mz_rt", "RAMClustR adduct")
annot_rcr_match <- annot_rcr_match[order(annot_rcr_match$X1, decreasing = F),]
cnr <- data.frame(raw_data_mz_rt$mz_time)
cnr <- cnr[order(cnr$raw_data_mz_rt.mz_time, decreasing = F),]
an_t_rcr <- an_t_rcr[order(an_t_rcr$`raw data name round`, decreasing = F),]
ind_rcr <- which(cnr %in% annot_rcr_match$X1)
an_t_rcr$`RAMClustR mz_rt`[ind_rcr] <- annot_rcr_match$X1
an_t_rcr$`RAMClustR adduct`[ind_rcr] <- annot_rcr_match$X2
an_t_rcr <- an_t_rcr[order(an_t_rcr$index, decreasing = F),]

# save
save(RC1, file = "all RumClust.RData")
fwrite(an_t_rcr, "xcms RAMClust.csv", row.names = T)

##############################################################################################################################################################
# xMSannotator
##############################################################################################################################################################

library(xMSannotator)
library(xcms)
library(stringr)
library(dplyr)

################
# BEFORE USE:
################
# load into environment R script (Ctrl+O): "fix for R 4.0 get_peak_blocks_modulesvhclust.R" (fix fun by yufree repo in GH: https://github.com/yufree/xMSannotator)
# load function (get_peak_blocks_modulesvhclust2) from script "fix for R 4.0 get_peak_blocks_modulesvhclust.R" into environment 
# then run script:
environment(get_peak_blocks_modulesvhclust2) <- asNamespace('xMSannotator')
assignInNamespace("get_peak_blocks_modulesvhclust", get_peak_blocks_modulesvhclust2, ns = "xMSannotator")

# load optimal params from IPO
load("IPO_optimiz_xcms_CDF_7QC.RData")
ppm_diff <- resultPeakpicking[["best_settings"]][["parameters"]][["ppm"]]
mz_diff <- 0.02  # resultPeakpicking[["best_settings"]][["parameters"]][["mzdiff"]]
rt_diff <- 15 # adjust to your data

# load peak table
load("xcms obj pk_fil.RData") # filename and the same folder of raw data as in section 2
# if data with MS2 scans:
# get only mslevel = 1 from your data
# data <- pk_fil@featureData@data
# cor_data <- data[which(data$msLevel == 1),]
# pk_fil@featureData@data <- cor_data
# pk_fil <- as(pk_fil, "xcmsSet")
pk_fil <- as(pk_fil, "xcmsSet") # convert XCMSnExp to xcmsSet object
data <- xcms::groupval(pk_fil, "medret", "into")
mz <- as.numeric(xcms::groups(pk_fil)[, 1])
time <- as.numeric(xcms::groups(pk_fil)[, 4])
data <- as.data.frame(cbind(mz, time, data)) 
data <- unique(data)
data[is.na(data)] <- 0
outloc <- paste0(getwd(), "/result xma/")

########################################## DATABASE SEARCH
# perform
num_nodes <- 3
xma <- HMDB.Annotation(dataA = data, max.mz.diff = ppm_diff, num_nodes = num_nodes, # adjust for your data
                       queryadductlist = c("M+2H", "M+H+NH4", "M+ACN+2H", 
                                           "M+2ACN+2H", "M+H", "M+NH4", "M+Na", "M+ACN+H",
                                           "M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na", "2M+ACN+H"), 
                       mode = "pos", outloc = outloc)

xma <- KEGG.Annotation(dataA = data, max.mz.diff = ppm_diff, num_nodes = num_nodes, # adjust for your data
                       queryadductlist = c("M+2H", "M+H+NH4", "M+ACN+2H", 
                                           "M+2ACN+2H", "M+H", "M+NH4", "M+Na", "M+ACN+H",
                                           "M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na", "2M+ACN+H"), 
                       mode = "pos", outloc = outloc)

# save
save(xma, file = "xma KEGG ipo ppm.RData")
fwrite(xma, "xcms xMSannotator.csv", row.names = T)

########################################## ANNOTATION
# perform
num_nodes <- 3 # adjust to your PC
xma <- multilevelannotation(dataA = data, max.mz.diff = ppm_diff, max.rt.diff = rt_diff, # adjust to your data
         num_nodes = num_nodes, queryadductlist = c("all"), mode = "pos", outloc = outloc, 
         db_name = c("HMDB"), mass_defect_window = mz_diff, mass_defect_mode = "pos")

# save
xma_df <-as.data.frame(fread(input = "result xma/Stage5.csv", header=T)) # load data
fwrite(xma_df, "xcms xMSannotator 5stage HMDB multilevel.csv", row.names = T)

########################################## combine results
# xma
xma_data_5 <-as.data.frame(fread(input = "result xma/Stage5.csv", header=T)) # load data
xma_data_5$mz_time <- paste(round(xma_data_5$mz, 3), round(xma_data_5$time, 1), sep = " / ") # round mz and rt
xma_data_5$mz <- round(xma_data_5$mz, 3) # round mz to 3
xma_data_5$time <- round(xma_data_5$time, 1) # round rt to 1
xma_data_5_un <- distinct(xma_data_5, mz_time, .keep_all = T)

# raw
raw_data <-as.data.frame(fread(input = "xcms after IPO MVI.csv", header=T)) # load data
rownames(raw_data) <- raw_data[,1] # load data
raw_data <- raw_data[,-1] # load data
cn <- colnames(raw_data)
raw_data_mz_rt <- data.frame(str_split_fixed(cn, pattern = " / ", 2))
raw_data_mz_rt$X1 <- round(as.numeric(raw_data_mz_rt$X1), 3) # round mz to 3
raw_data_mz_rt$X2 <- round(as.numeric(raw_data_mz_rt$X2), 1) # round rt to 1
raw_data_mz_rt$X3 <- paste(round(raw_data_mz_rt$X1, 3), round(raw_data_mz_rt$X2, 1), sep = " / ")
colnames(raw_data_mz_rt) <- c("mz", "time", "mz_time")

# match
xma_match <- filter(xma_data_5_un, xma_data_5_un$mz_time %in% raw_data_mz_rt$mz_time)
annot_xma_match <- xma_match[,c("chemical_ID", "mz_time", "Name", "Formula", "MonoisotopicMass", "Adduct")]

# Get annotation info
count_xma <- round(nrow(xma_match)/length(cn)*100, 0) # percent of annotated data
count_xma

# Get annotation table
cnr <- data.frame(raw_data_mz_rt$mz_time)
an_t_xma <- data.frame(rep(data.frame(cn),9))
an_t_xma[,3:9] <- NA
an_t_xma$cn <- c(1:nrow(an_t_xma))
colnames(an_t_xma) <- c("index", "raw data name", "raw data name round", "chemical_ID_xma", "mz_time xma", "Name_xma", "Formula_xma", "MonoisotopicMass_xma", "Adduct_xma")
an_t_xma[,3] <- cnr
annot_xma_match <- annot_xma_match[order(annot_xma_match$mz_time, decreasing = F),]
cnr <- data.frame(raw_data_mz_rt$mz_time)
cnr <- cnr[order(cnr$raw_data_mz_rt.mz_time, decreasing = F),]
an_t_xma <- an_t_xma[order(an_t_xma$`raw data name round`, decreasing = F),]
ind_xma <- which(cnr %in% annot_xma_match$mz_time)
an_t_xma$chemical_ID_xma[ind_xma] <- annot_xma_match$chemical_ID
an_t_xma$`mz_time xma`[ind_xma] <- annot_xma_match$mz_time
an_t_xma$Name_xma[ind_xma] <- annot_xma_match$Name
an_t_xma$Formula_xma[ind_xma] <- annot_xma_match$Formula
an_t_xma$MonoisotopicMass_xma[ind_xma] <- annot_xma_match$MonoisotopicMass
an_t_xma$Adduct_xma[ind_xma] <- annot_xma_match$Adduct
an_t_xma <- an_t_xma[order(an_t_xma$index, decreasing = F),]

# save
fwrite(an_t_xma, "xcms xMSannotator.csv", row.names = T)

##############################################################################################################################################################
# References
##############################################################################################################################################################

# 1. Kuhl, Carsten, et al. "CAMERA: an integrated strategy for compound spectra extraction and annotation of liquid chromatography/mass spectrometry data sets." Analytical chemistry 84.1 (2012): 283-289.
# 2. Broeckling, Corey David, et al. "RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data." Analytical chemistry 86.14 (2014): 6812-6817.
# 3. Fernández-Ochoa, Álvaro, et al. "A Case Report of Switching from Specific Vendor-Based to R-Based Pipelines for Untargeted LC-MS Metabolomics." Metabolites 10.1 (2020): 28.
# 4. Uppal, Karan, Douglas I. Walker, and Dean P. Jones. "xMSannotator: an R package for network-based annotation of high-resolution metabolomics data." Analytical chemistry 89.2 (2017): 1063-1067.