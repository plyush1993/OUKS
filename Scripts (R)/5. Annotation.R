##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                              TABLE OF CONTENTS                           ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Installation
# CAMERA
# RAMClustR
# xMSannotator
# mWISE
# metID
# MetaboAnnotation
# References

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Installation                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# setup environment
library(data.table)
library(stringr)
setwd("D:/...")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                   CAMERA                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(xcms)
library(CAMERA)

# if data with MS2 scans:
# get only mslevel = 1 from your data
# data <- pk_fil@featureData@data
# cor_data <- data[which(data$msLevel == 1),]
# pk_fil@featureData@data <- cor_data
# pk_fil <- as(pk_fil, "xcmsSet")

# load optimal params from IPO
load("IPO_optimiz_xcms_CDF_7QC.RData")
ppm_diff <- resultPeakpicking[["best_settings"]][["parameters"]][["ppm"]] # adjust to your data
mz_diff <- 0.02  # resultPeakpicking[["best_settings"]][["parameters"]][["mzdiff"]] # adjust to your data
rt_diff <- 15 # adjust to your data

# load xcms object after peak filling
load("xcms obj pk_fil.RData") # filename and the same folder of raw data as in section "XCMS with the best params (from IPO)" from script "2.Integration"
pk_fil <- as(pk_fil, "xcmsSet") # convert XCMSnExp to xcmsSet object

# perform
xsa <- xsAnnotate(pk_fil) # Create an xsAnnotate object 
xsaF <- groupFWHM(xsa, perfwhm=0.6) # Group after RT value of the xcms grouped peak, adjust to your data
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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  RAMClustR                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(xcms)
library(RAMClustR)
library(BiocParallel)
library(tidyverse)
library(dplyr)

# load optimal params from IPO
load("IPO_optimiz_xcms_CDF_7QC.RData")
ppm_diff <- resultPeakpicking[["best_settings"]][["parameters"]][["ppm"]] # adjust to your data
mz_diff <- 0.02  # resultPeakpicking[["best_settings"]][["parameters"]][["mzdiff"]] # adjust to your data
rt_diff <- 15 # adjust to your data

# load xcms object after peak filling
load("xcms obj pk_fil.RData") # filename and the same folder of raw data as in section "XCMS with the best params (from IPO)" from script "2.Integration"

# if data with MS2 scans:
# get only mslevel = 1 from your data
# data <- pk_fil@featureData@data
# cor_data <- data[which(data$msLevel == 1),]
# pk_fil@featureData@data <- cor_data
# pk_fil <- as(pk_fil, "xcmsSet")

# Function for optimization st parameter
getRamSt <- function(XObj) ({
  featInfo <- featureDefinitions(XObj)
  hist((featInfo$rtmax-featInfo$rtmin)/2)
  st <- round(median(featInfo$rtmax-featInfo$rtmin)/2, digits = 2)
  abline(v=st)
  return(st)})

# Define optimal value of the st parameter
st_val <- getRamSt(pk_fil)

# Function for optimization sr parameter
plotClust=function(ram,clustnr,xcmsData,samps,dtime=5,dmz=.05) ({ # adjust to your data
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
})

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
save(res_sr, file = "optimiz sr RamClustR.RData")

# perform
expDes=defineExperiment(force.skip = T)
RC1 <- ramclustR(ms='xcms after IPO MVI.csv', # you can use also ms2 data with argument "idmsms"
                 featdelim = " / ",
                 st = st_val, # you could also set it manually
                 sr = sr_val, # you could also set it manually
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
annot_rcr_match <- dplyr::filter(annot_rcr, annot_rcr$X1 %in% raw_data_mz_rt$mz_time)

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
save(RC1, file = "all RamClust.RData")
fwrite(an_t_rcr, "xcms RAMClust.csv", row.names = T)

# results
# see "spectra" folder for other info

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                xMSannotator                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(xcms)
library(xMSannotator)
library(stringr)
library(dplyr)
library(data.table)

# load optimal params from IPO
load("IPO_optimiz_xcms_CDF_7QC.RData")
ppm_diff <- resultPeakpicking[["best_settings"]][["parameters"]][["ppm"]] # adjust to your data
mz_diff <- 0.02  # resultPeakpicking[["best_settings"]][["parameters"]][["mzdiff"]] # adjust to your data
rt_diff <- 15 # adjust to your data

# load peak table
load("xcms obj pk_fil.RData") # filename and the same folder of raw data as in section "XCMS with the best params (from IPO)" from script "2.Integration"
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

#........................DATABASE SEARCH.........................
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
setwd("D:/...")
save(xma, file = "xma DB search.RData")
fwrite(xma, "xcms xMSannotator.csv", row.names = T)

#...........................ANNOTATION...........................
# perform
num_nodes <- 3 # adjust to your PC
xma <- multilevelannotation(dataA = data, max.mz.diff = ppm_diff, max.rt.diff = rt_diff, # adjust to your data
         num_nodes = num_nodes, queryadductlist = c("all"), mode = "pos", outloc = outloc, 
         db_name = c("HMDB"), mass_defect_window = mz_diff, mass_defect_mode = "pos")

# save
setwd("D:/...")
xma_df <-as.data.frame(fread(input = "result xma/Stage5.csv", header=T)) # load data "result xma/Stage5.csv" or "Stage5.csv"
fwrite(xma_df, "xcms xMSannotator 5stage HMDB multilevel.csv", row.names = T)

#........................combine results.........................
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
xma_match <- dplyr::filter(xma_data_5_un, xma_data_5_un$mz_time %in% raw_data_mz_rt$mz_time)
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

# results
# see "result xma" folder for other info

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    mWISE                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(mWISE)
library(stringr)
library(dplyr)
library(data.table)
library(parallel)
library(doParallel)

# stop parallel
stopCluster(cl)
stopImplicitCluster()

# start parallel processing
fc <- as.numeric(detectCores(logical = T))
cl <- makePSOCKcluster(fc-1)
registerDoParallel(cl)

# load optimal params from IPO
load("IPO_optimiz_xcms_CDF_7QC.RData")
ppm_diff <- resultPeakpicking[["best_settings"]][["parameters"]][["ppm"]]

# load peak table
raw_data <-as.data.frame(fread(input = "xcms after IPO MVI.csv", header=T)) # load data
rownames(raw_data) <- raw_data[,1] # load data
raw_data <- raw_data[,-1] # load data
cn <- colnames(raw_data)
raw_data_mz_rt <- data.frame(str_split_fixed(cn, pattern = " / ", 2))
raw_data_mz_rt$X1 <- round(as.numeric(raw_data_mz_rt$X1), 4) # round mz to 4 adjust to your data
raw_data_mz_rt$X2 <- round(as.numeric(raw_data_mz_rt$X2), 1) # round rt to 1 adjust to your data
raw_data_mz_rt$X3 <- paste(round(raw_data_mz_rt$X1, 3), round(raw_data_mz_rt$X2, 1), sep = " / ")
colnames(raw_data_mz_rt) <- c("mz", "rt", "mz_time")

# prepare matrix
data_m <- as.data.frame(t(raw_data))
data_m <- as.data.frame(cbind(rownames(data_m),raw_data_mz_rt[,1:2], data_m)) # adjust to your data
colnames(data_m)[1] <- "Peak.Id"

# perform
data("KeggDB") # load data
Cpd.Add <- CpdaddPreparation(KeggDB = KeggDB, do.Par = T, nClust = cl) # load data, adjust to your data
Intensity.idx <- c(4:ncol(data_m)) # adjust to your data
Annotated.List <- mWISE.annotation(Peak.List = data_m,
                                   polarity = "positive", # adjust to your data
                                   diffusion.input.type = "binary",
                                   score = "raw",
                                   ppm = ppm_diff, # adjust to your data
                                   Cpd.Add = Cpd.Add,
                                   Add.List = c("M+2H", "M+H+NH4", "M+ACN+2H", 
                                                "M+2ACN+2H", "M+H", "M+NH4", "M+Na", "M+ACN+H",
                                                "M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na", "2M+ACN+H") , # adjust to your data
                                   Unique.Annotation = TRUE,
                                   Intensity.idx = Intensity.idx,
                                   do.Par = T, nClust = cl) # adjust to your data

# save
nd <- dir.create("result mWise") 
save(Annotated.List, file = "result mWise/mWISE results.RData")
fwrite(Annotated.List$Annotated.Tab, "result mWise/xcms mWISE Annotated.Tab.csv", row.names = T)
fwrite(Annotated.List$MH.Tab, "result mWise/xcms mWISE MH.Tab.csv", row.names = T)
fwrite(Annotated.List$Diff.Tab, "result mWise/xcms mWISE Diff.Tab.csv", row.names = T)
fwrite(Annotated.List$Ranked.Tab, "result mWise/xcms mWISE Ranked.Tab.csv", row.names = T)

#........................combine results.........................
# mWISE
mwise_data <-as.data.frame(fread(input = "result mWise/xcms mWISE Annotated.Tab.csv", header=T)) # load data
mwise_data$mz <- round(as.numeric(mwise_data$mz), 4) # round mz to 4 adjust to your data
mwise_data$rt <- round(as.numeric(mwise_data$rt), 1) # round rt to 1 adjust to your data
mwise_data$mz_time <- paste(round(mwise_data$mz, 4), round(mwise_data$rt, 1), sep = " / ") # round mz and rt adjust to your data
mwise_data_un <- distinct(mwise_data, Peak.Id, .keep_all = T)

# raw
raw_data <-as.data.frame(fread(input = "xcms after IPO MVI.csv", header=T)) # load data
rownames(raw_data) <- raw_data[,1] # load data
raw_data <- raw_data[,-1] # load data
cn <- colnames(raw_data)
raw_data_mz_rt <- data.frame(str_split_fixed(cn, pattern = " / ", 2))
raw_data_mz_rt$X1 <- round(as.numeric(raw_data_mz_rt$X1), 4) # round mz to 4 adjust to your data
raw_data_mz_rt$X2 <- round(as.numeric(raw_data_mz_rt$X2), 1) # round rt to 1 adjust to your data
raw_data_mz_rt$X3 <- paste(round(raw_data_mz_rt$X1, 4), round(raw_data_mz_rt$X2, 1), sep = " / ") # adjust to your data
colnames(raw_data_mz_rt) <- c("mz", "time", "mz_time")

# match
mwise_match <- dplyr::filter(mwise_data_un, mwise_data_un$mz_time %in% raw_data_mz_rt$mz_time)
annot_mwise_match <- mwise_match[,c("mz_time", "Peak.Id", "mz", "rt", "Compound", "exact_mass", "Add.name")]

# Get annotation info
count_mwise <- round(nrow(mwise_match)/length(cn)*100, 0) # percent of annotated data
count_mwise

# Get annotation table
cnr <- data.frame(raw_data_mz_rt$mz_time)
an_t_mwise <- data.frame(rep(data.frame(cn),9))
an_t_mwise[,3:9] <- NA
an_t_mwise$cn <- c(1:nrow(an_t_mwise))
colnames(an_t_mwise) <- c("index", "raw name", "mz_time", "Peak.Id", "mz", "rt", "Compound", "exact_mass", "Add.name")
an_t_mwise[,3] <- cnr
annot_mwise_match <- annot_mwise_match[order(annot_mwise_match$mz_time, decreasing = F),]
cnr <- data.frame(raw_data_mz_rt$mz_time)
cnr <- cnr[order(cnr$raw_data_mz_rt.mz_time, decreasing = F),]
an_t_mwise <- an_t_mwise[order(an_t_mwise$mz_time, decreasing = F),]
ind_mwise <- which(cnr %in% annot_mwise_match$mz_time)
an_t_mwise$Peak.Id[ind_mwise] <- annot_mwise_match$Peak.Id
an_t_mwise$mz[ind_mwise] <- annot_mwise_match$mz
an_t_mwise$rt[ind_mwise] <- annot_mwise_match$rt
an_t_mwise$Compound[ind_mwise] <- annot_mwise_match$Compound
an_t_mwise$exact_mass[ind_mwise] <- annot_mwise_match$exact_mass
an_t_mwise$Add.name[ind_mwise] <- annot_mwise_match$Add.name
an_t_mwise <- an_t_mwise[order(an_t_mwise$index, decreasing = F),]

# save
fwrite(an_t_mwise, "xcms mWISE.csv", row.names = T)

# results
# see "result mWise" folder for other info

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    metID                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(metid)
library(tidyverse)
library(doParallel)
library(data.table)
library(stringr)
library(batchCorr)

# ppm error from IPO
load("IPO_optimiz_xcms_CDF_7QC.RData")
ppm <- resultPeakpicking[["best_settings"]][["parameters"]][["ppm"]]

# PEAK TABLE
# dataset with intensities and filenames column
ds <- as.data.frame(fread(input = "xcms after IPO MVI QC-XGB.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]

# extract m/z and rt from peak table
cn <- colnames(ds)
cn <- str_remove(cn, "X")
cn0 <- colnames(ds)
cn0 <- str_remove(cn0, "X")
cn <- gsub(pattern = "...", replacement = "_", x = cn, fixed = T)
cn2 <- t(data.frame(cn))
colnames(cn2) <- cn
peakIn <- peakInfo(PT = cn2, sep = "_", start = 1)
peakIn <- as.data.frame(cbind(name = cn0, peakIn))

# create a folder
path <- file.path(".", "for metID")
dir.create(path = path, showWarnings = FALSE)

# get MS1 peak table from metID
fwrite(peakIn, "for metID/for metID MS1.csv", row.names = F)

# get database from metID
database <- system.file("ms2_database", package = "metid")

#............................perform.............................
# identification by different databases
# copy databases (https://github.com/jaspershen/demoData/tree/master/inst/ms2_database) in "ms2_database" folder in metID library folder 

#..................identification by msDatabase..................
file.copy(from = file.path(database, "msDatabase_rplc0.0.2"), 
          to = path, overwrite = TRUE, recursive = TRUE)

cores <- detectCores()-1 # adjust to your data
annotate_result1 <-  identify_metabolites(ms1.data = "for metID MS1.csv",  
                     ms1.match.ppm = ppm, # adjust to your data
                     rt.match.tol = 1000000, # adjust to your data
                     polarity = "positive", # adjust to your data
                     column = "rp", # adjust to your data
                     path = path, 
                     candidate.num = 3, # adjust to your data
                     database = "msDatabase_rplc0.0.2", 
                     threads = cores)

# get results
table <- get_identification_table(annotate_result1, candidate.num = 3, type = "new") # adjust to your data

# save
fwrite(table, "for metID/metID MS1 msDatabase.csv", row.names = T)
# see "for metID" folder

#.....................identification by HMDB.....................
file.copy(from = file.path(database, "hmdbMS1Database0.0.1"), 
          to = path, overwrite = TRUE, recursive = TRUE)

cores <- detectCores()-1 # adjust to your data
annotate_result1 <-  identify_metabolites(ms1.data = "for metID MS1.csv",  
                                          ms1.match.ppm = ppm, # adjust to your data
                                          rt.match.tol = 1000000, # adjust to your data
                                          polarity = "positive", # adjust to your data
                                          column = "rp", # adjust to your data
                                          path = path, 
                                          candidate.num = 3, # adjust to your data
                                          database = "hmdbMS1Database0.0.1", 
                                          threads = cores)

# get results
table <- get_identification_table(annotate_result1, candidate.num = 3, type = "new") # adjust to your data

# save
fwrite(table, "for metID/metID MS1 HMDB.csv", row.names = T)
# see "for metID" folder

#.....................identification by MONA.....................
file.copy(from = file.path(database, "monaDatabase0.0.1"), 
          to = path, overwrite = TRUE, recursive = TRUE)

cores <- detectCores()-1 # adjust to your data
annotate_result1 <-  identify_metabolites(ms1.data = "for metID MS1.csv",  
                                          ms1.match.ppm = ppm, # adjust to your data
                                          rt.match.tol = 1000000, # adjust to your data
                                          polarity = "positive", # adjust to your data
                                          column = "rp", # adjust to your data
                                          path = path, 
                                          candidate.num = 3, # adjust to your data
                                          database = "monaDatabase0.0.1", 
                                          threads = cores)

# get results
table <- get_identification_table(annotate_result1, candidate.num = 3, type = "new") # adjust to your data

# save
fwrite(table, "for metID/metID MS1 MONA.csv", row.names = T)
# see "for metID" folder

#..................identification by Mass Bank...................
file.copy(from = file.path(database, "massbankDatabase0.0.2"), 
          to = path, overwrite = TRUE, recursive = TRUE)

cores <- detectCores()-1 # adjust to your data
annotate_result1 <-  identify_metabolites(ms1.data = "for metID MS1.csv",  
                                          ms1.match.ppm = ppm, # adjust to your data
                                          rt.match.tol = 1000000, # adjust to your data
                                          polarity = "positive", # adjust to your data
                                          column = "rp", # adjust to your data
                                          path = path, 
                                          candidate.num = 3, # adjust to your data
                                          database = "massbankDatabase0.0.2", 
                                          threads = cores)

# get results
table <- get_identification_table(annotate_result1, candidate.num = 3, type = "new") # adjust to your data

# save
fwrite(table, "for metID/metID MS1 Mass Bank.csv", row.names = T)
# see "for metID" folder

#...................identification by Orbitrap...................
file.copy(from = file.path(database, "orbitrapDatabase0.0.1"), 
          to = path, overwrite = TRUE, recursive = TRUE)

cores <- detectCores()-1 # adjust to your data
annotate_result1 <-  identify_metabolites(ms1.data = "for metID MS1.csv",  
                                          ms1.match.ppm = ppm, # adjust to your data
                                          rt.match.tol = 1000000, # adjust to your data
                                          polarity = "positive", # adjust to your data
                                          column = "rp", # adjust to your data
                                          path = path, 
                                          candidate.num = 3, # adjust to your data
                                          database = "orbitrapDatabase0.0.1", 
                                          threads = cores)

# get results
table <- get_identification_table(annotate_result1, candidate.num = 3, type = "new") # adjust to your data

# save
fwrite(table, "for metID/metID MS1 Orbitrap.csv", row.names = T)
# see "for metID" folder

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              MetaboAnnotation                            ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(MetaboAnnotation)
library(xcms)
library(stringr)
library(batchCorr)

# load library
library(curl)
dbname <- "CompDb.Hsapiens.HMDB.5.0.sqlite"
db_file <- file.path(tempdir(), dbname)
curl_download(
  paste0("https://github.com/plyush1993/OUKS/",
         "releases/download/v1.10.1/", dbname), destfile = db_file)

# load library
# other way
library(piggyback)
dbname <- "CompDb.Hsapiens.HMDB.5.0.sqlite"
pb_download(dbname,
            repo = "plyush1993/OUKS",
            tag = "v1.10.1",
            dest = tempdir())
db_file <- file.path(tempdir(), dbname)

# CompoundDb
library(CompoundDb)
cdb <- CompDb(db_file)
cdb

# Get peak table from XCMS
load("xcms obj pk_fil.RData") # filename and the same folder of raw data as in section "XCMS with the best params (from IPO)" from script "2.Integration"
ft_inf <- featureDefinitions(pk_fil) # features info
ft_inf <- as.data.frame(ft_inf)
ft_inf <- as.data.frame(cbind(peak_id = paste(ft_inf$mzmed, ft_inf$rtmed, sep = " / "), ft_inf))
colnames(ft_inf)[2] <- "mz"
                         
# Alternative way
cn <- colnames(ds)
cn <- str_remove(cn, "X")
cn0 <- colnames(ds)
cn0 <- str_remove(cn0, "X")
cn <- gsub(pattern = "...", replacement = "_", x = cn, fixed = T)
cn2 <- t(data.frame(cn))
colnames(cn2) <- cn
peakIn <- peakInfo(PT = cn2, sep = "_", start = 1)
peakIn <- as.data.frame(cbind(name = cn0, peakIn))
peakIn$mz <- as.numeric(peakIn$mz)

# Set parameters for the m/z-based annotation
load("IPO_optimiz_xcms_CDF_7QC.RData")
ppm_diff <- resultPeakpicking[["best_settings"]][["parameters"]][["ppm"]] # adjust to your data
mz_diff <- 0.02  # resultPeakpicking[["best_settings"]][["parameters"]][["mzdiff"]] # adjust to your data
param <- Mass2MzParam(adducts = c("[M+H]+", "[M+Na]+"), # adjust to your data; consider "Mass2MzRtParam" 
                      tolerance = mz_diff, ppm = ppm_diff) # mz column name: defaults to "exactmass" in reference and "mz" in query

# perform
pks_match <- matchMz(query = peakIn, compounds(cdb, c("compound_id", "exactmass", "formula", "name")), param = param) # query = ft_inf or query = peakIn; also "matchValues" function
pks_match
metan <- as.data.frame(matchedData(pks_match))
md <- na.omit(metan) # delete not annotated
annot_feat <- dplyr::distinct(md[,1:3]) # list of only annotated                         

# save
fwrite(metan, "xcms MetaboAnnotation.csv", row.names = T)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                 References                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. Kuhl, Carsten, et al. "CAMERA: an integrated strategy for compound spectra extraction and annotation of liquid chromatography/mass spectrometry data sets." Analytical chemistry 84.1 (2012): 283-289.
# 2. Broeckling, Corey David, et al. "RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data." Analytical chemistry 86.14 (2014): 6812-6817.
# 3. Fernández-Ochoa, Álvaro, et al. "A Case Report of Switching from Specific Vendor-Based to R-Based Pipelines for Untargeted LC-MS Metabolomics." Metabolites 10.1 (2020): 28.
# 4. Uppal, Karan, Douglas I. Walker, and Dean P. Jones. "xMSannotator: an R package for network-based annotation of high-resolution metabolomics data." Analytical chemistry 89.2 (2017): 1063-1067.
# 5. Barranco-Altirriba, Maria, et al. "mWISE: An Algorithm for Context-Based Annotation of Liquid Chromatography-Mass Spectrometry Features through Diffusion in Graphs." Analytical Chemistry (2021).
# 6. Shen, Xiaotao, et al. "metID: an R package for automatable compound annotation for LC2MS-based data." Bioinformatics (2021).
# 7. Rainer, Johannes, et al. "A Modular and Expandable Ecosystem for Metabolomics Data Annotation in R." Metabolites 12.2 (2022): 173.
