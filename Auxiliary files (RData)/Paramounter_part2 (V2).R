###############################################################
#This is the script to perform parameters estimation
#Jian Guo, Tao Huan 2021-07-28
#Copyright @ University of British Columbia
###############################################################

start_time <- Sys.time()
mzDiff <- c()
ppm <- c()
noiselevel <- c()
peakWidth <- c()
peakScans <- c()
SNRatio <- c()
peakHeight <- c()
massShiftALL <- c()
rtShiftALL <- c()
mzDiff2D <- as.data.frame(matrix(ncol = 3, nrow = 1))
colnames(mzDiff2D) <- c("mz", "rt", "mzdiff")
ppm2D <- as.data.frame(matrix(ncol = 3, nrow = 1))
colnames(ppm2D) <- c("mz", "rt", "ppm")
Shift <- as.data.frame(matrix(ncol = 3, nrow = 1))
colnames(Shift) <- c("mz", "rt", "Sample")

#smoothing function
peak_smooth <- function(x,level=smooth){
  n <- level
  if(length(x) < 2*n){
    return(x)
  } else if(length(unique(x))==1){
    return(x)
  } else{
    y <- vector(length=length(x))
    for(i in 1:n){
      y[i] <- sum(c((n-i+2):(n+1),n:1)*x[1:(i+n)])/sum(c((n-i+2):(n+1),n:1))
    }
    for(i in (n+1):(length(y)-n)){
      y[i] <-  sum(c(1:(n+1),n:1)*x[(i-n):(i+n)])/sum(c(1:(n+1),n:1))
    }
    for(i in (length(y)-n+1):length(y)){
      y[i] <- sum(c(1:n,(n+1):(n+i-length(x)+1))*x[(i-n):length(x)])/sum(c(1:n,(n+1):(n+i-length(x)+1)))
    }
    return(y)
  }
}
#########################################################################################################
for (q in 1:(length(filename))){
  # Parameter setting
  ms1data <- readMSData(files = filename[q], mode = "onDisk", msLevel. = 1)
  mzRange <- c(min(unlist(mz(ms1data))), max(unlist(mz(ms1data))))
  ROI <- seq(mzRange[1], mzRange[2], 0.05)
  mzData <- mz(ms1data)
  intData <- intensity(ms1data)
  rtime <- rtime(ms1data)
  snALL <- c()
  sALL <- c()
  peakWidthALL <- c()
  peakScansALL <- c()
  ShiftTable <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(ShiftTable) <- c("mz", "rt", "Sample")
  ppm2Ddist <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(ppm2Ddist) <- c("mz", "rt", "ppm")
  mzdiff2Ddist <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(mzdiff2Ddist) <- c("mz", "rt", "mzdiff")
  noiseALL <- c()
  
  # ROI detection and universal parameter estimation 
  for(i in 1:(length(ROI) - 1)) {
    # Obtain data lists in each m/z bin 
    currmzRange <- c(ROI[i], ROI[i+1])
    tmpMZdata <- mzData
    tmpINTdata <- intData
    for(j in 1:length(mzData)){
      index <- which(tmpMZdata[[j]] >= currmzRange[1] & tmpMZdata[[j]] < currmzRange[2])
      tmpMZdata[[j]] <- tmpMZdata[[j]][index]
      tmpINTdata[[j]] <- tmpINTdata[[j]][index]
    }
    # Extract the intensity vectors from each m/z bin 
    eicINTraw <- c()
    eicINT <- c()
    eicRT <- c()
    for(k in 1:length(mzData)){
      if(length(tmpINTdata[[k]]) > 0){
        eicINTraw[k] <- mean(tmpINTdata[[k]])
      }else{
        eicINTraw[k] <- 0
      }
      eicRT[k] <- rtime[k]
    }
    if(sum(eicINTraw != 0) == 0) next()
    # Sort the intensity vectors from each m/z bin, estimate the noise cut off and average
    eicINT <- peak_smooth(eicINTraw)
    eicNon0 <- sort(eicINT[eicINT > 0])
    if(length(eicNon0) > 10){
      for(x in seq(10,length(eicNon0), 10)){
        sd <- sd(eicNon0[1:x])
        blk <- sum(eicNon0[1:x])/x
        thres <- blk + 3*sd
        if(x+1 <= length(eicNon0)){
          if(eicNon0[x+1] >= thres) break()
        }
      }
      cutOFF <- eicNon0[x]
    }else{
      cutOFF <- max(eicNon0)
      sd <- 0
      blk <- 0
    }
    noiseALL[i] <- cutOFF
    
    
    # Find the Reference m/z in each ROI from each m/z bin
    aboveTHindex <- which(eicINT > cutOFF)
    if(length(aboveTHindex) == 0) next()
    candidateSegInd <- split(aboveTHindex, cumsum(c(1, diff(aboveTHindex) != 1)))
    peakInd <- c()
    for(x in 1:length(candidateSegInd)){
      peakInd[x] <- which(eicINT[candidateSegInd[[x]]] == max(eicINT[candidateSegInd[[x]]]))[1] + min(candidateSegInd[[x]]) - 1
    }
    refMZvec <- c()
    for(y in 1:length(peakInd)){
      highestINT <- which(tmpINTdata[[peakInd[y]]] == max(tmpINTdata[[peakInd[y]]]))[1]
      refMZvec[y] <- tmpMZdata[[peakInd[y]]][highestINT]
    }
    
    # Estimate the universal parameters (mass tolerance, peak height, and peak width) for each m/z bin
    for(z in 1:length(peakInd)){
      currPeakInd <- peakInd[z]
      currRefMz <- refMZvec[z]
      currSamePeakMass <- c()
      currSamePeakMass <- c(currSamePeakMass, currRefMz)
      leftInd <- currPeakInd-1
      rightInd <- currPeakInd+1
      if(leftInd > 0){
        while (length(tmpMZdata[[leftInd]]) > 0 & mean(tmpINTdata[[leftInd]]) >= cutOFF) {
          if (length(tmpMZdata[[leftInd]]) == 1){
            currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[leftInd]])
            if(eicINT[leftInd] > eicINT[leftInd+1] & length(currSamePeakMass) > 5){
              Q1 <- as.numeric(summary(currSamePeakMass)[2])
              Q3 <- as.numeric(summary(currSamePeakMass)[5])
              LB <- Q1 - 1.5 *(Q3 - Q1)
              RB <- Q3 + 1.5 *(Q3 - Q1)
              if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
            }
          } else {
            abvector <- abs(tmpMZdata[[leftInd]] - currRefMz)
            NearInd <- which(abvector == min(abvector))[1]
            currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[leftInd]][NearInd])
            if(eicINT[leftInd] > eicINT[leftInd+1] & length(currSamePeakMass) > 5){
              Q1 <- as.numeric(summary(currSamePeakMass)[2])
              Q3 <- as.numeric(summary(currSamePeakMass)[5])
              LB <- Q1 - 1.5 *(Q3 - Q1)
              RB <- Q3 + 1.5 *(Q3 - Q1)
              if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
            }
          }
          leftInd <- leftInd-1
          if(leftInd <= 0) break()
        }
      }
      if(rightInd <= length(tmpMZdata)){
        while (length(tmpMZdata[[rightInd]]) > 0 & mean(tmpINTdata[[rightInd]]) >= cutOFF) {
          if (length(tmpMZdata[[rightInd]]) == 1){
            currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[rightInd]])
            if(eicINT[rightInd] > eicINT[rightInd-1] & length(currSamePeakMass) > 5){
              Q1 <- as.numeric(summary(currSamePeakMass)[2])
              Q3 <- as.numeric(summary(currSamePeakMass)[5])
              LB <- Q1 - 1.5 *(Q3 - Q1)
              RB <- Q3 + 1.5 *(Q3 - Q1)
              if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
            }
          } else {
            abvector <- abs(tmpMZdata[[rightInd]] - currRefMz)
            NearInd <- which(abvector == min(abvector))[1]
            currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[rightInd]][NearInd])
            if(eicINT[rightInd] > eicINT[rightInd-1] & length(currSamePeakMass) > 5){
              Q1 <- as.numeric(summary(currSamePeakMass)[2])
              Q3 <- as.numeric(summary(currSamePeakMass)[5])
              LB <- Q1 - 1.5 *(Q3 - Q1)
              RB <- Q3 + 1.5 *(Q3 - Q1)
              if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
            }
          }
          rightInd <- rightInd+1
          if(rightInd > length(tmpMZdata)) break()
        }
      }
      if(sum(is.na(currSamePeakMass)) > 0) next()
      if(length(currSamePeakMass) > 1 && length(currSamePeakMass) < 200){
        ppmCheck <- (massSDrange*sd(currSamePeakMass))/currRefMz * 1e6
        if(ppmCheck < ppmCut){
          daDiff <- massSDrange*sd(currSamePeakMass)
          ppmDiff <- (massSDrange*sd(currSamePeakMass))/currRefMz * 1e6
          currPeakWidth <- rtime[[rightInd - 1]] - rtime[[leftInd + 1]]
          currPeakScans <- rightInd - leftInd
          ppm2Ddist <- rbind(ppm2Ddist, c(currRefMz, rtime[[peakInd[z]]], ppmDiff))
          mzdiff2Ddist <- rbind(mzdiff2Ddist, c(currRefMz, rtime[[peakInd[z]]], daDiff))
          peakWidthALL <- c(peakWidthALL, currPeakWidth)
          peakScansALL <- c(peakScansALL, currPeakScans)
          if (cutOFF > 0){
            snALL <- c(snALL, (eicINT[currPeakInd]-blk)/sd)
          } else {
            snALL <- c(snALL, eicINT[currPeakInd]) 
          }
          sALL <- c(sALL, eicINT[currPeakInd])
          if (z == 1) {
            if(is.na(peakInd[z+1])){
              ShiftTable <- rbind(ShiftTable, c(currRefMz, rtime[[currPeakInd]], q))
            }
            if(!is.na(peakInd[z+1])) {
              if(rtime[[currPeakInd]] < (rtime[[peakInd[z+1]]] - 300)) {
                ShiftTable <- rbind(ShiftTable, c(currRefMz, rtime[[currPeakInd]], q))
              }
            } 
          } 
          if (z > 1) { 
            if(is.na(peakInd[z+1]) & rtime[[currPeakInd]] > (rtime[[peakInd[z-1]]] + 300)){
              ShiftTable <- rbind(ShiftTable, c(currRefMz, rtime[[currPeakInd]], q))
            }
            if(!is.na(peakInd[z+1])) {
              if(rtime[[currPeakInd]] < (rtime[[peakInd[z+1]]] - 300)) {
                ShiftTable <- rbind(ShiftTable, c(currRefMz, rtime[[currPeakInd]], q))
              }
            }  
          }
        }  
      }
    }
  }
  mzDiff2D <- rbind(mzDiff2D, mzdiff2Ddist)
  ppm2D <- rbind(ppm2D, ppm2Ddist)
  peakWidth <- c(peakWidth, peakWidthALL)
  peakScans <- c(peakScans, peakScansALL)
  SNRatio <- c(SNRatio, snALL)
  peakHeight <- c(peakHeight, sALL)
  Shift <- rbind(Shift, ShiftTable)
  noiselevel <- c(noiselevel, noiseALL)
}

# Estimation of the instrumental shift parameters
if (length(filename) > 1) {
  Shift <- Shift[complete.cases(Shift),]
  Shiftlist <- split(Shift, Shift[,3]) 
  if (length(Shiftlist) < 3){
    AlignmentT1 <- data.frame(matrix(nrow = 0, ncol = (ncol(Shiftlist[[1]]) + ncol(Shiftlist[[2]]))))
    for(i in 1:nrow(Shiftlist[[1]])){
      mass.lower.limit <- as.numeric(Shiftlist[[1]]$mz[i]) - 0.015
      mass.upper.limit <- as.numeric(Shiftlist[[1]]$mz[i]) + 0.015
      rt.lower.limit <- as.numeric(Shiftlist[[1]]$rt[i]) - 30
      rt.upper.limit <- as.numeric(Shiftlist[[1]]$rt[i]) + 30
      temp <- Shiftlist[[2]][(as.numeric(Shiftlist[[2]]$mz) >= as.numeric(mass.lower.limit) & 
                                as.numeric(Shiftlist[[2]]$mz) <= as.numeric(mass.upper.limit) &
                                as.numeric(Shiftlist[[2]]$rt) >= as.numeric(rt.lower.limit) &
                                as.numeric(Shiftlist[[2]]$rt) <= as.numeric(rt.upper.limit)),]
      temp <- temp[complete.cases(temp),]
      if(nrow(temp) > 0) {
        AlignmentT1 <- rbind(AlignmentT1, c(Shiftlist[[1]][i,], temp[1,]))
      }  
    }
    AlignmentT1 <- AlignmentT1[complete.cases(AlignmentT1),]
    colnames(AlignmentT1) <- c(colnames(Shiftlist[[1]]), colnames(Shiftlist[[2]]))
  } else {
    AlignmentT1 <- data.frame(matrix(nrow = 0, ncol = (ncol(Shiftlist[[1]]) + ncol(Shiftlist[[2]]))))
    for(i in 1:nrow(Shiftlist[[1]])){
      mass.lower.limit <- as.numeric(Shiftlist[[1]]$mz[i]) - 0.015
      mass.upper.limit <- as.numeric(Shiftlist[[1]]$mz[i]) + 0.015
      rt.lower.limit <- as.numeric(Shiftlist[[1]]$rt[i]) - 30
      rt.upper.limit <- as.numeric(Shiftlist[[1]]$rt[i]) + 30
      temp <- Shiftlist[[2]][(as.numeric(Shiftlist[[2]]$mz) >= as.numeric(mass.lower.limit) & 
                                as.numeric(Shiftlist[[2]]$mz) <= as.numeric(mass.upper.limit) &
                                as.numeric(Shiftlist[[2]]$rt) >= as.numeric(rt.lower.limit) &
                                as.numeric(Shiftlist[[2]]$rt) <= as.numeric(rt.upper.limit)),]
      temp <- temp[complete.cases(temp),]
      if(nrow(temp) > 0) {
        AlignmentT1 <- rbind(AlignmentT1, c(Shiftlist[[1]][i,], temp[1,]))
      }  
    }
    AlignmentT1 <- AlignmentT1[complete.cases(AlignmentT1),]
    colnames(AlignmentT1) <- c(colnames(Shiftlist[[1]]), colnames(Shiftlist[[2]]))
    
    for(k in 3:length(Shiftlist)){
      output <- data.frame(matrix(nrow = 0, ncol = ncol(AlignmentT1) + ncol(Shiftlist[[k]])))
      for(j in 1:nrow(AlignmentT1)){
        mass.lower.limit <- as.numeric(AlignmentT1[j,1]) - 0.015
        mass.upper.limit <- as.numeric(AlignmentT1[j,1]) + 0.015
        rt.lower.limit <- as.numeric(AlignmentT1[j,2]) - 30
        rt.upper.limit <- as.numeric(AlignmentT1[j,2]) + 30
        temp <- Shiftlist[[k]][(as.numeric(Shiftlist[[k]]$mz) >= as.numeric(mass.lower.limit) & 
                                  as.numeric(Shiftlist[[k]]$mz) <= as.numeric(mass.upper.limit) &
                                  as.numeric(Shiftlist[[k]]$rt) >= as.numeric(rt.lower.limit) &
                                  as.numeric(Shiftlist[[k]]$rt) <= as.numeric(rt.upper.limit)),]
        temp <- temp[complete.cases(temp),]
        if(nrow(temp) > 0) {
          tmprow <- cbind(AlignmentT1[j,], temp[1,])
          colnames(tmprow) <- colnames(output)
          output <- rbind(output, tmprow)
        }
      }
      colnames(output) <- c(colnames(AlignmentT1), colnames(Shiftlist[[k]]))
      AlignmentT1 <- data.frame(matrix(nrow = 0, ncol = ncol(output)))
      AlignmentT1 <- output[complete.cases(output),]
    }
  }
  massShift <- as.data.frame(matrix(ncol = length(Shiftlist), nrow = 1))
  rtShift <- as.data.frame(matrix(ncol = length(Shiftlist), nrow = 1))
  massShift <- AlignmentT1[,seq(1, 3*length(Shiftlist)-2 ,3)]
  for (i in 1:nrow(massShift)){
    massShiftALL[i] <- max(massShift[i,]) - min(massShift[i,]) 
  }
  rtShift <- AlignmentT1[,seq(2, 3*length(Shiftlist)-1 ,3)]
  for (i in 1:nrow(rtShift)){
    rtShiftALL[i] <- max(rtShift[i,]) - min(rtShift[i,]) 
  }
}
# plot the distribution of the universal parameters
noiselevel <- noiselevel[!is.na(noiselevel)]
noiselevel <- noiselevel[order(noiselevel)]
noiselevel <- noiselevel[1:round(length(noiselevel)*0.97)]
ppm2D <- ppm2D[complete.cases(ppm2D),]
ppm2D <- ppm2D[order(ppm2D[,3]),]
mzDiff2D <- mzDiff2D[complete.cases(mzDiff2D),]
mzDiff2D <- mzDiff2D[order(mzDiff2D[,3]),]
mzDiff2D <- mzDiff2D[1:round(nrow(mzDiff2D)*0.97),]
peakWidth <- peakWidth[!is.na(peakWidth)]
peakWidth <- peakWidth[order(peakWidth)]
peakWidth <- peakWidth[1:round(length(peakWidth)*0.97)]
peakScans <- peakScans[!is.na(peakScans)]
peakScans <- peakScans[order(peakScans)]
peakScans <- peakScans[1:round(length(peakScans)*0.97)]
SNRatio <- SNRatio[!is.na(SNRatio)]
SNRatio <- SNRatio[order(-SNRatio)]
SNRatio <- SNRatio[1:round(length(SNRatio)*0.97)]
peakHeight <- peakHeight[!is.na(peakHeight)]
peakHeight <- peakHeight[order(-peakHeight)]
peakHeight <- peakHeight[1:round(length(peakHeight)*0.97)]
if (length(filename) > 1) {
  massShiftALL <- massShiftALL[!is.na(massShiftALL)]
  massShiftALL <- massShiftALL[order(massShiftALL)]
  massShiftALL <- massShiftALL[1:round(length(massShiftALL)*0.97)]
  rtShiftALL <- rtShiftALL[!is.na(rtShiftALL)]
  rtShiftALL <- rtShiftALL[order(rtShiftALL)]
  rtShiftALL <- rtShiftALL[1:round(length(rtShiftALL)*0.97)]
  maxmassshift <- max(massShiftALL)
  maxrtshift <- max(rtShiftALL)
}
maxppm <- max(ppm2D[,3])
maxppm <- ceiling(maxppm)
maxmzdiff <- max(mzDiff2D[,3])
maxmzdiff <- (ceiling(maxmzdiff*100))/100
minnoise <- min(noiselevel)
minnoise <- floor(minnoise)
W <- mean(peakWidth, trim=0.05, na.rm = TRUE)
H <- mean(peakHeight, trim=0.05, na.rm = TRUE)
ratio <- H/W
minpeakwidth <- min(peakWidth)
maxpeakwidth <- max(peakWidth)
if (maxpeakwidth > 35 & ratio > 515) {
  minpeakwidth <- 0
  maxpeakwidth <- (ceiling(maxpeakwidth)+7)/2
} else {
  minpeakwidth <- (ceiling(minpeakwidth))+4
  maxpeakwidth <- ceiling(maxpeakwidth)+5
}
minpeakscan <- min(peakScans)
minpeakscan <- floor(minpeakscan)
maxpeakscan <- max(peakScans)
maxpeakscan <- floor(maxpeakscan)
minSN <- min(SNRatio)
if (minSN < 3){
  minSN <- 3
}
minpeakheight <- min(peakHeight)
minpeakheight <- floor(minpeakheight)

if (Software  == "XCMS"){
  if (length(filename) > 1) {
    setwd(directory)
    dir.create("Parameters for XCMS")
    setwd("Parameters for XCMS")
    png(file="Parameters for XCMS.png",width=1500, height=1300)
    par(mfrow=c(3,2))
    hist(noiselevel, xlab = "noise", xlim = c(0, 10000), breaks = seq(0, 500000000, 200), xaxp  = c(0, 10000, 50),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "noise", col = "black")
    plot(ppm2D$mz, ppm2D$ppm, ylab = "ppm", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    hist(SNRatio[which(SNRatio < 10)], xlab = "S/N ratio", xlim = c(0, 10), breaks = seq(0,10, 0.5), xaxp  = c(0, 10, 20),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak height", col = "black")
    hist(peakWidth[which(peakWidth <= 300)], xlab = "peakwidth (seconds)", xlim = c(0,100), breaks = seq(0,300,5), 
         xaxp  = c(0, 300, 60), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width", col = "black")
    hist(massShiftALL[which(massShiftALL <= 0.03)], xlab = "massShift (Da)", xlim = c(0,0.03), breaks = seq(0,0.03,0.001), 
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "instrumental mass shift", col = "black")
    hist(rtShiftALL[which(rtShiftALL <= 60)], xlab = "rtShift (seconds)", xlim = c(0,60), breaks = seq(0,60,2), 
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "instrumental retention time shift", col = "black")
    dev.off()
    XCMSparameters <- as.data.frame(matrix(ncol = 2, nrow = 14))
    colnames(XCMSparameters) <- c("Parameters", "Value")
    P <- c("ppm", "minimum peakwidth", "maximum peakwidth", "signal/noise threshold", "mzdiff(Please refer to the user manual for more detailed explanation)", "Integration method", 
           "prefilter peaks", "prefilter intensity", "noise filter", "bw", "minfrac", "mzwid", "minsamp", "max")
    V <- c(maxppm, minpeakwidth, maxpeakwidth, minSN, -0.01, 1, minpeakscan, minnoise, minnoise, 5, 0.5, maxmassshift, 1, 100)
    XCMSparameters[,1] <- P
    XCMSparameters[,2] <- round(V, 3)
    pdf(file = "Parameters for XCMS.pdf", height=6, width=10.5)
    grid.table(XCMSparameters)
    dev.off()
  } else {
    setwd(directory)
    dir.create("Parameters for XCMS")
    setwd("Parameters for XCMS")
    png(file="Parameters for XCMS.png",width=1500, height=1300)
    par(mfrow=c(2,2))
    hist(noiselevel, xlab = "noise", xlim = c(0, 10000), breaks = seq(0, 500000000, 200), xaxp  = c(0, 10000, 50),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "noise", col = "black")
    plot(ppm2D$mz, ppm2D$ppm, ylab = "ppm", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    hist(SNRatio[which(SNRatio < 10)], xlab = "S/N ratio", xlim = c(0, 10), breaks = seq(0,10, 0.5), xaxp  = c(0, 10, 20),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak height", col = "black")
    hist(peakWidth[which(peakWidth <= 300)], xlab = "peakwidth (seconds)", xlim = c(0,100), breaks = seq(0,300,5), 
         xaxp  = c(0, 300, 60), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width", col = "black")
    dev.off()
    XCMSparameters <- as.data.frame(matrix(ncol = 2, nrow = 9))
    colnames(XCMSparameters) <- c("Parameters", "Value")
    P <- c("ppm", "minimum peakwidth", "maximum peakwidth", "signal/noise threshold", "mzdiff(Please refer to the user manual for more detailed explanation)", "Integration method", 
           "prefilter peaks", "prefilter intensity", "noise filter")
    V <- c(maxppm, minpeakwidth, maxpeakwidth, minSN, -0.01, 1, minpeakscan, minnoise, minnoise)
    XCMSparameters[,1] <- P
    XCMSparameters[,2] <- round(V, 3)
    pdf(file = "Parameters for XCMS.pdf", height=6, width=10.5)
    grid.table(XCMSparameters)
    dev.off()  
  }
}
if (Software == "MSDIAL"){
  if (length(filename) > 1) {
    setwd(directory)
    dir.create("Parameters for MSDIAL")
    setwd("Parameters for MSDIAL")
    png(file="Parameters for MSDIAL.png",width=1500, height=1300)
    par(mfrow=c(3,2))
    plot(mzDiff2D$mz, mzDiff2D$mzdiff, ylab = "mzDiff", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    hist(log10(peakHeight)[which(log10(peakHeight) < 5)], xlab = "log10(peakheight)", xlim = c(1,5), xaxp  = c(0, 5, 20),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak height", col = "black")
    hist(peakScans[which(peakScans < 25)], xlab = "peakscannumbers", xlim = c(0,25), breaks = seq(0,25,1), 
         xaxp  = c(0, 25, 25), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width", col = "black")
    hist(massShiftALL[which(massShiftALL <= 0.03)], xlab = "massShift (Da)", xlim = c(0,0.03), breaks = seq(0,0.03,0.001), 
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "instrumental mass shift", col = "black")
    hist(rtShiftALL[which(rtShiftALL <= 60)], xlab = "rtShift (seconds)", xlim = c(0,60), breaks = seq(0,60,2), 
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "instrumental retention time shift", col = "black")
    dev.off()
    MSDIALparameters <- as.data.frame(matrix(ncol = 2, nrow = 6))
    colnames(MSDIALparameters) <- c("Parameters", "Value")
    P <- c("mass accuracy : MS1 tolerance", "peak detection : minimum peak height(Please refer to the user manual for more detailed explanation)", "peak detection : mass slice width", 
           "peak detection : minimum peak width", "alignment : MS1 tolerance", "alignment : retention time tolerance")
    V <- c(maxmzdiff, minpeakheight, maxmzdiff, minpeakscan, maxmassshift, maxrtshift/60)
    MSDIALparameters[,1] <- P
    MSDIALparameters[,2] <- round(V, 3)
    pdf(file = "Parameters for MSDIAL.pdf", height=6, width=10.5)
    grid.table(MSDIALparameters)
    dev.off()
  } else {
    setwd(directory)
    dir.create("Parameters for MSDIAL")
    setwd("Parameters for MSDIAL")
    png(file="Parameters for MSDIAL.png",width=1500, height=1300)
    par(mfrow=c(2,2))
    plot(mzDiff2D$mz, mzDiff2D$mzdiff, ylab = "mzDiff", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    hist(log10(peakHeight)[which(log10(peakHeight) < 5)], xlab = "log10(peakheight)", xlim = c(1,5), xaxp  = c(0, 5, 20),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak height", col = "black")
    hist(peakScans[which(peakScans < 25)], xlab = "peakscannumbers", xlim = c(0,25), breaks = seq(0,25,1), 
         xaxp  = c(0, 25, 25), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width", col = "black")
    dev.off()
    MSDIALparameters <- as.data.frame(matrix(ncol = 2, nrow = 4))
    colnames(MSDIALparameters) <- c("Parameters", "Value")
    P <- c("mass accuracy : MS1 tolerance", "peak detection : minimum peak height(Please refer to the user manual for more detailed explanation)", "peak detection : mass slice width", 
           "peak detection : minimum peak width")
    V <- c(maxmzdiff, minpeakheight, maxmzdiff, minpeakscan)
    MSDIALparameters[,1] <- P
    MSDIALparameters[,2] <- round(V, 3)
    pdf(file = "Parameters for MSDIAL.pdf", height=6, width=10.5)
    grid.table(MSDIALparameters)
    dev.off()
  }
}
if (Software == "MZMINE2"){
  if (length(filename) > 1) {
    setwd(directory)
    dir.create("Parameters for MZMINE2")
    setwd("Parameters for MZMINE2")
    png(file="Parameters for MZMINE2.png",width=1500, height=1600)
    par(mfrow=c(4,2))
    hist(noiselevel, xlab = "noise", xlim = c(0, 4000), breaks = seq(0, 500000000, 200), xaxp  = c(0, 4000, 40),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "noise", col = "black")
    plot(ppm2D$mz, ppm2D$ppm, ylab = "ppm", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    hist(SNRatio[which(SNRatio < 10)], xlab = "S/N ratio", xlim = c(0, 10), breaks = seq(0,10, 0.5), xaxp  = c(0, 10, 20),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak height", col = "black")
    hist(peakScans[which(peakScans < 25)], xlab = "peakscannumbers", xlim = c(0,25), breaks = seq(0,25,1), 
         xaxp  = c(0, 25, 25), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width (scan numbers)", col = "black")
    hist(peakWidth[which(peakWidth <= 120)], xlab = "peakwidth (seconds)", xlim = c(0,100), breaks = seq(0,100,5), 
         xaxp  = c(0, 300, 60), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width", col = "black")
    hist(massShiftALL[which(massShiftALL <= 0.03)], xlab = "massShift (Da)", xlim = c(0,0.03), breaks = seq(0,0.03,0.001), 
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "instrumental mass shift", col = "black")
    hist(rtShiftALL[which(rtShiftALL <= 60)], xlab = "rtShift (seconds)", xlim = c(0,60), breaks = seq(0,60,2), 
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "instrumental retention time shift", col = "black")
    dev.off()
    MZMINE2parameters <- as.data.frame(matrix(ncol = 2, nrow = 9))
    colnames(MZMINE2parameters) <- c("Parameters", "Value")
    P <- c("mass detection : noise level", "ADAP chromatogram builder : min group size in # of scans", 
           "ADAP chromatogram builder : group intensity threshold", "ADAP chromatogram builder : min highest intensity", 
           "ADAP chromatogram builder : m/z tolerance", "chromatogram deconvolution : peak duration range left", 
           "chromatogram deconvolution : peak duration range right", "alignment : m/z tolerance", "alignment : RT tolerance")
    V <- c(minnoise, minpeakscan, minnoise, minpeakheight, maxppm, minpeakwidth/60, maxpeakwidth/60, maxmassshift, maxrtshift/60)
    
    MZMINE2parameters[,1] <- P
    MZMINE2parameters[,2] <- round(V, 3)
    pdf(file = "Parameters for MZMINE2.pdf", height=6, width=12)
    grid.table(MZMINE2parameters)
    dev.off()
  } else {
    setwd(directory)
    dir.create("Parameters for MZMINE2")
    setwd("Parameters for MZMINE2")
    png(file="Parameters for MZMINE2.png",width=1500, height=1600)
    par(mfrow=c(3,2))
    hist(noiselevel, xlab = "noise", xlim = c(0, 4000), breaks = seq(0, 500000000, 200), xaxp  = c(0, 4000, 40),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "noise", col = "black")
    plot(ppm2D$mz, ppm2D$ppm, ylab = "ppm", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    hist(SNRatio[which(SNRatio < 10)], xlab = "S/N ratio", xlim = c(0, 10), breaks = seq(0,10, 0.5), xaxp  = c(0, 10, 20),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak height", col = "black")
    hist(peakScans[which(peakScans < 25)], xlab = "peakscannumbers", xlim = c(0,25), breaks = seq(0,25,1), 
         xaxp  = c(0, 25, 25), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width (scan numbers)", col = "black")
    hist(peakWidth[which(peakWidth <= 120)], xlab = "peakwidth (seconds)", xlim = c(0,100), breaks = seq(0,100,5), 
         xaxp  = c(0, 300, 60), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width", col = "black")
    dev.off()
    MZMINE2parameters <- as.data.frame(matrix(ncol = 2, nrow = 7))
    colnames(MZMINE2parameters) <- c("Parameters", "Value")
    P <- c("mass detection : noise level", "ADAP chromatogram builder : min group size in # of scans", 
           "ADAP chromatogram builder : group intensity threshold", "ADAP chromatogram builder : min highest intensity", 
           "ADAP chromatogram builder : m/z tolerance", "chromatogram deconvolution : peak duration range left", 
           "chromatogram deconvolution : peak duration range right")
    V <- c(minnoise, minpeakscan, minnoise, minpeakheight, maxppm, minpeakwidth/60, maxpeakwidth/60)
    
    MZMINE2parameters[,1] <- P
    MZMINE2parameters[,2] <- round(V, 3)
    pdf(file = "Parameters for MZMINE2.pdf", height=6, width=12)
    grid.table(MZMINE2parameters)
    dev.off()  
  }
}
if (Software == "XCMS" | Software == "MSDIAL" | Software == "MZMINE2" | Software == "ALL" ){
  if (length(filename) > 1) {
    setwd(directory)
    dir.create("Universal Parameters")
    setwd("Universal Parameters")
    png(file="Universal Parameters.png",width=1500, height=1600)
    par(mfrow=c(4,2))
    hist(noiselevel, xlab = "noise", xlim = c(0, 4000), breaks = seq(0, 500000000, 200), xaxp  = c(0, 4000, 40),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "noise", col = "black")
    plot(ppm2D$mz, ppm2D$ppm, ylab = "ppm", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    plot(mzDiff2D$mz, mzDiff2D$mzdiff, ylab = "mzDiff", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    hist(SNRatio[which(SNRatio < 10)], xlab = "S/N ratio", xlim = c(0, 10), breaks = seq(0,10, 0.5), xaxp  = c(0, 10, 20),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak height", col = "black")
    hist(peakScans[which(peakScans < 25)], xlab = "peakscannumbers", xlim = c(0,25), breaks = seq(0,25,1), 
         xaxp  = c(0, 25, 25), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width (scan numbers)", col = "black")
    hist(peakWidth[which(peakWidth <= 120)], xlab = "peakwidth (seconds)", xlim = c(0,100), breaks = seq(0,100,5), 
         xaxp  = c(0, 300, 60), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width", col = "black")
    hist(massShiftALL[which(massShiftALL <= 0.03)], xlab = "massShift (Da)", xlim = c(0,0.03), breaks = seq(0,0.03,0.001), 
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "instrumental mass shift", col = "black")
    hist(rtShiftALL[which(rtShiftALL <= 60)], xlab = "rtShift (seconds)", xlim = c(0,60), breaks = seq(0,60,2), 
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "instrumental retention time shift", col = "black")
    dev.off()
    Parameters <- as.data.frame(matrix(ncol = 2, nrow = 10))
    colnames(Parameters) <- c("Parameters", "Value")
    P <- c("mass difference (ppm)", "mass difference (Da)", 
           "min peak width (seconds)", "max peak width (seconds)", "min peak width (scan number)", "max peak width (scan number)",
           "peak height (intensity)", "peak height (S/N ratio)", 
           "Instrumental shift (mass shift)", "Instrumental shift (RT shift)")
    V <- c(maxppm, maxmzdiff, minpeakwidth, maxpeakwidth, minpeakscan, maxpeakscan, minpeakheight, minSN, maxmassshift, maxrtshift)
    
    Parameters[,1] <- P
    Parameters[,2] <- round(V, 3)
    pdf(file = "Universal parameters.pdf", height=6, width=12)
    grid.table(Parameters)
    dev.off()
  } else {
    setwd(directory)
    dir.create("Universal Parameters")
    setwd("Universal Parameters")
    png(file="Universal Parameters.png",width=1500, height=1600)
    par(mfrow=c(3,2))
    hist(noiselevel, xlab = "noise", xlim = c(0, 4000), breaks = seq(0, 500000000, 200), xaxp  = c(0, 4000, 40),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "noise", col = "black")
    plot(ppm2D$mz, ppm2D$ppm, ylab = "ppm", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    plot(mzDiff2D$mz, mzDiff2D$mzdiff, ylab = "mzDiff", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    hist(SNRatio[which(SNRatio < 10)], xlab = "S/N ratio", xlim = c(0, 10), breaks = seq(0,10, 0.5), xaxp  = c(0, 10, 20),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak height", col = "black")
    hist(peakScans[which(peakScans < 25)], xlab = "peakscannumbers", xlim = c(0,25), breaks = seq(0,25,1), 
         xaxp  = c(0, 25, 25), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width (scan numbers)", col = "black")
    hist(peakWidth[which(peakWidth <= 120)], xlab = "peakwidth (seconds)", xlim = c(0,100), breaks = seq(0,100,5), 
         xaxp  = c(0, 300, 60), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width", col = "black")
    dev.off()
    Parameters <- as.data.frame(matrix(ncol = 2, nrow = 8))
    colnames(Parameters) <- c("Parameters", "Value")
    P <- c("mass difference (ppm)", "mass difference (Da)", 
           "min peak width (seconds)", "max peak width (seconds)", "min peak width (scan number)", "max peak width (scan number)",
           "peak height (intensity)", "peak height (S/N ratio)")
    V <- c(maxppm, maxmzdiff, minpeakwidth, maxpeakwidth, minpeakscan, maxpeakscan, minpeakheight, minSN)
    
    Parameters[,1] <- P
    Parameters[,2] <- round(V, 3)
    pdf(file = "Universal parameters.pdf", height=6, width=12)
    grid.table(Parameters)
    dev.off()
  }
}

if (Software == "ALL"){
  if (length(filename) > 1) {
    setwd(directory)
    dir.create("Parameters for XCMS")
    setwd("Parameters for XCMS")
    png(file="Parameters for XCMS.png",width=1500, height=1300)
    par(mfrow=c(3,2))
    hist(noiselevel, xlab = "noise", xlim = c(0, 10000), breaks = seq(0, 500000000, 200), xaxp  = c(0, 10000, 50),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "noise", col = "black")
    plot(ppm2D$mz, ppm2D$ppm, ylab = "ppm", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    hist(SNRatio[which(SNRatio < 10)], xlab = "S/N ratio", xlim = c(0, 10), breaks = seq(0,10, 0.5), xaxp  = c(0, 10, 20),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak height", col = "black")
    hist(peakWidth[which(peakWidth <= 300)], xlab = "peakwidth (seconds)", xlim = c(0,100), breaks = seq(0,300,5), 
         xaxp  = c(0, 300, 60), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width", col = "black")
    hist(massShiftALL[which(massShiftALL <= 0.03)], xlab = "massShift (Da)", xlim = c(0,0.03), breaks = seq(0,0.03,0.001), 
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "instrumental mass shift", col = "black")
    hist(rtShiftALL[which(rtShiftALL <= 60)], xlab = "rtShift (seconds)", xlim = c(0,60), breaks = seq(0,60,2), 
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "instrumental retention time shift", col = "black")
    dev.off()
    XCMSparameters <- as.data.frame(matrix(ncol = 2, nrow = 14))
    colnames(XCMSparameters) <- c("Parameters", "Value")
    P <- c("ppm", "minimum peakwidth", "maximum peakwidth", "signal/noise threshold", "mzdiff(Please refer to the user manual for more detailed explanation)", "Integration method", 
           "prefilter peaks", "prefilter intensity", "noise filter", "bw", "minfrac", "mzwid", "minsamp", "max")
    V <- c(maxppm, minpeakwidth, maxpeakwidth, minSN, -0.01, 1, minpeakscan, minnoise, minnoise, 5, 0.5, maxmassshift, 1, 100)
    XCMSparameters[,1] <- P
    XCMSparameters[,2] <- round(V, 3)
    pdf(file = "Parameters for XCMS.pdf", height=6, width=10.5)
    grid.table(XCMSparameters)
    dev.off()
    
    setwd(directory)
    dir.create("Parameters for MSDIAL")
    setwd("Parameters for MSDIAL")
    png(file="Parameters for MSDIAL.png",width=1500, height=1300)
    par(mfrow=c(3,2))
    plot(mzDiff2D$mz, mzDiff2D$mzdiff, ylab = "mzDiff", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    hist(log10(peakHeight)[which(log10(peakHeight) < 5)], xlab = "log10(peakheight)", xlim = c(1,5), xaxp  = c(0, 5, 20),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak height", col = "black")
    hist(peakScans[which(peakScans < 25)], xlab = "peakscannumbers", xlim = c(0,25), breaks = seq(0,25,1), 
         xaxp  = c(0, 25, 25), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width", col = "black")
    hist(massShiftALL[which(massShiftALL <= 0.03)], xlab = "massShift (Da)", xlim = c(0,0.03), breaks = seq(0,0.03,0.001), 
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "instrumental mass shift", col = "black")
    hist(rtShiftALL[which(rtShiftALL <= 60)], xlab = "rtShift (seconds)", xlim = c(0,60), breaks = seq(0,60,2), 
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "instrumental retention time shift", col = "black")
    dev.off()
    MSDIALparameters <- as.data.frame(matrix(ncol = 2, nrow = 6))
    colnames(MSDIALparameters) <- c("Parameters", "Value")
    P <- c("mass accuracy : MS1 tolerance", "peak detection : minimum peak height(Please refer to the user manual for more detailed explanation)", "peak detection : mass slice width", 
           "peak detection : minimum peak width", "alignment : MS1 tolerance", "alignment : retention time tolerance")
    V <- c(maxmzdiff, minpeakheight, maxmzdiff, minpeakscan, maxmassshift, maxrtshift/60)
    MSDIALparameters[,1] <- P
    MSDIALparameters[,2] <- round(V, 3)
    pdf(file = "Parameters for MSDIAL.pdf", height=6, width=10.5)
    grid.table(MSDIALparameters)
    dev.off()
    
    setwd(directory)
    dir.create("Parameters for MZMINE2")
    setwd("Parameters for MZMINE2")
    png(file="Parameters for MZMINE2.png",width=1500, height=1600)
    par(mfrow=c(4,2))
    hist(noiselevel, xlab = "noise", xlim = c(0, 4000), breaks = seq(0, 500000000, 200), xaxp  = c(0, 4000, 40),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "noise", col = "black")
    plot(ppm2D$mz, ppm2D$ppm, ylab = "ppm", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    hist(SNRatio[which(SNRatio < 10)], xlab = "S/N ratio", xlim = c(0, 10), breaks = seq(0,10, 0.5), xaxp  = c(0, 10, 20),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak height", col = "black")
    hist(peakScans[which(peakScans < 25)], xlab = "peakscannumbers", xlim = c(0,25), breaks = seq(0,25,1), 
         xaxp  = c(0, 25, 25), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width (scan numbers)", col = "black")
    hist(peakWidth[which(peakWidth <= 120)], xlab = "peakwidth (seconds)", xlim = c(0,100), breaks = seq(0,100,5), 
         xaxp  = c(0, 300, 60), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width", col = "black")
    hist(massShiftALL[which(massShiftALL <= 0.03)], xlab = "massShift (Da)", xlim = c(0,0.03), breaks = seq(0,0.03,0.001), 
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "instrumental mass shift", col = "black")
    hist(rtShiftALL[which(rtShiftALL <= 60)], xlab = "rtShift (seconds)", xlim = c(0,60), breaks = seq(0,60,2), 
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "instrumental retention time shift", col = "black")
    dev.off()
    MZMINE2parameters <- as.data.frame(matrix(ncol = 2, nrow = 9))
    colnames(MZMINE2parameters) <- c("Parameters", "Value")
    P <- c("mass detection : noise level", "ADAP chromatogram builder : min group size in # of scans", 
           "ADAP chromatogram builder : group intensity threshold", "ADAP chromatogram builder : min highest intensity", 
           "ADAP chromatogram builder : m/z tolerance", "chromatogram deconvolution : peak duration range left", 
           "chromatogram deconvolution : peak duration range right", "alignment : m/z tolerance", "alignment : RT tolerance")
    V <- c(minnoise, minpeakscan, minnoise, minpeakheight, maxppm, minpeakwidth/60, maxpeakwidth/60, maxmassshift, maxrtshift/60)
    
    MZMINE2parameters[,1] <- P
    MZMINE2parameters[,2] <- round(V, 3)
    pdf(file = "Parameters for MZMINE2.pdf", height=6, width=12)
    grid.table(MZMINE2parameters)
    dev.off()
  } else {
    setwd(directory)
    dir.create("Parameters for XCMS")
    setwd("Parameters for XCMS")
    png(file="Parameters for XCMS.png",width=1500, height=1300)
    par(mfrow=c(2,2))
    hist(noiselevel, xlab = "noise", xlim = c(0, 10000), breaks = seq(0, 500000000, 200), xaxp  = c(0, 10000, 50),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "noise", col = "black")
    plot(ppm2D$mz, ppm2D$ppm, ylab = "ppm", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    hist(SNRatio[which(SNRatio < 10)], xlab = "S/N ratio", xlim = c(0, 10), breaks = seq(0,10, 0.5), xaxp  = c(0, 10, 20),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak height", col = "black")
    hist(peakWidth[which(peakWidth <= 300)], xlab = "peakwidth (seconds)", xlim = c(0,100), breaks = seq(0,300,5), 
         xaxp  = c(0, 300, 60), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width", col = "black")
    dev.off()
    XCMSparameters <- as.data.frame(matrix(ncol = 2, nrow = 9))
    colnames(XCMSparameters) <- c("Parameters", "Value")
    P <- c("ppm", "minimum peakwidth", "maximum peakwidth", "signal/noise threshold", "mzdiff(Please refer to the user manual for more detailed explanation)", "Integration method", 
           "prefilter peaks", "prefilter intensity", "noise filter")
    V <- c(maxppm, minpeakwidth, maxpeakwidth, minSN, -0.01, 1, minpeakscan, minnoise, minnoise)
    XCMSparameters[,1] <- P
    XCMSparameters[,2] <- round(V, 3)
    pdf(file = "Parameters for XCMS.pdf", height=6, width=10.5)
    grid.table(XCMSparameters)
    dev.off()  
    
    setwd(directory)
    dir.create("Parameters for MSDIAL")
    setwd("Parameters for MSDIAL")
    png(file="Parameters for MSDIAL.png",width=1500, height=1300)
    par(mfrow=c(2,2))
    plot(mzDiff2D$mz, mzDiff2D$mzdiff, ylab = "mzDiff", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    hist(log10(peakHeight)[which(log10(peakHeight) < 5)], xlab = "log10(peakheight)", xlim = c(1,5), xaxp  = c(0, 5, 20),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak height", col = "black")
    hist(peakScans[which(peakScans < 25)], xlab = "peakscannumbers", xlim = c(0,25), breaks = seq(0,25,1), 
         xaxp  = c(0, 25, 25), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width", col = "black")
    dev.off()
    MSDIALparameters <- as.data.frame(matrix(ncol = 2, nrow = 4))
    colnames(MSDIALparameters) <- c("Parameters", "Value")
    P <- c("mass accuracy : MS1 tolerance", "peak detection : minimum peak height(Please refer to the user manual for more detailed explanation)", "peak detection : mass slice width", 
           "peak detection : minimum peak width")
    V <- c(maxmzdiff, minpeakheight, maxmzdiff, minpeakscan)
    MSDIALparameters[,1] <- P
    MSDIALparameters[,2] <- round(V, 3)
    pdf(file = "Parameters for MSDIAL.pdf", height=6, width=10.5)
    grid.table(MSDIALparameters)
    dev.off()
    
    setwd(directory)
    dir.create("Parameters for MZMINE2")
    setwd("Parameters for MZMINE2")
    png(file="Parameters for MZMINE2.png",width=1500, height=1600)
    par(mfrow=c(3,2))
    hist(noiselevel, xlab = "noise", xlim = c(0, 4000), breaks = seq(0, 500000000, 200), xaxp  = c(0, 4000, 40),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "noise", col = "black")
    plot(ppm2D$mz, ppm2D$ppm, ylab = "ppm", xlab = "m/z",
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "mass tolerance")
    hist(SNRatio[which(SNRatio < 10)], xlab = "S/N ratio", xlim = c(0, 10), breaks = seq(0,10, 0.5), xaxp  = c(0, 10, 20),
         cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak height", col = "black")
    hist(peakScans[which(peakScans < 25)], xlab = "peakscannumbers", xlim = c(0,25), breaks = seq(0,25,1), 
         xaxp  = c(0, 25, 25), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width (scan numbers)", col = "black")
    hist(peakWidth[which(peakWidth <= 120)], xlab = "peakwidth (seconds)", xlim = c(0,100), breaks = seq(0,100,5), 
         xaxp  = c(0, 300, 60), cex.main=4, cex.lab=1.7, cex.axis=2, main = "peak width", col = "black")
    dev.off()
    MZMINE2parameters <- as.data.frame(matrix(ncol = 2, nrow = 7))
    colnames(MZMINE2parameters) <- c("Parameters", "Value")
    P <- c("mass detection : noise level", "ADAP chromatogram builder : min group size in # of scans", 
           "ADAP chromatogram builder : group intensity threshold", "ADAP chromatogram builder : min highest intensity", 
           "ADAP chromatogram builder : m/z tolerance", "chromatogram deconvolution : peak duration range left", 
           "chromatogram deconvolution : peak duration range right")
    V <- c(minnoise, minpeakscan, minnoise, minpeakheight, maxppm, minpeakwidth/60, maxpeakwidth/60)
    
    MZMINE2parameters[,1] <- P
    MZMINE2parameters[,2] <- round(V, 3)
    pdf(file = "Parameters for MZMINE2.pdf", height=6, width=12)
    grid.table(MZMINE2parameters)
    dev.off() 
  }
}
message("The parameters estimation has been completed. Please check the pdf file in the corresponding folder.")

