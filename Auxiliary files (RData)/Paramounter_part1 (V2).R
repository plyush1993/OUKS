###############################################################
#This is the script to perform parameters estimation
#Jian Guo, Tao Huan 2021-07-28
#Copyright @ University of British Columbia
###############################################################

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
  ppm2Ddist <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(ppm2Ddist) <- c("mz", "rt", "ppm")
  mzdiff2Ddist <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(mzdiff2Ddist) <- c("mz", "rt", "mzdiff")
  
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
    }
    
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
    ppmDiff <- c()
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
      
      if(length(currSamePeakMass) > 1){
        ppmDiff[z] <- (massSDrange*sd(currSamePeakMass))/currRefMz * 1e6
        ppm2Ddist <- rbind(ppm2Ddist, c(currRefMz, rtime[[peakInd[z]]], ppmDiff[z]))
      }
    }
  }
  
  ppm2D <- rbind(ppm2D, ppm2Ddist)
}

ppm2D <- ppm2D[complete.cases(ppm2D),]
ppm2D <- ppm2D[order(ppm2D[,3]),]
ppm2D <- ppm2D[1:round(nrow(ppm2D)*0.97),]
plot(ppm2D$mz, ppm2D$ppm, ylab = "ppm", xlab = "m/z", pch=1, cex.main=4, cex.lab=1.7, cex.axis=2)
ppm2Ddash <- ppm2D[1:round(nrow(ppm2D)*cutoff),]
dashline <- max(ppm2Ddash[,3]) 
long <- length(ppm2D[,3])
cutoffvalue <- rep(dashline,long)
lines(ppm2D$mz, cutoffvalue, lty = "dashed", lwd = "3", col = "red")
print(Sys.time() - start_time)
message("Please find the cutoff line in the generated ppm distribution, and run Paramounter part 2 using the ppm cutoff")
