# no normalization.
none_norm = function(e,p,f){
  return(list(e=e,p=p,f=f))
}

# mTIC normalization
mTIC_norm= function(e,f,p, known = "KnownORUnknown", subset = rep(T,nrow(p))){

  if(sum(!unique(f[[known]])%in%c(T,F))>1){
    stop(paste0("Error in TIC Normalization: The column ",known," must contain 'TRUE' and 'FALSE' only."))
  }


  index = rep(T, nrow(f))
  index[f[[known]]==FALSE | f[[known]]=="FALSE"] = FALSE

  sums = apply(e[index,subset], 2, sum, na.rm=T)
  mean_sums = mean(sums)
  e_mTIC_norm = t(t(e)/(sums/mean_sums))
  rownames(e_mTIC_norm) = rownames(e[,subset])
  colnames(e_mTIC_norm) = colnames(e[,subset])
  return(list(e = e_mTIC_norm,p=p,f=f))
}

# loess normalization
loess_norm = function(e,f,p,
                      batch,
                      QC.index = which(p$Comments=='qc'), time = "Acq. Date-Time",
                      span.para = 'auto',loess.span.limit=0.25){
  e = data.matrix(e)
  norms = parSapply(cl, X = 1:nrow(e), function(i,e,f,p,QC.index,batch,time,remove_outlier,span.para,get_loess_para,
                                                loess.span.limit){
    # for(i in 1:nrow(e)){
    models = by(data.frame(v=e[i,QC.index],t=p[[time]][QC.index]),
                factor(batch[i,QC.index]),function(x){

                    # x = data.frame(v=e[i,QC.index],t=p[[time]][QC.index])[batch[i,QC.index]=="A",]
                    if(length(remove_outlier(x$v)[[2]])>0){# if outlier exists.
                      span = ifelse(span.para=='auto',
                                    get_loess_para(x=x$t[-remove_outlier(x$v)[[2]]],y=remove_outlier(x$v)[[1]],
                                                   loess.span.limit = loess.span.limit),span.para) # find a proper span.
                    }else{
                      span = ifelse(span.para=='auto',
                                    get_loess_para(x=x$t,y=x$v,
                                                   loess.span.limit = loess.span.limit),span.para) # find a proper span.

                    }
                    if(length(remove_outlier(x$v)[[2]])>0){
                      loess(v~t,data=x[-remove_outlier(x$v)[[2]],],span=span)
                    }else{
                      loess(v~t,data=x,span=span)
                    }


                })
    # }
    # predict using the models.
    norm = mapply(function(u,v){
      # o = tryCatch({
        predict(u,newdata = v)
      # },
      # error = function(e){
        # print(e)
        # v
      # })
    },models,by(p[[time]],batch[i,],function(x){x}))
    norm = unlist(norm)
    # replace NA with the closest value.
    if(length(which(is.na(norm)))>0){
      for(j in which(is.na(norm))){
        time_notNA = p[[time]][-which(is.na(norm))]
        closest_time = time_notNA[which.min(abs(time_notNA - p[[time]][j]))]
        norm[j] = norm[which(p[[time]]==closest_time)]
      }
    }
    return(norm)
  },e,f,p,QC.index,batch,time,remove_outlier,span.para,get_loess_para,loess.span.limit)
  norms = t(norms)
  e_norm = matrix(NA,nrow=nrow(f),ncol=nrow(p))

  for(i in 1:nrow(f)){
    e_norm[i,] = e[i,] / (norms[i,] / median(e[i,]))
  }

  rownames(e_norm) = rownames(e)
  colnames(e_norm) = colnames(e)
  return(list(e = e_norm,f=f,p=p,normalize_line = norms))
}

# SERRF normalization
SERRF_norm = function(e,f,p,
                      batch = define_batch(e,f,p),
                      QC.index, time = "Acq. Date-Time"){
  qc = rep(F, nrow(p))
  qc[QC.index] = T
  e. = e
  for(i in 1:nrow(e)){ # MAKE SURE THE QC AND SAMPLES ARE AT THE SAME LEVEL. This is critical for SERRF algorithm (and other tree-based machine learning algorithm) because when building each tree, the split on each leaf considers the level of the values. If the values are not consistant, then the RF models will be wrong and the RF will bias the intensity level after normalization (although the relative position won't change.)
    e.[i,qc] = unlist(by(data.frame(e.[i,],qc),batch[1,],function(x){# x = data.frame(e.[i,],qc)[batch[1,]=='A',]
      diff =  (median(x[x[,2],1]) - median(x[!x[,2],1]))
      x[x[,2],1] - diff
    }))
  }
  pred = parSapply(cl, X = 1:nrow(f), function(j,eData,batch,randomForest, QC.index, time){
    set.seed(1)
    data = data.frame(y = eData[j,], t(eData[-j,]), batch = batch[1,], time = time)
    colnames(data) = c("y", paste0("X",1:nrow(eData))[-j], "batch", "time")
    model = randomForest(y~., data = data,subset = QC.index, importance = F)
    newdata = data.frame(t(eData[-j,]), batch = batch[1,], time = time)
    colnames(newdata) =   c(paste0("X",1:nrow(eData))[-j], "batch", "time")
    new = (eData[j,]/predict(model,newdata = newdata))*median(eData[j,])
    return(new)
  }, e.,batch,randomForest, QC.index, p[[time]])

  e_SERRF_pred = t(pred)
  # put the QC level bach to where they were. Won't influence the value of samples.
  for(i in 1:nrow(e_SERRF_pred)){
    e_SERRF_pred[i,qc]  = e_SERRF_pred[i, qc] + (median(e[i,qc], na.rm = T) - median(e[i,!qc], na.rm = T))
    e_SERRF_pred[i,e_SERRF_pred[i,]<0] = .5 * min(e_SERRF_pred[i,e_SERRF_pred[i,]>0])
  }
  #
  return(list(e = e_SERRF_pred, p = p, f = f))
}
# SVM normalization (simplified from https://raw.githubusercontent.com/jaspershen/MetNormalizer/master/R/SXTsvrNor1.R)
SVM_norm = function(e,f,p, QC.index, multiple = 5, time = "Acq. Date-Time"){
  library(e1071)
  e_SVM_pred = e
  qc. = QC.index
  for(i in 1:nrow(f)){
    all.cor <- apply(e[,qc.], 1,function(x) {cor(e[1,qc.], x)}) # Can't involve external-validates when building model at all
    cor.peak <-match(sort(all.cor, decreasing = TRUE)[2:(as.numeric(multiple)+1)], all.cor)
    if (multiple != 1) {
      svr.reg <- svm(t(e[cor.peak,qc.]),e[i,qc.])
      newdata = t(e[cor.peak, ])
      pred = predict(svr.reg, newdata = newdata)
      e_SVM_pred[i,] = (e[i,]/pred)*median(e[i,qc.])
    } else{
      svr.reg <- svm(e[i,qc.] ~ p[[time]][qc.])
      pred = predict(svr.reg, newdata = p[[time]])
      e_SVM_pred[i,] = (e[i,]/pred)*median(e[i,qc.])
    }
  }
  return(list(e = e_SVM_pred, p = p, f = f))
}

# sum normalization
sum_norm = function(e,f,p){
  sums = colSums(e,na.rm = T)
  result = t(t(e)/sums*median(sums,na.rm = T))
  return(list(e=result,f=f,p=p))
}

# median normalization
median_norm = function(e,f,p){
  medians = apply(e, 2, median, na.rm = T)
  result = t(t(e)/medians*median(medians,na.rm = T))
  return(list(e=result,f=f,p=p))
}

# PQN normalization
PQN_norm = function(e,f,p) {
  reference <- apply(e, 1, median)
  reference[reference==0] = 1
  quotient <- e/reference
  quotient.median <- apply(quotient, 2, median)
  e.norm <- t(t(e)/quotient.median)
  return(list(e = e.norm,f=f,p=p))
}

# CONTRAST normalization
contrast_norm = function(e,f,p) {
  library(affy)
  threshold=1e-11
  e[e <= 0] <- threshold
  maffy.data <- maffy.normalize(data.matrix(e),
                                subset=1:nrow(e),
                                span=0.75,
                                verbose=FALSE,
                                family="gaussian",
                                log.it=FALSE)
  subtract <- function(x){
    t(t(x)-apply(x,2,quantile,0.1))
  }
  contrast.data <- subtract(maffy.data)
  rownames(contrast.data) <- rownames(e)
  colnames(contrast.data) <- colnames(e)
  return(list(e = contrast.data,f=f,p=p))
}


# Quantile normalization
quantile_norm = function(e,f,p) {
  library(affy)
  normalize.quantile <- get("normalize.quantiles",en=asNamespace("affy"))
  quantile.data <- normalize.quantile(data.matrix(e))
  rownames(quantile.data) <- rownames(e)
  colnames(quantile.data) <- colnames(e)
  return(list(e = quantile.data,f=f,p=p))
}

# linear normalization
linear_norm = function(e,f,p) {
  library(affy)
  linear.baseline <- apply(e, 1, median)
  baseline.mean <- mean(linear.baseline)
  sample.means <- apply(e, 2, mean)
  linear.scaling <- baseline.mean/sample.means
  linear.baseline.data <- t(t(e) * linear.scaling)

  return(list(e = linear.baseline.data,f=f,p=p))
}

# Li-Wong normalization
liwong_norm = function(e,f,p) {
  library(affy)
  #---First step: Find baseline sample
  average.intensity <- apply(e,2,mean)
  #R has an add way of rounding.
  median.number <- round(ncol(e)/2 + 0.1)
  #the additional 0.1 ensures that it rounds properly
  ordering <- order(average.intensity)
  median.sample.number <- ordering[median.number]
  median.sample <- e[,median.sample.number]
  #---Apply normalization
  liwong.data <- vector()
  for(i in 1:ncol(e)){
    liwong.model <- normalize.invariantset(data=e[,i],
                                           ref=median.sample,
                                           prd.td=c(0.003,0.007))
    #the threshold of the rank-invariant set might need to be adjusted from case to case
    liwong.sample <- predict(liwong.model$n.curve$fit, e[,i])
    liwong.data <- cbind(liwong.data,liwong.sample$y)
  }
  colnames(liwong.data) = colnames(e)
  return(list(e = liwong.data,f=f,p=p))
}

# cubic normalization
cubic_norm = function(e,f,p) {
  library(affy)
  spline.data <- normalize.qspline(e,samples=0.02,target=apply(e,1,mean),verbose = FALSE)
  rownames(spline.data) <- rownames(e)
  colnames(spline.data) <- colnames(e)
  return(list(e = spline.data,f=f,p=p))
}

# batch normalizatio
batchratio_norm = function(e = e, f=f, p=p,batch,QC.index){
  e_batch_norm = matrix(,nrow=nrow(e),ncol=ncol(e))
  for(i in 1:nrow(f)){
    means = by(as.numeric(e[i,QC.index]),batch[i,QC.index], mean, na.rm=T)
    mean_means = mean(means)
    e_batch_norm[i,] = as.numeric(e[i,])/(rep(means,times=table(batch[i,]))/mean_means)
  }
  rownames(e_batch_norm) = rownames(e)
  colnames(e_batch_norm) = colnames(e)
  return(list(e=e_batch_norm,f=f,p=p))
}

# cyclic normalization
cyclic_norm <- function(e, f, p) {
  library(affy)
  loess.data <- normalize.loess(e, subset = 1:nrow(e),
                                epsilon = 10^-2, maxit = 2,
                                log.it = FALSE, verbose = FALSE,
                                span = 0.75, family.loess = "gaussian")
  return(list(e=loess.data,f=f,p=p))
}



generate_PCA = function(e,f,p, batch , QC.index, method){
  pca = prcomp(t(e), center = T, scale. = T)
  variance = pca$sdev^2/sum(pca$sdev^2)
  pca.data = data.frame(pca$x,batch = batch[1,],order = 1:nrow(pca$x))
  batch.QC = batch[1,];batch.QC[QC.index] = "QC"
  qc = rep(F, ncol(e))
  qc[QC.index] = TRUE

  plot = ggplot(pca.data, aes(PC1, PC2, color = batch.QC,size = qc, order = order)) +
    geom_point(alpha = 3/4) +
    stat_ellipse( linetype = 2, size = 0.5) +
    labs(x = paste0("PC1: ",signif(variance[1]*100,3),"%"), y = paste0("PC2: ",signif(variance[2]*100,3),"%"),
         title = method)+
    theme.scatter
  return(plot)

}









