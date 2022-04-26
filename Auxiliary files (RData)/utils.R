# utils.
# read data function
readData = function(path =  "G:\\data\\D\\data project D.xlsx"){
  
  #check if it is csv of xlsx
  if(grepl("xlsx", path)){
    d <- openxlsx::read.xlsx(path, sheet = 1,colNames = FALSE)
  }else if(grepl("csv", path)){
    # file = "C:\\Users\\Sili Fan\\Downloads\\val (18).csv"
    d <- data.table::fread(path)
  }
  
  # make "" as NA
  d[d==""] <- NA
  
  #### fData
  fData <- d[!is.na(d[,1]),c(which(is.na(d[1,])),sum(is.na(d[1,]))+1)] # The first row and column is critical of formating the data.
  colnames(fData) = as.character(fData[1,]); fData = data.frame(fData[-1,],stringsAsFactors = F,check.names = FALSE);rownames(fData) = 1:nrow(fData);
  # following steps keeps the column type.
  fData.=lapply(fData,function(x){
    if(sum(!is.na(as.numeric(x))) == length(x)){
      as.numeric(x)
    }else{
      x
    }
  })
  fData. = do.call(cbind, lapply(fData., data.frame, stringsAsFactors=FALSE))
  colnames(fData.) = colnames(fData)
  fData = fData.
  
  fData = fData[,c(ncol(fData),2:ncol(fData)-1)]
  fData[[1]] = make.unique(make.names(fData[[1]]), sep = '_')
  
  #### pData
  pData <- d[c(which(is.na(d[,1])),max(which(is.na(d[,1])))+1) ,!is.na(d[1,])]
  pData <- t(pData); colnames(pData) = pData[1,]; pData = data.frame(pData[-1,],stringsAsFactors = F,check.names = FALSE)
  # following steps keeps the column type.
  pData.=lapply(pData,function(x){
    if(sum(!is.na(as.numeric(x))) == length(x)){
      as.numeric(x)
    }else{
      x
    }
  })
  pData. = do.call(cbind, lapply(pData., data.frame, stringsAsFactors=FALSE))
  colnames(pData.) = colnames(pData)
  pData = pData.
  
  pData = pData[,c(ncol(pData),2:ncol(pData)-1)]
  pData[[1]] = make.unique(make.names(pData[[1]]), sep = '_')
  
  #### eData
  eData <- d[!is.na(d[,1]),!is.na(d[1,])][-1,-1]
  eData <- sapply(eData, as.numeric)
  eData <- data.frame(eData,stringsAsFactors = F)
  colnames(eData) = pData[[1]]; rownames(eData) = fData[[1]]
  
  # # remove any unwanted character in columns of eData, fData and pData to _.
  # colnames(eData) = gsub("([_])|[[:punct:]]", "_", colnames(eData))
  # colnames(fData) = gsub("([_])|[[:punct:]]", "_", colnames(fData))
  # colnames(pData) = gsub("([_])|[[:punct:]]", "_", colnames(pData))
  
  # remove all the NA. And replace NA with "NA" Otherwise DataTables will give error.datatables warning requested unknown parameter
  # eData[is.na(eData)]="NA"
  # fData[is.na(fData)]="NA"
  # pData[is.na(pData)]="NA"
  
  # remove unwanted character in p.
  # for(i in 1:nrow(pData)){
  #   for(j in 1:ncol(pData)){
  #     pData[i,j] = gsub("\\+|~|-", " ", pData[i,j])
  #   }
  # }
  
  return(list(e = eData, f = fData, p = pData, original = d))
  
}



# remove outlier defined as >1.5 IQR.
remove_outlier = function(v){
  out = boxplot(v,plot=F)$out
  return(list(value = v[!v%in%out],index = which(v%in%out)))
}

# make sure there is no negative value after normalization and make sure each metabolite stays as the same level after normalization.
# can be used in case of negative value exsit. Not used in this paper.
shiftData = function(ori,norm){
              ori.min = apply(ori,1,min,na.rm=T)
              norm.min = apply(norm,1,min,na.rm=T)
              return(norm - c(norm.min - ori.min))
            }

# select best span parameter for loess.
get_loess_para = function(x,y,loess.span.limit = 0.5){ # use leave-one-out to select the best span parameter.
  j = 0
  error = rep(0, length(c(seq(loess.span.limit,1.5,0.1),1.75,2,2.25,2.5)))
  for(par in c(seq(loess.span.limit,1.25,0.1),1.5,2)){
    j = j+1
    for(i in 2:(length(y)-1)){ # if i from 1 or to length(y), the prediction would be NA
      o = loess(y[-i]~x[-i],span = par)
      if(sum(is.na(o))){
        error[j] = Inf
      }else{
        err = abs(predict(o,newdata = x[i])-y[i])
        error[j] = sum(error[j],err,na.rm = T)
      }
    }
  }
  return(c(seq(loess.span.limit,1.5,0.1),1.75,2,2.25,2.5)[which.min(error)])
}

# generate PCA plot.
generate_PCA = function(e, f, p, QC.index, batch,method){
  pca = prcomp(t(e), center = T, scale. = T)
  variance = pca$sdev^2/sum(pca$sdev^2)
  pca.data = data.frame(pca$x,batch = batch[1,],order = 1:nrow(pca$x))
  batch.QC = batch[1,];
  batch.QC[QC.index] = "QC"
  qc = rep(F, nrow(p))
  qc[QC.index] = TRUE 
  ggplot(pca.data, aes(PC1, PC2, color = batch.QC,size = qc, order = order)) +
    geom_point(alpha = 3/4) +
    stat_ellipse( linetype = 2, size = 0.5) +
    labs(x = paste0("PC1: ",signif(variance[1]*100,3),"%"), y = paste0("PC2: ",signif(variance[2]*100,3),"%"),
         title = method)+
    theme.scatter
}
