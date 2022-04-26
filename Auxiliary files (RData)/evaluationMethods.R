# evaluation methods.
# RSD
RSD = function(e,f,p,
         robust = T,
         cl){
  library(parallel)
  if(robust){
    result=parSapply(cl=cl,X=1:nrow(e),FUN = function(i,remove_outlier,e){
      x = remove_outlier(e[i,])[[1]]
      sd(x,na.rm=T)/mean(x,na.rm=T)
    },remove_outlier,e)
  }else{
    result=parSapply(cl=cl,X=1:nrow(e),FUN = function(i,e){
      x = e[i,]
      sd(x,na.rm=T)/mean(x,na.rm=T)
    },e)
  }
  
  
  return(result)
}