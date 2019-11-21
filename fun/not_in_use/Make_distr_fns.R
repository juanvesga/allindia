# Makes distribution functions to construct likelihood woith targetd data
#data : target list of data points
# W: 0=no weight, 1= weight likelihood  by number of points in category(for multi-target fitting)

Make_distr_fns<-function (data,weight){
  
  funclist <- list()
  ii<-0
  for (i in 1:length(data))
  {
    
    dat<-data[[i]]
    for (jj in 1:length(dat)){
      ii<-ii+1
      W<-length(dat)*weight + (1-weight)
      
      distr<-c()
      if (dat[jj]>=1) {
        distr<-"lognormal"
      }else if (dat[jj]<1){
        distr<-"beta"
      }
      prctles<-dat[jj]*c(0.8, 1,1.2)
      
      funclist[[ii]]<-get_distribution_fns(prctles,distr,W)
      
    }
  }
  return(funclist)
}
