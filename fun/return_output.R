return_output<-function(fx,x,ref,location,varargin){
  
  i <- ref$i; s <- ref$s
  xvals<-x
  nsims<-1
  runtype<-'mle'
  burn<-0; thin<-0
  if (!isempty(varargin)){
    tmp<-varargin; burn<-tmp[1];thin<-tmp[2]
    xvals<-x[burn:thin:length(x),]
    nsims<-size(xvals,1)
    runtype<-'mcmc';
  }
  # Pre allocate list of results
  res<-vector("list",19)
  for (jj in 1:8){
    res[[jj]]<-matrix(0,nsims,5)
  }
  res[[9]]<-matrix(0,nsims,3)
  res[[10]]<-matrix(0,nsims,3)
  
  res[[11]]<-matrix(0,nsims,i$nx)
  res[[12]]<-matrix(0,nsims,4)
  res[[13]]<-matrix(0,nsims,5)
  
  for (jj in 14:19){
    res[[jj]]<-matrix(0,nsims,1)
  }
  names(res)<-c('popu','inc_all','inc_hi','inc_mdr','prev_all','prev_hi','mort_tb',
                'notif_all','inc_remote','inc_recent',
                'sfin',
                'cas_ds','cas_mdr',
                'llk','pu','prpe','pr_mdr','prev_ratio','cas_prev')
  
  
  
  for (ii in 1:nsims){
    
    x0<-xvals[,ii]
    
    model <- fx(x0)
    
    t<-model$t
    sln<-model$soln
    
    tmp <- rowSums(sln[,1:i$nstates])
    pop <- (tmp[1:length(tmp)-1]+tmp[2:length(tmp)])/2
    
    tmp    <- rowSums(sln[,s$hi])
    pop_hi <- (tmp[1:length(tmp)-1]+tmp[2:length(tmp)])/2
    
    tmp_n <- apply(sln[,i$aux$inc],2,diff)
    inc   <- 1e5*tmp_n/cbind(pop,pop_hi,pop)
    
    inc_all<-inc[,1]
    inc_hi <-inc[,2]
    inc_mdr<-inc[,3]
    
    tmp <- rowSums(sln[,s$prevalent])
    pre <- (tmp[1:length(tmp)-1]+tmp[2:length(tmp)])/2;
    prev= 1e5*(pre/pop)
    
    tmp <- rowSums(sln[,intersect(s$hi,s$prevalent)])
    pre <- (tmp[1:length(tmp)-1]+tmp[2:length(tmp)])/2;
    prev_hi= 1e5*(pre/pop_hi)
    
    tmp_n <- diff(sln[,i$aux$mort])
    mort  <- 1e5*tmp_n/pop
    
    tmp_n <- diff(sln[,i$aux$noti[1]])
    not  <- 1e5*tmp_n/pop
    
    
    TBnoti<-model$TBnoti
    
    
    
    # Preventive cascade
    sfin <- sln[nrow(sln),]
   
    
    #Arrays
    
    res$popu[ii,]<-pop
    res$inc_all[ii,]<-inc_all
    res$inc_hi[ii,]<-inc_hi
    res$inc_mdr[ii,]<-inc_mdr
    res$prev_all[ii,]<-prev
    res$prev_hi[ii,]<-prev_hi
    res$mort_tb[ii,]<-mort
    res$notif_all[ii,]<-not
    res$sfin[ii,]<-sfin
    
    res$llk[ii]<-model$llk
    res$prpe[ii]<-model$TBepi[6]
    res$pr_mdr[ii]<-model$TBepi[7]
    
  }
  res$burn<-burn
  res$thin<-thin
  filename<-paste("res/","output","_",location,"_",runtype,sep="")
  saveRDS(res,filename)
  
  return(res)
}


