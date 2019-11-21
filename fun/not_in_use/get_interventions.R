get_interventions <- function (x,sfin,path,prm,ref,sel,agg,gps,location,varargin){
  
  r <- prm$r 
  p <- prm$p
  i <- ref$i 
  s <- ref$s
  i_all = ref$i_all 
  s_all <- ref$s_all
  xi <- ref$xi
  
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
  
  itv_labs<-c('+Improved TSR in FL','+Improved TSR in SL','+PSE','+ACF','+Preventive Therapy')
  
  
  res<-vector("list",9)
  res[[1]]<-array(0,c(nsims,(2025-2017),length(path)+1))
  res[[2]]<-array(0,c(nsims,(2025-2017),length(path)+1))
  
  res[[3]]<-matrix(0,nsims,length(path))
  res[[4]]<-matrix(0,nsims,length(path))
  res[[5]]<-matrix(0,nsims,length(path))
  
  res[[6]]<-matrix(0,nsims,length(path)+1)
  res[[7]]<-matrix(0,nsims,length(path)+1)
  
  res[[8]]<-array(0,c(nsims,4,length(path)+1))
  res[[9]]<-array(0,c(nsims,4,length(path)+1))
  
  names(res)<-c('inc' , 'mort','cas_prev', 'reduc_inc' , 'reduc_mrt','tsr_fl',
                'tsr_sl','cas_mdr', 'cas_ds')
  
  for (ii in 1:nsims){
    
    x0<-xvals[,ii]
    tmp <- allocate_parameters(x0,r,p,xi)
    r<-tmp$r ; p<-tmp$p
    
    
    
    
    # Baseline
    p0 <- p; r0<-r;
    M0 <- make_model(p0, r0, i, i_all, s, gps)
    Mset<-list() 
    Mset[[1]]<-M0
    intvn_seq<-0
    for (it in 2:(length(path)+1)){
      intvn_seq[it]<-path[it-1]
      Mset[[it]]<-make_intervention(p, r, i,i_all, s, gps,intvn_seq,M0)
    }
    
    #run models
    fx<-goveqs_basis
    init<-unlist(sfin)
    times  <- seq(2017, 2025, by=1)           # time scale
    intT   <- 3
    t.interv   <- c(times[2], times[2]+intT)
    allsol<-c() #array(0,c((2025-2017)+1, i$nx, length(path)+1))
    
    
    parms0<-list(agg=agg, sel=sel,s=s,p=p0, r=r0, i=i, i_all=i_all, M=M0)
    
    for (mi in 1:length(Mset)){
      Mfinal<-Mset[[mi]]
      fx_scale<-function(t, state, Mfinal) scale_up(t, state, Mfinal,parms0,t.interv,fx)
      out <- as.data.frame(lsoda(y = as.numeric(init), times = times, 
                                 func = fx_scale, parms = Mfinal))   #
      soln<-out[,2:ncol(out)]
      t   <-out[,1]
      
      if(mi==1){
        allsol<-soln
      }else{
        allsol<-  abind(x=allsol,y=soln,along=3)
        
      }
    }
    
    
    # Get Model outputs
    tmp<-allsol[,1:i$nstates,]
    tmp <-apply(tmp,c(1,3),sum)
    pop <- (tmp[1:nrow(tmp)-1,]+tmp[2:nrow(tmp),])/2
    
    Itmp1<-diff(allsol[,i$aux$inc[1],])
    tmp<-1-Itmp1[nrow(Itmp1),]/Itmp1[1,]
    reduc_inctmp=tmp[2:length(tmp)]
    inc<-Itmp1
    
    
    Itmp1<-diff(allsol[,i$aux$mort,])
    tmp<-1-Itmp1[nrow(Itmp1),]/Itmp1[1,]
    reduc_mrttmp=tmp[2:length(tmp)]
    mort<-Itmp1
    
    res$inc[ii,,]<-inc
    res$mort[ii,,]<-mort
    res$reduc_inc[ii,]<-reduc_inctmp 
    res$reduc_mrt[ii,]<-reduc_mrttmp 
  }
  
  return(res)
}

