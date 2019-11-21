get_objective <- function (x, prm, ref, sel, agg, gps, calfn, LLK=TRUE){
  
  r <- prm$r; p <- prm$p 
  i <- ref$i    ; s <- ref$s
  i_all <- ref$i_all; s_all <- ref$s_all 
  xi <- ref$xi
  
  tmp <- allocate_parameters(x,r,p,xi)
  
  r<-tmp$r ; p<-tmp$p
  
  
  
  
  # Check if parameters are within bounds
  tmp <- rbind(prm$bds[,1:length(x)], t(x))
  tmp <- diff(tmp[c(1,3,2),])
  # 
  
  if ((min(tmp)<0)==TRUE){
    #   
    llk<- -Inf
    return(llk)
    #   
  }
  
  #__________ Set up the necessary models 
  
  # PSE
  p3 <- p; r3<-r
  p3$pse<-p$pse_base
  M3    <- make_model(p3, r3, i,i_all, s, gps)

  # Full public
  p2 <- p; r2<-r
  M2 <- make_model(p2, r2, i,i_all, s, gps)
  
  
  # NTP
  p1 <- p; r1<-r
  p1$xpert  <- p$xpert*0
  p1$pu     <- p$pu*p$ntpcov
  M1 <- make_model(p1, r1, i, i_all, s, gps)
  
  
  # Equil
  p0 <- p; r0<-r
  p0$pu       <- 0
  r0$beta_mdr <- 0
  r0$MDR_acqu <- 0
  p0$xpert    <- p$xpert*0
  p0$growth   <- 0
  M0 <- make_model(p0, r0, i, i_all, s, gps)
  
  
  # Initial state
  init <- rep(0,i$nx)
  seed <- 1e-6
  init[i$U$hi]       <- p$hi*(1 - seed)
  init[i$U$lo]       <- (1-p$hi)*(1 - seed)
  init[i$Ia$hi$ds]   <- p$hi*seed
  init[i$Ia$lo$ds]   <- (1-p$hi)*seed
  # Set NTS population
  init[i_all$S$hi] <- p$resp_symptomatic*p$hi
  init[i_all$S$lo] <- p$resp_symptomatic*(1-p$hi)
  
  
  
  # --- Solve the models ----------------------------------------------------
  # Equilibrium Model
  times  <- seq(0, 200, by=1)          # time scale
  parms0<-list(agg=agg, sel=sel,s=s,p=p0, r=r0, i=i, i_all=i_all, M=M0)
  fx<-goveqs_basis
  out0 <- as.data.frame(lsoda(y = as.numeric(init), times = times, 
                              func = fx, parms = parms0))   
  soln0<-out0[,2:ncol(out0)]
  t0   <-out0[,1]
  
  
  
  #NTP roll out
  init<-out0[nrow(out0),2:ncol(out0)]
  times  <- seq(1970, 1997, by=1)           # time scale
  intT   <- 5
  t.interv   <- c(times[2], times[2]+intT)
  parms1<-list(agg=agg, sel=sel,s=s,p=p1, r=r1, i=i, i_all=i_all, M=M1)
  fx_scale<-function(t, state, M1) scale_up(t, state, M1,parms0,t.interv,fx)
  
  out1 <- as.data.frame(lsoda(y = as.numeric(init), times = times, 
                              func = fx_scale, parms = M1))   #
  
  soln1<-out1[,2:ncol(out1)]
  t1   <-out1[,1]
  
  #Full RNCTP
  init<-out1[nrow(out1),2:ncol(out1)]
  times  <- seq(1997, 2012, by=1)           # time scale
  intT   <- 5
  t.interv   <- c(times[2], times[2]+intT)
  parms2<-list(agg=agg, sel=sel,s=s,p=p2, r=r2, i=i, i_all=i_all, M=M2)
  fx_scale<-function(t, state, M2) scale_up(t, state, M2,parms1,t.interv,fx)
  
  out2 <- as.data.frame(lsoda(y = as.numeric(init), times = times, 
                              func = fx_scale, parms = M2))   #
  
  soln2<-out2[,2:ncol(out2)]
  t2   <-out2[,1]
  
  
  #PSE
  init<-out2[nrow(out2),2:ncol(out2)]
  times  <- seq(2012, 2017, by=1)           # time scale
  intT   <- 4
  t.interv   <- c(times[2], times[2]+intT)
  parms3<-list(agg=agg, sel=sel,s=s,p=p3, r=r3, i=i, i_all=i_all, M=M3)
  fx_scale<-function(t, state, M3) scale_up(t, state, M3,parms2,t.interv,fx)
  
  out3 <- as.data.frame(lsoda(y = as.numeric(init), times = times, 
                              func = fx_scale, parms = M3))   #
  
  soln3<-out3[,2:ncol(out3)]
  t3   <-out3[,1]
  
  
  soln = soln3;
  t=t3;
  
  # % --- Get the objectives --------------------------------------------------
  sfin<-soln[nrow(soln),]
  gs   <-function(cols) sum(sfin[cols])
  
  # Model output
  N           <-  gs(s$nstates) 
  Nhi         <-  gs(s$hi)
  
  #Prevalence
  prev_all<-(gs(s$prevalent)/N)*1e5
  prev_hi <-(gs(intersect(s$prevalent,s$hi))/Nhi)*1e5
  
  
  #Incidence
  tmp_n  <- apply(soln[,i$aux$inc],2,diff)
  tmp    <- 1e5*( tmp_n[nrow(tmp_n),]/c(N,Nhi,N))
  inc_all<- tmp[1];
  inc_hi <- tmp[2];

  
  #Notif
  tmp  <- apply(soln[,i$aux$notif],2,diff)
  tbnotif_n<-tmp[,1]
  tbnotif_n_risk <- tmp[,ncol(tmp)]
  n_notif  <- tmp[nrow(tmp),1]       
  notif    <-1e5*(n_notif/N )
  pr_notif_pe <- tmp[nrow(tmp),3]/tmp[nrow(tmp),1]
  pr_notif_mdr<- tmp[nrow(tmp),2]/tmp[nrow(tmp),1]
  
  
  #Mort
  tmp  <- diff(soln[,i$aux$mort])
  mort  <- 1e5*(tmp[nrow(tmp)]/N)
  
  #Model output to fit to data
  model<-c(prev_all, inc_all,prev_hi, inc_hi,notif,pr_notif_pe,pr_notif_mdr)
  
  llksum<-  sum( (lsq - model)^2 )
  
  # Uncomment for Likelihood method 
  # llk<-c(1:length(model))*0
  # for (l in 1:length(model)){
  #   
  #   llk[l]<-as.numeric(lhd[[l]](model[l]))
  #   if(is.nan( llk[l])){ llk[l]<-0}
  #   
  # }
  # 
  # 
  # llksum<-sum(llk)
  
  if (is.na(llksum)){
    llksum<- -Inf
  }
  
  
  
  
  if (LLK){ 
    
    return(as.numeric(llksum))
    
  }else{
    
    return( list(llk=as.numeric(llksum), soln0=soln0, t0=t0, soln=soln , t=t ,
                 TBepi=model))
  }
  
  
  
}
