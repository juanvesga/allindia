run_model <- function (ptpars,x,pt_in, prm, ref, sel, agg, gps, hivpoints, sfin){
  
  r <- prm$r; p <- prm$p ; i = ref$i; s <- ref$s; xi <- ref$xi
  
  tmp <- allocate_parameters(x,r,p,xi)
  
  r<-tmp$r ; p<-tmp$p
  
  
  # --- Set up the necessary models -----------------------------------------
  
  p1 <- p; r1<-r;
  p1$IPThiv <-ptpars[1]
  p1$IPT[,1]<-ptpars[2]
  p1$IPT[,2]<-ptpars[2]
  p1$IPT[,3]<-ptpars[2]
  #PT 
  tstart<-2017
  tmax<-2030
  p1$mono    <- 0                                            # is the regimen based on monotherapy
  r1$ptdur   <- c(pt_in[1],pt_in[1],pt_in[1])          # Duration
  p1$ptcure  <- pt_in[2]                                   # Steril 
  r1$pt_resist_acq<- 12/pt_in[1]*pt_in[3]/(1-pt_in[3]) #Rate of acquisition of resistance during PT
  p1$ptdef    <- 1-pt_in[4]   
  r1$ptdefault <- (12/r1$ptdur)*p1$ptdef/(1-p1$ptdef)        # estimated time lost of protection due to completion rates
  r1$outIPT    <- 12/r1$ptdur                                # 24 months protection Rangaka-Martins etc and 3HP
  
  M1 <- make_model(p1, r1, i, s, gps)
  
  # 2017
  p0 <- p; r0<-r;
  M0 <- make_model(p0, r0, i, s, gps)
  
  
  # --- Solve the models ----------------------------------------------------
  fx<-goveqs_basis
  
  times  = seq(tstart, tmax, by=1)           # time scale
  intT <-2
  t.interv   <- c(times[2], times[2]+intT)
  pars0<-list(agg=agg, sel=sel,s=s,p=p0, r=r0, i=i, hivpoints=hivpoints, M=M0)
  
  fx_scale<-function(t, state, M1) scale_up(t, state, M1,pars0,t.interv,fx)
  
  
  
  #Run the model
  out1 <- as.data.frame(lsoda(y = as.numeric(sfin), times = times, 
                            func = fx_scale, parms = M1))   #
  
  soln1<-out1[,2:ncol(out1)]
  t1   <-out1[,1]
  
  out<-list(soln=soln1 , t=t1)
  
  return(out)
  
}
