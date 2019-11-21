make_model<-function(p, r, i, i_all, s, gps){
  
  # --- Get the linear rates ------------------------------------------------
  # m <- rsparsematrix(i_all$nstates,i_all$nstates,0)
  # ms <- rsparsematrix(i_all$nstates,i_all$nstates,0)
  # mx <- rsparsematrix(i_all$nstates,i_all$nstates,0)
  # mac <- rsparsematrix(i_all$nstates,i_all$nstates,0)
  # mall <- rsparsematrix(i_all$nstates,i_all$nstates,0)
  # 
  m <- matrix(0,i_all$nstates,i_all$nstates)
  ms <- matrix(0,i_all$nstates,i_all$nstates)
  mx <- matrix(0,i_all$nstates,i_all$nstates)
  mxr <- matrix(0,i_all$nstates,i_all$nstates)
  mac <- matrix(0,i_all$nstates,i_all$nstates)

  
  sub2ind<-function(size,roi,coi)(coi-1)*size + roi
  
  
  
  for (ir in 1:length(gps$risk)){
    
    risk  <- gps$risk[ir]
    geti  <- function(st) i_all[[st]][[risk]]
    S     <- geti('S')
    SDxs  <- geti('SDx')
    STxs  <- geti('STx')
    
    #---------------ACF in resp symptomatic in public system
    
    acf_rate<- r$acf_sym[ir]
    
    acf_alg_xray <- p$xray_acf * (1-p$xray_spec) * 
      (  p$xpert_acf *(1-p$xpert_spec) +
           (1-p$xpert_acf)*(1-p$smear_spec))
    
    acf_alg_noxray <- (1-p$xray_acf)  * 
      (  p$xpert_acf *(1-p$xpert_spec) +
           (1-p$xpert_acf)*(1-p$smear_spec))
    
    
    source <- S; destin <- STxs$pu
    rate <- p$acf_k * acf_rate * (acf_alg_xray + acf_alg_noxray)
    m[destin,source] <- m[destin,source] + rate
    
    
    #-------------Count Dx
    acf_xprt_no <- acf_rate * (p$xray_acf * (1-p$xray_spec) * p$xpert_acf +
                            (1-p$xray_acf)* p$xpert_acf)
    
    acf_smr_no  <- acf_rate * (p$xray_acf * (1-p$xray_spec) * (1-p$xpert_acf) +
                            (1-p$xray_acf)* (1-p$xpert_acf))
    
    acf_xray_no <- acf_rate * p$xray_acf
    
    
    acf_no <- acf_rate
    
    
    
    source <- S; destin <- STxs$pu
    rate   <- acf_xprt_no 
    mx[destin,source] <- mx[destin,source] + rate
    
    rate <-  acf_smr_no 
    ms[destin,source] <- ms[destin,source] + rate
    
    rate <-  acf_xray_no 
    mxr[destin,source] <- mxr[destin,source] + rate
    
    rate <-  acf_no
    mac[destin,source] <- mac[destin,source] + rate
    
    #---------------Primary careseeking, in resp symptomatic in public system
    source <- S; destin <- SDxs$pu; rate <- r$careseek[ir]*p$pu
    m[destin, source] <- m[destin, source] + rate
    
    source <- S; destin <- SDxs$pr; rate <- r$careseek[ir]*(1-p$pu)*(1-p$pse)
    m[destin, source] <- m[destin, source] + rate
    
    source <- S; destin <- SDxs$pe; rate <- r$careseek[ir]*(1-p$pu)*p$pse
    m[destin, source] <- m[destin, source] + rate
    
    
    for (ip in 1:length(gps$sectors)){
      
      prov <- gps$sectors[ip]
      SDx <- SDxs[[prov]];  STx <- STxs[[prov]]
      
     
      ispu<-strcmp(prov, 'pu')
      
     
      
      
      
      #---------------Treatment of false positives
      
      DxTxrate <-r$Dx * p$Dx[ip] *p$Tx_init[ip]
      
      Dx_alg <-        p$xpert_upf[ip]  * (1- p$xpert_spec) + 
        ((1-ispu) * (1-p$xpert_upf[ip]) * (1- p$xray_spec)) +
        ((ispu)   * (1-p$xpert_upf[ip]) * (1- p$smear_spec))
      
      NoDx_alg <-      p$xpert_upf[ip]  * p$xpert_spec + 
        ((1-ispu) * (1-p$xpert_upf[ip]) * p$xray_spec) +
        ((ispu)   * (1-p$xpert_upf[ip]) * p$smear_spec)
      
      
      source <- SDx; destin<-STx
      rate<- DxTxrate * Dx_alg 
      m[destin,source] <- m[destin,source] + rate
      
      
      source <- SDx; destin<-S
      rate   <- DxTxrate * NoDx_alg 
      m[destin,source] <- m[destin,source] + rate
      
      source <- STx; destin<-S; rate <- r$Tx
      m[destin,source] <- m[destin,source] + rate
      
      
      # Count units in public/ pse
      if (ip !=2){
      
      pas_xprt_no <- r$Dx * p$Dx[ip] * p$xpert_upf[ip]
      pas_smr_no  <- r$Dx * p$Dx[ip] * (ispu) * (1-p$xpert_upf[ip])
      
      
      source <- SDx; destin<-STx
      rate<-  pas_xprt_no
      mx[destin,source] <- mx[destin,source] + rate
      
      rate <- pas_smr_no
      ms[destin,source] <- ms[destin,source] + rate
      }
      
    }
  }
  
  
  
  #---- TB cascade
  
  for (ir in 1:length(gps$risk)){
    
    
    for (istr in 1:length(gps$strain)){
      
      risk   <- gps$risk[ir]
      strain <- gps$strain[istr]
      
      gi <- function(st) i[[st]][[risk]][[strain]]
      Lf    <- gi('Lf')
      Ls    <- gi('Ls')
      Pt    <- gi('Pt')
      Ia    <- gi('Ia')
      Is    <- gi('Is')
      Dxs   <- gi('Dx')
      Txs   <- gi('Tx')
      Tx2s  <- gi('Tx2')
      E     <- gi('E')
      Rlo   <- gi('Rlo')
      Rhi   <- gi('Rhi')
      R     <- gi('R')
      
      ismdr<-strcmp(strain, 'mdr')
      
      
      #--- Reactivation
      
      source <- Lf; destin <- Ls; rate <- r$slow
      m[destin, source] <- m[destin, source] + rate
      
      source <- Lf; destin <- Ia; rate <- r$fast_react
      m[destin, source] <- m[destin, source] + rate
      
      source <- Ls; destin <- Ia; rate <- r$slow_react
      m[destin, source] <- m[destin, source] + rate
      
      
      #----- Preventive therapy
      source <- Lf; destin <- Pt; rate <- r$ptrate*(1-ismdr)
      m[destin, source] <- m[destin, source] + rate
      
      source <- Ls; destin <- Pt; rate <- r$ptrate*(1-ismdr)
      m[destin, source] <- m[destin, source] + rate
      
      source <- Pt; destin <- Ls; rate <- r$outpt
      m[destin, source] <- m[destin, source] + rate
      
      source <- Pt; destin <- Ls; rate <- r$ptdefault
      m[destin, source] <- m[destin, source] + rate
      
      
      #--- Symptom development
      source <- Ia; destin <- Is; rate <- r$symp_del
      m[destin, source] <- m[destin, source] + rate
      
     
      #Systematic screening in the community (Asymptomatic)
      #Note: all transitions in mac, mx, ms matrices are for
      #counting events for the costing model and have no effect in
      #the dynamics
      
      
      
      acf_rate_asy<- r$acf_asym[ir]  
      
      acf_rate_sy <- r$acf_sym[ir]
      
      
      acf_alg_xray <-p$xray_acf * p$xray_sens * 
                (  p$xpert_acf *  p$xpert_sens  +
                (1-p$xpert_acf)*  p$smear_sens)
      
      acf_alg_noxray <- (1-p$xray_acf)  * 
                (  p$xpert_acf  * p$xpert_sens  +
                (1-p$xpert_acf) * p$smear_sens)
      
      
      source <- Ia; destin <- Txs$pu
      rate <- p$acf_k   * acf_rate_asy * (1-ismdr) *(acf_alg_xray)
      m[destin, source] <- m[destin, source] + rate
      
      source <- Ia; destin <- Tx2s$pu
      rate  <-  p$acf_k   * acf_rate_asy * ismdr * 
                p$xray_acf * p$xray_sens * p$xpert_acf *  p$xpert_sens
      
      m[destin, source] <- m[destin, source] + rate
      
      
      
      
      #-------------Count Dx
      acf_xprt_no <- acf_rate_asy * (p$xray_acf * p$xray_sens * p$xpert_acf  +
                                  (1-p$xray_acf)* p$xpert_acf)
      
      acf_smr_no  <- acf_rate_asy * (p$xray_acf * p$xray_sens *  (1-p$xpert_acf) +
                              (1-p$xray_acf)* (1-p$xpert_acf))
     
      acf_no <- acf_rate_asy
      
      source <- Ia; destin <- Txs$pu
      rate <- acf_no*(1-ismdr)
      mac[destin, source] <- mac[destin, source] + rate
      
      rate <- (1-ismdr)* acf_xprt_no
      mx[destin, source] <- mx[destin, source] + rate
      
      rate <- (1-ismdr)* acf_smr_no
      ms[destin, source] <- ms[destin, source] + rate
      
      
      source <- Ia; destin <- Tx2s$pu
      rate <- acf_no*(ismdr)
      mac[destin, source] <- mac[destin, source] + rate
      
      rate <- (ismdr)* acf_xprt_no
      mx[destin, source] <- mx[destin, source] + rate
      

      
      
      #Systematic screening in the community (only symptomatic)
      source <- Is; destin <- Txs$pu
      rate <- p$acf_k   * acf_rate_sy * (1-ismdr) *(acf_alg_xray + acf_alg_noxray)
      m[destin, source] <- m[destin, source] + rate
      
      source <- Is; destin <- Tx2s$pu
      rate  <-  p$acf_k   * acf_rate_sy * ismdr * 
        (p$xray_acf * p$xray_sens * p$xpert_acf *  p$xpert_sens+
      (1-p$xray_acf)* p$xpert_acf *  p$xpert_sens)
      m[destin, source] <- m[destin, source] + rate
      
      
      source <- E; destin <- Txs$pu
      rate <- p$acf_k   * acf_rate_sy * (1-ismdr) *(acf_alg_xray + acf_alg_noxray)
      m[destin, source] <- m[destin, source] + rate
      
      source <- E; destin <- Tx2s$pu
      rate  <-  p$acf_k   * acf_rate_sy * ismdr * 
        (p$xray_acf * p$xray_sens * p$xpert_acf *  p$xpert_sens+
           (1-p$xray_acf)* p$xpert_acf *  p$xpert_sens)
      m[destin, source] <- m[destin, source] + rate
      
      
      #-------------Count Dx
      acf_xprt_no <- acf_rate_sy * (p$xray_acf * p$xray_sens * p$xpert_acf  +
                                       (1-p$xray_acf)* p$xpert_acf)
      
      acf_smr_no  <- acf_rate_sy * (p$xray_acf * p$xray_sens *  (1-p$xpert_acf) +
                                       (1-p$xray_acf)* (1-p$xpert_acf))
      
      
      acf_no <- acf_rate_asy
      
      source <- Is; destin <- Txs$pu
      rate <- acf_no*(1-ismdr)
      mac[destin, source] <- mac[destin, source] + rate
      
      rate <- (1-ismdr)* acf_xprt_no
      mx[destin, source] <- mx[destin, source] + rate
      
      rate <- (1-ismdr)* acf_smr_no
      ms[destin, source] <- ms[destin, source] + rate
      
      source <- E; destin <- Txs$pu
      rate <- acf_no*(1-ismdr)
      mac[destin, source] <- mac[destin, source] + rate
      
      rate <- (1-ismdr)* acf_xprt_no
      mx[destin, source] <- mx[destin, source] + rate
      
      rate <- (1-ismdr)* acf_smr_no
      ms[destin, source] <- ms[destin, source] + rate
      
            
      source <- Is; destin <- Tx2s$pu
      rate <- acf_no*(ismdr)
      mac[destin, source] <- mac[destin, source] + rate
      
      rate <- (ismdr)* acf_xprt_no
      mx[destin, source] <- mx[destin, source] + rate
      
      source <- E; destin <- Tx2s$pu
      rate <- acf_no*(ismdr)
      mac[destin, source] <- mac[destin, source] + rate
      
      rate <- (ismdr)* acf_xprt_no
      mx[destin, source] <- mx[destin, source] + rate
      
      
      #--- Primary careseeking, including access to public sector care
      source <- Is; destin <- Dxs$pu; rate <- r$careseek[ir]*p$pu
      m[destin, source] <- m[destin, source] + rate
      
      source <- Is; destin <- Dxs$pr; rate <- r$careseek[ir]*(1-p$pu)*(1-p$pse)
      m[destin, source] <- m[destin, source] + rate
      
      source <- Is; destin <- Dxs$pe; rate <- r$careseek[ir]*(1-p$pu)*p$pse
      m[destin, source] <- m[destin, source] + rate
      
      
      for (ip in 1:length(gps$sectors)){
        prov <- gps$sectors[ip]
        Dx <- Dxs[[prov]];  Tx <- Txs[[prov]]; Tx2 <- Tx2s[[prov]]
        
        DxTxrate <- p$Dx[ip] *p$Tx_init[ip]
        
        Dx_alg_fl <- DxTxrate * ((1-ismdr) * p$xpert_upf[ip]  * p$xpert_sens) + 
          ((1-ispu) * (1-p$xpert_upf[ip]) * p$xray_sens) +
          ((ispu)   * (1-p$xpert_upf[ip]) * p$smear_sens)
        
        Dx_alg_sl <- DxTxrate * (ismdr  * p$xpert_upf[ip]  * p$xpert_sens) + 
          (ismdr  * (1-ispu) * (1-p$xpert_upf[ip]) * p$xray_sens * p$xpert_upf[ip] * p$xpert_sens) +
          (ismdr  *   (ispu) * (1-p$xpert_upf[ip]) * p$smear_sens* p$xpert_upf[ip] * p$xpert_sens)
        
        
           
        p_ltfu  <- 1-(Dx_alg_fl+Dx_alg_sl)
        
        #--- Diagnosis
        source  <- Dx
        destins <-      c(Tx,        Tx2,        E)
        rates   <- r$Dx*c(Dx_alg_fl, Dx_alg_sl,  p_ltfu)
        m[destins, source] <- m[destins, source] + t(rates)
        
        #--- FL Treatment
        
        pFLcure <- p$cure[ip]
        rMDRacq <- r$MDR_acqu
        source  <- Tx
        destins <- c(Rlo,      E,             Rhi,             i[["Tx"]][[risk]][["mdr"]][[prov]])
        rates   <- c(r$Tx*pFLcure,  r$Tx*(1-pFLcure),  r$default[ip],   rMDRacq)
        m[destins, source] <- m[destins, source] + t(rates)
        
        if  (ismdr==1){
          pSLtran <- p$SL_trans[ip]
          
          source  <- Tx
          destins <- c(Tx2,             E)
          rates   <- c(r$Tx*pSLtran,   r$Tx*(1-pSLtran) + r$default[ip])
          m[destins, source] <- m[destins, source] + t(rates)
        }
        
        
        #--- SL Treatment
        source  <- Tx2
        destins <- c(Rlo,                 E)
        rates   <- c(r$Tx2*p$cure2[ip],  r$Tx2*(1-p$cure2[ip]) +  r$default2[ip])
        m[destins, source] <- m[destins, source] + t(rates)
        
        
        #--- Count Dx events by provider
        if (ip !=2){
          
          pas_xprt_flno <- (1-ismdr) * r$Dx * p$Dx[ip] * p$xpert_upf[ip]
          pas_xprt_slno<-    (ismdr  * r$Dx * p$Dx[ip] * p$xpert_upf[ip]) + 
            + (ismdr  *   (ispu) * (1-p$xpert_upf[ip]) * p$smear_sens* p$xpert_upf[ip]) 
          pas_smr_no    <-  r$Dx * p$Dx[ip] * (ispu) * (1-p$xpert_upf[ip])
          
          
          source  <- Dx; destins <- Tx
          rates   <- pas_xprt_flno
          mx[destins,source] <- mx[destins,source] + rates
          
          source  <- Dx; destins <- Tx2
          rates   <- pas_xprt_slno
          mx[destins,source] <- mx[destins,source] + t(rates)
          
          source  <- Dx; destins <- Tx
          rates   <- pas_smr_no
          ms[destins,source] <- ms[destins,source] + rates
        }
       
        
        
        
      }
      
      #--- Secondary careseeking
      source <- E; destin <- Dxs$pu; rate <- r$careseeking2[ir]*p$pu
      m[destin, source] <- m[destin, source] + rate
      
      source <- E; destin <- Dxs$pr; rate <- r$careseeking2[ir]*(1-p$pu)*(1-p$pse)
      m[destin, source] <- m[destin, source] + rate
      
      source <- E; destin <- Dxs$pe; rate <- r$careseeking2[ir]*(1-p$pu)*p$pse
      m[destin, source] <- m[destin, source] + rate
      
      #--- Relapse
      sources <- c(Rlo, Rhi, R)
      destin  <- Ia
      rates   <- r$relapse
      m[destin, sources] <- m[destin, sources] + rates
      
      sources <- c(Rlo, Rhi)
      destin  <- R
      rates   <- 0.5
      m[destin, sources] <- m[destin, sources] + rates
      
      
      # --- Self cure : To Rhi so Rlo reflects exclusively treatment
      sources <- c(Ia, Is,E,Dx)
      destin <- Rhi
      rate   <- r$selfcure
      m[destin, sources] <- m[destin, sources] + rate
    }
  }
  
  
  
  M<-list()
  M$lin <- m - diag(colSums(m))
  #-- Count Dx Matrices
  M$xp<-mx - diag(colSums(mx))
  M$sm<-ms - diag(colSums(ms))
  M$xr<-ms - diag(colSums(mxr))
  M$mac<-mac -diag(colSums(mac))

  # --- Get nonlinear rates of TB transmission
  M$nlin<-list()
  
  
  for (ir in 1:length(gps$risk)){
    risk <- gps$risk[ir]
    
    M$nlin[[risk]]<-list()
    for (istr in 1:length(gps$strain)){
      strain <- gps$strain[istr]
      
      m <- matrix(0,i_all$nstates,i_all$nstates)
      
      getso <-function(st) c(i[[st]][[risk]][['ds']] , i[[st]][[risk]][['mdr']])
      getdes <-function(st) i[[st]][[risk]][[strain]]
      
      U   <- i[['U']][[risk]]
      Lf  <- getdes('Lf')
      Ls  <- getso('Ls')
      Rlo <- getso('Rlo')
      Rhi <- getso('Rhi')
      R   <- getso('R')
      
      m[Lf, c(U,  Ls, Rlo, Rhi, R)] <- 1
      
      
      #immunity
      m[,c(s$Ls, s$Rlo,s$Rhi ,s$R)] <-   m[,c( s$Ls, s$Rlo,s$Rhi, s$R)]*p$imm           
      M$nlin[[risk]][[strain]] <- m-diag(colSums(m))
    }
  }
  
  
  # --- Getting force-of-infection for urban and rural settings
  tmp <-matrix(0,4,i_all$nstates)    # 1.DS/Urban, 2.MDR/Urban, 3.DS/Rural, 4.MDR/Rural
  int2 <-function(s1,s2) intersect(intersect(s[[s1]],s[[s2]]),s$infectious)
  tmp[1,int2('ds','lo')]  <- r$beta
  tmp[2,int2('mdr','lo')] <- r$beta_mdr
  tmp[3,int2('ds','hi')]  <- r$beta*r$bRRhi
  tmp[4,int2('mdr','hi')] <- r$beta_mdr*r$bRRhi
  
  # Now include the cross-geogr links, shifting perspective from infected to susceptible
  m <-matrix(0,4,i_all$nstates)    # 1.DS/Urban, 2.MDR/Urban, 3.DS/Rural, 4.MDR/Rural
  m[1,] <- tmp[1,] + p$crossg*tmp[3,]
  m[2,] <- tmp[2,] + p$crossg*tmp[4,]
  m[3,] <- tmp[3,] + p$crossg*tmp[1,]
  m[4,] <- tmp[4,] + p$crossg*tmp[2,]
  
  
  m[,setdiff(s$infectious,c(s$Ia,s$Is) )] <- m[,setdiff(s$infectious, c(s$Ia,s$Is))]*p$kappa
  M$lambda <- m
  
  # -- - Get the mortality rates
  #HIV negative
  m <- rep(r$mort,i_all$nstates)
  m[c(s$infectious,intersect(s$Tx,s$mdr))]<-r$mort_TB
  M$mortvec <- t(m)
  return(M)
}