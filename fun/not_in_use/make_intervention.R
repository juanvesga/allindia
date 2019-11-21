make_intervention<-function (p, r, i,i_all, s, gps,intvn_seq,M0){
  
  
  
  
  for (iseq in 1:length(intvn_seq)){
    intvn <- intvn_seq[iseq]
  }
  
  if (intvn == 0){
    
    Mset <-M0
    
  }else if(intvn == 1){
    # New adherence tools in FL (i.e. improved treatment completion)
    
    p$pdef[c(1,3)]<-1-0.95
    r$default <-  r$Tx*p$pdef/(1-p$pdef)
    Mset <- make_model(p, r, i,i_all, s, gps)
    
  }else if(intvn == 2){
    # New adherence tools in SL (i.e. improved treatment completion)
    
    p$pdef2[c(1, 3)]<-1-0.95
    r$default2 <- r$Tx2*p$pdef2/(1-p$pdef2)
    Mset <- make_model(p, r, i,i_all, s, gps)       
    
  }else if(intvn == 3){
    # Private sector engagement
    p$pse <- 0.95
    Mset <- make_model(p, r, i,i_all, s, gps)
    
  }else if(intvn == 4){
    # Active case-finding
    r$acf_sym[2] <- p$acfhi*p$acfsens[2]*p$acfconf_sens[2]*(1-p$acf_k)
    Mset <- make_model(p, r, i,i_all, s, gps)
    
    
  }else if(intvn == 5){
    
    # LTBI Treatment 
    
    r$ptrate <- p$ptreach[1]*p$hcscreen
    r$ptdefault <- r$outpt*p$ptdef/(1-p$ptdef)
    Mset <- make_model(p, r, i,i_all, s, gps)
    
  }
  M <- Mset
  
}

