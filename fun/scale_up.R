scale_up<- function (t, state, M1,pars0,t.interv,fx) {
  
  
  scale <- min((t-t.interv[1])/(t.interv[2]-t.interv[1]),1); 
  
  if (scale<0) 
    {
    scale<-0
    }
  
  Mt <- M1
  M0 <- pars0$M
  Mt$lin <- scale_matrix(M0$lin,M1$lin,scale)

  pars_scaled<-pars0
  pars_scaled$M<-Mt
     
  return(fx(t, state, pars_scaled))
   
}