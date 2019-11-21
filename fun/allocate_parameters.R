allocate_parameters<-function(x,r,p,xi){
  
  r$beta          <-  x[xi$beta]
  r$beta_mdr      <-  x[xi$beta_mdr]
  r$bRRhi         <-  x[xi$bRRhi]
  r$fast_react    <-  x[xi$fast]
  p$kappa         <-  x[xi$kappa]
  r$careseek[1]   <-  x[xi$careseek_lo]
  r$careseek[2]   <-  x[xi$careseek_hi]
  r$careseeking2  <-  r$careseek*x[xi$cs2]
  r$RRcs_rs       <-  x[xi$RRcs_rs]
  r$symp_del      <-  x[xi$symp_del]
  r$imm           <-  x[xi$imm]
  r$selfcure      <-  x[xi$selfcure]
  r$mort_TB       <-  x[xi$muTB]
  p$pu            <-  x[xi$pu]
  p$pse_base      <-  x[xi$pse_base]
  p$Dx[1]         <-  x[xi$Dxpu]
  p$Tx_init[1]    <-  x[xi$Txinitpu]
  p$ntpcov        <-  x[xi$ntpcov]
  p$crossg        <-  x[xi$cross]
  
  out<-list(p=p,r=r)
}