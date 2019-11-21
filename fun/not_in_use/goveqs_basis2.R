
goveqs_basis <- function (t, state, parameters) {
  
  
  
  with(as.list(c(state,parameters)),             
       {
         
         
         N      <- sum(state[1:Nco])# Total population
         morts  <-sum(c(state[1:Nco])*mu)  + sum(state[co[c("A","Ad")]]*mu_aids) + sum(state[co["At"]]*mu_aidstx)
         births <- morts * (birthrate==0) + birthrate*N;   # Births for stable population
         
         
         dx<-state*0
         
         for (ia in 1:length(gps$age)){
           age <- gps$age[ia]
           
           for (ig in 1:length(gps$hiv)){
             hiv <- gps$hiv[ig]
             
             gu <- function(st) i[[st]][[age]][[hiv]]
             U     <- gu('U');
             Pu    <- gu('Pu');
             
             for (istr in 1:length(gps$strain)){
               
               strain <- gps$strain[istr]
               
               gi <- function(st) i[[st]][[age]][[strain]][[hiv]]
               
               
               Lf    <- gi('Lf')
               Ls    <- gi('Ls')
               Pi    <- gi('Pi')
               Pl    <- gi('Pl')
               Ps    <- gi('Ps')
               Ia    <- gi('Ia')
               Is    <- gi('Is')
               Dx    <- gi('Dx')
               Tx    <- gi('Tx')
               Tx2   <- gi('Tx2')
               E     <- gi('E')
               Rlo   <- gi('Rlo')
               Rhi   <- gi('Rhi')
               R     <- gi('R');
               
               # Group Selectors
               ishiv <- ig > 1
               ishivcas <- ig > 2
               ishivart <- ig > 2
               ismdr   <- strcmp(strain, 'mdr')
               
               RRh <- max(1, r$progRRhiv*ishiv)*max(1-ishivart, r$ARTred)*r$RRprogage[ia] # Relative reduction for HIV
               
               
               # Susceptible  
               if (ia==1){
               dx[U] <- births + state[Pu]*r$outIPT[ig] + state[Pu]*r$ptdefault[ig] -
                 state[U] * (foi+(p$IPT[ia,ig]*p$housedist[1])+r$mort[ia]) 
               } else{
               dx[U] <- state[Pu]*r$outIPT[ig] + state[Pu]*r$ptdefault[ig] -
                 state[U] * (foi+(p$IPT[ia,ig]*p$housedist[1])+r$mort[ia]) 
               }
               
               dx[Pu]<- state[U] * p$IPT[ia,ig]*p$housedist[1] -
                 state[Pu]*(r$outIPT[ig] + r$ptdefault[ig] + +r$mort[ia])
               
               
               dx[Lf]<- state[U] * foi - 
                 state[Lf] * (p$IPT[ia,ig]*p$housedist[2]+r$slow[ia]+r$fast_react[ia]*RRh+r$mort[ia])
               
               
               dx[Ls]<- state[Lf] * r$slow[ia] - 
                 state[Ls] * (p$IPT[ia,ig]*p$housedist[3]+r.slow_react(ia)*RRh +r$mort[ia])
               
               dx[Pi]<- state[Lf]*p$IPT[ia,ig]*p$housedist[2]*(1-p$mono*ismdr) + 
                 state[Ls]*p$IPT[ia,ig]*p$housedist[3]*(1-p$mono*ismdr) -
                 state[Pi]*(r$outIPT[ig]+r$pt_resist_acq*(ismdr!=1)+
                              r$ptdefault[ig]+r$mort[ia])
               
               
               dx[Ps]<- state[Pi]*r$outIPT[ig]*p$ptcure -
                 state[Ps]*(foi+r$mort[ia])
               
               
               dx[Pl]<- state[Pi]*(r$outIPT[ig]*(1-p$ptcure) +r$ptdefault[ig])+
                 state[Lf]*p$IPT[ia,ig]*p$housedist[2]*p$mono*ismdr +
                 state[Ls]*p$IPT[ia,ig]*p$housedist[3]*p$mono*ismdr -
                 state[Pl]*(r$slow_react[ia]*RRh + foi + r$mort[ia])
               
               
               dx[Ia]<- state[c(Ls,Pl)]*r$slow_react[ia]*RRh+
                 state[Lf]*r$fast_react[ia]*RRh+
                 state[c(Rlo,Rhi,R)]r$relapse -
                 state[Ia]*(r$symp_del+ r$mort[ia]+r$muTB)
               
               csr<-r$careseek*max(1, r$csRRhiv*ishivcas)
               
               dx[Is]<- state[Ia]*r$symp_del -
                 state[Is]*(csr+ r$mort[ia]+r$muTB)
               
               
               #Dx Rates
               DxAttempt <- p$Dx[ig]*p$Dxage[ia]
               
               pFLinit<- DxAttempt*p$Tx_init*( ((1-ismdr)*p$xpert*p$xpert_sens) + ((1-p$xpert)*p$smear_sens) )
               
               pSLinit<- DxAttempt*p$Tx_init*ismdr*p$xpert*p$xpert_sens
               
               p_ltfu  <- 1-(DxAttempt*p$Tx_init* ((p$xpert*p$xpert_sens) + (1-p$xpert)*p$smear_sens))
               
               
               dx[Dx]<- state[Is]*r$careseek*max(1, r$csRRhiv*ishivcas)+
                 state[E]*(ishivcas*csr  + (1-ishivcas)*csr*r$cs2)-
                 state[Dx]*(r$Dx+ r$mort[ia]+r$muTB)
               
               dx[E]<- state[Dx]*p_ltfu*r$Dx + 
                 state[Tx]*r$Tx*(1-pFLcure)*(ismdr!=1) +
                 r$Tx*(1-pSLtran)*ismdr + 
                 state[Tx2]* (r$Tx2*(1-p$cure2)+  r$default2) -
                 state[E]*(ishivcas*csr  + (1-ishivcas)*csr*r$cs2) -
                 state[E]*(r$mort[ia]+r$muTB)
               
               
               dx[Tx]<- state[Dx]*r$Dx*pFLinit-
                 state[Tx]*(r$Tx + r$default+ r$mort[ia]+r$mort_TBtx+r$default*rMDRacq*(ismdr!=1) )
               
               dx[Tx2]<- state[Dx]*pSLinit + state[Tx]*r$Tx*pSLtran*ismdr-
                 state[Tx2]*( r$Tx2 +  r$default2 + r$mort[ia]+r$mort_TBtx)
               
               
               dx[Rlo]<- state[Tx]*r$Tx*pFLcure*(ismdr!=1) +
                 state[Tx2]*r$Tx2*p$cure2 + 
                 state[c(Ia, Is, E, Dx, i$Tx[[age]][[mdr]][[hiv]]) ] * r$selfcure -
                 state[Rlo]*(foi+ r$relapse+0.5+r$mort[ia])
               
               dx[Rhi]<- state[Tx]* r$default*(1-rMDRacq)*(ismdr!=1)-
                 state[Rhi]*(r$relapse+0.5+r$mort[ia])
               
               dx[R]<-state[c(Rlo, Rhi)]*0.5 -
                 state[R]*(r$relapse+r$mort[ia])
             }
           }
         }
             
         
         
         list(dx)
         
         
       }
  )
}