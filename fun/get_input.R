
get_input<-function(location){
  
  
  
  #__________________________________________________________________________
  #  Load data & Call distributions for lhd
  #__________________________________________________________________________
  
  # --Load data from file
  tmp      <-as.data.frame(read_excel("data/data.xlsx"))
  rownames(tmp)<-tmp[,1]  
  data_raw<-tmp[location,-1]

  # --Load data from file
  tmp      <-as.data.frame(read_excel("data/params.xlsx"))
  rownames(tmp)<-tmp[,1]  
  pars<-tmp[location,-1]
  
  
  #-- Load WHO TB Burden Data
  tb_burden<-getTBinR::get_tb_burden()
  
  data<-tb_burden[tb_burden$country==location,]
  data$e_inc_rr_100k<-1e5*(data$e_inc_rr_num/data$e_pop_num)
  data$e_inc_rr_100k_lo<-1e5*(data$e_inc_rr_num_lo/data$e_pop_num)
  data$e_inc_rr_100k_hi<-1e5*(data$e_inc_rr_num_hi/data$e_pop_num)
  
  
  nm<-c("year","e_inc_100k","e_inc_100k_lo","e_inc_100k_hi",
        "e_inc_tbhiv_100k","e_inc_tbhiv_100k_lo","e_inc_tbhiv_100k_hi",
        "e_inc_rr_100k","e_inc_rr_100k_lo","e_inc_rr_100k_hi",
        "e_mort_100k","e_mort_100k_lo","e_mort_100k_hi","c_newinc_100k")
  datapoints<-data[fields]
  datapoints$c_newinc_100k_lo<-NA
  datapoints$c_newinc_100k_hi<-NA
  

  # Create likelihood distributions for datapoints
  # Uncomment below for likelihood method
  # lhd<-Make_distr_fns(datapoints,1)
  
  lhd<-datapoints
  
  #__________________________________________________________________________
  #  Build model structure
  #__________________________________________________________________________
  
  # Model states
  gps<-list( sectors=c("pu","pr","pe"),
             strain =c("ds","mdr"),
             hiv    =c("hn","hu","ht"),
             dm     =c("dn","du","dt"))
  
  
  states1 = c('U')
  states2 = c('Lf','Ls','Pt','Ia','Is','E','Rlo','Rhi','R')
  states3 = c('Dx','Tx','Tx2')
  states4 = c('S')
  states5 = c('SDx','STx')
  
  
  
  i<-list()
  s<-list()
  groups<-list(states1,gps$hiv,gps$dm)
  ref_all <- get_addresses(groups, i, s, 0)
  groups<-list(states2,gps$hiv,gps$dm,gps$strain)
  ref_all<-  get_addresses(groups, ref_all$i, ref_all$s, ref_all$i$nstates)
  groups<-list(states3,gps$hiv,gps$dm,gps$strain, gps$sectors)
  ref_all<-  get_addresses(groups, ref_all$i, ref_all$s, ref_all$i$nstates)
  groups<-list(states4,gps$hiv,gps$dm)
  ref_all<-  get_addresses(groups, ref_all$i, ref_all$s, ref_all$i$nstates)
  groups<-list(states5,gps$hiv,gps$dm,gps$sectors)
  ref_all<-  get_addresses(groups, ref_all$i, ref_all$s, ref_all$i$nstates)
  
  #For indexing without NTS
  i<-list()
  s<-list()
  groups<-list(states1,gps$hiv,gps$dm)
  ref <- get_addresses(groups, i, s, 0)
  groups<-list(states2,gps$hiv,gps$dm,gps$strain)
  ref<-  get_addresses(groups, ref$i, ref$s, ref$i$nstates)
  groups<-list(states3,gps$hiv,gps$dm,gps$strain, gps$sectors)
  ref<-  get_addresses(groups, ref$i, ref$s, ref$i$nstates)
  
  
  # Include the auxiliaries
  
  auxnames <- c('inc', 'notif' , 'mort' ,'remo' ,'pt', 'txfp','dx','pmo','acf' )
  auxinds  <- c(  3      , 3    ,  1     ,3     ,1   ,  1    , 3  , 2   ,1)
  ref$i$aux<-list()
  lim_all<-ref_all$i$nstates
  for (ii in 1:length(auxnames)){
    inds <- lim_all + (1:auxinds[ii])
    ref$i$aux[[auxnames[ii]]] <- inds
    lim_all <- inds[length(inds)]
  }
  ref$i$nx <- lim_all
  
  
  #States
  ref$s$infectious <- c(ref$s$Ia, ref$s$Is, ref$s$E, ref$s$Dx)
  ref$s$prevalent  <- c(ref$s$infectious, ref$s$Tx , ref$s$Tx2)
  ref$s$TBmort     <- ref$s$infectious
  ref$s$sought     <- c(ref$s$Dx, ref$s$E )
  ref$s$treated_pr <- intersect(c(ref$s$Tx, ref$s$Tx2) ,ref$s$pr)
  ref$s$treated_pu <- intersect(c(ref$s$Tx, ref$s$Tx2) ,ref$s$pu)
  ref$s$treated_pe <- intersect(c(ref$s$Tx, ref$s$Tx2) ,ref$s$pe)
  ref$s$treated_pay <- intersect(c(ref$s$Tx, ref$s$Tx2) ,c(ref$s$pe,ref$s$pu))
  ref$s$postdx     <-c(ref$s$Dx, ref$s$Tx, ref$s$Tx2, ref$s$Rlo, ref$s$Rhi, ref$s$R)
  ref$s$tbcare     = c(ref$s$Dx, ref$s$Tx, ref$s$Tx2)
  ref$s$notbcare   = c(ref$s$Ia,ref$s$Is, ref$s$E, ref$s$Rhi, ref$s$Rlo, ref$s$R)
  ref$s$respSymp   <-c(ref_all$s$S, ref_all$s$SDx, ref_all$s$STx)
  ref$s$notbtreated_pay <-intersect(ref_all$s$STx ,c(ref_all$s$pe,ref_all$s$pu))
  ref$s$hivpositive <-intersect(ref$s$hu ,ref$s$ht)
  
  ref$s$nstates<-(1:ref$i$nstates)
  
  #__________________________________________________________________________
  #  Prepare selectors and aggregartors of model output
  #__________________________________________________________________________
  
  tmp <- matrix(1,ref_all$i$nstates,ref_all$i$nstates)
  tmp[ref_all$s$hn,ref_all$s$hu]<-0
  tmp[ref_all$s$hu,ref_all$s$hn]<-0
  tmp[ref_all$s$hn,ref_all$s$ht]<-0
  tmp[ref_all$s$ht,ref_all$s$hn]<-0
  tmp[ref_all$s$hu,ref_all$s$ht]<-0
  tmp[ref_all$s$ht,ref_all$s$hu]<-0
  tmp[ref_all$s$dn,ref_all$s$du]<-0
  tmp[ref_all$s$du,ref_all$s$dn]<-0
  tmp[ref_all$s$dn,ref_all$s$dt]<-0
  tmp[ref_all$s$dt,ref_all$s$dn]<-0
  tmp[ref_all$s$du,ref_all$s$dt]<-0
  tmp[ref_all$s$dt,ref_all$s$du]<-0
  
  check<-tmp - diag(diag(tmp))
  
  #  Incidence: 1.Total 2.TB/HIV 3. MDR
  tmp <- matrix(0,3,ref_all$i$nstates)
  tmp[1,ref$s$Ia] <- 1
  tmp[2,intersect(ref$s$Ia,ref$s$hivpositive) ] <- 1
  tmp[3,intersect(ref$s$Ia,ref$s$mdr)] <- 1
  agg <-list(inc= tmp)
  
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp[ref$s$Ia,] <- 1
  tmp<-tmp*check
  sel<-list(inc= tmp - diag(diag(tmp)))
  
  # DR acquisition selectors
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp[intersect(ref$s$Tx,ref$s$mdr),intersect(ref$s$Tx,ref$s$ds)] <- 1
  tmp<-tmp*check
  sel$acqu <- tmp - diag(diag(tmp))
  
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp [ref$s$Ia,ref$s$Ls] <- 1
  tmp <-tmp*check
  sel$remo <- tmp - diag(diag(tmp))
  
  # TB notification in RNTCP
  tmp <- matrix(0,3,ref_all$i$nstates)
  tmp[1,ref$s$treated_pay]   <- 1
  tmp[2,intersect(ref$s$mdr,ref$s$treated_pay)]  <- 1
  tmp[3,ref$s$treated_pe]  <- 1
  agg$notif <-tmp
  
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp[ref$s$treated_pay,c(ref$s$Ia, ref$s$Is, ref$s$E, ref$s$Dx)] <- 1
  tmp<-tmp*check
  sel$notif <- tmp - diag(diag(tmp))
  
  
  # Selectors mortality
  tmp <- matrix(0,1,ref_all$i$nstates)
  tmp[1,ref$s$infectious] <- 1
  agg$mort <-tmp
  
  # Preventive therapy
  tmp <- matrix(0,1,ref_all$i$nstates)
  tmp[1,intersect(ref$s$Pt,c(ref$s$pu, ref$s$pe))] <- 1
  agg$pt <-tmp
  
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp[ref$s$Pt,] <- 1
  tmp<-tmp*check
  sel$pt<-tmp - diag(diag(tmp))
  
  
  
  # Selectors for false positive treatment (Tx)
  tmp <- matrix(0,1,ref_all$i$nstates)
  tmp[1,ref$s$notbtreated_pay] <- 1
  agg$txfp <-tmp
  
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp[ref_all$s$STx,] <- 1
  tmp<-tmp*check
  sel$txfp<- tmp - diag(diag(tmp))
  
  #Aggregator for counting Dx : Will apply to matrices for smear, xpert and Xray
  tmp <- matrix(0,1,ref_all$i$nstates)
  tmp[1,c(ref$s$notbtreated_pay,ref$s$treated_pay)] <- 1
  agg$dx <-tmp
  
  
  #Selectors for Dx
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp[c(ref_all$s$STx,ref$s$Tx,ref$s$Tx2), ] <- 1
  tmp<-tmp*check
  sel$dx<- tmp - diag(diag(tmp))
  
  
  
  # Selector for ACF 
  tmp <- matrix(0,1,ref_all$i$nstates)
  tmp[1,c(ref_all$s$S, ref$s$Ia, ref$s$Is, ref$s$E)] <- 1
  agg$acf <-tmp
  
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp[c(ref$s$notbtreated_pay,ref$s$treated_pay), ]<- 1
  tmp<-tmp*check
  sel$acf <-tmp - diag(diag(tmp))
  
  
  # Aggregator for patient-month of tx
  tmp <- matrix(0,2,ref_all$i$nstates)
  tmp[1,intersect(c(ref$s$notbtreated_pay,ref$s$treated_pay),c(ref$s$Tx, ref_all$s$STx))] <- 1
  tmp[2,intersect(ref$s$treated_pay,ref$s$Tx2)] <- 1
  agg$pmo <-tmp
  
  #calib parameters
  xi <-list()
  
  xnames <-c('beta','beta_mdr','bRRhiv','bRRdm',"progRRHIV",'fast','kappa',
             'careseek_lo','careseek_hi','cs2','RRcs_rs','csRRhiv',
             'symp_del','imm','selfcure','muTB',
             'pse_base','Dxpu','Txinitpu','Dxhiv','IPThiv','ARTrec',
             'ntpcov','dmhzr','dmtx')
  
  xnums  <- rep(1,length(xnames))
  lim <- 0
  for (ii in 1:length(xnames)){
    inds <- lim + (1:xnums[ii])
    xi[[xnames[ii]]] <- inds
    lim <- inds[length(inds)]
  }
  
  
  xi$nx <- lim
  bds<-matrix(0,length(xnames),2)
  bds[xi$beta,]             <- c(2, 18)
  bds[xi$beta_mdr,]         <- c(2, 18)
  bds[xi$bRRhiv,]           <- c(0.75, 1)
  bds[xi$bRRdm,]            <- c(1, 6)
  bds[xi$fast,]             <- c(0.012, 4)
  bds[xi$progRRhiv,]        <- c(5, 24)
  bds[xi$kappa,]            <- c(0.2, 5)
  bds[xi$careseek_lo,]      <- c(0.5, 12)
  bds[xi$careseek_hi,]      <- c(0.5, 12)
  bds[xi$cs2,]              <- c(1, 5)
  bds[xi$RRcs_rs,]          <- c(0, 1)   #Factor reducing careseeking amongs resp symptomatic
  bds[xi$csRRhiv,]          <- c(1, 10) 
  bds[xi$symp_del,]         <- c(1.2, 12)
  bds[xi$imm,]              <- c(0.25, 0.75)
  bds[xi$selfcure,]         <- c(0.1, 0.22)
  bds[xi$muTB,]             <- c(0.1, 0.22)
  bds[xi$pse_base,]         <- c(0, 1)
  bds[xi$Dxpu,]             <- c(0.7, 1)
  bds[xi$Txinitpu,]         <- c(0.6, 0.9)
  bds[xi$Dxhiv,]            <- c(0.7, 1)
  bds[xi$IPThiv,]           <- c(0, 1)
  bds[xi$ARTrec,]           <- c(0, 10)
  bds[xi$ntpcov,]           <- c(0, 1)
  bds[xi$dmhzr,]            <- c(0, 1)
  bds[xi$dmtx,]            <- c(0, 1)
  
  
  
  #__________________________________________________________________________
  #  Prepare fixed parameters
  #__________________________________________________________________________
  
  
  r<-list()
  p<-list()
  
  #Epi
  
  r$fast_react    <- 0       # Ragonnet 2018_Epidemics
  r$slow_react    <- (0.0012 + 0.0012)/2  # Ragonnet 2018_Epidemics
  r$slow          <- 1.9710             # Ragonnet 2018_Epidemics
  r$careseek      <- c(0,0)
  r$careseeking2  <- c(0, 0, 0)
  r$access        <- 0
  p$kappa         <- 1
  r$selfcure  <- 0
  r$relapse   <- c(0.032, 0.14, 0.0015)
  r$mort_TB   <- 1/6
  p$Fast      <- c(0.035, 0.084, 0.087)
  p$imm       <- 0.5
  p$crossg    <- 0.3
  
  # Diagnosis stage
  p$strd_pr <-0
  if(strcmp(data_raw[location,'stdard_pr'],'good')){
    p$strd_pr<-1}
  else{
    p$strd_pr<-0.5
  }
  
  p$pu <- 0.5
  p$pse<- 0
  r$Dx <- 52
  p$Dx <- c(0.83,0.5,0.5)
  r$actC  <-1
  p$Tx_init <- c(0.88, 0.5, 0.5)
  p$smear_sens<-0.8 #Smear test speccificity (Swai F 2011, BMC resreach notes)
  p$smear_spec<-0.94#Smear test speccificity (Swai F 2011, BMC resreach notes)
  p$xpert_sens<-0.9 #Smear test speccificity (Swai F 2011, BMC resreach notes)
  p$xpert_spec<-0.99#Smear test speccificity (Swai F 2011, BMC resreach notes)
  p$xray_sens<-0.9 #Xray test speccificity (Swai F 2011, BMC resreach notes)
  p$xray_spec<-0.5 #Xray test speccificity (Swai F 2011, BMC resreach notes)
  
  p$smear<-c(0.93 , 0, 0)   #prob of microscopy by sector
  p$xray <-c(0,1,1)          #prob of Xray as Dx upfront
  p$xpert_upf<-c(0.07,0,0)  #prob of xpert upfront
  p$xpert_fup<-c(0.1,0,0)  #prob of xpert follow-up aas confirmaion of smear
  
  # Treatment stage
  r$Tx <- 2
  p$pdef <- (1-data_raw[location,'tsr_fl'])*c(1, 1.75, 1.75)
  r$default <- r$Tx*p$pdef/(1-p$pdef)
  p$cure      <- c(1, 1, 1)
  
  #MDR
  r$MDR_acqu <- 0.05 # MDR acquisition
  p$SL_trans <-c(0.8, 0, 0) # Transfer from FL
  r$Tx2      <-0.5
  p$newmol   <- 0 #Fraction recieving new molecuiles for their SL treatment
  p$pdef2 <- 0.15*c(1, 1.75, 1.75) 
  r$default2 <- r$Tx2*p$pdef2/(1-p$pdef2)
  p$cure2    <- c(0.5, 0, 0)
  
   # HIV cascades
  p$hivtest_notb = pars['hivtest_notb' ]#  care seeking 0.18pry as estimated by Olney et al 2016
  p$hivtest_tb   = pars['hivtest_tb']# WHO report % known status.. care seeking inside TB care system
  p$tbtest_hiv   = pars['tbtest_hiv']# Tb test among HIV pos
  p$art_notb     = pars['art_notb']# Coverage expected in 2015 * proportion estimated to adhreing and getting suppressed
  p$art_tb       = pars['art_tb']# Coverage expected in TB+ 2015 * proportion estimated to adhreing and getting suppressed
  r$art_dropout  = pars['art_dropout' ]#  Average UNGAS2016 + AMPATH/
  r$hivdecline= pars['hivdecline']# Annual Decline in HIV incidence after 2017
  r$reinfhiv=1 # Susceptibility to reinfection on HIV+
  r$OptimART <-1  
    
    # Preventive therapy
  # p$IPT<- matrix(0,3,3)
  # p$mono      <-0 # reg_profile{regimen,'mono'}   ;        % is the regimen based on monotherapy
  # p$ptdef     <-0.617# Alsdrurf 20161-reg_profile{regimen,'completion'} ;        % Churchyard et al 2014 NEJM
  # r$ptdur    = [reg_profile{regimen,'pt_dur'},...
  #               reg_profile{regimen,'pt_dur'},...
  #               reg_profile{regimen,'pt_durART'}];    % 6months
  # p$ptcure   = reg_profile{regimen,'sterile'};        % fraction LTBI getting cured by pt (see Nuermberger 2005, Zhang2009, Houben 2014)
  # r$pt_resist_acq= (12/r$ptdur(1)).*(1-0.95)./0.95; %Rate of acquisition of resistance during PT
  # r$pt_longeff= 0 ;          % Long lasting effect post PT  (months)
  # r$ptdefault =  (12./r$ptdur)*(1-p$ptdef)./p$ptdef;% estimated time lost of protection due to completion rates
  # r$outIPT    = 12./(r$ptdur+r$pt_longeff);      % 24 months protection Rangaka-Martins etc and 3HP
  # p$forg      =0.1;%forgiveness of regimen
  
  
  r$protection= 1/0.083#; % progression protection in Alaska studies (Comstock 1975)
  r$ARTred    = 0.35#;        % Reduction in progression given ART (Suthar2010PlosMed)
  p$d         = 0.63#;        % Reduction in progression given IPT (37% on top of ART reductions)
  p$Peff      = 0.63#        % Reduction to achieve with IPT;
  r$ARTrec    = 0#;           % Recruitment probability
  p$housedist =c(1, 1, 1)#; % rate control for household distr of TST

  
  # Demographic terms
  lex        <- pars['lex']    # Life expectancy
  r$mort     <- 1/lex #[(0.0253+	1.0000e-03)/2	(0.0065+	0.0933)/2];%[ 1/lex 1/(lex-15) 1/(lex-65)];    % Non-disease-related mortality
  p$growth   <- 0     #pars['growth'};
  p$frac_pop <- c(0.0910+0.1824	, 0.6643+0.0623) #	Population fraction <15 2016;
  # p$hi       <- data_raw[location,'size_hirisk'] # size of the vulnerable population
  p$resp_symptomatic<-0.045 # Fraction Respiratoty symnptomatic (From kenya survey (eligible only by symptoms)(TB survey)
  p$popN<-pars['pop']
  
  #intervs
  r$acf_asym<-c(0, 0)
  r$acf_sym <-c(0, 0)
  p$acfhi <-0.5
  p$xray_acf<-0
  p$xpert_acf<-1
  p$acf_k<-0 # Cumulative losses
  
  
  p$cfy_all<-0  #Switch parametyer to do contact trace
  p$cfy    <-0  #Switch parametyer to do contact trace
  r$nutri  <-0
  p$prevComm<-1 #TB prevalence in the community for CFY
  p$pteff  <-0.63
  p$ptreach<-0 # 0.1 for children, 1 for all population
  r$outpt <- 12/9 
  p$ptdef <- 1-0.75 
  r$ptdefault <- r$outpt*p$ptdef/(1-p$ptdef)
  r$ptrate<-0
  p$hcscreen<-1 
  
  
  # Diabetes
  p$dmtx<-0
  p$dmdef <- 1-0.85 
  r$dmdefault <- p$dmdef/(1-p$dmdef)
  
  
  # Unit costs (USD temporary costs )
  
  p$u_smear   <-10.87               #;cost.smear=rnu(u_smear,nruns);
  p$u_xpert   <-32.24               #cost.xpert=rnu(u_xpert,nruns);
  p$u_dst     <-20.03               #cost.dst=rnu(u_dst,nruns);
  p$u_xray    <-15               #cost.dst=rnu(u_dst,nruns);
  p$u_flmo    <- 33.91              #cost.flmoadu=rnu(u_flmoadu,nruns);
  p$u_slmo    <- 134.61             #cost.slmoadu=rnu(u_slmoadu,nruns);
  p$u_ipt     <-52.20               #cost.ipt=rnu(u_ipt,nruns);
  p$u_acf     <- ((30.69+1.93)) + 3 #cost.acf=rnu(u_acf,nruns);
  
  
  
  
  prm<-list(p=p, r=r, bds=t(bds),xnames=xnames )
  ref$i_all<-ref_all$i
  ref$s_all<-ref_all$s
  ref$xi<-xi
  
  return(list(prm=prm, ref=ref, lhd=lhd, sel=sel, agg=agg, gps=gps, datapoints=datapoints))
  
  
}
