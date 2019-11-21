plot_sim_intv<-function (itvs){
  
  
  x<-seq(2026-ncol(itvs$inc),2025,by=1)
  fields<-c('inc','mort')
  titles<-c('Incidence','Mortality','Incremental Cost')
  
  
  
  if (size(itvs$inc,1)>1){
    runtype<-'mcmc'
  }else{
    runtype<-'mle'
  }
  
  fs<-11
  
  red <- c(216, 23, 37)/ 255
  grey=c(0.9,0.9,0.9)
  lw<-1
  xlimits<-c(2018, 2026)
  
  
  

  plots<-vector("list",length(fields))
  for (ii in 1:length(fields)){
    
    df <- data.frame(year=rep(x,2),
                     y=as.numeric(itvs[[fields[ii]]]),
                     projection=c( rep('Baseline',length(x)), rep('Intervention',length(x))))
    
    p<- ggplot(df, aes(x=year, y=y, col =projection)) + 
      geom_line( size=lw) +
      ylim(0,max(df$y)*1.1)
    
    
    # Finished line plot
    p<-p+labs(title=titles[ii], x="Year", y = "Rate per 100K")+
      theme_classic() 
    plots[[ii]]<-p
    
  }


  
  df <- data.frame( Cost=c(itvs$icr_all, itvs$ic_fl, itvs$ic_sl,
                        itvs$ic_sm,itvs$ic_xp,itvs$ic_xr),
                   Item=c("Total","FL","SL","smear","xpert","xray"))
  
  p<-ggplot(data=df, aes(x=Item, y=Cost)) +
    geom_bar(stat = "identity",fill="steelblue")
  p<-p+labs(title="Incremental Cost", x="Item", y = "Cost (USD)")
  plots[[3]]<-p
  
  
  
  df <- data.frame( y=c(itvs$inc_av_pr, itvs$mort_av_pr),
                    x=c("Incidence","Mortality"))
  
  p<-ggplot(data=df, aes(x=x, y=y*100)) +
    geom_bar(stat = "identity",fill="steelblue")
  p<-p+labs(title="Relative cumulative reductions 2017-2025", x="Indicator", y = "Reduction(%)")
  plots[[4]]<-p
  
  
  
  
  grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],nrow=2)

  
}

