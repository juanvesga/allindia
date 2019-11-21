plot_targets<-function (runs,datapoints,location){
  
  
  x<-seq(2018-size(runs$popu,2),2017,by=1)
  labels<-names(datapoints) 
  fields<-c('prev_all','inc_all','prev_hi','inc_hi','mort_tb','notif_all')
  titles<-c('Prevalence','Incidence','Prevalence (slum)','Incidence (slum)','Mortality','Notified',
            'Notifications (2017)','NNS risk group')
  
  
  
  if (size(runs$popu,1)>1){
    runtype<-'mcmc'
  }else{
    runtype<-'mle'
  }
  
  fs<-11

  red <- c(216, 23, 37)/ 255
  grey=c(0.9,0.9,0.9)
  lw<-1
  xlimits<-c(2013, 2018)
  
  
  
  titletxt <- paste('Model fits,',location,runtype,sep=" ")
  
  plots<-vector("list",length(fields))
  for (ii in 1:length(fields)){
    
    df <- data.frame(year=x,
                     y=as.numeric(runs[[fields[ii]]]))
    
    fn<-function(x)strcmp(x,fields[ii])
    iddat<-which(lapply(labels,fn) ==TRUE)
    
    p<- ggplot(df, aes(x=year, y=y)) + 
      geom_line(color='red', size=lw) +
      ylim(0,max(df$y)*1.3)
    
    if (length(iddat)>0){ 
      
      df2<-data.frame(x=2017,y=as.numeric(datapoints[iddat]),sd=as.numeric(datapoints[iddat]*0.2))
      p <- p +
        geom_point(data=df2,aes(x=df2$x,y=y), size=3, color="blue")+
        geom_errorbar(data=df2, mapping=aes(x=df2$x,ymin=y-sd, ymax=y+sd),width=0.2, size=0.8, color="blue") 
    }
    # Finished line plot
    p<-p+labs(title=titles[ii], x="Year", y = "Rate per 100K")+
      theme_classic() 
    plots[[ii]]<-p
    
  }
  filename<-paste("plots/","fits","_",location,"_","mle",".jpeg",sep="")

  jpeg(file=filename)
  
  grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],nrow=3)
  
  dev.off()
  
  }

