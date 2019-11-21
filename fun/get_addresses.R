get_addresses<- function (groups, i, s, lim){
  
  
  # Initiate any sets not so far covered by s
  if  (!length(s)){
    fnames = {}
  } else {
    fnames <- names(s)
  }
  
  # for (ig in  1:length(groups)){
  #   gp <- groups[[ig]]
  #   for (ig2 in 1:length(gp)){
  #     if ( !(gp[ig2] %in% names(fnames))){  
  #       s[[gp[ig2]]]<-{}
  #     }
  #   }
  # }
  
  
  if (length(groups) == 1){
    gp1 <- groups[[1]]
    for (ig1 in 1:length(gp1)){
      lim <- lim+1
      i[[gp1[ig1]]] <- lim
      s[[gp1[ig1]]] <- c(s[[gp1[ig1]]], lim)
      
    }
  }
  
  
  if (length(groups) == 2){
    gp1 = groups[[1]] 
    gp2 = groups[[2]]
    for (ig1 in 1:length(gp1)){
      i[[gp1[ig1]]]<-list()
      for (ig2 in 1:length(gp2)){
        lim <- lim+1
        i[[gp1[ig1]]][[gp2[ig2]]]<-lim
        s[[gp1[ig1]]] <- c(s[[gp1[ig1]]], lim)
        s[[gp2[ig2]]] <- c(s[[gp2[ig2]]], lim)
        
      }
    }
  }
  
  if (length(groups) == 3){
    gp1 = groups[[1]] 
    gp2 = groups[[2]]
    gp3 = groups[[3]]
    for (ig1 in 1:length(gp1)){
      i[[gp1[ig1]]]<-list()
      for (ig2 in 1:length(gp2)){
        i[[gp1[ig1]]][[gp2[ig2]]]<-list()
        for (ig3 in 1:length(gp3)){
          lim <- lim+1
          i[[gp1[ig1]]][[gp2[ig2]]][[gp3[ig3]]]<-lim
          
          s[[gp1[ig1]]] <- c(s[[gp1[ig1]]], lim)
          s[[gp2[ig2]]] <- c(s[[gp2[ig2]]], lim)
          s[[gp3[ig3]]] <- c(s[[gp3[ig3]]], lim)
          
        }
      }
    }
  }
  
  if (length(groups) == 4){
    gp1 = groups[[1]] 
    gp2 = groups[[2]]
    gp3 = groups[[3]]
    gp4 = groups[[4]]
    for (ig1 in 1:length(gp1)){
      i[[gp1[ig1]]]<-list()
      for (ig2 in 1:length(gp2)){
        i[[gp1[ig1]]][[gp2[ig2]]]<-list()
        for (ig3 in 1:length(gp3)){
          i[[gp1[ig1]]][[gp2[ig2]]][[gp3[ig3]]]<-list()
          for (ig4 in 1:length(gp4)){
            lim <- lim+1
            i[[gp1[ig1]]][[gp2[ig2]]][[gp3[ig3]]][[gp4[ig4]]]<-lim
            
            s[[gp1[ig1]]] <- c(s[[gp1[ig1]]], lim)
            s[[gp2[ig2]]] <- c(s[[gp2[ig2]]], lim)
            s[[gp3[ig3]]] <- c(s[[gp3[ig3]]], lim)
            s[[gp4[ig4]]] <- c(s[[gp4[ig4]]], lim)
            
          }
        }
      }
    }
  }
  
 

  i$nstates<-lim
  return(list(s=s, i=i))
}
