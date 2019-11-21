# %  Created by juan fernando vesga on 09/02/2018.
# %  Copyright © 2018 juan fernando vesga. All rights reserved.

get_foi_hiv <- function (t,r,points) {
  
  
  hiv <- 0
  yrstart <-1980 
  hivIR <- 0
  if (t >= yrstart){
    hiv <- 1
    it <- round(t) # integer of t to used as index in pre-def vectors.
    ii <- (it - yrstart)+1;
    if (t > 2017){
      ii <- (2016 - yrstart)+1;
    }
    iii <- ((it - yrstart) + 1)+1;
    if (t > 2017){
      iii <- ((2016 - yrstart) + 1)+1;
    }
    x0 <- it;
    x1 <- it + 1;
    
    # Built in
    hivIR <- (points[ii]*(x1-t) + points[iii]*(t-x0))/(x1-x0)*r$turnoffHIV
    # Faster performance
    #hivIR <- (hivInc[ii]*(x1-tt) + hivInc[iii]*(tt-x0))/(x1-x0)
    
    # -- When year is > 2017 use decline as stated 
    if (t > 2017){
      incdecline <- (1 + (t - 2017)) * 0.04
      hivIR <- hivIR * (1 - incdecline);
    }
    
  }
  FOI<-as.numeric(hivIR*hiv)
  
  return(FOI)
  
}