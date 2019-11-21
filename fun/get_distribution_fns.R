get_distribution_fns<-function(prctiles, distribution,W)
{
  dat  <- sort(prctiles)
  opts = optimset(TolX =1e-12,TolFun=1e-12, MaxIter=1e6, MaxFunEvals=1e6)
  #Estimators to help setting up initial guesses
  mn <- dat[2]; var <- (dat[3]-dat[1])^2/4
  
  #::::::::::::::::::: poiss
  
  if (strcmp(distribution, 'poisson')){
    logfn <-function(x) sum(dpois(x = dat[2],lambda = x,log = TRUE))  # ??
    
  #::::::::::::::::::: Log Normal
  
 }else if (strcmp(distribution, 'lognormal')){
    ff  <-function (x) c(plnorm(dat[1],x[1],x[2]), plnorm(dat[2],x[1],x[2]), plnorm(dat[3],x[1],x[2]))
    #Set up initial guess
    sig <- sqrt(log(var/mn^2+1)); mu  <- log(mn)-sig^2/2; init = c(mu, sig)
    #Do the optimisation
    obj <-function(x) sum((ff(x)/(c(2.5, 50, 97.5)/100) - 1)^2);
    res <- neldermead::fminsearch(obj,init,opts)
    out<-res$simplexopt$x[1,]
    val<-res$simplexopt$fv[1]
    if (val > 1e-2){
      stop('Calibration setup not converged')
    }
    #Get the log-pdf
    mu <- out[1]; sig <- out[2]
    logfn <-function(x) (-(log(x)-mu)^2/(2*sig^2) - log(x*sig*sqrt(2*pi)))/W
    
  #::::::::::::::::::: Beta 
  }else if (strcmp(distribution, 'beta')){
    ff <-function (x) c(pbeta(dat[1],x[1],x[2]), pbeta(dat[2],x[1],x[2]), pbeta(dat[3],x[1],x[2]))
    #Set up initial guess
    tmp <- mn*(1-mn)/var-1; a <- tmp*mn; b <- tmp*(1-mn); init <- c(a, b);
    #Do the optimisation
    obj<-function(x) sum((ff(x)/(c(2.5, 50, 97.5)/100) - 1)^2)
    res <- neldermead::fminsearch(obj,init,opts)
    out<-res$simplexopt$x[1,]
    val<-res$simplexopt$fv[1]
    if (val > 1e-2){
      stop('Calibration setup not converged');
    }
    #Get the log-pdf
    a <- out[1]; b <- out[2]
    logfn <-function(x) ((a-1)*log(x) + (b-1)*log(1-x) - lbeta(a,b))/W
  }
  #::::::::::::::::::: Beta 
  else if (strcmp(distribution, 'uniform')){
    ff <-function (x) c(punif(dat[1],x[1],x[2]), punif(dat[2],x[1],x[2]), punif(dat[3],x[1],x[2]))
    #Set up initial guess
    tmp = mn*(1-mn)/var-1; a = tmp*mn; b = tmp*(1-mn); init = c(a, b)
    #Do the optimisation
    obj<-function(x) sum((ff(x)/(c(2.5, 50, 97.5)/100) - 1)^2)
    res <- neldermead::fminsearch(obj,init,opts)
    out<-res$simplexopt$x[1,]
    val<-res$simplexopt$fv[1]
    if (val > 1e-2){
      stop('Calibration setup not converged');
    }
    #Get the log-pdf
    a <- out[1]; b <- out[2]
    logfn <-function(x)  (1/(x*(log(b)-log(a))))/W
    
  }
  
  # aux.sim = ff(out);                                                         #Simulated values of CDF at given percentile points
  # aux.val = val;                                                             #The final objective function
  
  return(logfn)
  
}