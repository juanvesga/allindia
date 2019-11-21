# % Adaptive MCMC, using Haario et al:
#   % https://link.springer.com/article/10.1007/s11222-008-9110-y
# % and
# % http://probability.ca/jeff/ftpdir/adaptex.pdf
# 
# % F:          Function giving log-posterior density for a parameter set x
# % x0:         Initial value of parameter set x
# % n:          Number of iterations
# % cov0:       Initial covariance matrix
# % fac:        Scaling factor for covariance matrix. Set fac = 1 for default
# % displ:      Structure with display options. Set displ = true to show progress

MCMC_adaptive<-function(fn, x0, n, sigma, displ=TRUE){
  
  d = length(x0); b = 0.05; sd = sigma*2.4^2/d;
  inds <- c(); vals <- c();
  
  # Checks on the initial covariance matrix
  cov0 <- eye(d); cov0[inds,] <- 0; cov0[,inds] <- 0;
  
  # Initiate the output matrices
  xsto <- matrix(0,d,n); outsto <- c(1:n)*0;
  history <- matrix(0,d+1,n)# zeros(d+1,n);                                                  
  

  xsto[,1] <- x0; xbar = xsto;
  FX <- fn(x0); outsto[1] <- FX;
  acc <- 0;
  
  #Cholesky like decomposition 
  cholcov<-function(T){
    E<-eigen(T, symmetric = TRUE)
    V <- E$vectors
    D <- sqrt(E$values)  ## root eigen values
    rankM <-length(D)-sum(is.na(D))
    V1 <- V[, 1:rankM]
    D1 <- D[1:rankM]
    R <- D1 * t(V1)
    return(R)
  }
  isPosDef <- function(M) { if ( all(M == t(M) ) ) {
    # first test  symmetric-ity
    if (  all(eigen(M)$values > 0) ) {TRUE}
    else {FALSE} } #
    else {FALSE}  # not symmetric

  }
  
  
  # --- Start the MCMC loop -------------------------------------------------
  for (t in 2:n){
    
    X = xsto[ ,t-1]
    
    # --- Make a proposal from the distribution
    Y0 <- mvrnorm(1,X,0.1^2*cov0*sigma/d);
    if (t < 3*d){
      #Y = max(Y0,0); 
      Y <- Y0;
      Y[inds] <- vals;
    }else{
      
      ind0 <- 1; ind1 <- t-1
      covmat <- cov(t(xsto[ ,ind0:ind1]))
      covmat[inds,] <- 0; covmat[,inds] = 0;
      #covmat[1:blockind,(blockind+1:end)) = 0; covmat((blockind+1):end,1:blockind) = 0;
      # Need to modify this bit so it goes recursively - faster
      
      # Precaution to make sure values don't get negative
      # Y = max((1-b)*mvnrnd(X,sd*covmat) + b*Y0,0);
      
      ##Here
       # T<- cholcov(sd*covmat)
       if (isposdef(sd*covmat)==FALSE){
           print("stats:mvnrnd:BadCovariance2DSymPos")
       }
      Y <- (1-b)*mvrnorm(1,X,sd*covmat) + b*Y0;
      
      Y[inds] = vals;
    }
    history[1:d,t] <- Y
    
    # --- Decide whether to accept or not
    FY = fn(Y)
    
    if (runif(1) < exp(FY-FX)){
      # Accept
      xsel<- Y
      FX  <- FY
      acc <- acc+1
      history[nrow(history),t] = 1
    }else{
      # Reject
      xsel <- xsto[,t-1]
    }
    xsto[,t] <- xsel
    outsto[t] <- FX
    xbar[,t] <- (xbar[ ,t-1]*(t-1) + xsel)/t
    
    # % Display options
    if (displ==TRUE){
      # lines(t(outsto[1:t-1]),xlim=c(0, n))
           # lines(t(outsto[1:t-1]),xlim=c(0, n))
      print(c(t,FX))
    }
    # if displ && (mod(t,round(n/25))==0); fprintf('%0.5g ', t/n*25); end
    # if displ && (mod(t,200)==0)
    #     plot(outsto(1:t-1)'); xlim([0 n]); drawnow; 
    #          %         save mcmc_progress;
    #          end
  }
  results<-list(accept_rate = acc/n, xsto = t(xsto), history=t(history), outsto=t(outsto))  
  
  return(results)
}
