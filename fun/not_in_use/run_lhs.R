run_lhs <- function (fx,bds,n) {
  

  
  
  # Create grid
  pars<-nrow(bds)
  
  grid <- matrix(
    (1:(pars*n))*0,
    nrow = n
  )
  
  for (ii in 1:n){
    grid[ii,]<- runif(pars, min=bds[,1], max=bds[,2])
  }
  
  
  # run 
  
  
  
  numCores <- detectCores()
  cl <- makeCluster(numCores[1])
  doParallel::registerDoParallel(cl)
  finalLLK <- foreach(i=1:n, .combine=cbind, .packages = c("deSolve","splines","data.table"), .export = ls(globalenv())) %dopar% {
    tempLLK = fx(grid[i,]) #calling a function
    #do other things if you want

    tempLLK #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
  }
  #stop cluster
  stopCluster(cl)
  
  # finalLLK<-rep(0,n)
  # for (ii in 1:n){
  #   x<-grid[ii,]
  #   finalLLK[ii]<-fx(x)
  # 
  # }
  
  id_best<-which.min(finalLLK)
  best_set<-grid[id_best,]
  
  out<-list(id_best=id_best, best_set=best_set, grid=grid, llk=finalLLK)
  return(out)
  
}