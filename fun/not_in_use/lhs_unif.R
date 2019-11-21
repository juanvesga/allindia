lhs_unif <- function (bds,n) {
  
  pars<-nrow(bds)
  
  grid <- matrix(
    (1:(pars*n))*0,
    nrow = n
  )
  
  for (ii in 1:n){
    grid[ii,]<- runif(pars, min=bds[,1], max=bds[,2])
  }
  
  return(grid)
}