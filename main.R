

rm(list = ls())
# All India model   : created by Juan F vesga 6/11/2019

#install.packages("Rcpp","RcppArmadillo","inline","rbenchmark")
#install.packages("RcppArmadillo")
#__________________________________________________________________________
#
#  LOAD LIBRARIES AND SOURCE FUNCTIONS
#__________________________________________________________________________
library(Matrix)
library(readxl)
library(deSolve)
library(pracma)
library(neldermead)
library(profvis)
library(Rcpp)
library(ggplot2)
library(gridExtra)
library(abind)
library(getTBinR)

source("fun/get_addresses.R")
source("fun/get_distribution_fns.R")
source("fun/get_objective.R")
source("fun/goveqs_basis.R")
source("fun/scale_up.R")
source("fun/allocate_parameters.R")
source("fun/make_model.R")
source("fun/return_output.R")
source("fun/plot_targets.R")
source("fun/get_input.R")
source("fun/sim_pack_intervention.R")
source("fun/plot_sim_intv.R")

sourceCpp("fun/compute_dx.cpp")
sourceCpp("fun/scale_matrix.cpp")

#__________________________________________________________________________
#  Select setting & action
#  Action options: "simplex" , "plot_fits","sim_intv"
#__________________________________________________________________________
location <- "India"
action  <- "sim_intv"

nint <- 5     #Number of interventions for pathway method (ITV_MLE)

#__________________________________________________________________________
#  Select Intervention parameters
#__________________________________________________________________________

itv <- list()
#---------- 1)  Improve first-line treatment outcomes
# Increase completion rates in FL public sector

itv$pdef <- 1 - 0.95

#---------- 2) Improve second-line treatment outcomes
# Increase completion rates  in SL public

itv$pdef2 <- 1 - 0.75

#---------- 3) Universal DST in routine, public TB services
# Increase probability of being diagnosed with xpert up front
# and getting a DST (xpert) test done

itv$xpert_uf <- 0.5
itv$xpert_fu <- 0.95



#---------- 4) Private sector engagement
# Engae a fraction of the private providers
# and improve probability of Dx
itv$pse <- 0.95
itv$pDx_pe <- 0.8
itv$pDx_pe <- 0.8
itv$xpert_fu_pe <- 0.8

# ----------5)Preventive therapy amongst household contacts

itv$ptreach <- 1   # 0.1 for childen, 1 for al population
itv$hcscreen <- 4   #number of household contacts to e screened
itv$ptdef <- 1 - 0.6 #Expected completion rates


# ----------5)Active case finding (ACF)

itv$acfhi    <- 0.2 # Prportion of slum population reached
itv$xray_acf <- 0   # ACF screening by Xray =1 ; Verbally =0;
itv$xpert_acf <- 1   # Confirmation by Xpert =1 ; smear =0;
itv$acf_loss <- 0.5 # Operational losses on expected yield



input <- get_input(location)
prm <- input$prm
ref <- input$ref

prm$itv <- itv

sel <- input$sel
agg <- input$agg

gps <- input$gps
datapoints <- input$datapoints
lsq <- input$lhd
#__________________________________________________________________________
#  Create instances of model functions
#__________________________________________________________________________

obj       <-  function(x)
  get_objective(x, prm, ref, sel, agg, gps, lsq, FALSE)
obj_spx   <-  function(x)
  - get_objective(x, prm, ref, sel, agg, gps, lsq)
obj_mcmc  <-  function(x)
  get_objective(x, prm, ref, sel, agg, gps, lsq)


#Random draw of parameters
xu = rep(0, length(input$prm$xnames))
for (u in 1:length(input$prm$xnames)) {
  a <- input$prm$bds[1, u]
  b <- input$prm$bds[2, u]
  xu[u] <- (b - a) * runif(1) + a
}

#__________________________________________________________________________
#  Execute procedure according to initial selection
#__________________________________________________________________________


if (strcmp(action, "simplex")) {
  #### Simplex: Finds a single best set of parameters to fit data (MLE)
  f <- paste("res/", "bestset", "_", location, "_", "mle", sep = "")
  if (file.exists(f)) {
    x = readRDS(f)
  } else{
    x = xu
  }
  fminsearch <- neldermead::fminsearch
  opt <- optimset(TolX = 1.e-2,  Display = "iter")
  x1 <- fminsearch(fun = obj_spx,                   x0 = x,
                   options = opt)
  x  <- x1$optbase$xopt
  #
  filename <- paste("res/", "bestset", "_", location, "_", "mle", 
                    sep ="")
  saveRDS(x, filename)
  runs <- return_output(obj, x, ref, location, c())
  
}


if (strcmp(action, "plot_mle")) {
  ######### PLots the results of best fit (model vs data)
  
  f <- paste("res/", "output", "_", location, "_", "mle", sep = "")
  if (file.exists(f)) {
    runs <- readRDS(f)
    plot_targets(runs, datapoints, location)
    
  } else{
    message('No best fits runs saved')
    
  }
}

if (strcmp(action, "sim_intv")) {
  ######### Simulates and plots intervention on best fitted model
  # To explore all result numbers see object "itvs"
  
  f <- paste("res/", "output", "_", location, "_", "mle", sep = "")
  ff <- paste("res/", "bestset", "_", location, "_", "mle", sep = "")
  if (file.exists(f)) {
    runs <- readRDS(f)
    x1 <- readRDS(ff)
    itvs <-
      sim_pack_intervention(x1, runs$sfin, prm, ref, sel, agg, gps, location, c())
    plot_sim_intv(itvs)
    
  } else{
    message('No best fit runs saved')
    
  }
}





#
