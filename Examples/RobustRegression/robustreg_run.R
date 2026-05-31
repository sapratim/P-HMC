############################################################################
################## Robust compressed sensing ###############################
############################################################################

set.seed(209)
library(Rcpp)
library(mcmcse)
library(foreach)
library(doParallel)
source("robustreg_data.R")          # loads model settings as well
sourceCpp("pre_robustreg_functions.cpp")
load("marginal_vars.RData")

###############   Obtain MAP estimate
MAP <- map_estimate(B, y, alpha, nu, sigma, w_truth)

# start a little off the MAP to avoid zero-gradient issues 
w_start <- MAP + rnorm(length(MAP), 0, 0.01)

# post_var_diag from loading marginal variances
precond_diag <- post_var_diag

iter <- 1e5
lambda_prox <- .002
L_px        <- 10    # leapfrog steps for pHMC / MALA
L_guo       <- 10    # leapfrog steps for Guo-HMC

parallel::detectCores()
num_cores <- 10
doParallel::registerDoParallel(cores = num_cores)
reps <- 100

output_rreg <- foreach(b = 1:reps) %dopar% 
  {
    ## Run samplers
    print(b)
    eps_p    <- 0.045
    phmc_time <- system.time(phmc_run <- phmc_cpp(B, y,
                                                  lambda = lambda_prox, alpha = alpha, sigma = sigma,
                                                  iter   = iter, eps_hmc = eps_p, L = L_px, nu = nu,
                                                  start  = w_start, precond = precond_diag, blather = TRUE))[3]
    
    cat("\n--- Guo-HMC ---\n")
    eps_guo    <- 0.0008
    guohmc_time <- system.time(guohmc_run <- guohmc_cpp(B, y,
                                                        lambda = lambda_prox, alpha = alpha, sigma = sigma,
                                                        iter   = iter, eps_hmc = eps_guo, L = L_guo, nu = nu,
                                                        start  = w_start, precond = precond_diag, blather = TRUE))[3]
    
    cat("\n--- MALA (pHMC with L = 1, lambda = eps/2) ---\n")
    eps        <- 0.05
    mymala_time <- system.time(mymala_run <- phmc_cpp(B, y,
                                                      lambda = eps / 2, alpha = alpha, sigma = sigma,
                                                      iter   = iter, eps_hmc = eps, L = 1, nu = nu,
                                                      start  = w_start, precond = precond_diag, blather = TRUE))[3]
    
    cat("\n--- RWM ---\n")
    rwm_time <- system.time(rwm_run <- rwm_cpp(B, y,
                                               iter  = iter, h = 0.02,
                                               start = w_start, alpha = alpha, sigma = sigma, nu = nu,
                                               precond = precond_diag, blather = TRUE))[3]
  
    # Means
    all_means <- cbind(colMeans(rwm_run[[1]]), colMeans(phmc_run[[1]]), 
                       colMeans(mymala_run[[1]]),colMeans(guohmc_run[[1]]))
    colnames(all_means) <- c("RWM", "pHMC", "myMALA", "guoHMC")
    
    # ESS
    all_ess <- cbind(ess(rwm_run[[1]]),ess(phmc_run[[1]]), 
                     ess(mymala_run[[1]]), ess(guohmc_run[[1]]))
    colnames(all_ess) <- c("RWM", "pHMC", "myMALA", "guoHMC")
    
    all_time <- c(rwm_time, phmc_time, 
                  mymala_time, guohmc_time)
    names(all_time) <- c("RWM", "pHMC", "myMALA", "guoHMC")
    
    list(all_means, all_ess, all_time)
  }

save(output_rreg, file = "Output/outputrreg.Rdata")

