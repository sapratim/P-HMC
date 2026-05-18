############################################################################
################## Robust regression (L2E criterion) #######################
############################################################################

library(mcmcse)
library(Matrix)
library(ks)
library(stats)
library(fasta)
library(glmnet)
library(Rcpp)
library(doParallel)
source("robustreg_data.R")
sourceCpp("robustreg_functions.cpp")

w_start <- w_truth + 0.5

L_px <- 10
L_guo <- 10
iter <- 1e3
eps_p <- 0.085
eps_guo <- 0.003

parallel::detectCores()
num_cores <- 4
doParallel::registerDoParallel(cores = num_cores)
reps <- 4

output_rreg <- foreach(b = 1:reps) %dopar% 
  {
    ## Run samplers
    print(b)
    L_px <- ifelse(runif(1) <= 0.05, 1, 10)
    L_guo <- ifelse(runif(1) <= 0.05, 1, 10)
    
    phmc_time <- system.time(phmc_run <- phmc_cpp(Phi_mat, y, lambda = .01, iter = iter, eps_hmc = eps_p, 
                                      L = L_px, nu = 3, start = w_start, alpha = alpha, blather = T))
    
    guohmc_time <- system.time(guohmc_run <- guohmc_cpp(Phi_mat, y, lambda = .01, iter = iter, eps_hmc = eps_guo, 
                                        L = L_guo, nu = nu, start = w_start, alpha = alpha, blather = T))
    
    eps <-  0.0005
    mymala_time <- system.time(mymala_run <- phmc_cpp(Phi_mat, y, lambda = eps/2, iter = iter, eps_hmc = eps, 
                                    L = 1, nu = 3, start = w_start, alpha = alpha, blather = T))
    
    rwm_time <- system.time(rwm_run <- rwm_cpp(Phi_mat, y, iter = iter, h = .075 ,
                   start = w_start, alpha = alpha, nu = nu, blather = T))
    
    # Means
    all_means <- cbind(colMeans(rwm_run[[1]]), colMeans(phmc_run[[1]]), colMeans(mymala_run[[1]]),
                       colMeans(nshmc_run[[1]]), colMeans(pmala_run[[1]]), colMeans(guohmc_run[[1]]))
    colnames(all_means) <- c("RWM", "pHMC", "myMALA", "nsHMC", "pMALA", "guoHMC")
    
    # ESS
    all_ess <- cbind(ess(rwm_run[[1]]),ess(phmc_run[[1]]), ess(mymala_run[[1]]), 
                     ess(nshmc_run[[1]]), ess(pmala_run[[1]]), ess(guohmc_run[[1]]))
    colnames(all_ess) <- c("RWM", "pHMC", "myMALA", "nsHMC", "pMALA", "guoHMC")
    
    all_time <- c(rwm_time[3], phmc_time[3], mymala_time[3], 
                  nshmc_time[3], pmala_time[3], guohmc_time[3])
    names(all_time) <- c("RWM", "pHMC", "myMALA", "nsHMC", "pMALA", "guoHMC")
    
    list(all_means, all_ess, all_time)
  }

save(output_rreg, file = "Output/outputrreg.Rdata")

