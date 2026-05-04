######################################################################
################# Nuclear Norm regression   ##########################
######################################################################

library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(mcmcse)

set.seed(1)
source("nn_data.R")
sourceCpp("nn_functions.cpp")
load("warmup_chain.Rdata")

iter <- 1e5
lamb_coeff <- 1e-4
eps_ns <-  .00008
eps_px <- 0.0075
eps_guo <- 0.00008
L <- 10


output_poisson <- list()
parallel::detectCores()
num_cores <- 10
doParallel::registerDoParallel(cores = num_cores)
reps <- 100

blat <- TRUE
output_nn <- foreach(b = 1:reps) %dopar% 
{
## Run samplers
  print(b)
  
  L <- ifelse(runif(1) <= 0.05, 1, 10)
  nshmc_time <- system.time(nshmc_run <- nshmc_cpp(y=y, alpha = alpha_hat, lambda = 1, sigma2 = sigma2_hat, 
                                                   iter = iter, eps_hmc = eps_ns, L = L, start = warmup_end_iter, blather = blat))
  eps <-  0.01
  pmala_time <- system.time(pmala_run <- nshmc_cpp(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat, 
                                                   iter = iter, eps_hmc = eps, L = 1, start = warmup_end_iter, blather = blat))
  
  phmc_time <- system.time(phmc_run <- phmc_cpp(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                                                iter = iter, eps_hmc = eps_px, L=L, start = warmup_end_iter, blather = blat))
  
  guohmc_time <- system.time(guohmc_run <- guohmc_cpp(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                                                iter = iter, eps_hmc = eps_guo, L=L, start = warmup_end_iter, blather = blat))
  eps <-  0.011
  mymala_time <- system.time(mymala_run <- phmc_cpp(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                                                    iter = iter, eps_hmc = eps, L = 1, start = warmup_end_iter, blather = blat) )
  
  rwm_time <- system.time(rwm_run <- rwm_cpp(y=y, alpha = alpha_hat, sigma2 = sigma2_hat, 
                                                    iter = iter, proposal_sd = .002, start = warmup_end_iter, blather = blat))

  # Means
  all_means <- cbind(colMeans(rwm_run[[1]]), colMeans(phmc_run[[1]]), colMeans(mymala_run[[1]]),
                        colMeans(nshmc_run[[1]]), colMeans(pmala_run[[1]]), colMeans(guohmc_run[[1]]))
  colnames(all_means) <- c("RWM", "pHMC", "myMALA", "nsHMC", "pMALA", "guoHMC")
  
  # ESS
  all_ess <- cbind(ess(rwm_run[[1]]),ess(phmc_run[[1]]), ess(mymala_run[[1]]), 
                   ess(nshmc_run[[1]]), ess(pmala_run[[1]]), ess(guohmc_run[[1]]))
  colnames(all_ess) <- c("RWM", "pHMC", "myMALA", "nsHMC", "pMALA", "guoHMC")
  
  all_time <- c(rwm_time[3], phmc_time[3], mymala_time[3], nshmc_time[3], pmala_time[3], guohmc_time[3])
  names(all_time) <- c("RWM", "pHMC", "myMALA", "nsHMC", "pMALA", "guoHMC")
  
  list(all_means, all_ess, all_time)
}

save(output_nn, file = "outputnn.Rdata")
