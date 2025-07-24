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
L <- 10


output_poisson <- list()
parallel::detectCores()
num_cores <- 10
doParallel::registerDoParallel(cores = num_cores)
reps <- 100

blat <- TRUE
output_slog <- foreach(b = 1:reps) %dopar% 
{
## Run samplers
  print(b)
  nshmc_time <- system.time(nshmc_run <- nshmc_cpp(y=y, alpha = alpha_hat, lambda = 1, sigma2 = sigma2_hat, 
                                                   iter = iter, eps_hmc = eps_ns, L = L, start = warmup_end_iter, blather = blat))
  
  eps <-  0.01
  pmala_time <- system.time(pmala_run <- nshmc_cpp(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat, 
                                                   iter = iter, eps_hmc = eps, L = 1, start = warmup_end_iter, blather = blat))
  
  phmc_time <- system.time(phmc_run <- phmc_cpp(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                                                iter = iter, eps_hmc = eps_px, L=L, start = warmup_end_iter, blather = blat) )
  
  eps <-  0.011
  mymala_time <- system.time(mymala_run <- phmc_cpp(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                                                    iter = iter, eps_hmc = eps, L = 1, start = warmup_end_iter, blather = blat) )
  
  rwm_time <- system.time(rwm_run <- rwm_cpp(y=y, alpha = alpha_hat, sigma2 = sigma2_hat, 
                                                    iter = iter, proposal_sd = .002, start = warmup_end_iter, blather = blat))

  # ESS
  all_ess <- cbind(ess(rwm_run[[1]]),ess(phmc_run[[1]]), 
               ess(mymala_run[[1]]), ess(nshmc_run[[1]]), ess(pmala_run[[1]]))
  colnames(all_ess) <- c("RWM", "pHMC", "myMALA", "nsHMC", "pMALA")
  
  all_time <- c(rwm_time[3], phmc_time[3], mymala_time[3], nshmc_time[3], pmala_time[3])
  names(all_time) <- c("RWM", "pHMC", "myMALA", "nsHMC", "pMALA")
  
  list(all_ess, all_time)
}

save(output_slog, file = "outputnn.Rdata")
