######################################################################
################# Sparse logistic regression   #######################
######################################################################

library(mcmcse)
library(Matrix)
library(ks)
library(stats)
library(fasta)
library(glmnet)
library(Rcpp)
library(doParallel)
source("slogistic_data.R")
sourceCpp("slogistic_functions.cpp")

# starting values
logistic_fit <- glmnet(x, y, family = "binomial",
                       alpha = 1, lambda = alpha/length(y), nlambda = 1,
                       standardize = FALSE, intercept = FALSE)$beta

beta_start <- as.matrix(unname(logistic_fit ))

L_ns <- 10
L_px <- 10
iter <- 1e6
eps_ns <-  0.00015


output_poisson <- list()
parallel::detectCores()
num_cores <- 10
doParallel::registerDoParallel(cores = num_cores)
reps <- 100

output_slog <- foreach(b = 1:reps) %dopar% 
{
## Run samplers
  print(b)
  nshmc_time <- system.time(nshmc_run <- nshmc_cpp(x, y, lambda = 1e-4, alpha = alpha, iter = iter,
                                     eps_hmc = eps_ns, L = L_ns, start = beta_start, blather = FALSE) )
  
  eps <-  0.0018
  pmala_time <- system.time(pmala_run <- nshmc_cpp(x, y, lambda = eps/2, alpha = alpha, iter = iter,
                                     eps_hmc = eps, L = 1, start = beta_start, blather = FALSE) )
  
  phmc_time <- system.time(phmc_run <- phmc_cpp(x, y, lambda = .01, iter = iter, eps_hmc = 0.002, 
                                                L = L_px, start = beta_start, alpha = alpha, blather = FALSE) )
  
  eps <-  0.002
  mymala_time <- system.time(mymala_run <- phmc_cpp(x, y, lambda = eps/2, alpha = alpha, iter = iter,
                                     eps_hmc = eps, L = 1, start = beta_start, blather = FALSE) )
  
  rwm_time <- system.time(rwm_run <- rwm_cpp(x, y, iter = iter, 
                                             h = .005 ,start = beta_start, alpha = alpha, blather = FALSE) )

  # ESS
  all_ess <- cbind(ess(rwm_run[[1]]),ess(phmc_run[[1]]), 
               ess(mymala_run[[1]]), ess(nshmc_run[[1]]), ess(pmala_run[[1]]))
  colnames(all_ess) <- c("RWM", "pHMC", "myMALA", "nsHMC", "pMALA")
  
  all_time <- c(rwm_time[3], phmc_time[3], mymala_time[3], nshmc_time[3], pmala_time[3])
  names(all_time) <- c("RWM", "pHMC", "myMALA", "nsHMC", "pMALA")
  
  list(all_ess, all_time)
  # list(asymp_covmat_ism, asymp_covmat_pxm, asymp_covmat_isb, asymp_covmat_pxb, asymp_covmat_trubark,
  #      asymp_covmat_ishmc, asymp_covmat_pxhmc, n_eff_mala, n_eff_bark, n_eff_hmc)
}

save(output_slog, file = "outputslog.Rdata")
