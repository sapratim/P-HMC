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
L_guo <- 10
iter <- 1e5
eps_ns <- 0.00012
eps_p <-  0.0019
eps_guo <- 0.00012
blat <- FALSE

parallel::detectCores()
num_cores <- 10
doParallel::registerDoParallel(cores = num_cores)
reps <- 100


output_slog <- foreach(b = 1:reps) %dopar% 
{
## Run samplers
  print(b)
  nshmc_time <- system.time(nshmc_run <- nshmc_cpp(x, y, lambda = 1, alpha = alpha, iter = iter,
                             eps_hmc = eps_ns, L_val = L_ns, start = beta_start, blather = blat))
  
  eps <-  0.0016
  pmala_time <- system.time(pmala_run <- nshmc_cpp(x, y, lambda = eps/2, alpha = alpha, iter = iter,
                            eps_hmc = eps, L_val = 1, start = beta_start, blather = blat))
  
  phmc_time <- system.time(phmc_run <- phmc_cpp(x, y, lambda = .01, iter = iter, eps_hmc = eps_p, 
                          L_val = L_px, start = beta_start, alpha = alpha, blather = blat))
  
  guohmc_time <- system.time(guohmc_run <- guohmc_cpp(x, y, lambda = .01, iter = iter, eps_hmc = eps_guo, 
                          L_val = L_guo, start = beta_start, alpha = alpha, blather = blat))
  
  eps <-  0.0019
  mymala_time <- system.time(mymala_run <- phmc_cpp(x, y, lambda = eps/2, alpha = alpha, iter = iter,
                            eps_hmc = eps, L_val = 1, start = beta_start, blather = blat))
  
  rwm_time <- system.time(rwm_run <- rwm_cpp(x, y, iter = iter, 
                             h = .0045 ,start = beta_start, alpha = alpha, blather = blat))

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

save(output_slog, file = "Output/outputslog.Rdata")
