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
eps_px <- 0.0075
L <- 10


parallel::detectCores()
num_cores <- 10
doParallel::registerDoParallel(cores = num_cores)
reps <- 100

blat <- TRUE
output_nn_true <- foreach(b = 1:reps) %dopar% 
  {
    ## Run samplers
    print(b)
    
    phmc_time <- system.time(phmc_run <- phmc_cpp(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                                                iter = iter, eps_hmc = eps_px, L_val = L, start = warmup_end_iter, blather = blat))
  
    post_mean <- colMeans(phmc_run[[1]])
    list(post_mean)
  }

save(output_nn_true, file = "Output/outputnn_true.Rdata")
