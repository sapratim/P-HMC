
##########################################
# Single run of the Nuclear Norm Example
##########################################
library(Rcpp)
library(RcppArmadillo)

set.seed(1)
source("nn_data.R")
#source("nuclear_norm_functions.R")
sourceCpp("nn_functions.cpp")
load("warmup_chain.Rdata")


iter <- 1e6
lamb_coeff <- 1e-4
eps_dur <- 0.008
eps_px <-  0.004
eps <- .004
L <- 10


system.time(result_pxhmc <- pxhmc_cpp(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = eps_px, L=L, start = warmup_end_iter) )

system.time(result_durhmc <- durhmc_cpp(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = eps_dur, L=L, start = warmup_end_iter) )

output_single_hmc <- list(result_durhmc[[1]], result_pxhmc[[1]])

save(output_single_hmc, eps_dur, eps_px, L, lamb_coeff, file = "single_run.Rdata")

