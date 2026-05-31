
############################################################################################
############################# Evaluate optimal lambda ######################################
############################################################################################

library(Rcpp)
library(RcppArmadillo)

set.seed(1)
source("nn_data.R")
sourceCpp("nn_functions.cpp")
load("warmup_chain.Rdata")

iter <- 5e3
eps_ns <-  .0045
lambda_grid <- seq(1e-6, 0.0001, length = 100)
L <- 10
acf_vec_min <- numeric(length = length(lambda_grid))
acf_vec_mean <- numeric(length = length(lambda_grid))
blat <- TRUE

for (i in 1:length(lambda_grid)) {
  print(i)
  nshmc_time <- system.time(nshmc_run <- nshmc_cpp(y=y, alpha = alpha_hat, lambda = lambda_grid[i], 
                                  sigma2 = sigma2_hat, iter = iter, eps_hmc = eps_ns, 
                                  L_val = L, start = warmup_end_iter, blather = blat))
  samp_mat <- nshmc_run[[1]]
  acf_lag1 <- sapply(seq_len(ncol(samp_mat)), function(j) {
    cor(samp_mat[-1, j], samp_mat[-nrow(samp_mat), j])})
  acf_vec_min[i] <- min(acf_lag1)
  acf_vec_mean[i] <- mean(acf_lag1)
}
acf_output <- list(acf_vec_min, acf_vec_mean)

#####opt_lambda <- lambda_grid[which.min(acf_vec)]

save(acf_output, file = "Output/opt_lambda.Rdata")
