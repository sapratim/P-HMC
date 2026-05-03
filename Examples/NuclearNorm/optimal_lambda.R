
############################################################################################
############################# Evaluate optimal lambda ######################################
############################################################################################

library(Rcpp)
library(RcppArmadillo)

set.seed(1)
source("nn_data.R")
sourceCpp("nn_functions.cpp")
load("warmup_chain.Rdata")

iter <- 1e3
eps_ns <-  .00008
lambda_grid <- seq(1e-5, 2, length = 100)
acf_vec <- numeric(length = length(lambda_grid))
blat <- TRUE

for (i in 1:length(lambda_grid)) {
  print(i)
  L <- ifelse(runif(1) <= 0.05, 1, 10)
  nshmc_time <- system.time(nshmc_run <- nshmc_cpp(y=y, alpha = alpha_hat, lambda = lambda_grid[i], 
                                  sigma2 = sigma2_hat, iter = iter, eps_hmc = eps_ns, 
                                  L = L, start = warmup_end_iter, blather = blat))
  samp_mat <- nshmc_run[[1]]
  acf_lag1 <- sapply(seq_len(ncol(samp_mat)), function(j) {
    cor(samp_mat[-1, j], samp_mat[-nrow(samp_mat), j])})
  acf_vec[i] <- mean(acf_lag1)
}
opt_lambda <- lambda_grid[which.min(acf_vec)]

save(opt_lambda, file = "Output/opt_lambda.Rdata")
