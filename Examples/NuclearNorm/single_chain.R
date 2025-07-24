##########################################
# Single run of the Nuclear Norm Example
##########################################
library(Rcpp)
library(RcppArmadillo)

set.seed(1)
source("nn_data.R")
sourceCpp("nn_functions.cpp")
load("warmup_chain.Rdata")

## Map calculator
nuclear_prox_sigma <- function(Y, alpha, sigma2) {
  svd_Y <- svd(Y)
  lambda <- sigma2 * alpha
  d_thresh <- pmax(svd_Y$d - lambda, 0)
  return(svd_Y$u %*% diag(d_thresh) %*% t(svd_Y$v))
}
# 
MAP1 <- nuclear_prox_sigma(image_mat, alpha_hat, sigma2_hat)

iter <- 1e5
lamb_coeff <- 1e-4
eps_ns <-  .00008
eps_px <- 0.0075
L <- 10


system.time(result_nshmc <- nshmc_cpp(y=y, alpha = alpha_hat, lambda = 1, sigma2 = sigma2_hat,
                      iter = iter, eps_hmc = eps_ns, L=L, start = warmup_end_iter) )

system.time(result_phmc <- phmc_cpp(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = eps_px, L = L, start = warmup_end_iter) )

lag.max <- 100
ns_acfs <- matrix(0, nrow = lag.max + 1, ncol = 4096)
p_acfs <- matrix(0, nrow = lag.max + 1, ncol = 4096)
for(i in 1:4096)
{
  ns_acfs[,i] <- acf(result_nshmc[[1]][ ,i], plot = FALSE, lag.max = lag.max)$acf
  p_acfs[, i] <- acf(result_phmc[[1]][ ,i], plot = FALSE, lag.max = lag.max)$acf
}


upper_quantile <- apply(result_phmc[[1]], 2, quantile, .95)
lower_quantile <- apply(result_phmc[[1]], 2, quantile, .05)
cred_interval <- upper_quantile - lower_quantile

post_mean <- apply(result_phmc[[1]], 2, mean)
post_mat <- matrix(post_mean, nrow = n, ncol = n)

save(ns_acfs, p_acfs, checker, image_mat, MAP1, cred_interval, file = "nn_single_image.Rdata")
