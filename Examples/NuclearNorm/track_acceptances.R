library(Rcpp)
sourceCpp("nn_functions.cpp")
source("nuclear_norm_functions.R")
load("warmup_chain.Rdata")
iter <- 1e2
lamb_coeff <- 1e-4
sigma2_hat <- 0.01
alpha_hat <- 1.15/sigma2_hat

L <- 10


lam_vec <- seq(1e-6, 1e-2, length = 10)
eps_vec <- seq(1e-6, 1e-2, length = 10)
acc_dur <- matrix(0, nrow =  length(eps_vec), ncol = length(lam_vec))
acc_px <- matrix(0, nrow =  length(eps_vec), ncol = length(lam_vec))

pt <- proc.time()
for(j in 1:length(lam_vec))
{
  lambda_coeff <- lam_vec[j]
  print(j)
  for(k in 1:length(eps_vec))
  {
    eps <- eps_vec[k]
    result_pxhmc <- pxhmc_cpp(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat, 
                          iter = iter, eps_hmc = eps, L=L, start = warmup_end_iter)
    
    result_durhmc <- durhmc_cpp(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                            iter = iter, eps_hmc = eps, L=L, start = warmup_end_iter)
    
    acc_px[k, j] <- result_pxhmc[[2]]
    acc_dur[k, j] <- result_durhmc[[2]]
  }
}
proc.time() - pt

save(acc_px, acc_dur, lam_vec, eps_vec, file = "acceptances.Rdata")
