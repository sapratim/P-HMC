
rm(list = ls())
source("nuclear_norm_functions.R")
iter_chain <- 1e5
lamb_coeff <- 1e-5
sigma2_hat <- 0.01
alpha_hat <- 1.15/sigma2_hat
delta_mym <- 0.0001 

warmup_chain <- mymala(y,alpha_hat,lamb_coeff,sigma2_hat,iter_chain,delta_mym,y)
warmup_end_iter <- warmup_chain[[1]][iter_chain,]

save(warmup_end_iter, file = "warmup_chain.Rdata")
