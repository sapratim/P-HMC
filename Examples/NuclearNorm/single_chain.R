
source("nuclear_norm_functions.R")
load("warmup_chain.Rdata")
iter <- 1e2
lamb_coeff <- 1e-4
sigma2_hat <- 0.01
alpha_hat <- 1.15/sigma2_hat
step_ismala <- 0.00012
step_pxmala <- 0.0001
step_isb <- 0.0001
step_pxb <- 0.0001
eps_dur <- 0.008
eps_px <-  0.004
eps <- .004
L <- 10


system.time(result_pxhmc <- pxhmc(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = eps_px, L=L, start = warmup_end_iter) )

system.time(result_durhmc <- durhmc(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                      iter = iter, eps_hmc = eps_dur, L=L, start = warmup_end_iter) )



library(Rcpp)
sourceCpp("nn_functions.cpp")
load("warmup_chain.Rdata")
iter <- 1e3
lamb_coeff <- 1e-4
sigma2_hat <- 0.01
alpha_hat <- 1.15/sigma2_hat
step_ismala <- 0.00012
step_pxmala <- 0.0001
step_isb <- 0.0001
step_pxb <- 0.0001
eps_dur <- 0.008
eps_px <-  0.004
eps <- .004
L <- 10


system.time(result_pxhmc_cpp <- pxhmc_cpp(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat, 
                                  iter = iter, eps_hmc = eps_px, L=L, start = warmup_end_iter) )

system.time(result_durhmc <- durhmc(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                                    iter = iter, eps_hmc = eps_dur, L=L, start = warmup_end_iter) )



lam_vec <- seq(1e-6, 1e-2, length = 20)
eps_vec <- seq(1e-6, 1e-2, length = 20)
acc_dur <- matrix(0, nrow =  length(eps_vec), ncol = length(lam_vec))
acc_px <- matrix(0, nrow =  length(eps_vec), ncol = length(lam_vec))
for(j in 1:length(lam_vec))
{
  print(j)
  for(k in 1:length(eps_vec))
  {
    eps <- eps_vec[k]
    result_pxhmc <- pxhmc(y=y, alpha = alpha_hat, lambda = lamb_coeff, sigma2 = sigma2_hat, 
                                      iter = iter, eps_hmc = eps, L=L, start = warmup_end_iter)
    
    result_durhmc <- durhmc(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                                        iter = iter, eps_hmc = eps, L=L, start = warmup_end_iter)
    
    acc_px[k, j] <- result_pxhmc[[2]]
    acc_dur[k, j] <- result_durhmc[[2]]
  }
}
output_single_hmc <- list(result_durhmc[[1]], result_pxhmc[[1]])


rand <- 1:length(y) #sample(c(1:length(y)), subset)
dim <- length(y)
lag.max <- 100
acf_dur_hmc <- acf(output_single_hmc[[1]][,rand[1]], plot = FALSE, lag.max = lag.max)$acf
acf_pxhmc <- acf(output_single_hmc[[2]][,rand[1]], plot = FALSE, lag.max = lag.max)$acf

diff.acf <- matrix(0, ncol = dim, nrow = lag.max + 1)
diff.acf[,1] <- acf_is_hmc - acf_pxhmc


acf(result_durhmc[[1]][, 2])
acf(result_pxhmc[[1]][, 2])

for (i in 1:dim) 
{
  if(i %% 1000 == 0) print(i)
  acf_dur_hmc <- acf(output_single_hmc[[1]][,i], plot = FALSE, lag.max = lag.max)$acf
  acf_pxhmc <- acf(output_single_hmc[[2]][,i], plot = FALSE, lag.max = lag.max)$acf
  diff.acf[,i] <- acf_dur_hmc - acf_pxhmc
}


# Make boxplot of ACFs
boxplot(t(diff.acf),
        xlab = "Lags", col = "pink",
        ylab = "Difference in ACFs of HMCs",ylim = range(diff.acf),
        names = 0:lag.max, show.names = TRUE, range = 3)
#dev.off()



# save(output_single_mala, file = "output_single_chain_mala.Rdata")
# save(output_single_bark, file = "output_single_chain_bark.Rdata")
# save(output_single_hmc, file = "output_single_chain_hmc.Rdata")

