
################  Check if it works

library(Rcpp)
source("robustreg_data.R")
sourceCpp("robustreg_functions.cpp")
L_px <- 10
L_guo <- 10
eps_p <- 0.075
eps_guo <- 0.003
w_start <- w_truth + 0.02
alpha <- 1
iter <- 1e4
nu <- 5

phmc_run <- phmc_cpp(Phi_mat, y, lambda = .01, iter = iter, eps_hmc = eps_p, 
                     L = L_px, nu = nu, start = w_start, alpha = alpha, blather = T)

guohmc_run <- guohmc_cpp(Phi_mat, y, lambda = .01, iter = iter, eps_hmc = eps_guo, 
           L = L_guo, nu = nu, start = w_start, alpha = alpha, blather = T)

eps <-  0.0804
mymala_run <- phmc_cpp(Phi_mat, y, lambda = eps/2, iter = iter, eps_hmc = eps, 
                       L = 1, nu = nu, start = w_start, alpha = alpha, blather = T)

rwm_run <- rwm_cpp(Phi_mat, y, iter = iter, h = .042 ,
      start = w_start, alpha = alpha, nu = nu, blather = T)

################ Plots and check  ################

#####  P-HMC
plot(w_truth, type = 'l')
lines(colMeans(phmc_run[[1]]), type = 'l', col = "red")
samp_nonzero <- phmc_run[[1]][,t_index]
(CI <- apply(samp_nonzero, 2, function(x) quantile(x, probs = c(0.025, 0.975))))

##### Guo-HMC
plot(w_truth, type = 'l')
lines(colMeans(guohmc_run[[1]]), type = 'l', col = "red")
samp_nonzero <- guohmc_run[[1]][,t_index]
(CI <- apply(samp_nonzero, 2, function(x) quantile(x, probs = c(0.025, 0.975))))

###### MyMALA
plot(w_truth, type = 'l')
lines(colMeans(mymala_run[[1]]), type = 'l', col = "red")
samp_nonzero <- mymala_run[[1]][,t_index]
(CI <- apply(samp_nonzero, 2, function(x) quantile(x, probs = c(0.025, 0.975))))

###### RWM
plot(w_truth, type = 'l')
lines(colMeans(rwm_run[[1]]), type = 'l', col = "red")
samp_nonzero <- rwm_run[[1]][,t_index]
(CI <- apply(samp_nonzero, 2, function(x) quantile(x, probs = c(0.025, 0.975))))
