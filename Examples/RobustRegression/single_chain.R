################################################################
######## Single chain for acf and credible interval#############
################################################################

set.seed(209)
library(Rcpp)
library(mcmcse)
library(foreach)
library(doParallel)
source("robustreg_data.R")          # loads model settings as well
sourceCpp("pre_robustreg_functions.cpp")
load("marginal_vars.RData")

###############   Obtain MAP estimate
MAP <- map_estimate(B, y, alpha, nu, sigma, w_truth)

# start a little off the MAP to avoid zero-gradient issues 
w_start <- MAP + rnorm(length(MAP), 0, 0.01)

# post_var_diag from loading marginal variances
precond_diag <- post_var_diag

iter <- 1e5
lambda_prox <- .002
L_px        <- 20    # leapfrog steps for pHMC / MALA

eps_p    <- 0.04
phmc_time <- system.time(phmc_run <- phmc_cpp(B, y,
                                     lambda = lambda_prox, alpha = alpha, sigma = sigma,
                                     iter   = iter, eps_hmc = eps_p, L = L_px, nu = nu,
                                     start  = w_start, precond = precond_diag, blather = TRUE))[3]

save(phmc_run, file = "Output/single_phmc_chain.Rdata")

# CI_mat <- matrix(0, nrow = 2, ncol = length(t_index))
# CI_mat[1,] <- upper_quantiles[t_index]
# CI_mat[2,] <- lower_quantiles[t_index] 

# 
# # Lower and upper credible limits
# lower <- CI_mat[2, ]
# upper <- CI_mat[1, ]
# 
# # Point estimates (optional)
# est <- (lower + upper)/2
# 
# pdf(file = "Output/Credible_Intervals.pdf", height = 6, width = 12)
# # Plot signal
# plot(w_truth,
#      type = "l",
#      lwd  = 1)
# 
# # # Vertical credible interval bars
# # segments(
# #   x0  = t_index,
# #   y0  = lower,
# #   x1  = t_index,
# #   y1  = upper,
# #   lwd = 1,
# #   col = "red"
# # )
# 
# # Width of horizontal caps
# cap <- 1.5
# 
# # Lower caps
# segments(
#   x0 = t_index - cap,
#   y0 = lower,
#   x1 = t_index + cap,
#   y1 = lower,
#   lwd = 2,
#   col = "red"
# )
# 
# # Upper caps
# segments(
#   x0 = t_index - cap,
#   y0 = upper,
#   x1 = t_index + cap,
#   y1 = upper,
#   lwd = 2,
#   col = "red"
# )
# 
# # Optional: estimated points
# points(t_index, est,
#        pch = 16,
#        col = "blue")
# 
# dev.off()
# 
# plot(w_truth, type = 'l')
# lines(colMeans(phmc_run[[1]]), type = 'l', col = "red")
# lines(colMeans(rwm_run[[1]]), type = 'l', col = "red")

# dens_values <- apply(phmc_run[[1]], 1, function(x) log_pi(x, y, Phi_mat, nu, alpha, sigma))
# MAP_sample <- phmc_run[[1]][which.max(dens_values),]

pdf(file = "Output/MAP_plot.pdf", width = 12, height = 6)
par(mar = c(5, 6, 4, 2))  # bottom, left, top, right
plot(w_truth, type = 'l', lwd = 1, cex.lab = 2.5, cex.axis = 3, ylab = "True signal")
lines(MAP, type = 'l', col = "red", lwd = 1)
dev.off()

pdf(file = "Output/noisy_output.pdf", width = 12, height = 6)
par(mar = c(5, 6, 4, 2))  # bottom, left, top, right
plot(y, type = 'l', lwd = 1, cex.lab = 2.5, cex.axis = 3)
dev.off()

pdf(file = "Output/truth_plot.pdf", width = 12, height = 6)
par(mar = c(5, 6, 4, 2))  # bottom, left, top, right
plot(w_truth, type = 'l', lwd = 1, cex.lab = 2.5, cex.axis = 3, ylab = "True signal")
dev.off()


