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


