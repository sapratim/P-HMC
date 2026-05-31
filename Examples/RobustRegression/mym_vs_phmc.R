#########################################################
############### Code for myMALA vs pHMC acf #############
#########################################################

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
eps_p  <- 0.045
eps_mym  <- 0.05
L_px   <- 10    # leapfrog steps for pHMC / MALA
L_guo  <- 10    # leapfrog steps for Guo-HMC


phmc_time <- system.time(phmc_run <- phmc_cpp(B, y,
                                              lambda = lambda_prox, alpha = alpha, sigma = sigma,
                                              iter   = iter, eps_hmc = eps_p, L = L_px, nu = nu,
                                              start  = w_start, precond = precond_diag, blather = TRUE))[3]


eps        <- 0.05
mymala_time <- system.time(mymala_run <- phmc_cpp(B, y,
                                                  lambda = eps_mym / 2, alpha = alpha, sigma = sigma,
                                                  iter   = iter, eps_hmc = eps_mym, L = 1, nu = nu,
                                                  start  = w_start, precond = precond_diag, blather = TRUE))[3]


save(mymala_run, phmc_run, file = "Output/acf_mymvsphmc_samp.Rdata")

#############  Plots  #############

load("Output/acf_mymvsphmc_samp.Rdata")
pdf("Output/mym_vs_phmc_rreg.pdf", height = 3, width = 4.5)
par(mar = c(5, 4, 2, 2))

lag.max <- 100
plot(0:lag.max, rep(1, lag.max + 1), type = "n", ylim = c(-.02, 1),
     ylab = "Estimated autocorrelations", xlab = "Lags")
for(i in 1:7)
{
  mym_acfs <- acf(mymala_run[[1]][ ,i], plot = FALSE, lag.max = lag.max)$acf
  p_acfs <- acf(phmc_run[[1]][ ,i], plot = FALSE, lag.max = lag.max)$acf
  lines(0:lag.max, mym_acfs, lwd = 1.5, col = "red")
  lines(0:lag.max, p_acfs, lwd = 1.5, col = "purple", lty = 2)
}
# Add a legend on top of the plot
legend("top",
       legend = c("myMALA", "p-HMC"),
       col = c("red", "purple"),
       lty = c(1,2),
       lwd = 1.5,
       horiz = TRUE,
       bty = "n",
       inset = c(0, -0.19),  # pushes legend into the top margin
       xpd = TRUE)               # allow drawing outside plot region
dev.off()
