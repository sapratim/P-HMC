######################################################################
################# Nuclear Norm regression   ##########################
######################################################################

library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(mcmcse)

set.seed(1)
source("nn_data.R")
sourceCpp("nn_functions.cpp")
load("warmup_chain.Rdata")

iter <- 1e5
lamb_coeff <- 1e-4
eps_ns <-  0.00008
eps_px <- 0.0075
L <- 10
blat <- TRUE

nshmc_time <- system.time(nshmc_run <- nshmc_cpp(y=y, alpha = alpha_hat, lambda = 1, sigma2 = sigma2_hat, 
                                                 iter = iter, eps_hmc = eps_ns, L_val = L, start = warmup_end_iter, blather = blat))

phmc_time <- system.time(phmc_run <- phmc_cpp(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                                              iter = iter, eps_hmc = eps_px, L_val=L, start = warmup_end_iter, blather = blat))

lag.max <- 100
ns_acfs <- matrix(0, nrow = lag.max + 1, ncol = 4096)
p_acfs <- matrix(0, nrow = lag.max + 1, ncol = 4096)
for(i in 1:4096)
{
  ns_acfs[,i] <- acf(nshmc_run[[1]][ ,i], plot = FALSE, lag.max = lag.max)$acf
  p_acfs[, i] <- acf(phmc_run[[1]][ ,i], plot = FALSE, lag.max = lag.max)$acf
}

save(ns_acfs, p_acfs, file = "Output/nn_acf_old.Rdata")

###########  Plots

# load("Output/nn_acf_old.Rdata")
# pdf("Output/nn_acf_old.pdf", height = 3, width = 4.5)
# par(mar = c(5, 4, 2, 2))
# 
# lag.max <- 100
# plot(0:lag.max, rep(1, lag.max + 1), type = "n", ylim = c(-.02, 1),
#      ylab = "Estimated autocorrelations", xlab = "Lags")
# for(i in 1:4096)
# {
#   lines(0:lag.max, ns_acfs[,i][1:(lag.max + 1)], lwd = 1.5, col = "orange")
#   lines(0:lag.max, p_acfs[,i][1:(lag.max + 1)], lwd = 1.5, col = "purple")
# }
# # Add a legend on top of the plot
# legend("top",
#        legend = c("ns-HMC", "p-HMC"),
#        col = c("orange", "purple"),
#        lty = c(1,2),
#        lwd = 1.5,
#        horiz = TRUE,
#        bty = "n",
#        inset = c(0, -0.19),  # pushes legend into the top margin
#        xpd = TRUE)               # allow drawing outside plot region
# dev.off()
