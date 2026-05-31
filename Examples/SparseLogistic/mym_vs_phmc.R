#########################################################
############### Code for myMALA vs pHMC acf #############
#########################################################

library(mcmcse)
library(Matrix)
library(ks)
library(stats)
library(fasta)
library(glmnet)
library(Rcpp)
library(doParallel)
source("slogistic_data.R")
sourceCpp("slogistic_functions.cpp")

# starting values
logistic_fit <- glmnet(x, y, family = "binomial",
                       alpha = 1, lambda = alpha/length(y), nlambda = 1,
                       standardize = FALSE, intercept = FALSE)$beta

beta_start <- as.matrix(unname(logistic_fit ))

L_px <- 10
eps_p <-  0.0019
iter <- 1e5
blat <- TRUE

phmc_time <- system.time(phmc_run <- phmc_cpp(x, y, lambda = .01, iter = iter, eps_hmc = eps_p, 
                                   L_val = L_px, start = beta_start, alpha = alpha, blather = blat))
eps_mym <-  0.0019
mymala_time <- system.time(mymala_run <- phmc_cpp(x, y, lambda = eps_mym/2, alpha = alpha, iter = iter,
                                   eps_hmc = eps_mym, L_val = 1, start = beta_start, blather = blat))

save(mymala_run, phmc_run, file = "Output/acf_mymvsphmc_samp.Rdata")

#############  Plots  #############

load("Output/acf_mymvsphmc_samp.Rdata")
pdf("Output/mym_vs_phmc_slr.pdf", height = 3, width = 4.5)
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
       legend = c("my-MALA", "p-HMC"),
       col = c("red", "purple"),
       lty = c(1,2),
       lwd = 1.5,
       horiz = TRUE,
       bty = "n",
       inset = c(0, -0.19),  # pushes legend into the top margin
       xpd = TRUE)               # allow drawing outside plot region
dev.off()
