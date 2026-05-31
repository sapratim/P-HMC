#########################################################
############### Code for myMALA vs pHMC acf #############
#########################################################
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
eps_px <- 0.0075
eps_mym <-  0.0038
L <- 10
blat <- TRUE

phmc_time <- system.time(phmc_run <- phmc_cpp(y=y, alpha = alpha_hat,lambda = lamb_coeff, sigma2 = sigma2_hat, 
                  iter = iter, eps_hmc = eps_px, L_val=L, start = warmup_end_iter, blather = blat))


mymala_time <- system.time(mymala_run <- phmc_cpp(y=y, alpha = alpha_hat,lambda = eps_mym/2, sigma2 = sigma2_hat, 
                 iter = iter, eps_hmc = eps_mym, L_val = 1, start = warmup_end_iter, blather = blat) )

lag.max <- 100
mym_acfs <- matrix(0, nrow = lag.max, ncol = 4096)
p_acfs <- matrix(0, nrow = lag.max, ncol = 4096)
for(i in 1:4096)
{
  mym_acfs[,i] <- acf(mymala_run[[1]][ ,i], plot = FALSE, lag.max = lag.max)$acf
  p_acfs[,i] <- acf(phmc_run[[1]][ ,i], plot = FALSE, lag.max = lag.max)$acf
}
acf_list <- list(mym_acfs, p_acfs)

save(acf_list, file = "Output/acf_mymvsphmc_samp.Rdata")

#############  Plots  #############

load("Output/acf_mymvsphmc_samp.Rdata")
pdf("Output/mym_vs_phmc_nn.pdf", height = 3, width = 4.5)
par(mar = c(5, 4, 2, 2))

lag.max <- 100
plot(0:lag.max, rep(1, lag.max + 1), type = "n", ylim = c(-.02, 1),
     ylab = "Estimated autocorrelations", xlab = "Lags")
for(i in 1:4096)
{
  mym_acfs <- acf_list[[1]]  #acf(mymala_run[[1]][ ,i], plot = FALSE, lag.max = lag.max)$acf
  p_acfs <- acf_list[[2]]                  #acf(phmc_run[[1]][ ,i], plot = FALSE, lag.max = lag.max)$acf
  lines(0:lag.max, mym_acfs[,i], lwd = 1.5, col = "red")
  lines(0:lag.max, p_acfs[,i], lwd = 1.5, col = "purple", lty = 2)
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
