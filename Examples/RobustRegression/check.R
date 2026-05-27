################################################################################
##  Robust EEG
################################################################################

set.seed(209)
library(Rcpp)
library(mcmcse)
# loads model settings as well
source("robustreg_data.R")
sourceCpp("pre_robustreg_functions.cpp")
load("marginal_vars.RData")

## Sampler / problem settings
iter        <- 5e4
lambda_prox <- .002
L_px        <- 10    # leapfrog steps for pHMC / MALA
L_guo       <- 10    # leapfrog steps for Guo-HMC

## ---- 2. MAP and sanity checks ------------------------------------------------

MAP <- map_estimate(B, y, alpha, nu, sigma, rep(0, length(w_truth)))

cbind(summary(MAP), summary(w_truth))
log_pi(w_truth, y, B, nu, alpha, sigma)
log_pi(MAP,     y, B, nu, alpha, sigma)

# start a little off the MAP to avoid zero-gradient issues 
w_start <- MAP + rnorm(length(MAP), 0, 0.01)

# post_var_diag from loading marginal variances
precond_diag <- post_var_diag


## ---- 4. Production runs with the estimated preconditioner --------------------
cat("--- pHMC ---\n")
eps_p    <- 0.048
phmc_time <- system.time(phmc_run <- phmc_cpp(B, y,
                                              lambda = lambda_prox, alpha = alpha, sigma = sigma,
                                              iter   = iter, eps_hmc = eps_p, L = L_px, nu = nu,
                                              start  = w_start, precond = precond_diag, blather = TRUE))[3]


cat("\n--- Guo-HMC ---\n")
eps_guo    <- 0.0008
guohmc_time <- system.time(guohmc_run <- guohmc_cpp(B, y,
                                                    lambda = lambda_prox, alpha = alpha, sigma = sigma,
                                                    iter   = iter, eps_hmc = eps_guo, L = L_guo, nu = nu,
                                                    start  = w_start, precond = precond_diag, blather = TRUE))[3]

cat("\n--- MALA (pHMC with L = 1, lambda = eps/2) ---\n")
eps        <- 0.05
mymala_time <- system.time(mymala_run <- phmc_cpp(B, y,
                                                  lambda = eps / 2, alpha = alpha, sigma = sigma,
                                                  iter   = iter, eps_hmc = eps, L = 1, nu = nu,
                                                  start  = w_start, precond = precond_diag, blather = TRUE))[3]

cat("\n--- RWM ---\n")
rwm_time <- system.time(rwm_run <- rwm_cpp(B, y,
                                           iter  = iter, h = 0.02,
                                           start = w_start, alpha = alpha, sigma = sigma, nu = nu,
                                           precond = precond_diag, blather = TRUE))[3]

cat("\nAcceptance rates:\n")
print(round(c(pHMC   = phmc_run$accept_rate,
              guoHMC = guohmc_run$accept_rate,
              myMALA = mymala_run$accept_rate,
              RWM    = rwm_run$accept_rate), 3))

## ---- 5. Diagnostics & credible-band summary ----------------------------------

acf_all <- function(i)
{
  rw <- acf(rwm_run[[1]][,i], plot = FALSE)$acf
  phmc <- acf(phmc_run[[1]][,i], plot = FALSE)$acf
  mymala <- acf(mymala_run[[1]][,i], plot = FALSE)$acf
  guohmc <- acf(guohmc_run[[1]][,i], plot = FALSE)$acf
  
  plot(rw, col = "black", type = 'l',
       ylim = c(0,1))
  lines(phmc, col = "blue")
  lines(mymala, col = "darkgreen")
  lines(guohmc, col = "red")
  legend("bottomleft", legend = c("RWM", "pHMC", "myMALA", "Guo-HMC"),
         col = c("black", "blue", "darkgreen", "red"), lty = 1, cex = 0.4)
}

den_all <- function(i)
{
  
  plot(density(rwm_run[[1]][,i]), col = "black", 
       main = paste("Density coord", i), 
       xlim = range(c(density(rwm_run[[1]][,i])$x, density(phmc_run[[1]][,i])$x,
                      density(mymala_run[[1]][,i])$x, density(guohmc_run[[1]][,i])$x)))
  lines(density(phmc_run[[1]][,i]), col = "blue")
  lines(density(mymala_run[[1]][,i]), col = "darkgreen")
  lines(density(guohmc_run[[1]][,i]), col = "red")
  legend("topright", legend = c("RWM", "pHMC", "myMALA", "Guo-HMC"),
         col = c("black", "blue", "darkgreen", "red"), lty = 1, cex = 0.4)
}

foo <- phmc_run[[1]]
sample_logp <- apply(foo, 1, function(x) log_pi(x, y, B, nu, alpha, sigma))
plot.ts(sample_logp)
cbind(summary(MAP), summary(w_truth), summary(colMeans(foo)))
## Per-coordinate trace / density spot-checks
i <- 470
acf_all(i)
den_all(i)
i <- i+1

# calculating ess
ess_phmc <- ess(phmc_run$samples)
ess_guohmc <- ess(guohmc_run$samples)
ess_mymala <- ess(mymala_run$samples)
ess_rwm <- ess(rwm_run$samples)

summarize_ess <- function(ess_vec)
{
  c(mean = mean(ess_vec), median = median(ess_vec), min = min(ess_vec), max = max(ess_vec))
}

cat("\nEffective sample size (ESS) summaries:\n")
print(round(rbind(pHMC = summarize_ess(ess_phmc),
                  GuoHMC = summarize_ess(ess_guohmc),
                  myMALA = summarize_ess(ess_mymala),
                  RWM = summarize_ess(ess_rwm)), 2))

# ESS per time
cat("\nESS per second:\n")
print(round(rbind(pHMC = summarize_ess(ess_phmc / phmc_time),
                  GuoHMC = summarize_ess(ess_guohmc / guohmc_time),
                  myMALA = summarize_ess(ess_mymala / mymala_time),
                  RWM = summarize_ess(ess_rwm / rwm_time)), 4))



######################  Check alpha for credible intervals  #####################

set.seed(209)
library(Rcpp)
library(mcmcse)
# loads model settings as well
source("robustreg_data.R")
sourceCpp("pre_robustreg_functions.cpp")
load("marginal_vars.RData")

## Sampler / problem settings
iter        <- 1e4
lambda_prox <- .002
L_px        <- 20    # leapfrog steps for pHMC / MALA

MAP <- map_estimate(B, y, alpha, nu, sigma, rep(0, length(MAP)), n_init = 10)
# start a little off the MAP to avoid zero-gradient issues 
w_start <- MAP + rnorm(length(MAP), 0, 0.01)

# post_var_diag from loading marginal variances
precond_diag <- post_var_diag

alpha <- 25

eps_p    <- 0.0012
phmc_time <- system.time(phmc_run <- phmc_cpp(B, y,
                                              lambda = lambda_prox, alpha = alpha, sigma = sigma,
                                              iter   = iter, eps_hmc = eps_p, L = L_px, nu = nu,
                                              start  = w_start, precond = precond_diag, blather = TRUE))[3]

upper_quantiles <- apply(phmc_run[[1]], 2, function(x) quantile(x, 0.975))
lower_quantiles <- apply(phmc_run[[1]], 2, function(x) quantile(x, 0.025))

CI_mat <- matrix(0, nrow = 2, ncol = length(t_index))
CI_mat[1,] <- upper_quantiles[t_index]
CI_mat[2,] <- lower_quantiles[t_index]


inside <- d2 <= threshold
X_inside <- X[inside, ]
joint_intervals <- apply(X_inside, 2, range)
#pdf(file = "Output/Credible_Intervals.pdf", height = 6, width = 12)

plot(w_truth, type = "l",lwd  = 1)
# segments(x0 = t_index, y0 = CI_mat[2,],        # Vertical credible interval bars
#          x1 = t_index, y1 = CI_mat[1,], lwd = 1, col = "red")
cap <- 1.5      # Width of horizontal caps
segments(x0 = t_index - cap, y0 = CI_mat[2,],   # Lower caps
         x1 = t_index + cap, y1 = CI_mat[2,], lwd = 2,col = "red")
segments(x0 = t_index - cap, y0 = CI_mat[1,],   # Upper caps
         x1 = t_index + cap, y1 = CI_mat[1,], lwd = 2, col = "red")
# Optional: estimated points
points(t_index, MAP[t_index],
       pch = 16,
       col = "blue")

#dev.off()



