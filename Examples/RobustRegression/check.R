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
L_px        <- 20    # leapfrog steps for pHMC / MALA
L_guo       <- 20    # leapfrog steps for Guo-HMC

## ---- 2. MAP and sanity checks ------------------------------------------------

MAP <- map_estimate(B, y, alpha, nu, sigma, w_truth)

cbind(summary(MAP), summary(w_truth))
log_pi(w_truth, y, B, nu, alpha, sigma)
log_pi(MAP,     y, B, nu, alpha, sigma)

# start a little off the MAP to avoid zero-gradient issues 
w_start <- MAP + rnorm(length(MAP), 0, 0.01)

# post_var_diag from loading marginal variances
precond_diag <- post_var_diag


## ---- 4. Production runs with the estimated preconditioner --------------------
cat("--- pHMC ---\n")
eps_p    <- 0.04
phmc_time <- system.time(phmc_run <- phmc_cpp(B, y,
                                              lambda = lambda_prox, alpha = alpha, sigma = sigma,
                                              iter   = iter, eps_hmc = eps_p, L = L_px, nu = nu,
                                              start  = w_start, precond = precond_diag, blather = TRUE))[3]

cat("\n--- Guo-HMC ---\n")
eps_guo    <- 0.0004
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
         col = c("black", "blue", "darkgreen", "red"), lty = 1)
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
         col = c("black", "blue", "darkgreen", "red"), lty = 1)
}
## Per-coordinate trace / density spot-checks
i <- 1
acf_all(i)
den_all(i)


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




