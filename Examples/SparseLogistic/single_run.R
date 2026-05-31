######################################################################
########## Sparse logistic regression Single run  ####################
######################################################################
library(mcmcse)
library(Matrix)
library(ks)
library(stats)
library(fasta)
library(glmnet)
library(Rcpp)
source("slogistic_data.R")
sourceCpp("slogistic_functions.cpp")

set.seed(1)
# starting values
logistic_fit <- glmnet(x, y, family = "binomial",
                      alpha = 1, lambda = alpha/length(y), nlambda = 1,
                      standardize = FALSE, intercept = FALSE)$beta
  
beta_start <- as.matrix(unname(logistic_fit ))
freq_mode <<- beta_start
blat <- TRUE

L_ns <- 10
L_px <- 10
L_guo <- 10
iter <- 1e5
eps_ns <- 0.00012
eps_p <-  0.0019
eps_guo <- 0.00012


nshmc_time <- system.time(nshmc_run <- nshmc_cpp(x, y, lambda = 1, alpha = alpha, iter = iter,
                                                 eps_hmc = eps_ns, L_val = L_ns, start = beta_start, blather = blat))

eps <-  0.0016
pmala_time <- system.time(pmala_run <- nshmc_cpp(x, y, lambda = eps/2, alpha = alpha, iter = iter,
                                                 eps_hmc = eps, L_val = 1, start = beta_start, blather = blat))

phmc_time <- system.time(phmc_run <- phmc_cpp(x, y, lambda = .01, iter = iter, eps_hmc = eps_p, 
                                              L_val = L_px, start = beta_start, alpha = alpha, blather = blat))

guohmc_time <- system.time(guohmc_run <- guohmc_cpp(x, y, lambda = .01, iter = iter, eps_hmc = eps_guo, 
                                                    L_val = L_guo, start = beta_start, alpha = alpha, blather = blat))

eps <-  0.0019
mymala_time <- system.time(mymala_run <- phmc_cpp(x, y, lambda = eps/2, alpha = alpha, iter = iter,
                                                  eps_hmc = eps, L_val = 1, start = beta_start, blather = blat))

rwm_time <- system.time(rwm_run <- rwm_cpp(x, y, iter = iter, 
                                           h = .0045 ,start = beta_start, alpha = alpha, blather = blat))


save(nshmc_run, pmala_run, phmc_run, mymala_run, rwm_run, guohmc_run,
     nshmc_time, pmala_time, mymala_time, rwm_time, guohmc_time, file = "Output/slog_single.Rdata")


