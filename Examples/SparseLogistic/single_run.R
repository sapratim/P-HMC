######################################################################
########## Sparse logistic regression Single run  ####################
######################################################################
library(mcmcse)
library(Matrix)
library(ks)
library(stats)
library(fasta)
library(glmnet)
source("slogistic_data.R")
source("slogistic_functions.R")


# starting values
logistic_fit <- glmnet(x, y, family = "binomial",
                      alpha = 1, lambda = alpha/length(y), nlambda = 1,
                      standardize = FALSE, intercept = FALSE)$beta
  
beta_start <- as.matrix(unname(logistic_fit ))
freq_mode <<- beta_start


L_ns <- 10
L_px <- 10

iter <- 1e5
lamb_coeff <- 1e-4
eps_px <-  0.0019
eps_ns <- 0.00014

tau <- 5


nshmc_run <- nshmc(x, y, lambda = lamb_coeff, alpha = alpha, iter = iter,
	eps_hmc = eps_px, L = L_px, start = beta_start,
	fasta_start = beta_start, fasta_step_start = tau)


phmc_time <- system.time(phmc_run <- phmc(x, y, lambda = lamb_coeff, iter = iter, 
	eps_hmc = eps_px, L = L_px, start = beta_start) )

rwm_time <- system.time(rwm_run <- rwm(x, y, iter = iter, h = .005 ,start = beta_start) )

library(SimTools)
plot(as.Smcmc(rwm_run[[1]][,5:7]) )
plot(as.Smcmc(phmc_run[[1]][,5:7]) )

ess(rwm_run[[1]])
ess(phmc_run[[1]])

cbind(ess(rwm_run[[1]])/rwm_time[3], ess(phmc_run[[1]])/phmc_time[3] )




