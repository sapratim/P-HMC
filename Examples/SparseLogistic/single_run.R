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

# starting values
logistic_fit <- glmnet(x, y, family = "binomial",
                      alpha = 1, lambda = alpha/length(y), nlambda = 1,
                      standardize = FALSE, intercept = FALSE)$beta
  
beta_start <- as.matrix(unname(logistic_fit ))
freq_mode <<- beta_start


L_ns <- 10
L_px <- 10

iter <- 1e4
eps_px <-  0.00015
eps_ns <- 0.00014



system.time(nshmc_run <- nshmc_cpp(x, y, lambda = 1e-4, alpha = alpha, iter = iter,
	eps_hmc = eps_px, L = L_px, start = beta_start) )

eps <-  0.0018
system.time(pmala_run <- nshmc_cpp(x, y, lambda = eps/2, alpha = alpha, iter = iter,
                                   eps_hmc = eps, L = 1, start = beta_start) )

phmc_time<- system.time(phmc_run <- phmc_cpp(x, y, lambda = .01, iter = iter, 
	eps_hmc = 0.002, L = L_px, start = beta_start, alpha = alpha) )

eps <-  0.002
system.time(mymala_run <- phmc_cpp(x, y, lambda = eps/2, alpha = alpha, iter = iter,
                                   eps_hmc = eps, L = 1, start = beta_start) )

rwm_time <- system.time(rwm_run <- rwm_cpp(x, y, iter = iter, 
    h = .005 ,start = beta_start, alpha = alpha) )


ESS <- round(cbind(ess(rwm_run[[1]]), ess(phmc_run[[1]]), 
                   ess(mymala_run[[1]]), ess(nshmc_run[[1]]), 
                   ess(pmala_run[[1]])), 0)
colnames(ESS) <- c("RWM", "pHMC", "MYMala", "nsHMC", "PMala")
ESS


library(SimTools)
plot(as.Smcmc(nshmc_run[[1]][ ,5:7]) )
plot(as.Smcmc(phmc_run[[1]][,5:7]) )
plot(as.Smcmc(rwm_run[[1]][,5:7]) )

cbind(colMeans(nshmc_run[[1]]), colMeans(phmc_run[[1]]), colMeans(rwm_run[[1]]) )
ess(rwm_run[[1]])
ess(phmc_run[[1]])


cbind(ess(rwm_run[[1]])/rwm_time[3], ess(phmc_run[[1]])/phmc_time[3] )




