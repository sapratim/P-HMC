######################################################################
################# Sparse logistic regression   #######################
######################################################################

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
fit <- glmnet::glmnet(x, y, family = "gaussian", alpha = 1, intercept = FALSE,
                      standardize = FALSE, lambda = alpha/length(y), nlambda = 1)$beta
beta_start <- as.matrix(unname(fit))

L_px <- 10
iter <- 1e6
eps_p <-  0.0002
sigma2 <- 0.1

parallel::detectCores()
num_cores <- 4
doParallel::registerDoParallel(cores = num_cores)
reps <- 4

output_rreg_truth <- foreach(b = 1:reps) %dopar% 
  {
    ## Run samplers
    print(b)
    
    if(runif(1) < 0.05) {L_px <- 1}
    phmc_time <- system.time(phmc_run <- phmc_cpp(x, y, lambda = .01, iter = iter, eps_hmc = eps_p, 
                   L = L_px, sigma2 = sigma2, start = beta_start, alpha = alpha, blather = FALSE))
    
    post_mean <- colMeans(phmc_run[[1]])
    list(post_mean)
  }

save(output_rreg_truth, file = "Output/outputrregtruth.Rdata")
