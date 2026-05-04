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
logistic_fit <- glmnet(x, y, family = "binomial",
                       alpha = 1, lambda = alpha/length(y), nlambda = 1,
                       standardize = FALSE, intercept = FALSE)$beta

beta_start <- as.matrix(unname(logistic_fit ))

L_px <- 10
iter <- 1e6
eps_p <-  0.00192


#output_poisson <- list()
parallel::detectCores()
num_cores <- 4
doParallel::registerDoParallel(cores = num_cores)
reps <- 8

output_slog_truth <- foreach(b = 1:reps) %dopar% 
  {
    ## Run samplers
    print(b)
    
    if(runif(1) < 0.05) {L_px <- 1}
    phmc_time <- system.time(phmc_run <- phmc_cpp(x, y, lambda = .01, iter = iter, eps_hmc = eps_p, 
                                                  L = L_px, start = beta_start, alpha = alpha, blather = FALSE) )
    
    post_mean <- colMeans(phmc_run[[1]])

    list(post_mean)
  }

save(output_slog_truth, file = "Output/outputslogtruth.Rdata")
