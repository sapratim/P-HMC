############################################################################################
############################# Evaluate optimal lambda ######################################
############################################################################################

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

lambda_grid <- seq(1e-5, 2, length = 100)

L_ns <- 10
iter <- 2e3
eps_ns <- 0.00012
lag_max <- 1
acf_vec <- numeric(length = length(lambda_grid))

for (i in 1:length(lambda_grid)) {
  print(i)
  nshmc_time <- system.time(nshmc_run <- nshmc_cpp(x, y, lambda = lambda_grid[i], 
                    alpha = alpha, iter = iter, eps_hmc = eps_ns, L_val = L_ns, start = beta_start))
  samp_mat <- nshmc_run[[1]]
  acf_lag1 <- sapply(seq_len(ncol(samp_mat)), function(j) {
                          cor(samp_mat[-1, j], samp_mat[-nrow(samp_mat), j])})
  acf_vec[i] <- min(acf_lag1)
}
#opt_lambda <- lambda_grid[which.min(acf_vec)]

save(acf_vec, file = "Output/opt_lambda.Rdata")


###  Benchmark for faster simulation

# function_apply <- function()
# {
# for (i in 1:length(lambda_grid)) {
#   print(i)
#   nshmc_time <- system.time(nshmc_run <- nshmc_cpp(x, y, lambda = lambda_grid[i], 
#                   alpha = alpha, iter = iter, eps_hmc = eps_ns, L = L_ns, start = beta_start))
#   
#   acf_val <- acf(nshmc_run[[1]], plot = FALSE, lag.max = lag_max)$acf
#   term_acf <- apply(acf_val, 1, diag)
#   acf_vec[i] <- colMeans(term_acf)[lag_max+1]
# }
# min_lambda <- which.min(acf_vec)}
# 
# 
# function_means <- function()
# {
# for (i in 1:length(lambda_grid)) {
#  print(i)
#   nshmc_time <- system.time(nshmc_run <- nshmc_cpp(x, y, lambda = lambda_grid[i],
#                    alpha = alpha, iter = iter, eps_hmc = eps_ns, L = L_ns, start = beta_start))
# 
#   acf_val <- lapply(seq_len(ncol(nshmc_run[[1]])), function(j) {
#     acf(nshmc_run[[1]][, j], plot = FALSE, lag.max = lag_max)$acf[lag_max+1,1,1]
#   })
#   acf_vec[i] <- mean(unlist(acf_val))
# }
# min_lambda <- which.min(acf_vec1)}
# # 
# 
# function_lag1 <- function()
# {for (i in 1:length(lambda_grid)) {
#   print(i)
#   nshmc_time <- system.time(nshmc_run <- nshmc_cpp(x, y, lambda = lambda_grid[i], 
#                                                    alpha = alpha, iter = iter, eps_hmc = eps_ns, L = L_ns, start = beta_start))
#   samp_mat <- nshmc_run[[1]]
#   acf_lag1 <- sapply(seq_len(ncol(samp_mat)), function(j) {
#     cor(samp_mat[-1, j], samp_mat[-nrow(samp_mat), j])})
#   acf_vec[i] <- mean(acf_lag1)
# }
#   min_lambda <- which.min(acf_vec)}
# 
#  rbenchmark::benchmark(function_lag1(), function_means(), replications = 5)
# # 
