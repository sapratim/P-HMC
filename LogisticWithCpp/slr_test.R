#################################################
# File contains main functions
# for Sparse logistic regression Example
#################################################

library(mcmcse)
library(Matrix)
library(ks)
library(stats)
library(fasta)
library(glmnet)
library(Rcpp)
sourceCpp("fasta.cpp")
##############------Functions------##################################

log_pi <- function(x,y,beta)
{
  f <- sum(log(1+exp(x %*% beta)) - y*(x%*%beta))
  g <- alpha*sum(abs(beta))
  return(-(f+g))
}

#########  Soft threshold function

softthreshold <- function(u, pen) {       ####  u is a vector
  return(sign(u)*sapply(u, FUN=function(x) {max(abs(x)-pen,0)}))
}

######### Proximity mapping functions for Chaari

# f <- function(z) {colSums(log(1 + exp(x%*%z)) - y*(x%*%z)) + sum((beta_point - z)^2)/(2*lamb)}
# 
# gradf <- function(z) {colSums(c(1/(1+exp(-x%*%z)) - y)*x) + (beta_point - z)/lamb}
# 
# g <- function(z) {alpha*sum(abs(z))}
# 
# proxg <- function(z, tau_fasta) {softthreshold(z, alpha*tau_fasta)}

######### gradient of log target

grad_logpiLam <- function(x, y, beta, lambda, alpha, fasta_start, fasta_step_start)  
{
  temp_mode <- rcpp_fasta(x = x, y = y, x0 = freq_mode,
                          beta_point = beta, alpha = alpha,
                          lamb = lambda, tau1 = fasta_step_start)
  temp_curr <- rcpp_fasta(x, y, c(fasta_start), beta, alpha, lambda, fasta_step_start)
  temp_curr$objective <- temp_curr$objective[temp_curr$objective!= 0]
  temp_mode$objective <- temp_mode$objective[temp_mode$objective!= 0]
  beta_prox <- `if`(min(temp_curr$objective) <= min(temp_mode$objective), temp_curr$x, temp_mode$x)
  ans <-  (beta-beta_prox)/lambda
  return(-ans)
}

######### Proximity mapping for Durmus

prox_func_dur <- function(beta, lambda) {   #### input x and y as a vector
  proxval <- softthreshold(beta, alpha*lambda)
  return(proxval)
}

grad_logpiLam_dur <- function(x, y, beta, lambda)  # gradient of log target for Durmus
{
  beta_prox <- prox_func_dur(beta, lambda)
  grad_f <- gradf_dur(beta, x, y)
  ans <-  grad_f + (beta-beta_prox)/lambda
  return(-ans)
}

##############  Chaari P-HMC  #################################### 

pxhmc_chaari <- function(x, y, lambda, alpha, iter, 
                         eps_hmc, L, start, fasta_start, fasta_step_start)
{
  nvar <- length(start)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  # lamb <<- lambda
  
  # starting value computations
  samp <- start
  samp.hmc[1,] <- samp
  
  # For HMC
  mom_mat <- matrix(rnorm(iter*nvar), nrow = iter, ncol = nvar)
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    beta_point <<- samp
    U_samp <- -grad_logpiLam(x, y, samp, lambda, alpha, beta_point, fasta_step_start)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      beta_point <<- samp
      U_samp <- -grad_logpiLam(x, y, samp, lambda, alpha, beta_point, fasta_step_start)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- -log_pi(x, y, q_current)
    U_prop <- -log_pi(x, y, samp)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- samp
      accept <- accept + 1
    }else
    {
      samp.hmc[i,] <- q_current
      samp <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, acc_rate)
  return(object)
}


######################################## Durmus #########################################


pxhmc_dur <- function(x, y, lambda, iter, eps_hmc, L, start)
{
  nvar <- length(start)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  
  # starting value computations
  samp <- start
  #print(start)
  samp.hmc[1,] <- samp
  
  # For HMC
  mom_mat <- matrix(rnorm(iter*nvar), nrow = iter, ncol = nvar)
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_samp <- -grad_logpiLam_dur(x, y, samp, lambda)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -grad_logpiLam_dur(x, y, samp, lambda)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- -log_pi(x, y, q_current)
    U_prop <- -log_pi(x, y, samp)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- samp
      accept <- accept + 1
    }else
    {
      samp.hmc[i,] <- q_current
      samp <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, acc_rate)
  return(object)
}


#######################################################################################
########################## Sparse logistic regression run #############################
#######################################################################################

################### PIMA Indian diabetes dataset ###################

data <- MASS::Pima.tr
x <- as.matrix(data[,c(1:7)])
y <- as.matrix(ifelse(data$type == "Yes", 1, 0))
colnames(x) <- NULL
colnames(y) <- NULL
alpha <- 2

colnames(data) = c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "y")
logistic_fit <- glmnet(x, y, family = "binomial",
                       alpha = 1, lambda = alpha/length(y), nlambda = 1,
                       standardize = FALSE, intercept = FALSE)$beta

beta <- logistic_fit #c(unlist(logistic_fit$coefficients[-1]))
beta_start <- as.matrix(unname(beta))
freq_mode <<- beta_start



L_pxch <- 10
iter <- 1e4
lamb_coeff <- 10^seq(-7, 20, by = 1)
eps_px_dur <-  0.0019
L_pxdur <- 10
tau <- 5

mcse_multi_chaari <- numeric(length = length(lamb_coeff))
mcse_multi_dur <- numeric(length = length(lamb_coeff))
ess_mat_chaari <- matrix(0, nrow = length(lamb_coeff), ncol = length(beta_start))
ess_mat_dur <- matrix(0, nrow = length(lamb_coeff), ncol = length(beta_start))
acc_chaari <- numeric(length = length(lamb_coeff))
acc_dur <- numeric(length = length(lamb_coeff))

# for(i in 1:length(lamb_coeff)){
i <- 1


####################################
# checking the hamiltonian
####################################
## Durmus
Leap_Durmus <- function(samp, p_prop, eps_hmc, L, lambda)
{
  U_samp <- -grad_logpiLam_dur(x, y, samp, lambda)
  p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
  q_current <- samp
  for (j in 1:L)
  {
    samp <- samp + eps_hmc*p_current   # full step for position
    U_samp <- -grad_logpiLam_dur(x, y, samp, lambda)
    if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
  }
  p_current <- p_current - eps_hmc*U_samp/2
  p_current <- - p_current  # negation to make proposal symmetric
  
  potential <- log_pi(x, y, samp)
  return(potential)
}

## Chaari
Leap_Chaari <- function(samp, p_prop, eps_hmc, L, lambda)
{
  beta_point <<- samp
  U_samp <- -grad_logpiLam(x, y, samp, lambda, alpha, beta_point, fasta_step_start)
  p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
  q_current <- samp
  for (j in 1:L)
  {
    samp <- samp + eps_hmc*p_current   # full step for position
    beta_point <<- samp
    U_samp <- -grad_logpiLam(x, y, samp, lambda, alpha, beta_point, fasta_step_start)
    if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
  }
  p_current <- p_current - eps_hmc*U_samp/2
  p_current <- - p_current  # negation to make proposal symmetric
  
  potential <- log_pi(x, y, samp)
  return(potential)
}

p_prop <- rnorm(length(samp))


lambda <- 100
samp <- rnorm(dim(x)[2])
log_pi(x, y, samp)
Leap_Durmus(samp, p_prop, eps_hmc = 1e-4, L = 10, lambda = 1e-6)
Leap_Durmus(samp, p_prop, eps_hmc = 1e-4, L = 10, lambda = 1e1)
Leap_Chaari(samp, p_prop, eps_hmc = 1e-4, L = 10, lambda = 1e-6)
Leap_Chaari(samp, p_prop, eps_hmc = 1e-4, L = 10, lambda = 2)


eps_px_chaari <- 1e-5
system.time(pxhmc_chaari_run <- pxhmc_chaari(x, y, lambda = .001, alpha = alpha, iter = iter,
                                             eps_hmc = eps_px_chaari, L=L_pxch, start = beta_start,
                                             fasta_start = beta_start, fasta_step_start = tau))

eps_px_chaari <- 1e-5
system.time(pxhmc_chaari_run_alt <- pxhmc_chaari_alt(x, y, lambda = .001, alpha = alpha, iter = 1e2,
                                             eps_hmc = eps_px_chaari, L=L_pxch, start = beta_start,
                                             fasta_start = beta_start, fasta_step_start = tau))

system.time(pxhmc_dur_run <- pxhmc_dur(x, y, lambda = lamb_coeff[i], iter = iter, 
                                       eps_hmc = eps_px_dur, L=L_pxdur, start = beta_start))

mcse_multi_chaari[i] <- multiESS(pxhmc_chaari_run[[1]])
mcse_multi_dur[i] <- multiESS(pxhmc_dur_run[[1]])
acc_chaari[i] <- pxhmc_chaari_run[[2]]
ess_mat_chaari[i,] <- ess(pxhmc_chaari_run[[1]])
ess_mat_dur[i,] <- ess(pxhmc_dur_run[[1]])
acc_dur[i] <- pxhmc_dur_run[[2]]
# }


out <- list(mcse_multi_chaari, mcse_multi_dur, 
            acc_chaari, acc_dur, ess_mat_chaari, ess_mat_dur)

save(out, file = "slr_test_output.Rdata")

# pdf(file = "slr_chaari.pdf")
# par(mfrow = c(2,1))
# plot(out[[1]], type = "o", xlab = "lambda index", ylab = "MultiESS")
# plot(out[[3]], type = "o", xlab = "lambda index", ylab = "acceptance")
# dev.off()
# 
# pdf(file = "slr_dur.pdf")
# par(mfrow = c(2,1))
# plot(out[[2]], type = "o", xlab = "lambda index", ylab = "MultiESS")
# plot(out[[4]], type = "o", xlab = "lambda index", ylab = "acceptance")
# dev.off()



# par(mfrow = c(1,2))
# plot(rowMeans(ess_mat), type = "o")
# plot(mcse_multi, type = "o")
# 
# 
# 
# acf_mat <- matrix(0, nrow = length(L), ncol = length(beta_start))
# lag <- 100
# for (i in 1:length(L)) {
#   pxhmc_dur_run <- pxhmc_dur(y, lambda = 10 , iter = 1e4, eps_hmc = eps_vals[i], L = L[i],
#                              start = freq_mode$x)
#   for (j in 1:length(rand)) {
#     acf_dur <- acf(pxhmc_dur_run[[1]][,rand[j]], plot = FALSE, lag.max = lag)$acf
#     acf_mat[i,j] <- acf_dur[-c(1:lag)]
#   }
# }
# 
# 
# par(mfrow = c(3,3))
# lag <- 100
# 
# for (i in 1:length(beta_start)) {
#  # acf_chaari_hmc <- acf(output[[1]][,i], plot = FALSE, lag.max = lag)$acf
#   acf_dur_hmc <- acf(pxhmc_dur_run[[1]][,i], plot = FALSE, lag.max = lag)$acf
#   plot(1:length(acf_dur_hmc), acf_dur_hmc, col = "blue", type = 'l', ylim = c(-0.1, 1),
#        xlab = "Lag", ylab = "Autocorrelation", main = paste("ACF plot for component",i))
#   #lines(1:length(acf_dur_hmc), acf_dur_hmc, col = "blue", type = 'l')
#   #legend("bottomleft", c("Chaari", "Durmus"), lty = 1,
#      #    col = c("red", "blue"), cex = 0.7, bty = "n")
# }
# 


# i <- 1
# plot(density(pxhmc_dur_run[[1]][,i]))
# abline(v=colMeans(pxhmc_dur_run[[1]])[i], col = "red")
# i <- i+1
# 
# library(SimTools)
# cbind(colMeans(pxhmc_chaari_run[[1]]), colMeans(pxhmc_dur_run[[1]]))
# plot(as.Smcmc(pxhmc_chaari_run[[1]]), which = 1:4)
# plot(as.Smcmc(pxhmc_dur_run[[1]]), which = 5:7)
# 
# dim <- length(beta_start)
# rand <- 1:dim
# 
# pdf("sparse_log_reg_acf.pdf", height = 6, width = 6)
# 
# lag.max <- 100
# acf_chaari_hmc <- acf(pxhmc_chaari_run[[1]][,rand[1]], plot = FALSE, lag.max = lag.max)$acf
# acf_dur_hmc <- acf(pxhmc_dur_run[[1]][,rand[1]], plot = FALSE, lag.max = lag.max)$acf
# 
# diff.acf <- matrix(0, ncol = dim, nrow = lag.max + 1)
# diff.acf[,1] <- acf_dur_hmc - acf_chaari_hmc
# 
# for (i in 2:dim) 
# {
#   acf_chaari_hmc <- acf(pxhmc_chaari_run[[1]][,rand[i]], plot = FALSE, lag.max = lag.max)$acf
#   acf_dur_hmc <- acf(pxhmc_dur_run[[1]][,rand[i]], plot = FALSE, lag.max = lag.max)$acf
#   diff.acf[,i] <- acf_dur_hmc - acf_chaari_hmc
# }
# 
# # Make boxplot of ACFs
# boxplot(t(diff.acf),
#         xlab = "Lags", col = "pink",
#         ylab = "Difference in ACFs of HMCs",ylim = range(diff.acf),
#         names = 0:lag.max, show.names = TRUE, range = 3)
# dev.off()
