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
  if(min(temp_curr$objective) <= min(temp_mode$objective)){
    beta_prox <- temp_curr$x
  }
  else{
    beta_prox <- temp_mode$x
  }
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

##### Chaari P-HMC

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


rwm_alg <- function(y, start, iter, h)
{
  nvar <- length(start)
  samp.rwm <- matrix(0, nrow = iter, ncol = nvar)
  
  # starting value computations
  samp <- start
  samp.rwm[1,] <- samp
  
  # count number of acceptances
  accept <- 0        
  gaussian_mat <- matrix(rnorm((iter-1)*nvar), nrow = (iter-1), ncol = nvar)
  
  for (i in 2:iter) 
  {
    propval <- samp.rwm[i-1,] + h*gaussian_mat[i-1,]
    log_ratio <- log_pi(x, y, propval) - log_pi(x, y, samp.rwm[i-1,])
    if(log(runif(1)) < log_ratio)
    {
      samp.rwm[i,] <- propval
      accept <- accept + 1
    } else{
      samp.rwm[i,] <- samp.rwm[i-1,]
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))}
  }
  print(accept/iter)
  return(samp.rwm)
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

iter <- 1e3
lamb_coeff <- 1e-1
eps_px_chaari <- 13e-5
eps_px_dur <-  0.0019
L_pxch <- 10
L_pxdur <- 10
tau <- 5

system.time(pxhmc_chaari_run <- pxhmc_chaari(x, y, lambda = lamb_coeff, alpha = alpha, iter = iter,
                                             eps_hmc = eps_px_chaari, L=L_pxch, start = beta_start, 
                                             fasta_start = beta_start, fasta_step_start = tau))

system.time(pxhmc_dur_run <- pxhmc_dur(x, y, lambda = lamb_coeff, iter = iter, 
                                       eps_hmc = eps_px_dur, L=L_pxdur, start = beta_start))

output <- list(pxhmc_chaari_run[[1]], pxhmc_dur_run[[1]])

save(output, file = "chains.Rdata")

