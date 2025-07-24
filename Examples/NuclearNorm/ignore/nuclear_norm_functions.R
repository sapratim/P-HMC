library(mcmcse)
library(Matrix)
library(expm)
library(foreach)
library(doParallel)
library(ks)

## There are also Rcpp versions of these functions that we eventually used

##############------Functions------##################################

nucl_norm <- function(vect)    ## vector input
{
  A <- matrix(vect, nrow = n, ncol = n)
  norm_val <- sum(svd(A)$d)
  return(norm_val)
}

log_pi <- function(x,y,sigma2,alpha)
{
  n_norm <- nucl_norm(x)
  dens_val <- alpha*n_norm + sum((y - x)^2)/(2*sigma2)
  return(-dens_val)
}

log_pilambda <- function(eta,x,lambda,y,sigma2,alpha)   # log target function of pi^lambda
{
  n_norm <- nucl_norm(eta)
  dens_val <- alpha*n_norm + sum((eta-x)^2)/(2*lambda) + 
    sum((y-eta)^2)/(2*sigma2)
  return(-dens_val)
}

#########  Soft threshold function

softthreshold <- function(u, lambda) {       ####  u is a vector
  return(sign(u)*sapply(u, FUN=function(x) {max(abs(x)-lambda,0)}))
}

######### Proximity mapping

prox_func <- function(x,lambda,y,sigma2,alpha) {   #### input x and y as a vector
  num <- lambda*y + sigma2*x
  denom <- lambda + sigma2
  mat <- matrix(num/denom, nrow = n, ncol = n)
  svdsol <- svd(mat)
  s <- softthreshold(svdsol$d, (alpha*sigma2*lambda)/denom)
  output <- svdsol$u %*% (s*t(svdsol$v))  # Multiply each row of V’ by singular values
  return(vec(output))
}

dur_prox_func <- function(x,lambda, alpha) {   #### input x and y as a vector
  mat <- matrix(x, nrow = n, ncol = n)
  svdsol <- svd(mat)
  #s <- softthreshold(svdsol$d, (alpha*sigma2*lambda)/denom)
  output <- svdsol$u %*% (pmax((svdsol$d - lambda*alpha), 0)*t(svdsol$v))  # Multiply each row of V’ by singular values
  return(vec(output))
}

grad_logpiLam <- function(x,lambda,y,sigma2,alpha)  # gradient of log target
{
  x_prox <- prox_func(x,lambda,y,sigma2,alpha)
  ans <-  (x-x_prox)/lambda
  return(-ans)
}

grad_log_durpiLam <- function(x,lambda,y,sigma2,alpha)  # gradient of log target
{
  x_prox <- dur_prox_func(x,lambda, alpha)
  term2 <-  -(x-x_prox)/lambda
  term1 <- (y - x)/(sigma2)
  return(term1 + term2)
}

##  MYHMC samples

durhmc <- function(y, alpha, lambda, sigma2, iter, eps_hmc, L, start)
{
  nvar <- length(y)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  
  # starting value computations
  samp <- start
  samp.hmc[1,] <- samp
  
  # For HMC
  mom_mat <- matrix(rnorm(iter*length(y)), nrow = iter, ncol = length(y))
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_samp <- -grad_log_durpiLam(samp, lambda,y,sigma2,alpha)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -grad_log_durpiLam(samp, lambda,y,sigma2,alpha)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- sum((y - q_current)^2)/(2*sigma2) + alpha*nucl_norm(q_current)
    U_prop <- sum((y - samp)^2)/(2*sigma2) + alpha*nucl_norm(samp)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- samp
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i,] <- q_current
      samp <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/i
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, acc_rate)
  return(object)
}

## pxhmc samples

pxhmc <- function(y, alpha, lambda, sigma2, iter, eps_hmc, L, start)
{
  nvar <- length(y)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  
  # starting value computations
  samp <- start
  samp.hmc[1,] <- samp
  
  # For HMC
  mom_mat <- matrix(rnorm(iter*length(y)), nrow = iter, ncol = length(y))
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_samp <- -grad_logpiLam(samp, lambda,y,sigma2,alpha)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -grad_logpiLam(samp, lambda,y,sigma2,alpha)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- sum((y - q_current)^2)/(2*sigma2) + alpha*nucl_norm(q_current)
    U_prop <- sum((y - samp)^2)/(2*sigma2) + alpha*nucl_norm(samp)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- samp
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i,] <- q_current
      samp <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/i
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, acc_rate)
  return(object)
}