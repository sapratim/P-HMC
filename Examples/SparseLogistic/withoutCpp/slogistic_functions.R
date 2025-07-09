#################################################
# File contains main functions
# for Sparse logistic regression Example
#################################################


##############------Functions------##################################

# log of the target
log_pi <- function(x, y, beta)
{
  f <- sum(log(1+exp(x %*% beta)) - y*(x%*%beta))
  g <- alpha*sum(abs(beta))
  return(-(f+g))
}

#########  Soft threshold function
softthreshold <- function(u, pen) {       ####  u is a vector
  return(sign(u)*sapply(u, FUN=function(x) {max(abs(x)-pen,0)}))
}

######### Proximity mapping functions for Chaari et. al. 
## needed for fasta function
f <- function(z) 
{
  colSums(log(1 + exp(x%*%z)) - y*(x%*%z)) + sum((beta_point - z)^2)/(2*lamb)
}

gradf <- function(z) 
{
  colSums(c(1/(1 + exp(-x%*%z)) - y)*x) + (beta_point - z)/lamb
}

g <- function(z) alpha*sum(abs(z))
proxg <- function(z, tau_fasta) softthreshold(z, alpha*tau_fasta)



######### gradient of log target_lambda (Chaari)
grad_logpiLam <- function(beta, lambda, f, gradf, g, proxg, fasta_start, fasta_step_start)  
{
  temp <- fasta(f, gradf, g, proxg, c(beta), fasta_step_start, stepsizeShrink = .1, max_iters = 50)
 # print(length(temp$objective))
  beta_prox <- temp$x
  ans <-  (beta-beta_prox)/lambda
  return(-ans)
}

######### Proximity mapping for Partial Proximal
prox_func_dur <- function(beta, lambda) 
{   
  proxval <- softthreshold(beta, alpha*lambda)
  return(proxval)
}

# gradient of log target U^lambda_g
grad_logpiLam_g <- function(x, y, beta, lambda)  
{
  beta_prox <- prox_func_dur(beta, lambda)
  grad_f <- colSums(c(1/(1 + exp(-x%*%beta)) - y)*x)
  ans <-  grad_f + (beta - beta_prox)/lambda
  return(-ans)
}

##### Chaari ns-HMC

nshmc_chaari <- function(x, y, lambda, iter, eps_hmc, L, start, fasta_start, fasta_step_start)
{
  nvar <- length(start)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  lamb <<- lambda
  
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
    U_samp <- -grad_logpiLam(samp, lambda, f, gradf, g, proxg, fasta_start, fasta_step_start)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      beta_point <<- samp
      U_samp <- -grad_logpiLam(samp, lambda, f, gradf, g, proxg, fasta_start, fasta_step_start)
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
    U_samp <- -grad_logpiLam_g(x, y, samp, lambda)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -grad_logpiLam_g(x, y, samp, lambda)
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



