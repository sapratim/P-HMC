#################################################
# File contains main functions
# for Sparse logistic regression Example
#################################################

##############------Functions------##################################

# log of the target
log_pi <- function(x,y,beta)
{
  eta <- x %*% beta
  f <- sum(y * (log1p(exp(-eta))) + (1 - y) * (eta + log1p(exp(-eta))))
  g <- alpha*sum(abs(beta))
  return(-(f+g))
}

#########  Soft threshold function

softthreshold <- function(u, pen) {       ####  u is a vector
  return(sign(u)*sapply(u, FUN=function(x) {max(abs(x)-pen,0)}))
}


######### gradient of log target

grad_logpiLam <- function(x, y, beta, lambda, alpha)  
{
  fista <- calc_prox_fista(x = beta, X = x, y = y, lambda = lambda, alpha = alpha,
                   step_size = .1)
  beta_prox <- fista
  ans <-  -(beta-beta_prox)/lambda
  return(ans)
}

######### Proximity mapping for Durmus

prox_phmc <- function(beta, lambda) {   #### input x and y as a vector
  proxval <- softthreshold(beta, alpha*lambda)
  return(proxval)
}

# gradient of log target with partial proximal approximation
grad_logpiLam_dur <- function(x, y, beta, lambda)  
{
  beta_prox <- prox_phmc(beta, lambda)
  grad_f <- gradf_dur(beta, x, y)
  ans <-  grad_f + (beta-beta_prox)/lambda
  return(-ans)
}


##############  Chaari ns-HMC  #################################### 

nshmc <- function(x, y, lambda, alpha, iter, 
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
    U_samp <- -grad_logpiLam(x, y, samp, lambda, alpha)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -grad_logpiLam(x, y, samp, lambda, alpha)
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



#################  Proposed Proximal HMC  ###################
phmc <- function(x, y, lambda, iter, eps_hmc, L, start)
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


#################  Random Walk  ###################

rwm <- function(x, y, iter, h, start)
{
  nvar <- length(start)
  samp.rwm <- matrix(0, nrow = iter, ncol = nvar)
  
  # starting value computations
  samp <- start
  #print(start)
  samp.rwm[1,] <- samp
  
  # For HMC
  accept <- 0
  
  for (i in 2:iter) 
  {
    prop <- samp.rwm[i-1, ] + rnorm(nvar, 0, sd = h)
    log.ratio <- log_pi(x, y, prop) - log_pi(x, y, samp.rwm[i-1, ])
    
    if(log(runif(1)) <= log.ratio)
    {
      samp.rwm[i,] <- prop
      accept <- accept + 1
    }else
    {
      samp.rwm[i,] <- samp.rwm[i-1, ]
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.rwm, acc_rate)
  return(object)
}


