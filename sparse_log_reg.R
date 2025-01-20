###############################
# File contains main functions
# for Nuclear Norm Example
###############################

library(mcmcse)
library(Matrix)
library(expm)
library(ks)

##############------Functions------##################################

log_pi <- function(x,y,beta)
{
  f <- sum(log(1+exp(x %*% beta)) - y*(x%*%beta))
  g <- sum(abs(beta))
  return(-(f+g))
}

#########  Soft threshold function

softthreshold <- function(u, lambda) {       ####  u is a vector
  return(sign(u)*sapply(u, FUN=function(x) {max(abs(x)-lambda,0)}))
}

######### Proximity mapping for Durmus

prox_func_dur <- function(beta, lambda) {   #### input x and y as a vector
  proxval <- softthreshold(beta, lambda)
  return(proxval)
}

grad_logpiLam_dur <- function(x, y, beta, lambda)  # gradient of log target for Durmus
{
  beta_prox <- prox_func_dur(beta, lambda)
  exp_term <- exp(x%*%beta)
  grad_f_coeff <- exp_term/(1+exp_term) - y
  grad_f <- t(colSums(grad_f_coeff*x))
  ans <-  grad_f + (beta-beta_prox)/lambda
  return(-ans)
}

############### Durmus 


pxhmc_dur <- function(x, y, lambda, iter, eps_hmc, L, start)
{
  nvar <- length(start)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  
  # starting value computations
  samp <- start
  print(start)
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
    }
    else
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
########################## Nuclear norm run ###########################################
#######################################################################################

data <- MASS::Pima.tr
x <- data[,c(1:7)]
y <- ifelse(data$type == "Yes", 1, 0)
colnames(x) <- NULL
colnames(y) <- NULL

iter <- 1e2
lamb_coeff <- 1e-4
eps_px_dur <-  0.008
L <- 10

system.time(pxhmc_dur_run <- pxhmc_dur(as.matrix(x), as.matrix(y), lambda = lamb_coeff, iter = iter, 
                                       eps_hmc = eps_px_dur, L=L, start = colMeans(x)))

dim <- length(y)
rand <- 1:dim

pdf("nn_acf_HMC_PervsDur.pdf", height = 6, width = 6)

lag.max <- 100
acf_per_hmc <- acf(pxhmc_per_run[[1]][,rand[1]], plot = FALSE, lag.max = lag.max)$acf
acf_dur_hmc <- acf(pxhmc_dur_run[[1]][,rand[1]], plot = FALSE, lag.max = lag.max)$acf

diff.acf <- matrix(0, ncol = dim, nrow = lag.max + 1)
diff.acf[,1] <- acf_dur_hmc - acf_per_hmc

for (i in 2:dim) 
{
  acf_per_hmc <- acf(pxhmc_per_run[[1]][,rand[i]], plot = FALSE, lag.max = lag.max)$acf
  acf_dur_hmc <- acf(pxhmc_dur_run[[1]][,rand[i]], plot = FALSE, lag.max = lag.max)$acf
  diff.acf[,i] <- acf_dur_hmc - acf_per_hmc
}

# Make boxplot of ACFs
boxplot(t(diff.acf),
        xlab = "Lags", col = "pink",
        ylab = "Difference in ACFs of HMCs",ylim = range(diff.acf),
        names = 0:lag.max, show.names = TRUE, range = 3)
dev.off()
