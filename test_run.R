############################################
## Main functions for trendfiltering example
## Plus data generation
############################################

# loading libraries needed
library(mcmcse)
library(glmgen)
library(Matrix)
library(expm)

# generating data
set.seed(12345)
alpha_hat <- 5   # obtained from the first dataset
sigma2_hat <- 9  # obtained from the first dataset
x <- seq(1,100,len=100)
f <- Vectorize(function(x){if(x<=35){x} else if(x<=70){70-x} else{0.5*x-35}})
fx_linear <- f(x)
y <- fx_linear + rnorm(length(x), sd = 3)

# setting parameters for the proximal mappings
tol <- 1e-6
max_iter <- 200L

# function calculates the D matrix of the penalty term
getD <- function(k, n, x=NULL){
  if(is.null(x)){
    x <- 1:n
  }
  diags <- list(rep(-1,n),rep(1,n))
  D <- Matrix::bandSparse(n-1,n,k=c(0,1),diag=diags,symm=F)
  if(k>=1){
    for(i in 1:k){
      leftD <- Matrix::bandSparse(n-i-1,n-i,k=c(0,1),diag=diags,symm=F)
      xdiag <- Matrix::Diagonal(n-i,i/diff(x,lag=i))
      D <- leftD %*% xdiag %*% D
    }
  }
  return(D)
}

log_pi <- function(beta,y,sigma2,alpha)
{
  dens_val <- alpha*(sum(abs(D_mat%*%beta))) + sum((y-beta)^2)/(2*sigma2)
  return(-dens_val)
}

# function calculates the value of the proximal function
prox_func <- function(beta,lambda,alpha,sigma2,k,grid)
{
  betaval <- (beta*sigma2+(lambda*y))/(sigma2+lambda)
  lambdaval <- (alpha)/ ((lambda + sigma2)/ (lambda*sigma2) )
  temp = trendfilter(grid,betaval, k=k,lambda = lambdaval,
                     control = trendfilter.control.list(obj_tol = tol, max_iter = max_iter))$beta
  return(as.vector(temp))
}

# function calculates the value of the proximal function Durmus style
prox_func_dur <- function(beta,lambda,alpha,sigma2,k,grid)
{
  lambdaval <- alpha*lambda
  temp = trendfilter(grid,beta, k=k,lambda = lambdaval,
                     control = trendfilter.control.list(obj_tol = tol, max_iter = max_iter))$beta
  return(as.vector(temp))
}

# gradient of log target (pi-lambda)
grad_logpiLam <- function(beta,lambda,y,sigma2,alpha,k,grid)  
{
  beta_prox <- prox_func(beta,lambda,alpha,sigma2,k,grid)
  ans <-  (beta-beta_prox)/lambda
  return(-ans)
}

# gradient of log target (pi-lambda) Durmus setup
grad_logpiLam_dur <- function(beta,lambda,y,sigma2,alpha,k,grid)  
{
  beta_prox <-prox_func_dur(beta,lambda,alpha,sigma2,k,grid)
  ans <-  - (y - beta)/sigma2 + (beta-beta_prox)/lambda
  return(-ans)
}

## Pereyra's P-HMC samples
pxhmc <- function(y, alpha, sigma2, k, grid, iter, eps_hmc, L, start)
{
  nvar <- length(y)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  lambda <- lamb_coeff
  
  # starting value computations
  beta <- start
  samp.hmc[1,] <- beta
  
  # For HMC
  mom_mat <- matrix(rnorm(iter*nvar), nrow = iter, ncol = nvar)
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_beta <- -grad_logpiLam(beta, lambda,y,sigma2,alpha,k,grid)
    p_current <- p_prop - eps_hmc*U_beta /2  # half step for momentum
    q_current <- beta
    for (j in 1:L)
    {
      beta <- beta + eps_hmc*p_current   # full step for position
      U_beta <- -grad_logpiLam(beta, lambda,y,sigma2,alpha,k,grid)
      if(j!=L) p_current <- p_current - eps_hmc*U_beta  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_beta/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- - log_pi(q_current, y, sigma2, alpha)
    U_prop <- - log_pi(beta, y, sigma2, alpha)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- beta
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i,] <- q_current
      beta <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, acc_rate)
  return(object)
}

## Durmus style P-HMC samples

pxhmc_dur <- function(y, alpha, sigma2, k, grid, iter, eps_hmc, L, start)
{
  nvar <- length(y)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  lambda <- lamb_coeff
  
  # starting value computations
  beta <- start
  samp.hmc[1,] <- beta
  
  # For HMC
  mom_mat <- matrix(rnorm(iter*nvar), nrow = iter, ncol = nvar)
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_beta <- -grad_logpiLam_dur(beta, lambda,y,sigma2,alpha,k,grid)
    p_current <- p_prop - eps_hmc*U_beta /2  # half step for momentum
    q_current <- beta
    for (j in 1:L)
    {
      beta <- beta + eps_hmc*p_current   # full step for position
      U_beta <- -grad_logpiLam_dur(beta, lambda,y,sigma2,alpha,k,grid)
      if(j!=L) p_current <- p_current - eps_hmc*U_beta  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_beta/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- - log_pi(q_current, y, sigma2, alpha)
    U_prop <- - log_pi(beta, y, sigma2, alpha)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- beta
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i,] <- q_current
      beta <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, acc_rate)
  return(object)
}

####################################################################################
############################# Trendfiltering run ###################################
####################################################################################

load("warmup.Rdata")
iter_hmc <- 5e3
lamb_coeff <- 0.001
D_mat <- getD(k=1, n=1e2, x)   #  D matrix
lag.max <- 30
dim <- 100

system.time(px.hmc <- pxhmc(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_hmc,
                eps_hmc = 0.0003, L = 100, start = warmup_end_iter)) 

system.time(px.hmc_dur <- pxhmc_dur(y, alpha_hat,sigma2_hat,k=1, grid=x,iter = iter_hmc,
                eps_hmc = 0.0003, L = 100, start = warmup_end_iter)) 


pdf("tf_acf_HMC_PervsDur.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
acf_hmc_per <- acf(px.hmc[[1]][,1], plot = FALSE, lag.max = lag.max)$acf
acf_hmc_dur <- acf(px.hmc_dur[[1]][,1], plot = FALSE, lag.max = lag.max)$acf

diff.acf <- matrix(0, ncol = dim, nrow = lag.max + 1)
diff.acf[,1] <- acf_hmc_dur - acf_hmc_per 

for (i in 2:100) 
{
  acf_hmc_per <- acf(px.hmc[[1]][,i], plot = FALSE, lag.max = lag.max)$acf
  acf_hmc_dur <- acf(px.hmc_dur[[1]][,i], plot = FALSE, lag.max = lag.max)$acf
  diff.acf[,i] <- acf_hmc_dur - acf_hmc_per 
}

# Make boxplot of ACFs
boxplot(t(diff.acf),
        xlab = "Lags", col = "pink",
        ylab = "Difference in ACF of HMCs",ylim = range(diff.acf),
        names = 0:lag.max, show.names = TRUE)

dev.off()
