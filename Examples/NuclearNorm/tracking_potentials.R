######################################
## Tracking Hamiltonians
######################################
library(Rcpp)
library(RcppArmadillo)

set.seed(1)
source("nn_data.R")
sourceCpp("nn_functions.cpp")
load("warmup_chain.Rdata")

# 
# 
nuclear_prox_sigma <- function(Y, alpha, sigma2) {
  svd_Y <- svd(Y)
  lambda <- sigma2 * alpha
  d_thresh <- pmax(svd_Y$d - lambda, 0)
  return(svd_Y$u %*% diag(d_thresh) %*% t(svd_Y$v))
}

MAP <- nuclear_prox_sigma(Y = image_mat, alpha = alpha_hat, sigma2 = sigma2_hat)
vec.MAP <- vec(MAP)
# plot(as.cimg(image_mat))
# plot(as.cimg(MAP))


####################################
# checking the hamiltonian
####################################
## Our method
Leap_pHMC <- function(samp, p_prop, eps_hmc, L, lambda)
{
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
  
  potential <- sum((y - samp)^2)/(2*sigma2) + alpha*nucl_norm(samp) - sum(dnorm(p_current, log = TRUE))
  return(potential)
}

## Chaari
Leap_ns <- function(samp, p_prop, eps_hmc, L, lambda)
{
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
  
  potential <- sum((y - samp)^2)/(2*sigma2) + alpha*nucl_norm(samp) - sum(dnorm(p_current, log = TRUE))
  return(potential)
}

####################################
# tracking potentials

# random place
set.seed(3)
p_prop <- rnorm(length(warmup_end_iter))
samp <- vec.MAP 

sigma2 = sigma2_hat
alpha = alpha_hat
############################
# Choosing lambda
# x = 1
lambda.seq <- seq(1e-6, 3e-4, length = 1e2)

phmc_ham <- numeric(length = length(lambda.seq))
ns_ham <- numeric(length = length(lambda.seq))

potential <- sum((y - samp)^2)/(2*sigma2) + alpha*nucl_norm(samp) - sum(dnorm(p_prop, log = TRUE))


for(i in 1:length(lambda.seq))
{
  print(i)
  phmc_state <- Leap_pHMC(samp, p_prop, eps_hmc = 1e-7, L = 1, lambda = lambda.seq[i])
  ns_state <- Leap_ns(samp, p_prop, eps_hmc = 1e-7, L= 1, lambda = lambda.seq[i])

  phmc_ham[i] <- abs(potential - phmc_state)/abs(potential)
  ns_ham[i] <- abs(potential - ns_state)/abs(potential)
}

pdf("Output/lambda_nn.pdf", height = 2.5, width = 4)
par(
  mar = c(3.2, 3.5, 1.0, 1.5),  # bottom, left, top, right margins
  mgp = c(2.0, 1.0, 0)          # axis title, axis labels, axis line
)

plot(
  lambda.seq, phmc_ham,
  type = "l", lwd = 2,
  xlab = expression(lambda[g]),
  ylab = expression(R^{lambda[g]}),
  main = ""
)

abline(v = 1e-4, lwd = 2, lty = 2)

legend(
  "bottomright",
  legend = expression("Choice of " * lambda[g]),
  lty = 2, lwd = 2,
  bty = "n"
)

dev.off()

