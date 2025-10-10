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

pdf("Output/lambda_nn.pdf", height = 3.5, width = 4.2)
plot(lambda.seq, phmc_ham, type = 'l', lwd = 2,
     xlab = expression(lambda[g]), ylab = expression(R[lambda[g]]), 
     main = "")
abline(v = 1e-4, lwd = 2, lty = 2)
legend("bottomright", legend = expression("Choice of " * lambda[g]), 
       col = c( "black"), lty = c(2), lwd = 2, bty = "n")
dev.off()


L <- 10
lamb.vec <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10)
eps.vec <-  seq(1e-5, 1e-3, length = 20)

truth <- log_pi(x, y, samp) + sum(dnorm(p_prop, log = TRUE))
ham_dur <- matrix(0, nrow = length(eps.vec), ncol = length(lamb.vec))
ham_cha <- matrix(0, nrow = length(eps.vec), ncol = length(lamb.vec))
for(i in 1:length(eps.vec))
{
  for(j in 1:length(lamb.vec))
  {
    ham_dur[i, j] <- Leap_pHMC(samp, p_prop, eps_hmc = eps.vec[i], L = L, lambda = lamb.vec[j])
    ham_cha[i, j] <- Leap_ns(samp, p_prop, eps_hmc = eps.vec[i], L = L, lambda = lamb.vec[j])
  }
}

par(mfrow = c(2,3))
for(i in 1:min(6,length(eps.vec)))
{
  plot(log(lamb.vec), ham_dur[i, ], type = 'l', 
       ylim = range(c(truth, ham_dur[i,], ham_cha[i,])),
       main = paste("Eps = ", eps.vec[i]), col = "blue")
  lines(log(lamb.vec), ham_cha[i, ], col = "red")
  abline(h = truth, col = "black")
}

par(mfrow = c(2,3))
for(i in 1:length(lamb.vec))
{
  plot(log(eps.vec), ham_dur[,i], type = 'l', 
       ylim = range(c(truth, ham_dur[,i], ham_cha[,i])),
       main = paste("Lambda = ", lamb.vec[i]), col = "blue")
  lines(log(eps.vec), ham_cha[,i], col = "red")
  abline(h = truth, col = "black")
}

