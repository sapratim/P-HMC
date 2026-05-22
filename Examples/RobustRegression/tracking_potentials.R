######################################
## Tracking Hamiltonians
######################################
set.seed(1)
library(Rcpp)
sourceCpp("pre_robustreg_functions.cpp")
source("robustreg_data.R")

##############  Give some starting value of sampler

####################################
# checking the hamiltonian
####################################
## Durmus
Leap_pHMC <- function(samp, B, p_prop, eps_hmc, L, lambda, sigma)
{
  U_samp <- -grad_logpiLamg(B, y, samp, lambda, alpha = alpha, sigma = sigma, nu)
  p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
  q_current <- samp
  for (j in 1:L)
  {
    samp <- samp + eps_hmc*p_current   # full step for position
    U_samp <- -grad_logpiLamg(B, y, samp, lambda, alpha = alpha, sigma = sigma, nu)
    if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
  }
  p_current <- p_current - eps_hmc*U_samp/2
  p_current <- - p_current  # negation to make proposal symmetric
  
  potential <- log_pi(samp, y, B, nu, alpha, sigma) + sum(dnorm(p_current, log = TRUE))
  return(potential)
}

####################################
# tracking potentials

# random place
set.seed(3)
p_prop <- rnorm(dim(Phi_mat)[2])
samp <- w_truth 

############################
# Choosing lambda
# x = 1
lambda.seq <- seq(2e-6, .005, length = 1e3)

phmc_ham <- numeric(length = length(lambda.seq))
ns_ham <- numeric(length = length(lambda.seq))

potential <- log_pi(samp, y, Phi_mat, nu, alpha, sigma) + sum(dnorm(p_prop, log = TRUE))


for(i in 1:length(lambda.seq))
{
  phmc_state <- Leap_pHMC(samp, Phi_mat, p_prop, eps_hmc = 1e-7, L = 1, lambda = lambda.seq[i], sigma = sigma)
  phmc_ham[i] <- abs(potential - phmc_state)/abs(potential)
}

# plot(lambda.seq, phmc_ham, type = 'l', lwd = 2,
#      xlab = expression(lambda[g]), ylab = expression(R[lambda[g]]),
#      main = "")
# legend("bottomright", legend = expression("Choice of " * lambda[g]),
#        col = c( "black"), lty = c(2), lwd = 2, bty = "n")
# abline(v = .002, lty = 2)

pdf("Output/lambda_robustreg.pdf", height = 3.5, width = 4.8)
plot(lambda.seq, phmc_ham, type = 'l', lwd = 2,
     xlab = expression(lambda[g]), ylab = expression(R[lambda[g]]),
     main = "")
abline(v = .002, lwd = 2, lty = 2)
legend("bottomright", legend = expression("Choice of " * lambda[g]),
       col = c( "black"), lty = c(2), lwd = 2, bty = "n")
dev.off()

