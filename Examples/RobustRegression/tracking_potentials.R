######################################
## Tracking Hamiltonians
######################################
set.seed(1)
library(Rcpp)
sourceCpp("robustreg_functions.cpp")
source("robustreg_data.R")

##############  Give some starting value of sampler

####################################
# checking the hamiltonian
####################################
## Durmus
Leap_pHMC <- function(samp, B, p_prop, eps_hmc, L, lambda)
{
  U_samp <- -grad_logpiLamg(B, y, samp, lambda, alpha = alpha, nu)
  p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
  q_current <- samp
  for (j in 1:L)
  {
    samp <- samp + eps_hmc*p_current   # full step for position
    U_samp <- -grad_logpiLamg(B, y, samp, lambda, alpha = alpha, nu)
    if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
  }
  p_current <- p_current - eps_hmc*U_samp/2
  p_current <- - p_current  # negation to make proposal symmetric
  
  potential <- log_pi(samp, y, B, nu, alpha) + sum(dnorm(p_current, log = TRUE))
  return(potential)
}

####################################
# tracking potentials

# random place
set.seed(3)
p_prop <- rnorm(dim(Phi_mat)[2])
samp <- w_truth 
alpha <- 100

############################
# Choosing lambda
# x = 1
lambda.seq <- seq(1e-3, 2e-2, length = 1e3)

phmc_ham <- numeric(length = length(lambda.seq))
ns_ham <- numeric(length = length(lambda.seq))

potential <- log_pi(samp, y, Phi_mat, nu, alpha) + sum(dnorm(p_prop, log = TRUE))


for(i in 1:length(lambda.seq))
{
  phmc_state <- Leap_pHMC(samp, Phi_mat, p_prop, eps_hmc = 1e-7, L = 1, lambda = lambda.seq[i])
  phmc_ham[i] <- abs(potential - phmc_state)/abs(potential)
}

plot(lambda.seq, phmc_ham, type = 'l', lwd = 2,
     xlab = expression(lambda[g]), ylab = expression(R[lambda[g]]),
     main = "")
#abline(v = lambda.seq[which.min(phmc_ham)], lwd = 2, lty = 2)
legend("bottomright", legend = expression("Choice of " * lambda[g]),
       col = c( "black"), lty = c(2), lwd = 2, bty = "n")


pdf("Output/lambda_robustreg.pdf", height = 3.5, width = 4.8)
plot(lambda.seq, phmc_ham, type = 'l', lwd = 2,
     xlab = expression(lambda[g]), ylab = expression(R[lambda[g]]),
     main = "")
abline(v = lambda.seq[which.min(phmc_ham)], lwd = 2, lty = 2)
legend("bottomright", legend = expression("Choice of " * lambda[g]),
       col = c( "black"), lty = c(2), lwd = 2, bty = "n")
dev.off()


# ignore 
# 
# L <- 10
# lamb.vec <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10)
# eps.vec <-  seq(1e-5, 1e-3, length = 20)
# 
# truth <- log_pi(x, y, samp) + sum(dnorm(p_prop, log = TRUE))
# ham_dur <- matrix(0, nrow = length(eps.vec), ncol = length(lamb.vec))
# ham_cha <- matrix(0, nrow = length(eps.vec), ncol = length(lamb.vec))
# for(i in 1:length(eps.vec))
# {
#   for(j in 1:length(lamb.vec))
#   {
#     ham_dur[i, j] <- Leap_pHMC(samp, p_prop, eps_hmc = eps.vec[i], L = L, lambda = lamb.vec[j])
#     ham_cha[i, j] <- Leap_Chaari(samp, p_prop, eps_hmc = eps.vec[i], L = L, lambda = lamb.vec[j])
#   }
# }
# 
# par(mfrow = c(2,3))
# for(i in 1:min(6,length(eps.vec)))
# {
#   plot(log(lamb.vec), ham_dur[i, ], type = 'l', 
#        ylim = range(c(truth, ham_dur[i,], ham_cha[i,])),
#        main = paste("Eps = ", eps.vec[i]), col = "blue")
#   lines(log(lamb.vec), ham_cha[i, ], col = "red")
#   abline(h = truth, col = "black")
# }
# 
# par(mfrow = c(2,3))
# for(i in 1:length(lamb.vec))
# {
#   plot(log(eps.vec), ham_dur[,i], type = 'l', 
#        ylim = range(c(truth, ham_dur[,i], ham_cha[,i])),
#        main = paste("Lambda = ", lamb.vec[i]), col = "blue")
#   lines(log(eps.vec), ham_cha[,i], col = "red")
#   abline(h = truth, col = "black")
# }
# 
