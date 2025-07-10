######################################
## Tracking Hamiltonians
######################################
set.seed(1)
library(glmnet)
source("slogistic_functions.R")
source("slogistic_data.R")

### data stuff
logistic_fit <- glmnet(x, y, family = "binomial",
                       alpha = 1, lambda = alpha/length(y), nlambda = 1,
                       standardize = FALSE, intercept = FALSE)$beta

beta <- logistic_fit #c(unlist(logistic_fit$coefficients[-1]))
beta_start <- as.matrix(unname(beta))
freq_mode <<- beta_start


iter <- 1e4
lamb_coeff <- 10^seq(-7, 20, by = 1)
eps_px_dur <-  0.0019
tau <- 5


####################################
# checking the hamiltonian
####################################
## Durmus
Leap_Durmus <- function(samp, p_prop, eps_hmc, L, lambda)
{
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
  
  potential <- log_pi(x, y, samp) + sum(dnorm(p_current, log = TRUE))
  return(potential)
}

## Chaari
Leap_Chaari <- function(samp, p_prop, eps_hmc, L, lambda)
{
  beta_point <<- samp
  U_samp <- -grad_logpiLam(x, y, samp, lambda, alpha, beta_point, tau)
  p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
  q_current <- samp
  for (j in 1:L)
  {
    samp <- samp + eps_hmc*p_current   # full step for position
    beta_point <<- samp
    U_samp <- -grad_logpiLam(x, y, samp, lambda, alpha, beta_point, tau)
    if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
  }
  p_current <- p_current - eps_hmc*U_samp/2
  p_current <- - p_current  # negation to make proposal symmetric
  
  potential <- log_pi(x, y, samp) + sum(dnorm(p_current, log = TRUE))
  return(potential)
}

####################################
# tracking potentials

# random place
set.seed(3)
p_prop <- rnorm(dim(x)[2])
samp <- beta_start #+ rnorm(dim(x)[2], 0, sd = .001)


############################
# Choosing lambda
# x = 1
lambda.seq <- seq(1e-5, .1, length = 1e2)

dur_ham <- numeric(length = length(lambda.seq))
px_ham <- numeric(length = length(lambda.seq))

potential <- log_pi(x, y, samp) + sum(dnorm(p_prop, log = TRUE))


for(i in 1:length(lambda.seq))
{
  dur_state <- Leap_Durmus(samp, p_prop, eps_hmc = 1e-7, L = 10, lambda = lambda.seq[i])
  px_state <- Leap_Chaari(samp, p_prop, eps_hmc = 1e-7, L = 10, lambda = lambda.seq[i])

  dur_ham[i] <- abs(potential - dur_state)/abs(potential)
  px_ham[i] <- abs(potential - px_state)/abs(potential)
}

plot(lambda.seq, px_ham, type = 'l', lwd = 2)
plot(lambda.seq, dur_ham, col = "blue", lwd = 2, type = 'l')



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
    ham_dur[i, j] <- Leap_Durmus(samp, p_prop, eps_hmc = eps.vec[i], L = L, lambda = lamb.vec[j])
    ham_cha[i, j] <- Leap_Chaari(samp, p_prop, eps_hmc = eps.vec[i], L = L, lambda = lamb.vec[j])
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

