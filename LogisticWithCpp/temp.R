library(Rcpp)
library(RcppArmadillo)
library(fasta)
library(glmnet)
sourceCpp("fasta.cpp")

### R functions   ###
softthreshold <- function(u, pen) {       ####  u is a vector
  return(sign(u)*sapply(u, FUN=function(x) {max(abs(x)-pen,0)}))
}
f <- function(z) {colSums(log(1 + exp(x%*%z)) - y*(x%*%z)) + sum((beta_point - z)^2)/(2*lamb)}
gradf <- function(z) {colSums(c(1/(1+exp(-x%*%z)) - y)*x) + (beta_point - z)/lamb}
g <- function(z) {alpha*sum(abs(z))}
proxg <- function(z, tau_fasta) {softthreshold(z, alpha*tau_fasta)}

###########
# Dataset and params
data <- MASS::Pima.tr
x <- as.matrix(data[,c(1:7)])
y <- as.matrix(ifelse(data$type == "Yes", 1, 0))
colnames(x) <- NULL
colnames(y) <- NULL
alpha <- 1
lamb <- 2
###########

tau <- seq(0.5, 20, length = 40)

alpha <- 3

beta_point <- rnorm(7)
initial_step <- sample(tau, 1)
shrink_step <- sample(seq(0.1, 0.9, length = 20), 1)
dum <- fasta(f, gradf, g, proxg, x0 = rep(0,7), tau1 = initial_step, stepsizeShrink = shrink_step)
foo <- rcpp_fasta(x = x, y = y, x0 = rep(0,7), tau1 = initial_step,
           beta_point = beta_point, alpha = alpha,
           lamb = lamb, stepsizeShrink = shrink_step)

# checking solution
cbind(dum$x, foo$x, dum$x-foo$x)



# checking how fast
library(rbenchmark)
benchmark(fasta(f, gradf, g, proxg, x0 = rep(0,7), tau1 = 1, stepsizeShrink = .1),
          rcpp_fasta(x = x, y = y, x0 = rep(0,7), tau1 = 1,
                     beta_point = beta_point, alpha = alpha,
                     lamb = lamb, stepsizeShrink = .1),
          replications = 100)


