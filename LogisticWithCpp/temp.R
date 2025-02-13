library(Rcpp)
library(RcppArmadillo)
library(fasta)
library(glmnet)
sourceCpp("fasta.cpp")
load("chains.Rdata")

### R functions   ###
softthreshold <- function(u, pen) {       ####  u is a vector
  return(sign(u)*sapply(u, FUN=function(x) {max(abs(x)-pen,0)}))
}
f <- function(z) {colSums(log(1 + exp(x%*%z)) - y*(x%*%z)) + sum((beta_point - z)^2)/(2*lamb)}
gradf <- function(z) {colSums(c(1/(1+exp(-x%*%z)) - y)*x) + (beta_point - z)/lamb}
g <- function(z) {alpha*sum(abs(z))}
proxg <- function(z, tau_fasta) {softthreshold(z, alpha*tau_fasta)}



##########################################################################################
##########################################################################################
##########################################################################################

data <- MASS::Pima.tr
x <- as.matrix(data[,c(1:7)])
y <- as.matrix(ifelse(data$type == "Yes", 1, 0))
colnames(x) <- NULL
colnames(y) <- NULL
alpha <- 2

colnames(data) = c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "y")
logistic_fit <- glmnet(x, y, family = "binomial",
                       alpha = 1, lambda = alpha/length(y), nlambda = 1,
                       standardize = FALSE, intercept = FALSE)$beta

beta <- logistic_fit #c(unlist(logistic_fit$coefficients[-1]))
beta_start <- as.matrix(unname(beta))


lamb <- 1e-2
rand <- sample(c(1:5e4), 1)
beta_point <- output[[1]][rand,]
alpha <- 2
stepshrink <- 0.1

dum <- fasta(f, gradf, g, proxg, x0 = beta_point, tau1 = 5, stepsizeShrink = stepshrink)
foo <- rcpp_fasta(x = x, y = y, x0 = beta_point,
           beta_point = beta_point, alpha = alpha,
           lamb = lamb, tau1 = 5, stepsizeShrink = stepshrink)

# checking solution
cbind(dum$x, foo$x, dum$x-foo$x, beta_start)

##########################################################################################
##########################################################################################
##########################################################################################


# checking how fast
library(rbenchmark)
benchmark(fasta(f, gradf, g, proxg, x0 = rep(0,7), tau1 = 1, stepsizeShrink = .1),
          rcpp_fasta(x = x, y = y, x0 = rep(0,7), tau1 = 1,
                     beta_point = beta_point, alpha = alpha,
                     lamb = lamb, stepsizeShrink = .1),
          replications = 100)


