set.seed(1)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("hmc_cpp.cpp")
source("cox_functions.R")
load("estimated-cov.RData")
load("cox-data.RData")

ns <- xn[[1]]
x <- xn[[2]]
c <- xn[[3]]
t <- xn[[4]]
cov <- xn[[5]]
m <- xn[[6]]
delta_m <- xn[[7]]
N0 <- xn[[8]]
mu <- xn[[9]]


#The posterior
cov.svd <- svd(cov)
sqrt.cov <- cov.svd$u %*% diag(cov.svd$d^(1/2), m) %*% t(cov.svd$v)
max(cov - sqrt.cov %*% sqrt.cov)
inv.cov <- qr.solve(cov)

#Running MCMC using exact proposal
N <- 1e5

eta_bf <- 0.003 #step size
bf_chain <- cox_bf_cpp(N, init = rep(1,  m), ns = ns, x, c, t, cov, eta_bf, sqrt.cov, m, delta_m)
bf_chain[[3]]

# eta_mh <- 0.005
# mh_chain <- cox_mh(N, init = rep(1, m), ns = ns, x, c, t, cov, eta_mh)
# mh_chain[[2]]



#hmc_chain <- cox_hmc(N, init = rep(1, m), ns = ns, x, c, t, cov, lambda = .001, eps_hmc = .002, L = 10)
hmc_chain <- cox_hmc_cpp(1e5, init = bf_chain[[1]][N, ], ns = ns, x, c, t, sqrt_cov = sqrt.cov, 
	inv_cov = inv.cov,  lambda = .000001, eps_hmc = .005, L = 10, delta_m = delta_m)

#hmc_chain2 <- cox_hmc_cpp
#mh_chain[[2]]
plot.ts(hmc_chain[[1]][,100])
acf(hmc_chain[[1]][,100])
acf(bf_chain[[1]][,100])
library(mcmcse)
cbind(ess(bf_chain[[1]]), ess(hmc_chain[[1]]))


i <- 100
plot(density((bf_chain[[1]][, i]) ) )
lines(density((mh_chain[[1]][, i]) ), col = "red")
lines(density((hmc_chain[[1]][, i]) ), col = "blue")	

#save(bf_chain, mh_chain, file = "output_cox_single_run.RData")


