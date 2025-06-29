set.seed(1)
source("cox_functions.R")
load("estimated-cov.RData")
load("cox-data.RData")
library(foreach) 
library(doParallel)
library(mcmcse)

cov <- xn[[5]]
m <- xn[[6]]

# #The posterior
cov.svd <- svd(cov)
sqrt.cov <- cov.svd$u %*% diag(cov.svd$d^(1/2), m) %*% t(cov.svd$v)
max(cov - sqrt.cov %*% sqrt.cov)
inv.cov <- qr.solve(cov)

n <- 2e5 #Length of the chain
eta_bf <- 0.004
eta_mh <- 0.005

#extracting data
ns <- xn[[1]]
x <- xn[[2]]
c <- xn[[3]]
t <- xn[[4]]
cov <- xn[[5]]
m <- xn[[6]]
delta_m <- xn[[7]]
N0 <- xn[[8]]
mu <- xn[[9]]



output_ram <- list()
num_cores <- 50
doParallel::registerDoParallel(cores = num_cores)

#Number of repetitions
reps <- 100

output_cox <- foreach(b = 1:reps) %dopar% {
  bf_time <- system.time(bf <- cox_bf(n, init = rep(1, m), ns = ns, x, c, t, cov, eta_bf))
  mh_time <- system.time(mh <- cox_mh(n, init = rep(1, m), ns = ns, x, c, t, cov, eta_mh))
  bf_chain <- bf[[1]]
  mh_chain <- mh[[1]]
  
  bf_loops_avg <- mean(bf[[2]], round = 2)
  bf_loops_max <- max(bf[[2]], round = 2)
  
  bf_multi_ess <- multiESS(bf_chain)
  mh_multi_ess <- multiESS(mh_chain)
  
  bf_time <- bf_time[3]
  mh_time <- mh_time[3]
  
  bf_ess <- ess(bf_chain)
  mh_ess <- ess(mh_chain)
  
  print(paste('Replication:', b))
  
  list(bf_time, mh_time, bf_loops_avg, bf_loops_max, bf_multi_ess, mh_multi_ess, bf_ess, mh_ess)
}


save(output_cox, file = "output_Cox.RData")

