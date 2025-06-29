##########################################
## Code for reproducing the results from
## the Cox processes example
##########################################
load("cox-data.RData")
load("output_Cox.RData")


ns <- xn[[1]]
xs <- xn[[2]] #generated data
m <- xn[[6]] #number of observations
t <- xn[[4]] #grid


N <- 1e6
reps <- length(output_cox)


#To store outputs
bf_loops_avg <- numeric(reps)
bf_loops_max <- numeric(reps)
bf_MultiESS <- numeric(reps)
bf_time <- numeric(reps)
bf_ess <- matrix(nrow = reps, ncol = m)

mh_MultiESS <- numeric(reps)
mh_time <- numeric(reps)
mh_ess <- matrix(nrow = reps, ncol = m)


for(i in 1:reps)
{
  foo <- output_cox[[i]]
  
  bf_time[i] <- foo[[1]]
  mh_time[i] <- foo[[2]]
  bf_loops_avg[i] <- foo[[3]]
  bf_loops_max[i] <- foo[[4]]
  bf_MultiESS[i] <- foo[[5]]
  mh_MultiESS[i] <- foo[[6]]
  bf_ess[i, ] <- foo[[7]]
  mh_ess[i, ] <- foo[[8]]
}

print(paste("Average number of mean loops: ", round(mean(bf_loops_avg), 4)))
print(paste("Average number of max loops for: ", round(mean(bf_loops_max), 4)))

bf_mESS_per_unit_time <- bf_MultiESS/bf_time
mh_mESS_per_unit_time <- mh_MultiESS/mh_time

Multi_ESS_df <- data.frame(
  Method = c("Bernoulli Factory", "Metropolis Hastings"),
  MultiESS = c(round(mean(bf_MultiESS), 0), round(mean(mh_MultiESS), 0)),
  MultiESS_by_time = c(round(mean(bf_mESS_per_unit_time), 4), round(mean(mh_mESS_per_unit_time), 4)),
  Avg_compute_time = c(round(mean(bf_time), 0), round(mean(mh_time), 0))
)

print(Multi_ESS_df)

avg_ess_bf <- round(colMeans(bf_ess), 0)
avg_ess_mh <- round(colMeans(mh_ess), 0)

ESS_df <- data.frame(
  Method = c("Component 1", "Component 2", "Component 3", "Component 4", "Component 5", "Component 6", 
             "Component 7", "Component 8", "Component 9", "Component 10"),
  Bernoulli_ESS = avg_ess_bf,
  Auxiliary_ESS = avg_ess_mh
)

print(ESS_df)

##########################################
# Single chain output (plots)
##########################################

source("cox_functions.R")
load("output_cox_single_run.RData")
#PLOTS
delta_m <- 50/(m-1)
grid <- seq(0, 50, length = 100)
est_fun1 <- numeric(length = length(grid))
est_fun2 <- numeric(length = length(grid))
bf_samps <- bf_chain[[1]]
mh_samps <- mh_chain[[1]]

#log posterior
log_post_bf <- bf_chain[[4]]
log_post_mh <- mh_chain[[3]]


#True density
temp <- seq(0, 50, length = 1e4)
y_temp <- (lam1(temp))



pdf("plots/cox-component-density.pdf")
j <- 100
plot(density(bf_samps[-c(1:1000), j]), col = "blue", ylab = "Estimated Density", xlab = "x", main = "")
lines(density(mh_samps[-c(1:1000), j]), col = "red")
legend("topright", legend = c("Bernoulli factory MCMC", "Inexact Metropolis-Hastings"), col = c("blue", "red"), cex = 1.2, lty = 1, lwd = 2, bty = "n")
dev.off()
