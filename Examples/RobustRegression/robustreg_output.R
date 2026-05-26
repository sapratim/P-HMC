##############################################################################
########################## Robust regression output ##########################
##############################################################################

load("Output/outputrreg.Rdata")
load("Output/outputrregtruth.Rdata")
source("robustreg_data.R")

output_rreg[[1]]

all_ess <- lapply(output_rreg, function(t) t[[2]])
avg_ess <- Reduce("+", all_ess)/length(output_rreg)
avg_ess <- round(avg_ess, 5)

all_time <- lapply(output_rreg, function(t) t[[3]])
avg_time <- Reduce("+", all_time)/length(output_rreg)
avg_time 

# ESS/sec
ESS_time <- t(round(apply(avg_ess, 1, function(t) t/avg_time), 3))
min <- apply(ESS_time, 2, min)
median <- apply(ESS_time, 2, median)
max <- apply(ESS_time, 2, max)

# print this table
cbind(min, median, max)


#############################################
### Single run ACF          
#############################################
# load("Output/slog_single.Rdata")
# lag.max <- 100
# 
# pdf("Output/slog_acf.pdf", height = 3, width = 4.5)
# par(mar = c(5, 4, 2, 2))
# 
# plot(0:lag.max, rep(1, lag.max + 1), type = "n", ylim = c(-.02, 1),
#      ylab = "Estimated autocorrelations", xlab = "Lags")
# for(i in 1:7)
# {
#   ns_acfs <- acf(nshmc_run[[1]][ ,i], plot = FALSE, lag.max = lag.max)$acf
#   p_acfs <- acf(phmc_run[[1]][ ,i], plot = FALSE, lag.max = lag.max)$acf
#   lines(0:lag.max, ns_acfs, lwd = 1.5, col = "orange")
#   lines(0:lag.max, p_acfs, lwd = 1.5, col = "purple", lty = 2)
# }
# # Add a legend on top of the plot
# legend("top",
#        legend = c("ns-HMC", "p-HMC"),
#        col = c("orange", "purple"),
#        lty = c(1,2),
#        lwd = 1.5,
#        horiz = TRUE, 
#        bty = "n", 
#        inset = c(0, -0.19),  # pushes legend into the top margin
#        xpd = TRUE)               # allow drawing outside plot region
# dev.off()

########################  Function for metrics  ########################

alg_means <- lapply(output_rreg, '[[', 1)    #### list of posterior means for all reps
true_means <- do.call(cbind, lapply(output_rreg_truth, '[[', 1))  #### actual means (dim x reps)
truth_estimate <- mean(true_means)

names <- c("RWM", "pHMC", "myMALA", "guoHMC")

metric_fun <- function(name)# enter one of ("RWM","pHMC","myMALA","nsHMC","pMALA","guoHMC")
{
  i <- which(names == name)
  
  # posterior mean matrix (dim x reps)
  mat_means <- sapply(alg_means, function(x) x[,i])
  
  # mean deviations for different reps
  mean_devs_mat <- mat_means - truth_estimate
  norm_vec <- colSums(mean_devs_mat^2)
  # errors
  mean_error <- mean(sqrt(norm_vec))  # mean
  avg_mse_mat <- mean(colSums(mean_devs_mat^2))
  stan_error <- sd(norm_vec)/length(norm_vec)
  
  output <- list(mean_error, avg_mse_mat, stan_error)
  return(output)
}

rwm_output <- metric_fun("RWM")
pHMC_output <- metric_fun("pHMC")
myMALA_output <- metric_fun("myMALA")
guoHMC_output <- metric_fun("guoHMC")

output_mat <- matrix(0, nrow = length(names), ncol = 3)

output_mat[1,] <- sapply(rwm_output, rbind)
output_mat[2,] <- sapply(pHMC_output, rbind)
output_mat[3,] <- sapply(myMALA_output, rbind)
output_mat[4,] <- sapply(guoHMC_output, rbind)

rownames(output_mat) <- c("RWM", "pHMC", "myMALA", "guoHMC")
colnames(output_mat) <- c("Mean Error", "Avg MSE", "Standard Error")

round(output_mat, 4)




#################################################################################
#################################################################################
#################################################################################
#################################################################################
###### New output


# acf_all <- function(i)
# {
#   rw <- acf(rwm_run[[1]][,i], plot = FALSE)$acf
#   phmc <- acf(phmc_run[[1]][,i], plot = FALSE)$acf
#   mymala <- acf(mymala_run[[1]][,i], plot = FALSE)$acf
#   guohmc <- acf(guohmc_run[[1]][,i], plot = FALSE)$acf
#   
#   plot(rw, col = "black", type = 'l',
#        ylim = c(0,1))
#   lines(phmc, col = "blue")
#   lines(mymala, col = "darkgreen")
#   lines(guohmc, col = "red")
#   legend("bottomleft", legend = c("RWM", "pHMC", "myMALA", "Guo-HMC"),
#          col = c("black", "blue", "darkgreen", "red"), lty = 1, cex = 0.4)
# }
# 
# den_all <- function(i)
# {
#   
#   plot(density(rwm_run[[1]][,i]), col = "black", 
#        main = paste("Density coord", i), 
#        xlim = range(c(density(rwm_run[[1]][,i])$x, density(phmc_run[[1]][,i])$x,
#                       density(mymala_run[[1]][,i])$x, density(guohmc_run[[1]][,i])$x)))
#   lines(density(phmc_run[[1]][,i]), col = "blue")
#   lines(density(mymala_run[[1]][,i]), col = "darkgreen")
#   lines(density(guohmc_run[[1]][,i]), col = "red")
#   legend("topright", legend = c("RWM", "pHMC", "myMALA", "Guo-HMC"),
#          col = c("black", "blue", "darkgreen", "red"), lty = 1, cex = 0.4)
# }
# ## Per-coordinate trace / density spot-checks
# i <- 1
# acf_all(i)
# den_all(i)
# i <- i+1
# 
# # calculating ess
# ess_phmc <- ess(phmc_run$samples)
# ess_guohmc <- ess(guohmc_run$samples)
# ess_mymala <- ess(mymala_run$samples)
# ess_rwm <- ess(rwm_run$samples)
# 
# summarize_ess <- function(ess_vec)
# {
#   c(mean = mean(ess_vec), median = median(ess_vec), min = min(ess_vec), max = max(ess_vec))
# }
# 
# cat("\nEffective sample size (ESS) summaries:\n")
# print(round(rbind(pHMC = summarize_ess(ess_phmc),
#                   GuoHMC = summarize_ess(ess_guohmc),
#                   myMALA = summarize_ess(ess_mymala),
#                   RWM = summarize_ess(ess_rwm)), 2))
# 
# # ESS per time
# cat("\nESS per second:\n")
# print(round(rbind(pHMC = summarize_ess(ess_phmc / phmc_time),
#                   GuoHMC = summarize_ess(ess_guohmc / guohmc_time),
#                   myMALA = summarize_ess(ess_mymala / mymala_time),
#                   RWM = summarize_ess(ess_rwm / rwm_time)), 4))
# 
# 
