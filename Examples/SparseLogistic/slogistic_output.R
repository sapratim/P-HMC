############################################################################
###################### Output visualisation for ############################
##################### sparse logistic regression ###########################
############################################################################

load("Output/outputslog.Rdata")
load("Output/outputslogtruth.Rdata")

output_slog[[1]]

all_ess <- lapply(output_slog, function(t) t[[1]])
avg_ess <- Reduce("+", all_ess)/length(output_slog)
avg_ess <- round(avg_ess, 5)

all_time <- lapply(output_slog, function(t) t[[2]])
avg_time <- Reduce("+", all_time)/length(output_slog)
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
load("Output/slog_single.Rdata")
lag.max <- 100

pdf("Output/slog_acf.pdf", height = 3, width = 4.5)
par(mar = c(5, 4, 2, 2))

plot(0:lag.max, rep(1, lag.max + 1), type = "n", ylim = c(-.02, 1),
     ylab = "Estimated autocorrelations", xlab = "Lags")
for(i in 1:7)
{
  ns_acfs <- acf(nshmc_run[[1]][ ,i], plot = FALSE, lag.max = lag.max)$acf
  p_acfs <- acf(phmc_run[[1]][ ,i], plot = FALSE, lag.max = lag.max)$acf
  lines(0:lag.max, ns_acfs, lwd = 1.5, col = "orange")
  lines(0:lag.max, p_acfs, lwd = 1.5, col = "purple", lty = 2)
}
# Add a legend on top of the plot
legend("top",
       legend = c("ns-HMC", "p-HMC"),
       col = c("orange", "purple"),
       lty = c(1,2),
       lwd = 1.5,
       horiz = TRUE, 
       bty = "n", 
       inset = c(0, -0.19),  # pushes legend into the top margin
       xpd = TRUE)               # allow drawing outside plot region
dev.off()

########################  Function for metrics  ########################

alg_means <- lapply(output_slog, '[[', 1)    #### list of posterior means for all reps
true_means <- do.call(cbind, lapply(output_slog_truth, '[[', 1))  #### actual means (dim x reps)
truth_estimate <- rowMeans(true_means)

names <- c("RWM", "pHMC", "myMALA", "nsHMC", "pMALA", "guoHMC")

metric_fun <- function(name)# enter one of ("RWM","pHMC","myMALA","nsHMC","pMALA","guoHMC")
{
  i <- which(names == name)
  
  # posterior mean matrix (dim x reps)
  mat_means <- sapply(alg_means, function(x) x[,i])
  
  # mean deviations for different reps
  mean_devs_mat <- mat_means - truth_estimate
  
  # relative error
  mean_re <- mean(sqrt(colSums(mean_devs_mat^2)) / sqrt(sum(truth_estimate^2))) # mean
  max_re <-  max(sqrt(colSums(mean_devs_mat^2)) / sqrt(sum(truth_estimate^2))) # maximum
  pooled_re <- sqrt((rowMeans(mat_means) - truth_estimate)^2)/sqrt(sum(truth_estimate^2)) # pooled
  
  # average mean square error
  avg_mse_mat <- mean(colSums(mean_devs_mat^2))
  
  output <- list(mean_re, max_re, pooled_re, avg_mse_mat)
  return(output)
}

rwm_output <- metric_fun("RWM")
pHMC_output <- metric_fun("pHMC")
myMALA_output <- metric_fun("myMALA")
nsHMC_output <- metric_fun("nsHMC")
pMALA_output <- metric_fun("pMALA")
guoHMC_output <- metric_fun("guoHMC")

