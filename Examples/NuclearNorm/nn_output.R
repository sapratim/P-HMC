############################################################################
###################### Output visualisation for ############################
##################### nuclear norm matrix estimation #######################
############################################################################

load("Output/outputnn.Rdata")
load("Output/outputnn_true.Rdata")

output_nnorm[[1]]

all_ess <- lapply(output_nnorm, function(t) t[[2]])
avg_ess <- Reduce("+", all_ess)/length(output_nnorm)
avg_ess <- round(avg_ess, 5)

all_time <- lapply(output_nnorm, function(t) t[[3]])
avg_time <- Reduce("+", all_time)/length(output_nnorm)
#avg_time 

# ESS/sec
ESS_time <- t(round(apply(avg_ess, 1, function(t) t/avg_time), 3))


min <- apply(ESS_time, 2, min)
median <- apply(ESS_time, 2, median)
max <- apply(ESS_time, 2, max)

cbind(min, median, max)


######################################
### Single run
######################################

load("Output/nn_single_image.Rdata")
library(imager)

cred_interval_trunc <- 0*(cred_interval < .22) + 1*(cred_interval > .22)
#cred_mat <- matrix(cred_interval_trunc, nrow = 64, ncol = 64)
cred_matrix <- matrix(cred_interval, nrow = 64, ncol = 64)

pdf("Output/nn_image.pdf", height = 2, width = 7)
par(mfrow = c(1, 4),
    mar = c(1, 1, 1, 1),   # bottom, left, top, right
    oma = c(0, 0, 0, 0)) 
par(mfrow = c(1,4))
plot(as.cimg(checker), axes = FALSE)
plot(as.cimg(image_mat), axes = FALSE)
plot(as.cimg(MAP1), axes = FALSE)
#plot(as.cimg(cred_mat), axes = FALSE)
plot(as.cimg(cred_matrix), axes = FALSE)
# image(cred_matrix,
#       col = gray(seq(0, 1, length = 256)),
#       axes = FALSE,
# )
dev.off()


############ acfs #########

lag.max <- 100
pdf("Output/nn_acf.pdf", height = 3, width = 4.5)
par(mar = c(5, 4, 2, 2))

plot(0:lag.max, rep(1, lag.max + 1), type = "n", ylim = c(-.02, 1),
     ylab = "Estimated autocorrelations", xlab = "Lags")
for(i in 1:4096)
{
  lines(0:lag.max, ns_acfs[,i], lwd = 1.5, col = "orange")
  lines(0:lag.max, p_acfs[,i], lwd = 1.5, col = "purple", lty = 2)
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

alg_means <- lapply(output_nnorm, '[[', 1)    #### list of posterior means for all reps
true_means <- do.call(cbind, lapply(output_nn_true, '[[', 1))  #### actual means (dim x reps)
truth_estimate <- rowMeans(true_means)

names <- c("RWM", "pHMC", "myMALA", "nsHMC", "pMALA", "guoHMC")

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
nsHMC_output <- metric_fun("nsHMC")
pMALA_output <- metric_fun("pMALA")
guoHMC_output <- metric_fun("guoHMC")

output_mat <- matrix(0, nrow = length(names), ncol = 3)

output_mat[1,] <- sapply(rwm_output, rbind)
output_mat[2,] <- sapply(pHMC_output, rbind)
output_mat[3,] <- sapply(myMALA_output, rbind)
output_mat[4,] <- sapply(nsHMC_output, rbind)
output_mat[5,] <- sapply(pMALA_output, rbind)
output_mat[6,] <- sapply(guoHMC_output, rbind)

rownames(output_mat) <- c("RWM", "pHMC", "myMALA", "nsHMC", "pMALA", "guoHMC")
colnames(output_mat) <- c("Mean Error", "Avg MSE", "Standard Error")

round(output_mat, 4)



# ###################### Comparison metrics  ######################
# 
# alg_means <- lapply(output_nn, '[[', 1)    #### list of posterior means for all reps
# true_means <- do.call(cbind, lapply(output_nn, '[[', 1))  #### actual means (dim x reps)
# truth_estimate <- rowMeans(true_means)

# rwm_means <- sapply(alg_means, function(x) x[,1])
# pHMC_means <- sapply(alg_means, function(x) x[,2])
# myMALA_means <- sapply(alg_means, function(x) x[,3])
# nsHMC_means <- sapply(alg_means, function(x) x[,4])
# pMALA_means <- sapply(alg_means, function(x) x[,5])
# 
# 
# ############################  Relative error  ############################
# 
# mean_devs_rwm <- rwm_means - truth_estimate
# mean_devs_phmc  <- pHMC_means - truth_estimate  
# mean_devs_myMALA  <- myMALA_means - truth_estimate
# mean_devs_nsHMC  <- nsHMC_means - truth_estimate
# mean_devs_pMALA  <- pMALA_means - truth_estimate
#  
# ########### Mean relative error ############
#   
# mean_re_rwm <-  mean(sqrt(colSums(mean_devs_rwm^2)) / sqrt(sum(truth_estimate^2)))
# mean_re_pHMC <-  mean(sqrt(colSums(mean_devs_pHMC^2)) / sqrt(sum(truth_estimate^2)))
# mean_re_myMALA <-  mean(sqrt(colSums(mean_devs_myMALA^2)) / sqrt(sum(truth_estimate^2)))
# mean_re_nsHMC <-  mean(sqrt(colSums(mean_devs_nsHMC^2)) / sqrt(sum(truth_estimate^2)))
# mean_re_pMALA <-  mean(sqrt(colSums(mean_devs_pMALA^2)) / sqrt(sum(truth_estimate^2)))
# 
# ########### Maximum relative error ############
# 
# max_re_rwm <-  max(sqrt(colSums(mean_devs_rwm^2)) / sqrt(sum(truth_estimate^2)))
# max_re_pHMC <-  max(sqrt(colSums(mean_devs_pHMC^2)) / sqrt(sum(truth_estimate^2)))
# max_re_myMALA <-  max(sqrt(colSums(mean_devs_myMALA^2)) / sqrt(sum(truth_estimate^2)))
# max_re_nsHMC <-  max(sqrt(colSums(mean_devs_nsHMC^2)) / sqrt(sum(truth_estimate^2)))
# max_re_pMALA <-  max(sqrt(colSums(mean_devs_pMALA^2)) / sqrt(sum(truth_estimate^2)))
# 
# ############  Pooled relative error ############
# 
# pooled_re_rwm <- sqrt((rowMeans(rwm_means) - truth_estimate)^2)/sqrt(sum(truth_estimate^2))
# pooled_re_pHMC <- sqrt((rowMeans(pHMC_means) - truth_estimate)^2)/sqrt(sum(truth_estimate^2))
# pooled_re_myMALA <- sqrt((rowMeans(myMALA_means) - truth_estimate)^2)/sqrt(sum(truth_estimate^2))
# pooled_re_nsHMC <- sqrt((rowMeans(nsHMC_means) - truth_estimate)^2)/sqrt(sum(truth_estimate^2))
# pooled_re_pMALA <- sqrt((rowMeans(pMALA_means) - truth_estimate)^2)/sqrt(sum(truth_estimate^2))
# 
# ########################## Average Mean Squared error ###########################
# 
# avg_mse_rwm <- mean(colSums(mean_devs_rwm^2))
# avg_mse_pHMC <- mean(colSums(mean_devs_pHMC^2))
# avg_mse_myMALA <- mean(colSums(mean_devs_myMALA^2))
# avg_mse_nsHMC <- mean(colSums(mean_devs_nsHMC^2))
# avg_mse_pMALA <- mean(colSums(mean_devs_pMALA^2))

