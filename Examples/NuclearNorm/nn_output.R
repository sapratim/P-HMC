############################################################################
###################### Output visualisation for ############################
##################### sparse logistic regression ###########################
############################################################################

load("Output/outputnn.Rdata")


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

cbind(min, median, max)


######################################
### Single run
######################################

load("Output/nn_single_image.Rdata")
library(imager)

cred_interval_trunc <- 0*(cred_interval < .22) + 1*(cred_interval > .22)
cred_mat <- matrix(cred_interval_trunc, nrow = 64, ncol = 64)


pdf("Output/nn_image.pdf", height = 2, width = 7)
par(mfrow = c(1, 4),
    mar = c(1, 1, 1, 1),   # bottom, left, top, right
    oma = c(0, 0, 0, 0)) 
par(mfrow = c(1,4))
plot(as.cimg(checker), axes = FALSE)
plot(as.cimg(image_mat), axes = FALSE)
plot(as.cimg(MAP1), axes = FALSE)
plot(as.cimg(cred_mat), axes = FALSE)
dev.off()


############ acfs #########


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
       legend = c("nsHMC", "pHMC"),
       col = c("orange", "purple"),
       lty = c(1,2),
       lwd = 1.5,
       horiz = TRUE, 
       bty = "n", 
       inset = c(0, -0.19),  # pushes legend into the top margin
       xpd = TRUE)               # allow drawing outside plot region
dev.off()

# library(mcmcse)
# library(glmnet)
# load("chains.Rdata")
# load("rwm.Rdata")
# 
# data <- MASS::Pima.tr
# x <- as.matrix(data[,c(1:7)])
# y <- as.matrix(ifelse(data$type == "Yes", 1, 0))
# colnames(x) <- NULL
# colnames(y) <- NULL
# alpha <- 2
# 
# colnames(data) = c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "y")
# logistic_fit <- glmnet(x, y, family = "binomial",
#                        alpha = 1, lambda = alpha/length(y), nlambda = 1,
#                        standardize = FALSE, intercept = FALSE)$beta
# 
# beta <- logistic_fit 
# beta_start <- as.matrix(unname(beta))  ## starting value
# 
# ##################  Density plots  #####################
# 
# pdf("density_plot.pdf", width = 12, height = 10)
# 
# par(mfrow = c(3,3))
# for (i in 1:length(beta_start)) {
#   if (i != 6) {
#     plot(density(output[[1]][,i]), type = 'l', col = "red", main = paste("Component", i))
#     lines(density(output[[2]][,i]), type = 'l', col = "blue")
#     lines(density(rwm_run[,i]), type = 'l', col = "green")
#     abline(v = beta_start[i], col = "black")
#     legend("topright", c("Chaari", "Durmus", "RWM"), lty = 1,
#            col = c("red", "blue", "green"), cex = 1, bty = "n")
#   } else {
#     plot(density(output[[1]][,i]), type = 'l', col = "red", main = paste("Component", i),
#          xlim = c(-0.6, 2))
#     lines(density(output[[2]][,i]), type = 'l', col = "blue")
#     lines(density(rwm_run[,i]), type = 'l', col = "green")
#     abline(v = beta_start[i], col = "black")
#     legend("topright", c("Chaari", "Durmus", "RWM"), lty = 1,
#            col = c("red", "blue", "green"), cex = 1, bty = "n")
#   }
# }
# dev.off()
# 
# #################  Acf plots  #################
# 
# pdf("acf_plot.pdf", width = 12, height = 10)
# 
# par(mfrow = c(3,3))
# lag <- 100
# 
# for (i in 1:length(beta_start)) {
#   acf_chaari_hmc <- acf(output[[1]][,i], plot = FALSE, lag.max = lag)$acf
#   acf_dur_hmc <- acf(output[[2]][,i], plot = FALSE, lag.max = lag)$acf
#   acf_rwm <- acf(rwm_run[,i], plot = FALSE, lag.max = lag)$acf
#   plot(1:length(acf_chaari_hmc), acf_chaari_hmc, col = "red", type = 'l', ylim = c(-0.1, 1),
#        xlab = "Lag", ylab = "Autocorrelation", main = paste("ACF plot for component",i))
#   lines(1:length(acf_dur_hmc), acf_dur_hmc, col = "blue", type = 'l')
#   lines(1:length(acf_rwm), acf_rwm, col = "green", type = 'l')
#   legend("bottomleft", c("Chaari", "Durmus", "RWM"), lty = 1,
#          col = c("red", "blue", "green"), cex = 0.7, bty = "n")
# }
# dev.off()
# 
# ##################  ESS evaluation  #####################
# 
# ess_uni <- cbind(ess(output[[1]]), ess(output[[2]]))
# ess_multi <- cbind(multiESS(output[[1]]), multiESS(output[[2]]))
# ess_output <- list(ess_uni, ess_multi)
# mcse_output <- cbind(mcse.mat(output[[1]]), mcse.mat(output[[2]]))
# save(ess_output, file = "slr_ess.Rdata")
