############################################################################
###################### Output visualisation for ############################
##################### sparse logistic regression ###########################
############################################################################

load("Output/outputslog.Rdata")


output_slog[[1]]

all_ess <- lapply(output_slog, function(t) t[[1]])
avg_ess <- Reduce("+", all_ess)/length(output_slog)
avg_ess <- round(avg_ess, 0)

all_time <- lapply(output_slog, function(t) t[[2]])
avg_time <- Reduce("+", all_time)/length(output_slog)
avg_time 

# ESS/sec
t(round(apply(avg_ess, 1, function(t) t/avg_time), 2))


# 
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
