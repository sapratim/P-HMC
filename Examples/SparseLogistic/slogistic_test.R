#######################################################################################
########################## Sparse logistic regression run #############################
#######################################################################################

logistic_fit <- glmnet(x, y, family = "binomial",
                       alpha = 1, lambda = alpha/length(y), nlambda = 1,
                       standardize = FALSE, intercept = FALSE)$beta

beta <- logistic_fit #c(unlist(logistic_fit$coefficients[-1]))
beta_start <- as.matrix(unname(beta))
freq_mode <<- beta_start



L_pxch <- 10
iter <- 1e4
lamb_coeff <- 10^seq(-7, 20, by = 1)
eps_px_dur <-  0.0019
L_pxdur <- 10
tau <- 5

mcse_multi_chaari <- numeric(length = length(lamb_coeff))
mcse_multi_dur <- numeric(length = length(lamb_coeff))
ess_mat_chaari <- matrix(0, nrow = length(lamb_coeff), ncol = length(beta_start))
ess_mat_dur <- matrix(0, nrow = length(lamb_coeff), ncol = length(beta_start))
acc_chaari <- numeric(length = length(lamb_coeff))
acc_dur <- numeric(length = length(lamb_coeff))

# for(i in 1:length(lamb_coeff)){
i <- 1


eps_px_chaari <- 1e-5
system.time(pxhmc_chaari_run <- pxhmc_chaari(x, y, lambda = .001, alpha = alpha, iter = iter,
                                             eps_hmc = eps_px_chaari, L=L_pxch, start = beta_start,
                                             fasta_start = beta_start, fasta_step_start = tau))

eps_px_chaari <- 1e-5
system.time(pxhmc_chaari_run_alt <- pxhmc_chaari_alt(x, y, lambda = .001, alpha = alpha, iter = 1e2,
                                             eps_hmc = eps_px_chaari, L=L_pxch, start = beta_start,
                                             fasta_start = beta_start, fasta_step_start = tau))

system.time(pxhmc_dur_run <- pxhmc_dur(x, y, lambda = lamb_coeff[i], iter = iter, 
                                       eps_hmc = eps_px_dur, L=L_pxdur, start = beta_start))

mcse_multi_chaari[i] <- multiESS(pxhmc_chaari_run[[1]])
mcse_multi_dur[i] <- multiESS(pxhmc_dur_run[[1]])
acc_chaari[i] <- pxhmc_chaari_run[[2]]
ess_mat_chaari[i,] <- ess(pxhmc_chaari_run[[1]])
ess_mat_dur[i,] <- ess(pxhmc_dur_run[[1]])
acc_dur[i] <- pxhmc_dur_run[[2]]
# }


out <- list(mcse_multi_chaari, mcse_multi_dur, 
            acc_chaari, acc_dur, ess_mat_chaari, ess_mat_dur)

save(out, file = "slr_test_output.Rdata")

# pdf(file = "slr_chaari.pdf")
# par(mfrow = c(2,1))
# plot(out[[1]], type = "o", xlab = "lambda index", ylab = "MultiESS")
# plot(out[[3]], type = "o", xlab = "lambda index", ylab = "acceptance")
# dev.off()
# 
# pdf(file = "slr_dur.pdf")
# par(mfrow = c(2,1))
# plot(out[[2]], type = "o", xlab = "lambda index", ylab = "MultiESS")
# plot(out[[4]], type = "o", xlab = "lambda index", ylab = "acceptance")
# dev.off()



# par(mfrow = c(1,2))
# plot(rowMeans(ess_mat), type = "o")
# plot(mcse_multi, type = "o")
# 
# 
# 
# acf_mat <- matrix(0, nrow = length(L), ncol = length(beta_start))
# lag <- 100
# for (i in 1:length(L)) {
#   pxhmc_dur_run <- pxhmc_dur(y, lambda = 10 , iter = 1e4, eps_hmc = eps_vals[i], L = L[i],
#                              start = freq_mode$x)
#   for (j in 1:length(rand)) {
#     acf_dur <- acf(pxhmc_dur_run[[1]][,rand[j]], plot = FALSE, lag.max = lag)$acf
#     acf_mat[i,j] <- acf_dur[-c(1:lag)]
#   }
# }
# 
# 
# par(mfrow = c(3,3))
# lag <- 100
# 
# for (i in 1:length(beta_start)) {
#  # acf_chaari_hmc <- acf(output[[1]][,i], plot = FALSE, lag.max = lag)$acf
#   acf_dur_hmc <- acf(pxhmc_dur_run[[1]][,i], plot = FALSE, lag.max = lag)$acf
#   plot(1:length(acf_dur_hmc), acf_dur_hmc, col = "blue", type = 'l', ylim = c(-0.1, 1),
#        xlab = "Lag", ylab = "Autocorrelation", main = paste("ACF plot for component",i))
#   #lines(1:length(acf_dur_hmc), acf_dur_hmc, col = "blue", type = 'l')
#   #legend("bottomleft", c("Chaari", "Durmus"), lty = 1,
#      #    col = c("red", "blue"), cex = 0.7, bty = "n")
# }
# 


# i <- 1
# plot(density(pxhmc_dur_run[[1]][,i]))
# abline(v=colMeans(pxhmc_dur_run[[1]])[i], col = "red")
# i <- i+1
# 
# library(SimTools)
# cbind(colMeans(pxhmc_chaari_run[[1]]), colMeans(pxhmc_dur_run[[1]]))
# plot(as.Smcmc(pxhmc_chaari_run[[1]]), which = 1:4)
# plot(as.Smcmc(pxhmc_dur_run[[1]]), which = 5:7)
# 
# dim <- length(beta_start)
# rand <- 1:dim
# 
# pdf("sparse_log_reg_acf.pdf", height = 6, width = 6)
# 
# lag.max <- 100
# acf_chaari_hmc <- acf(pxhmc_chaari_run[[1]][,rand[1]], plot = FALSE, lag.max = lag.max)$acf
# acf_dur_hmc <- acf(pxhmc_dur_run[[1]][,rand[1]], plot = FALSE, lag.max = lag.max)$acf
# 
# diff.acf <- matrix(0, ncol = dim, nrow = lag.max + 1)
# diff.acf[,1] <- acf_dur_hmc - acf_chaari_hmc
# 
# for (i in 2:dim) 
# {
#   acf_chaari_hmc <- acf(pxhmc_chaari_run[[1]][,rand[i]], plot = FALSE, lag.max = lag.max)$acf
#   acf_dur_hmc <- acf(pxhmc_dur_run[[1]][,rand[i]], plot = FALSE, lag.max = lag.max)$acf
#   diff.acf[,i] <- acf_dur_hmc - acf_chaari_hmc
# }
# 
# # Make boxplot of ACFs
# boxplot(t(diff.acf),
#         xlab = "Lags", col = "pink",
#         ylab = "Difference in ACFs of HMCs",ylim = range(diff.acf),
#         names = 0:lag.max, show.names = TRUE, range = 3)
# dev.off()
