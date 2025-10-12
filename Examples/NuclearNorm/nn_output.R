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
