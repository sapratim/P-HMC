
load("single_run.Rdata")

durhmc <- output_single_hmc[[1]]
pxhmc <- output_single_hmc[[2]]

library(mcmcse)
ess_dur <- ess(durhmc)
ess_px <- ess(pxhmc)

hist(ess_dur, col = adjustcolor("purple", alpha.f = .3), xlim = c(0, max(ess_dur)))
hist(ess_px, add = TRUE, col = adjustcolor("red", alpha.f = .3))

plot(ess_dur, ess_px, ylim = c(0, max(ess_dur)), xlim = c(0, max(ess_dur)))
abline(a = 0, b = 1)



rand <- 1:length(y) #sample(c(1:length(y)), subset)
dim <- length(y)
lag.max <- 100
acf_dur_hmc <- acf(output_single_hmc[[1]][,rand[1]], plot = FALSE, lag.max = lag.max)$acf
acf_pxhmc <- acf(output_single_hmc[[2]][,rand[1]], plot = FALSE, lag.max = lag.max)$acf

diff.acf <- matrix(0, ncol = dim, nrow = lag.max + 1)
diff.acf[,1] <- acf_is_hmc - acf_pxhmc


acf(result_durhmc[[1]][, 2])
acf(result_pxhmc[[1]][, 2])

for (i in 1:dim) 
{
  if(i %% 1000 == 0) print(i)
  acf_dur_hmc <- acf(output_single_hmc[[1]][,i], plot = FALSE, lag.max = lag.max)$acf
  acf_pxhmc <- acf(output_single_hmc[[2]][,i], plot = FALSE, lag.max = lag.max)$acf
  diff.acf[,i] <- acf_dur_hmc - acf_pxhmc
}


# Make boxplot of ACFs
boxplot(t(diff.acf),
        xlab = "Lags", col = "pink",
        ylab = "Difference in ACFs of HMCs",ylim = range(diff.acf),
        names = 0:lag.max, show.names = TRUE, range = 3)
#dev.off()



# save(output_single_mala, file = "output_single_chain_mala.Rdata")
# save(output_single_bark, file = "output_single_chain_bark.Rdata")
# save(output_single_hmc, file = "output_single_chain_hmc.Rdata")

