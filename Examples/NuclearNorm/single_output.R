



# library(mcmcse)
# 
# durhmc <- output_single_hmc[[1]]
# pxhmc <- output_single_hmc[[2]]
# 
# 
# ess_dur <- ess(durhmc)
# ess_px <- ess(pxhmc)
# 
# pdf("nn_hist_ess.pdf", height = 5, width = 7)
# hist(ess_dur, col = adjustcolor("purple", alpha.f = .5), xlim = c(0, max(ess_dur)),
#      xlab = "Effective sample size for all components", main = "")
# hist(ess_px, add = TRUE, col = adjustcolor("orange", alpha.f = .7))
# legend("top", fill = c(adjustcolor("purple", .5), adjustcolor("orange", .7)), 
#        legend = c("Partial Proximal", "Full Proximal"))
# dev.off()