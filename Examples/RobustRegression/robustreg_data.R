################### PIMA Indian diabetes dataset ###################
data <- MASS::Boston 
x <- as.matrix(data[,c(1:13)])
y <- as.matrix(data[,-c(1:13)])
colnames(x) <- NULL
colnames(y) <- NULL
alpha <- 2

# ## Lipschitz coefficient for the smooth part
# xtx <- t(x) %*% x
# 
# C_f <- eigen(xtx)$values[1]/4
# 
# 1/C_f
#
# 
# library(tidyverse)
# 
# # If your data 'df' has columns: Subject, Time, Channel, Amplitude
# eeg_wide <- eegdata %>%
#   pivot_wider(names_from = Channel, values_from = Amplitude)
# 
# # Now your data looks like this:
# # Subject | Time | Fz  | Cz  | Pz 
# # 1       | 100  | -2.5| 1.2 | 4.8
# 
# # Now you can run your Multiple Regression
# model <- lm(TargetVariable ~ Fz + Cz + Pz, data = df_wide)
# 
# ###############
# wide_eeg <- eegdata %>%
#   pivot_wider(
#     id_cols = c(subject, group, condition, trial, time),
#     names_from = channel,
#     values_from = voltage,
#     values_fn = mean
#   )


