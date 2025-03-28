
source("image_deconv_functions.R")
library(RcppArmadillo)
load("pixel_mat.Rdata")

true_pixel_vec <- c(pixel_mat)
dimen <- sqrt(length(true_pixel_vec))
H <- blur_func(blur_size = 5)    ######  Blur matrix
y <- convolve_image(true_pixel_vec, dimen, dimen, H)
beta_pen <- 0.02   ###### penalty parameter
sigma2 <- 1.76     ###### noise variance
#haar_level <<- 3

############################  DWT/IDWT test  ############################

h_pass <- wave.filter("haar")$hpf     ###### high pass filter
l_pass <- wave.filter("haar")$lpf     ###### low pass filter

img_rand <- rnorm(dimen^2, y, sd = 8)
img_mat <- matrix(img_rand, dimen, dimen)
wave_r <- dwt.2d(img_mat, "haar", J=3)
wave_cpp <- two_D_dwt_multi(img_mat, dimen, dimen, 2, h_pass, l_pass, haar_level)

LL_diff <- max(abs(wave_r$LL3) - abs(wave_cpp$LL))
HL_diff <- cbind(max(wave_r$HL1 - wave_cpp$HL_1), max(wave_r$HL2 - wave_cpp$HL_2),
                 max(wave_r$HL3 - wave_cpp$HL_3))
LH_diff <- cbind(max(wave_r$LH1 - wave_cpp$LH_1), max(wave_r$LH2 - wave_cpp$LH_2),
                 max(wave_r$LH3 - wave_cpp$LH_3))
HH_diff <- cbind(max(wave_r$HH1 - wave_cpp$HH_1), max(wave_r$HH2 - wave_cpp$HH_2),
                 max(wave_r$HH3 - wave_cpp$HH_3))

print(LL_diff)
print(HL_diff)
print(LH_diff)
print(HH_diff)

##########################  Inverse transform  ##########################

inv_trans_R <- idwt.2d(wave_r)
inv_trans_cpp <- two_D_idwt_multi(wave_cpp, haar_level, 2, h_pass, l_pass)
inv_trans_cpp - inv_trans_R
inv_trans_R - img_mat
inv_trans_cpp - img_mat

################  Wavelet function check  ################


wavelet_l1 <- function(image_vec, nlev = 3){
    image_mat <- matrix(image_vec, dimen, dimen)
    trans <- dwt.2d(image_mat, wf = "haar", nlev)
    wave_sum <- sum(abs(unlist(trans)))
    return(wave_sum)
  }

img_rand <- rnorm(dimen^2, y, sd = 8)
wavelet_l1_cpp(img_rand, dimen, h_pass, l_pass) - wavelet_l1(img_rand)   ### small

###################### fasta functions check ######################

f_R <- function(z) {
  H.z <- convolve_image(z, dimen, dimen, H)
  t <- sum((H.z - y)^2)/(2*sigma2) + sum((x_true - z)^2)/(2*lamb)
  return(t)
}

gradf_R <- function(z) {
  H.z <- convolve_image(z, dimen, dimen, H)
  t1 <- convolve_image((H.z - y), dimen, dimen, t(H))/sigma2
  t2 <- (x_true - z)/lamb
  return(t1+t2)
}

g_R <- function(z) {
  Psi.z <- wavelet_l1(z)
  t <- beta_pen*Psi.z
  return(t)
}

proxg_R <- function(z, tau_fasta) {
  z <- matrix(z, dimen, dimen)
  wave_trans <- dwt.2d(z, wf = "haar", 3)
  for (i in 1:(length(wave_trans)-1)) {
    wave_trans[[i]] <- matrix(softthreshold(wave_trans[[i]],
                                            beta_pen*tau_fasta), dimen, dimen)
  }
  proxval <- idwt.2d(wave_trans)
  return(vec(proxval))
}

####################  functions for mode  ####################

f_freq <- function(z) {
  H.z <- convolve_image(z, dimen, dimen, H)
  t <- sum((y - H.z)^2)/(2*sigma2)
  return(t)
}

gradf_freq <- function(z) {
  H.z <- convolve_image(z, dimen, dimen, H)
  t1 <- convolve_image((H.z - y), dimen, dimen, t(H))/sigma2
  return(t1)
}

freq_mode <- fasta(f_freq, gradf_freq, g_R, proxg_R, true_pixel_vec, tau1 = 5, stepsizeShrink = .1)

lamb <- .5
img_rand <- rnorm(dimen^2, y, sd = 8)
x_true <- freq_mode$x

(f_R(img_rand) - f(img_rand, dimen, H, y, sigma2, x_true, lamb))
(gradf_R(img_rand) - gradf(img_rand, dimen, H, y, sigma2, x_true, lamb))
(g_R(img_rand) - g(img_rand, dimen, h_pass, l_pass, beta_pen))
(proxg_R(y,lamb) - proxg(y, dimen, h_pass, l_pass, lamb, beta_pen))


#######################  Fasta check  #######################

# lamb <- 100
# x_true <- freq_mode$x
# start <- rnorm(dimen^2, x_true, sd = .2)
# prox_mode_R <- fasta(f_R, gradf_R, g_R, proxg_R, start, 10, max_iters = 100)
# prox_mode_cpp <- fasta_cpp(start, 10, H, y, h_pass, l_pass, 
#                       sigma2, x_true, lamb, beta_pen, dimen)
# prox_mode_R$x - prox_mode_cpp$x
# plot.ts(prox_mode_R$objective, xlim = c(0, max(length(prox_mode_R$objective), length(prox_mode_cpp$objective))))
# lines(prox_mode_cpp$objective, col = "red")
# 
# 
# prox_mode_R$objective
# prox_mode_cpp$objective

####################### efficacy of fasta_cpp #######################

# lamb <- 150
# x_true <- rnorm(dimen^2, y, sd = 10)
# start <- rnorm(dimen^2, x_true, sd = .1)
# u <- fasta(f_R, gradf_R, g_R, proxg_R, start, 5, max_iters = 100)
# v <- fasta_cpp(start, 5, H, y, h_pass, l_pass,
#                sigma2, x_true, lamb, beta_pen, dimen, max_iters = 100)
# hist(u$x - v$x)
# (psi_R <- -log_pi(y, u$x))
# (psi_cpp <- -log_pi(y, v$x))
# 
# plot.ts(u$objective, xlim = c(0, max(length(u$objective), length(v$objective))))
# lines(v$objective, col = "red")
# 
# u$totalBacktracks
# v$totalBacktracks
# 
# u$objective
# v$objective
# 
# u$taus
# v$taus
# 
# u$residual
# v$residual
######################  Benchmark results  ######################
# 
# library(rbenchmark)
# lamb <- .15
# x_true <- rnorm(dimen^2, y, sd = 10)
# start <- rnorm(dimen^2, x_true, sd = .1)
# benchmark(fasta(f_R, gradf_R, g_R, proxg_R, start, 5), fasta_cpp(start, 5, H, y, h_pass, l_pass, sigma2, x_true, lamb, beta_pen, dimen), replications = 20)

lamb <- 1.5
x_true <- rnorm(dimen^2, y, sd = 10)
start <- rnorm(dimen^2, x_true, sd = .01)
u <- fasta(f_R, gradf_R, g_R, proxg_R, start, 5, recordIterates = TRUE)
v <- fasta_cpp(start, 5, H, y, h_pass, l_pass,
               sigma2, x_true, lamb, beta_pen, dimen, recordIterates = TRUE)
hist(u$x - v$x)
(psi_R <- -log_pi(y, u$x))
(psi_cpp <- -log_pi(y, v$x))

plot.ts(u$objective, xlim = c(0, max(length(u$objective), length(v$objective))))
lines(v$objective, col = "red")

u$objective
v$objective

u$x - vec(v$x)

u$totalBacktracks
v$totalBacktracks


u$taus
v$taus

