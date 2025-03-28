
########################################################################
###########################  Problem run  ##############################
########################################################################

library(magick)
library(fasta)
library(wavelets)
library(waveslim)
source("image_deconv_functions.R")
# load("pixel_mat.Rdata")

################  Define global variables here ################

img_input <- image_read("Boat_gray.png")
img <- image_resize(img_input, "16x16")
pixel_mat <- as.matrix(as.integer(image_data(img)))
#save(pixel_mat, file = "pixel_mat.Rdata")
true_pixel_vec <- c(pixel_mat)
dimen <- sqrt(length(true_pixel_vec))
H <- blur_func(blur_size = 5)    ######  Blur matrix
y <- convolve_image(true_pixel_vec, dimen, dimen, H)
beta_pen <- 0.02   ###### penalty parameter
sigma2 <- 1.76     ###### noise variance

# #######################   Frequentist mode evaluation   #######################

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

g_R <- function(z) {
  Psi.z <- wavelet_l1_cpp(z, dimen, h_pass, l_pass)
  t <- beta_pen*Psi.z
  return(t)
}

proxg_R <- function(z, tau_fasta) {
  z <- matrix(z, dimen, dimen)
  wave_trans <- dwt.2d(z, "haar", 3)
  for (i in 1:(length(wave_trans)-1)) {
    wave_trans[[i]] <- matrix(softthreshold(wave_trans[[i]],
                                               beta_pen*tau_fasta), dimen, dimen)
  }
  proxval <- idwt.2d(wave_trans)
  return(vec(proxval))
}

# #######################  Optimal tau and stepshrink #######################
# 
taus <- c(1:10)
steps <- seq(0.1, 0.8, by = .1)
dens_values <- matrix(0, length(taus), length(steps))
for (i in 1:length(taus)) {
  for (j in 1:length(steps)) {
    val <- fasta(f_freq, gradf_freq, g, proxg_R, y, tau1 = taus[i], stepsizeShrink = steps[j])
    dens_values[i, j] = log_pi(y, val$x)
  }
}
# 
# 
# 
 freq_mode <- fasta(f_freq, gradf_freq, g_R, proxg_R, y, tau1 = 1, stepsizeShrink = .3, 
                          max_iters = 100)


########################  Chaari/Durmus runs  ########################

system.time(pxhmc_chaari_run <- pxhmc_chaari(y = y, lambda <- 10, 
                           iter = 1e2, eps_hmc <- 0.01, L = 10, freq_mode$x,freq_mode$x, 3))


system.time(pxhmc_dur_run <- pxhmc_dur(y, lambda = 1 , iter = 1e4, eps_hmc = 0.13, L = 10,
                           start = freq_mode$x))

output <- list(pxhmc_chaari_run[[1]], pxhmc_dur_run[[1]])

save(output, file = "samples.Rdata")

###############################  Mode matching  ###############################
###############################   Validation    ###############################
##############################  Run sequentially ##############################

######  Chaari ######

# i <- 1
# plot(density(pxhmc_chaari_run[[1]][,i]), type = 'l')
# abline(v = freq_mode$x[i], col = "red")
# i <- i+1
# 
# # ######  Durmus ###### 
# # 
# i <- 1
# plot(density(pxhmc_dur_run[[1]][,i]), type = 'l', col = "blue")
# lines(density(pxhmc_chaari_run[[1]][,i]), type = 'l', col = "red")
# abline(v = freq_mode$x[i], col = "black")
# i <- i+1
# 


# dens_val <- numeric(length = 1e4)
# for(i in 1:1e4)
# {
#  dens_val[i] <- log_pi(y, output[[1]][i,]) 
# }
# 
# map_est <- freq_mode_cpp$x
#   #output[[1]][which.max(dens_val),]
# 
# 
# mat <- matrix(map_est, dimen, dimen)
#  mat_blur <- matrix(y, dimen, dimen)
#  img_blur <- image(t(mat_blur)[, nrow(mat_blur):1], col = gray.colors(256), axes = FALSE)
#  img_map <- image(t(mat)[, nrow(mat):1], col = gray.colors(256), axes = FALSE)
#  
# 
# 
# 
#  
#  
 # 
 # dens_vals <- numeric(length = nrow(output[[1]]))
 # 
 # for(i in 1:nrow(output[[1]])){
 #   dens_vals[i] <- log_pi(y, output[[1]][i,])
 # }
 # 
 # plot(dens_vals[1:1e4], type = "l")
 # 
 # 
 # 
 # dens_vals_dur <- numeric(length = nrow(output[[1]]))
 # 
 # for(i in 1:nrow(output[[2]])){
 #   dens_vals_dur[i] <- log_pi(y, output[[2]][i,])
 # }
 # 
 # plot(dens_vals_dur[1:1e4], type = "l")
 # 
 
 
 
 
 
 
#  load("samples.Rdata")
#  
#  grad_logpiLam <- function(x, fasta_start, fasta_step_start, H, y, sigma2, 
#                            lambda, beta_pen, dimen)  
#  {
#    temp <- fasta_cpp(x, fasta_step_start, H, y, sigma2, x, lambda, beta_pen, dimen)
#    # temp <- fasta(f, gradf, g, proxg, fasta_start, fasta_step_start,
#    #              stepsizeShrink = .1, max_iters = 100)
#    x_prox <- temp$x
#    ans <-  (x-x_prox)/lambda
#    plot.ts(temp$objective)
#    return(-ans)
#  }
#  
#  
# # grad_logpiLam <- function(x, fasta_start, fasta_step_start, H, y, sigma2, 
#                    #        lambda, beta_pen, dimen)  
#    
#  
#  wavelet_l1 <- function(image_vec, nlev = 3){
#    image_mat <- matrix(image_vec, dimen, dimen)
#    trans <- dwt.2d(image_mat, wf = "haar", nlev)
#    wave_sum <- sum(abs(unlist(trans)))
#    return(wave_sum)
#  }
#  
#  
#  foo <- output[[1]][1,]
#  
#  
#  lambda <- 10
#  
#  cpp_val <- -grad_logpiLam(foo, freq_mode$x, 3, H, y,
#                            sigma2, lambda, beta_pen, dimen)
#  
#  
#  
#  x_true <- foo
#  lamb <- 10
#    
#   softthreshold <- function(u, pen) {       ####  u is a vector
#     return(sign(u)*sapply(u, FUN=function(x) {max(abs(x)-pen,0)}))
#   }
#  
#   f <- function(z) {
#     H.z <- convolve_image(z, dimen, dimen, H)
#     t <- sum((H.z - y)^2)/(2*sigma2) + sum((x_true - z)^2)/(2*lamb)
#     return(t)
#   }
#  
#   gradf <- function(z) {
#     H.z <- convolve_image(z, dimen, dimen, H)
#     t1 <- convolve_image((H.z - y), dimen, dimen, t(H))/sigma2
#     t2 <- (x_true - z)/lamb
#     return(t1+t2)
#   }
#  
#   g <- function(z) {
#     Psi.z <- wavelet_l1(z)
#     t <- beta_pen*Psi.z
#     return(t)
#   }
#  
#   proxg <- function(z, tau_fasta) {
#     z <- matrix(z, dimen, dimen)
#     wave_trans <- dwt.2d(z, wf = "haar", 3)
#     for (i in 1:(length(wave_trans)-1)) {
#       wave_trans[[i]] <- matrix(softthreshold(wave_trans[[i]],
#                                               beta_pen*tau_fasta), dimen, dimen)
#     }
#     proxval <- idwt.2d(wave_trans)
#     return(vec(proxval))
#   }
#  
#     
#   R_val <- fasta(f, gradf, g, proxg, foo, tau1 = 3, stepsizeShrink = .1, 
#                      max_iters = 100)
#  
#   x_prox <- R_val$x
#   ans <-  (foo-x_prox)/lambda
#   r_grad <- -ans
# hist(cpp_val - r_grad) 
#  
