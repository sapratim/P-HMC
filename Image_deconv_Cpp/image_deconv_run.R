
########################################################################
###########################  Problem run  ##############################
########################################################################

library(magick)
library(fasta)
library(wavelets)
library(waveslim)
source("image_deconv_functions.R")

################  Define global variables here ################

img_input <- image_read("Boat_gray.png") 
img <- image_resize(img_input, "16x16")
true_pixel_vec <- as.matrix(as.integer(image_data(img)))
dimen <- sqrt(nrow(true_pixel_vec))
H <- blur_func(blur_size = 5)    ######  Blur matrix
Psi <- haar_matrix(8)   #####  Haar frame
y <- convolve_image(true_pixel_vec, dimen, dimen, H)
beta_pen <- 0.02   ###### penalty parameter
sigma2 <- 1.76     ###### noise variance

#######################   Frequentist mode evaluation   #######################

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

freq_mode <- fasta(f_freq, gradf_freq, g, proxg, y, tau1 = 5)

########################  Chaari/Durmus runs  ########################

system.time(pxhmc_chaari_run <- pxhmc_chaari(y = y, lambda <- 10, f, gradf, g, proxg, 
                           iter = 1e3, eps_hmc <- 0.04, L = 10, freq_mode$x,freq_mode$x, 3))

system.time(pxhmc_dur_run <- pxhmc_dur(y, lambda = 10, iter = 1e4, eps_hmc = 0.03, L = 10,
                           start = freq_mode$x))

###############################  Mode matching  ###############################
###############################   Validation    ###############################
##############################  Run sequentially ##############################

######  Chaari ######

i <- 1
plot(density(pxhmc_chaari_run[[1]][,i]), type = 'l')
abline(v = freq_mode$x[i], col = "red")
i <- i+1

######  Durmus ###### 

i <- 1
plot(density(pxhmc_dur_run[[1]][,i]), type = 'l')
abline(v = freq_mode$x[i], col = "red")
i <- i+1

