
########################################################################
###########################  Problem run  ##############################
########################################################################

library(magick)
library(fasta)
library(wavelets)
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
  t1 <- convolve_image((y - H.z), dimen, dimen, t(H))/sigma2
  return(t1)  
}

freq_mode <- fasta(f_freq, gradf_freq, g, proxg, true_pixel_vec,
                   tau1 = 1, stepsizeShrink = 0.1)

########################  Chaari/Durmus runs  ########################

pxhmc_chaari_run <- pxhmc_chaari(y = y, lambda <- 10, f, gradf, g, proxg, iter = 1e2,
                                 eps_hmc <- 0.015, L = 10, y,
                                 blur_pixel_vec, 3)

pxhmc_dur_run <- pxhmc_dur(y, lambda = 10, iter = 1e4, eps_hmc = 0.003, L = 10,
                           start = y)

 i <- 1
 plot(density(pxhmc_chaari_run[[1]][,i]), type = 'l')
 abline(v = freq_mode$x[i], col = "red")
 i <- i+1

# Create Haar matrix for 128x128 image

######  For plotting the image from matrix: image_read(as.raster(noisy_mat))

