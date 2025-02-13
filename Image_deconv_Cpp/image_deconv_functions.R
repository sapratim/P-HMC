
###################################################################################### 
################################  Image deconvolution ################################ 
################################  using wavelet frame ################################ 
###################################################################################### 

library(Matrix)
library(expm)
library(ks)
library(magick)
library(Rcpp)
library(dplyr)
sourceCpp("fasta_image_prob.cpp")

##############################  Model is (y = Hx + w) ##############################  
# y = Noisy image
# x = True image
# H = a uniform blur operator
# w = Gaussian noise N(0, sigma2)
# sigma2 to be specified globally
####################################################################################

###################### L_1 norm Wavelet function  ######################

wavelet_l1 <- function(image_vec, nlev = 3){
  image_mat <- matrix(image_vec, dimen, dimen)
  trans <- dwt.2d(image_mat, wf = "haar", nlev)
  wave_sum <- sum(abs(unlist(trans)))
  return(wave_sum)
}

###################### Log of target posterior #######################

log_pi <- function(y, x)   ##### input x as vector
{
  Hx <- convolve_image(x, dimen, dimen, H) 
  Psi.x <- wavelet_l1(x)
  f <- (norm(y - Hx,"2")^2)/(2*sigma2)
  g <- beta_pen*Psi.x
  return(-(f+g))
}

#########  Soft threshold function

softthreshold <- function(u, pen) {       ####  u is a vector
  return(sign(u)*sapply(u, FUN=function(x) {max(abs(x)-pen,0)}))
}

f <- function(z) {
  H.z <- convolve_image(z, dimen, dimen, H)
  t <- sum((y - H.z)^2)/(2*sigma2) + sum((x_true - z)^2)/(2*lamb)
  return(t)  
}

gradf <- function(z) {
  H.z <- convolve_image(z, dimen, dimen, H)
  t1 <- convolve_image((y - H.z), dimen, dimen, t(H))/sigma2
  t2 <- (x_true - z)/lamb
  return(t1+t2)  
}

g <- function(z) {
  Psi.z <- wavelet_l1(z)
  t <- beta_pen*Psi.z
  return(t)
}

proxg <- function(z, tau_fasta) {
  z <- matrix(z, dimen, dimen)
  wave_trans <- dwt.2d(z, wf = "haar", 3)
  for (i in 1:(length(wave_trans)-1)) {
    wave_trans[[i]] <- matrix(softthreshold(wave_trans[[i]],
                                            beta_pen*tau_fasta), dimen, dimen)
  }
  proxval <- idwt.2d(wave_trans)
  return(vec(proxval))
}


######### Proximity mapping for Chaari

grad_logpiLam <- function(x, f, gradf, g, proxg, fasta_start, fasta_step_start)  
{
  temp <- fasta(f, gradf, g, proxg, fasta_start, fasta_step_start,
                stepsizeShrink = .1, max_iters = 100)
  x_prox <- temp$x
  ans <-  (x-x_prox)/lambda
  return(-ans)
}

###################### Proximity mapping for Durmus ######################

prox_func_dur <- function(x, lambda) {   #### input x 
  x <- matrix(x, dimen, dimen)
  wave_trans <- dwt.2d(x, wf = "haar", 3)
  for (i in 1:(length(wave_trans)-1)) {
    wave_trans[[i]] <- matrix(softthreshold(wave_trans[[i]],
                                            beta_pen*lambda), dimen, dimen)
  }
  proxval <- idwt.2d(wave_trans)
  return(vec(proxval))
}

grad_logpiLam_dur <- function(x, y, lambda)  # gradient of log target for Durmus
{
  x_prox <- prox_func_dur(x, lambda)
  H.x <- convolve_image(x, dimen, dimen, H) 
  ratio_term <- pmax((y - H.x)/sigma2, 0)
  grad_f <- convolve_image(ratio_term, dimen, dimen, H) 
  ans <-  grad_f + (x-x_prox)/lambda
  return(-ans)
}

blur_func <- function(blur_size){
  mat <- matrix(rep(1/blur_size^2, blur_size^2), nrow = blur_size, ncol = blur_size)
  return(mat)
}

# Function to generate the Haar matrix for a given size (e.g., 128x128)
haar_matrix <- function(n) {
  # Initialize the Haar matrix as an identity matrix
  H <- diag(1, n)
  
  # Apply Haar transform recursively (log2(n) levels)
  for (level in 1:log2(n)) {
    # Size of the block
    block_size <- 2^level
    half_block_size <- block_size / 2
    
    # Loop through the matrix
    for (i in seq(1, n, by = block_size)) {
      for (j in i:(i + half_block_size - 1)) {
        # Low-pass (1) and high-pass (-1) coefficients
        H[j, i:(i + block_size - 1)] <- c(rep(1/sqrt(2), half_block_size), rep(-1/sqrt(2), half_block_size))
      }
    }
  }
  return(H)
}

######################################## Chaari #########################################

pxhmc_chaari <- function(y, lambda, f, gradf, g, proxg, iter, 
                         eps_hmc, L, start, fasta_start, fasta_step_start)
{
  nvar <- length(start)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  lamb <<- lambda
  
  # starting value computations
  samp <- start
  samp.hmc[1,] <- samp
  
  # For HMC
  mom_mat <- matrix(rnorm(iter*nvar), nrow = iter, ncol = nvar)
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    x_true <<- samp
    U_samp <- -grad_logpiLam(samp, f, gradf, g, proxg, x_true, fasta_step_start)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      x_true <<- samp
      U_samp <- -grad_logpiLam(samp, f, gradf, g, proxg, x_true, fasta_step_start)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- -log_pi(y, q_current)
    U_prop <- -log_pi(y, samp)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- samp
      accept <- accept + 1
    }else
    {
      samp.hmc[i,] <- q_current
      samp <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, acc_rate)
  return(object)
}


######################################## Durmus #########################################


pxhmc_dur <- function(y, lambda, iter, eps_hmc, L, start)
{
  nvar <- length(start)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  
  # starting value computations
  samp <- start
  samp.hmc[1,] <- samp
  
  # For HMC
  mom_mat <- matrix(rnorm(iter*nvar), nrow = iter, ncol = nvar)
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i,]
    U_samp <- -grad_logpiLam_dur(samp, y, lambda)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -grad_logpiLam_dur(samp, y, lambda)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- -log_pi(y, q_current)
    U_prop <- -log_pi(y, samp)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i,] <- samp
      accept <- accept + 1
    }else
    {
      samp.hmc[i,] <- q_current
      samp <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, acc_rate)
  return(object)
}
