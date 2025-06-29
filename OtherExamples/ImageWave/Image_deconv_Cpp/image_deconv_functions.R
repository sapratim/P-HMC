
###################################################################################### 
################################  Image deconvolution ################################ 
################################  using wavelet frame ################################ 
###################################################################################### 

library(Matrix)
library(expm)
library(ks)
library(Rcpp)
library(dplyr)
library(fasta)
library(waveslim)
sourceCpp("image_functions_copy.cpp")
haar_level <<- 3
h_pass <- wave.filter("haar")$hpf
l_pass <- wave.filter("haar")$lpf

##############################  Model is (y = Hx + w) ##############################  
# y = Noisy image
# x = True image
# H = a uniform blur operator
# w = Gaussian noise N(0, sigma2)
# sigma2 to be specified globally
####################################################################################

###################### L_1 norm Wavelet function  ######################

# wavelet_l1 <- function(image_vec, nlev = 3){
#   image_mat <- matrix(image_vec, dimen, dimen)
#   trans <- dwt.2d(image_mat, wf = "haar", nlev)
#   wave_sum <- sum(abs(unlist(trans)))
#   return(wave_sum)
# }

###################### Log of target posterior #######################

log_pi <- function(y, x)   ##### input x as vector
{
  Hx <- convolve_image(x, dimen, dimen, H) 
  Psi.x <- wavelet_l1_cpp(x, dimen, h_pass, l_pass)
  f <- (norm(y - Hx,"2")^2)/(2*sigma2)
  g <- beta_pen*Psi.x
  return(-(f+g))
}

#########  Soft threshold function

# softthreshold <- function(u, pen) {       ####  u is a vector
#   return(sign(u)*sapply(u, FUN=function(x) {max(abs(x)-pen,0)}))
# }

# f <- function(z) {
#   H.z <- convolve_image(z, dimen, dimen, H)
#   t <- sum((H.z - y)^2)/(2*sigma2) + sum((x_true - z)^2)/(2*lamb)
#   return(t)
# }
#
# gradf <- function(z) {
#   H.z <- convolve_image(z, dimen, dimen, H)
#   t1 <- convolve_image((H.z - y), dimen, dimen, t(H))/sigma2
#   t2 <- (x_true - z)/lamb
#   return(t1+t2)
# }
#
# g <- function(z) {
#   Psi.z <- wavelet_l1(z)
#   t <- beta_pen*Psi.z
#   return(t)
# }
#
# proxg <- function(z, tau_fasta) {
#   z <- matrix(z, dimen, dimen)
#   wave_trans <- dwt.2d(z, wf = "haar", 3)
#   for (i in 1:(length(wave_trans)-1)) {
#     wave_trans[[i]] <- matrix(softthreshold(wave_trans[[i]],
#                                             beta_pen*tau_fasta), dimen, dimen)
#   }
#   proxval <- idwt.2d(wave_trans)
#   return(vec(proxval))
# }
#

######### Proximity mapping for Chaari

grad_logpiLam <- function(fasta_start, fasta_step_start, dimen, H, y, sigma2, 
                          x, lambda, beta_pen)  
{
  #temp <- fasta_cpp(fasta_start, fasta_step_start, H, y, sigma2, x, lambda, beta_pen, dimen)
  temp <- fasta_cpp(fasta_start, fasta_step_start, dimen, h_pass, l_pass, H, y, sigma2, 
                    x, lambda, beta_pen, stepsizeShrink = .1)
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
  ratio_term <- (H.x - y)
  grad_f <- convolve_image(ratio_term, dimen, dimen, H)/sigma2 
  ans <-  grad_f + (x-x_prox)/lambda
  return(-ans)
}

############################  Blurring matrix ############################

blur_func <- function(blur_size){
  mat <- matrix(rep(1/blur_size^2, blur_size^2), nrow = blur_size, ncol = blur_size)
  return(mat)
}

######################################## Chaari #########################################

pxhmc_chaari <- function(y, lambda, iter, 
                         eps_hmc, L, start, fasta_start, fasta_step_start)
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
   # x_true <<- samp
    U_samp <- -grad_logpiLam(fasta_start = samp, fasta_step_start, dimen, H, y, sigma2, 
                             samp, lambda, beta_pen)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      #x_true <<- samp
      U_samp <- -grad_logpiLam(fasta_start = samp, fasta_step_start, dimen, H, y, sigma2, 
                               samp, lambda, beta_pen)
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


rwm_alg <- function(y, start, iter, h)
{
  nvar <- length(start)
  samp.rwm <- matrix(0, nrow = iter, ncol = nvar)
  
  # starting value computations
  samp <- start
  samp.rwm[1,] <- samp

  # count number of acceptances
  accept <- 0        
  gaussian_mat <- matrix(rnorm((iter-1)*nvar), nrow = (iter-1), ncol = nvar)
  
  for (i in 2:iter) 
  {
    propval <- samp.rwm[i-1,] + h*gaussian_mat[i-1,]
    log_ratio <- log_pi(y, propval) - log_pi(y, samp.rwm[i-1,])
    if(log(runif(1)) < log_ratio)
    {
      samp.rwm[i,] <- propval
      accept <- accept + 1
    } else{
    samp.rwm[i,] <- samp.rwm[i-1,]
    }
    if(i %% (iter/10) == 0){
      j <- accept/iter
      print(cat(i, j))}
  }
  print(accept/iter)
  return(samp.rwm)
}



