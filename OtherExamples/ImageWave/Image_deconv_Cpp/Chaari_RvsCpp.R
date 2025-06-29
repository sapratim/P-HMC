############# Cpp vs R #############


###################################################################################### 
################################  Image deconvolution ################################ 
################################  using wavelet frame ################################ 
###################################################################################### 

library(Matrix)
library(expm)
library(ks)
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
#library(fasta)
library(waveslim)
sourceCpp("image_functions_copy.cpp")
haar_level <<- 3
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
  Psi.x <- wavelet_l1(x, dimen)
  f <- (norm(y - Hx,"2")^2)/(2*sigma2)
  g <- beta_pen*Psi.x
  return(-(f+g))
}

#########  Soft threshold function

# softthreshold <- function(u, pen) {       ####  u is a vector
#   return(sign(u)*sapply(u, FUN=function(x) {max(abs(x)-pen,0)}))
# }

f <- function(z) {
  H.z <- convolve_image(z, dimen, dimen, H)
  t <- sum((H.z - y)^2)/(2*sigma2) + sum((x_true - z)^2)/(2*lamb)
  return(t)
}

gradf <- function(z) {
  H.z <- convolve_image(z, dimen, dimen, H)
  t1 <- convolve_image((H.z - y), dimen, dimen, t(H))/sigma2
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


###################################  Fasta_R version ###################################

my_fasta <- function(f, gradf, g, proxg, x0, tau1, max_iters = 100, w = 10, 
                      backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                      eps_n = 1e-15) 
{
  residual <- double(max_iters)
  normalizedResid <- double(max_iters)
  taus <- double(max_iters)
  fVals <- double(max_iters)
  objective <- double(max_iters + 1)
  totalBacktracks <- 0
  backtrackCount <- 0
  x1 <- x0
  d1 <- x1
  f1 <- f(d1)
  fVals[1] <- f1
  gradf1 <- gradf(d1)
  if (recordIterates) {
    iterates <- matrix(0, length(x0), max_iters + 1)
    iterates[, 1] <- x1
  }
  else {
    iterates <- NULL
  }
  maxResidual <- -Inf
  minObjectiveValue <- Inf
  objective[1] <- f1 + g(x0)
  for (i in 1:max_iters) {
    x0 <- x1
    gradf0 <- matrix(gradf1)
    tau0 <- tau1
    x1hat <- x0 - tau0 * c(gradf0)
    x1 <- proxg(x1hat, tau0)
    Dx <- matrix(x1 - x0)
    d1 <- x1
    f1 <- f(d1)
    if (backtrack) {
      M <- max(fVals[max(i - w, 1):max(i - 1, 1)])
      backtrackCount <- 0
      prop <- (f1 - 1e-12 > M + t(Dx) %*% gradf0 + 0.5 * 
                 (norm(Dx, "f")^2)/tau0) && (backtrackCount < 
                                               20)
      while (prop) {
        tau0 <- tau0 * stepsizeShrink
        x1hat <- x0 - tau0 * c(gradf0)
        x1 <- proxg(x1hat, tau0)
        d1 <- x1
        f1 <- f(d1)
        Dx <- matrix(x1 - x0)
        backtrackCount <- backtrackCount + 1
        prop <- (f1 - 1e-12 > M + t(Dx) %*% gradf0 + 
                   0.5 * (norm(Dx, "f")^2)/tau0) && (backtrackCount < 
                                                       20)
      }
      totalBacktracks <- totalBacktracks + backtrackCount
    }
    taus[i] <- tau0
    residual[i] <- norm(Dx, "f")/tau0
    maxResidual <- max(maxResidual, residual[i])
    normalizer <- max(norm(gradf0, "f"), norm(as.matrix(x1 - 
                                                          x1hat), "f")/tau0) + eps_n
    normalizedResid[i] <- residual[i]/normalizer
    fVals[i] <- f1
    objective[i + 1] <- f1 + g(x1)
    newObjectiveValue <- objective[i + 1]
    if (recordIterates) {
      iterates[, i + 1] <- x1
    }
    if (newObjectiveValue < minObjectiveValue) {
      bestObjectiveIterate <- x1
      minObjectiveValue <- min(minObjectiveValue, newObjectiveValue)
    }
    gradf1 <- gradf(d1)
    Dg <- matrix(gradf1 + (x1hat - x0)/tau0)
    dotprod <- t(Dx) %*% Dg
    tau_s <- norm(Dx, "f")^2/dotprod
    tau_m <- dotprod/norm(Dg, "f")^2
    tau_m <- max(tau_m, 0)
    if (abs(dotprod) < 1e-15) 
      break
    if (2 * tau_m > tau_s) {
      tau1 <- tau_m
    }
    else {
      tau1 <- tau_s - 0.5 * tau_m
    }
    if ((tau1 <= 0) || is.infinite(tau1) || is.nan(tau1)) {
      tau1 <- tau0 * 1.5
    }
  }
  if (recordIterates) {
    iterates <- iterates[, 1:(i + 1), drop = FALSE]
  }
  return(list(x = bestObjectiveIterate, objective = objective[1:(i + 1)], fVals = fVals[1:i], 
              totalBacktracks = totalBacktracks, residual = residual[1:i], taus = taus[1:i],
              iterates = iterates, iternumber = i))
}

########################### Proximity mapping for Chaari ###########################

grad_logpiLam_R <- function(x, y, f, gradf, g, proxg, fasta_start, fasta_step_start)  
{
  temp <- my_fasta(f, gradf, g, proxg, fasta_start, fasta_step_start, eps_n = 1e-10)
  x_prox <- temp$x
  ans <-  (x-x_prox)/lambda
  iters <- temp$iternumber
  out <- list(-ans, iters)
  return(out)
}

grad_logpiLam_cpp <- function(x, fasta_start, fasta_step_start, H, y, h_fil, l_fil, sigma2, 
                              lambda, beta_pen, dimen)  
{
  temp <- fasta_cpp(fasta_start, fasta_step_start, H, y, h_fil, l_fil, sigma2, x_true, 
                    lamb, beta_pen, dimen,  eps_n = 1e-10)
  x_prox <- temp$x
  ans <-  (x-x_prox)/lambda
  iters <- temp$iternumber
  out <- list(-ans, iters)
  return(out)
}

blur_func <- function(blur_size){
  mat <- matrix(rep(1/blur_size^2, blur_size^2), nrow = blur_size, ncol = blur_size)
  return(mat)
}

######################################## Chaari ########################################

pxhmc_chaari_R <- function(y, lambda, iter, 
                         eps_hmc, L, start, fasta_start, fasta_step_start)
{
  nvar <- length(start)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  lamb <<- lambda
  iter_num_vec <- numeric(length = iter - 1)
  
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
    grad_list <- grad_logpiLam_R(samp, y, f, gradf, g, proxg, fasta_start = samp, fasta_step_start)
    U_samp <- -grad_list[[1]]
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    sum_iter <- 0
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      x_true <<- samp
      grad_list <- grad_logpiLam_R(samp, y, f, gradf, g, proxg, fasta_start = samp, fasta_step_start)
      U_samp <- -grad_list[[1]]
      sum_iter <- sum_iter + grad_list[[2]]
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    iter_num_vec[i-1] <- sum_iter
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
  #print(acc_rate <- accept/iter)
  object <- list(samp.hmc, iter_num_vec)
  return(object)
}


######################################## Chaari #########################################

pxhmc_chaari_cpp <- function(y, lambda, iter, 
                         eps_hmc, h_fil, l_fil, L, start, fasta_start, fasta_step_start)
{
  nvar <- length(start)
  samp.hmc <- matrix(0, nrow = iter, ncol = nvar)
  iter_num_vec <- numeric(length = iter - 1)
  
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
    grad_list <- grad_logpiLam_cpp(samp, fasta_start = samp, fasta_step_start, H, y, h_fil, l_fil,
                             sigma2, lambda, beta_pen, dimen)
    U_samp <- -grad_list[[1]]
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    sum_iter <- 0
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      x_true <<- samp
      grad_list <- grad_logpiLam_cpp(samp, fasta_start = samp, fasta_step_start, H, y, h_fil, l_fil,
                               sigma2, lambda, beta_pen, dimen)
      U_samp <- -grad_list[[1]]
      sum_iter <- sum_iter + grad_list[[2]]
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    iter_num_vec[i-1] <- sum_iter
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
  #print(acc_rate <- accept/iter)
  object <- list(samp.hmc, iter_num_vec)
  return(object)
}




######################################  Run file ######################################

load("pixel_mat.Rdata")

true_pixel_vec <- c(pixel_mat)
dimen <- sqrt(length(true_pixel_vec))
H <- blur_func(blur_size = 5)    ######  Blur matrix
y <- convolve_image(true_pixel_vec, dimen, dimen, H)
beta_pen <<- 0.02   ###### penalty parameter
sigma2 <- 1.76     ###### noise variance
h_pass <- wave.filter("haar")$hpf     ###### high pass filter
l_pass <- wave.filter("haar")$lpf     ###### low pass filter

##################  frequentist mode evaluation ##################

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


freq_mode <- my_fasta(f_freq, gradf_freq, g, proxg, y, tau1 = 5, stepsizeShrink = .5, 
                   max_iters = 100)

system.time(pxhmc_chaari_run.R <- pxhmc_chaari_R(y = y, lambda <- 10, 
                                             iter = 1e3, eps_hmc <- 0.045, L = 10, freq_mode$x,freq_mode$x, 5))

system.time(pxhmc_chaari_run.cpp <- pxhmc_chaari_cpp(y = y, lambda <- 10, 
                                             iter = 1e3, eps_hmc <- 0.038, h_pass, l_pass, L = 10, freq_mode$x,freq_mode$x, 5))

#samples_chaari_R <- pxhmc_chaari_run[[1]]

#save(samples_chaari_R, file = "chaari_R.Rdata")

