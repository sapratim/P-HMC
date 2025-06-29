########################################
# One-dimensional example
########################################
set.seed(1)
n <- 1e2
y <- rnorm(n, mean = 1, sd = .5)
################################

# U(x) = sum(y - x)^2/(2) + |x|
log_target <- function(x)  -sum((y - x)^2)/2 - abs(x)

grad_log_pi_dur <- function(x, lambda)
{
  prox <- sign(x)*pmax(abs(x) - lambda, 0)
  sum(y - x) - (x - prox)/lambda
}

grad_log_pi_px <- function(x, lambda)
{
  foo <- lambda/(1 + n*lambda)
  z <-  (x + lambda*sum(y))/(1 + n*lambda)
  prox <- sign(z)* pmax(abs(z) - foo, 0)
  -(x - prox)/lambda
}

##################################
# plotting target contour

log_kinetic <- function(p) -p^2/2
log_joint <- function(q, p)
{
  log_target(q) + log_kinetic(p)
}
p_vec <- seq(-3, 3, length = 1e2)
q_vec <- seq(-2, 4, length = 1e2)
z <- matrix(0, nrow = length(q_vec), ncol = length(p_vec))
for(i in 1:length(q_vec))
{
  for(j in 1:length(p_vec))
  {
    z[i,j] <- exp(log_joint(q_vec[i], p_vec[j]))
  }
}

dur_path <- function(q, p, lambda, eps, L)
{
  store <- matrix(0, nrow = L + 1, ncol = 2)
  q_current <- q
  p_current <- p
  
  store[1, ] <- c(q,p)
  for (j in 2:(L+1))
  {
    # half step
    U_samp <- -grad_log_pi_dur(q_current, lambda)
    p_current <- p_current - eps*U_samp /2  
    q_current <- q_current + eps*p_current   # full step for position
    U_samp <- -grad_log_pi_dur(q_current, lambda)
    
    p_current <- p_current - eps*U_samp/2 # another half step
    
    store[j, ] <- c(q_current, p_current)
    # samp <- q_current
    # p_prop <- p_current
  }
  return(store)
}

px_path <- function(q, p, lambda, eps, L)
{
  store <- matrix(0, nrow = L + 1, ncol = 2)
  q_current <- q
  p_current <- p
  
  store[1, ] <- c(q,p)
  for (j in 2:(L+1))
  {
    # half step
    U_samp <- -grad_log_pi_px(q_current, lambda)
    p_current <- p_current - eps*U_samp /2  
    q_current <- q_current + eps*p_current   # full step for position
    U_samp <- -grad_log_pi_px(q_current, lambda)
    
    p_current <- p_current - eps*U_samp/2 # another half step
    
    store[j, ] <- c(q_current, p_current)
    # samp <- q_current
    # p_prop <- p_current
  }
  return(store)
}

draw_contour <- function(q, p, lambda, eps = .02, L = 10)
{
  path_dur <- dur_path(q, p, lambda = lambda, eps = eps, L = L)
  path_px <- px_path(q, p, lambda = lambda, eps = eps, L = L)
  
  contour(q_vec, p_vec, z, nlevels = 2, 
          xlim = c(.93, 1.15), ylim = c(-1.2, 1.2),
          main = bquote(lambda == .(lambda)), 
          xlab = bquote(x), ylab = bquote(p), lwd = 2)
  points(path_dur, col = "purple", pch = 19)
  points(path_px, col = "orange", pch = 19)
}

pdf("toy_contours.pdf", height = 6, width = 10)
layout(matrix(c(1:8, rep(9, 4)), nrow = 3, byrow = TRUE),
       heights = c(1, 1, 0.3))  # Last row is for legend
# Set tight margins for each plot (bottom, left, top, right)
par(mar = c(2, 2, 2, 1), oma = c(4, 4, 2, 1))  # oma allows outer text space

q <- .947
p <- 0
exp(log_joint(q, p))
draw_contour(q, p, lambda = .001, eps = .01, L = 20)
draw_contour(q, p, lambda = .01, eps = .01, L = 20)
draw_contour(q, p, lambda = 1, eps = .01, L = 20)
draw_contour(q, p, lambda = 10, eps = .01, L = 20)

q <- 1.1
p <- -.8
exp(log_joint(q, p))

draw_contour(q, p, lambda = .001, eps = .01, L = 20)
draw_contour(q, p, lambda = .01, eps = .01, L = 20)
draw_contour(q, p, lambda = 1, eps = .01, L = 20)
draw_contour(q, p, lambda = 10, eps = .01, L = 20)

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = c("Partial Proximal", "Full Proximal"),
       col = c("purple", "orange"), pch = 19, horiz = TRUE, bty = "n", cex = 1.5)


# Add shared axis labels using outer margin text
mtext("Potential component: x", side = 1, line = 2.2, outer = TRUE, cex = 1.2)
mtext("Momentum component: p", side = 2, line = 2.2, outer = TRUE, cex = 1.2)

dev.off()


q <- 1
p <- -.816110003
exp(log_joint(q, p))

draw_contour(q, p, lambda = .001, eps = .01, L = 20)
draw_contour(q, p, lambda = .01, eps = .01, L = 20)
draw_contour(q, p, lambda = 1, eps = .01, L = 20)
draw_contour(q, p, lambda = 10, eps = .01, L = 0)

###########################
# HMC with Durmus split
durhmc <- function(lambda, iter, eps_hmc, L, start)
{
  samp.hmc <- numeric(length = iter)
  
  # starting value computations
  samp <- start
  samp.hmc[1] <- samp
  
  # For HMC
  mom_mat <- rnorm(iter)
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i]
    U_samp <- -grad_log_pi_dur(samp, lambda)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -grad_log_pi_dur(samp, lambda)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- -log_target(q_current)
    U_prop <- -log_target(samp)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i] <- samp
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i] <- q_current
      samp <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/i
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, acc_rate)
  return(object)
}


###########################
# HMC on whole proximal
pxhmc <- function(lambda, iter, eps_hmc, L, start)
{
  samp.hmc <- numeric(length = iter)
  
  # starting value computations
  samp <- start
  samp.hmc[1] <- samp
  
  # For HMC
  mom_mat <- rnorm(iter)
  accept <- 0
  
  for (i in 2:iter) 
  {
    p_prop <- mom_mat[i]
    U_samp <- -grad_log_pi_px(samp, lambda)
    p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
    q_current <- samp
    for (j in 1:L)
    {
      samp <- samp + eps_hmc*p_current   # full step for position
      U_samp <- -grad_log_pi_px(samp, lambda)
      if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
    }
    p_current <- p_current - eps_hmc*U_samp/2
    p_current <- - p_current  # negation to make proposal symmetric
    
    U_curr <- -log_target(q_current)
    U_prop <- -log_target(samp)
    K_curr <-  sum((p_prop^2)/2)
    K_prop <-  sum((p_current^2)/2)
    
    log_acc_prob = U_curr - U_prop + K_curr - K_prop
    
    if(log(runif(1)) <= log_acc_prob )
    {
      samp.hmc[i] <- samp
      accept <- accept + 1
    }
    else
    {
      samp.hmc[i] <- q_current
      samp <- q_current
    }
    if(i %% (iter/10) == 0){
      j <- accept/i
      print(cat(i, j))}
  } 
  print(acc_rate <- accept/iter)
  object <- list(samp.hmc, acc_rate)
  return(object)
}


chain_dur <- durhmc(lambda = .01, iter = 1e4, eps_hmc = .1, L = 10, start = 2)
chain_px <- pxhmc(lambda = .01, iter = 1e4, eps_hmc = .1, L = 10, start = 2)

plot(density(chain_dur[[1]]))
plot(density(chain_px[[1]]))
acf(chain_dur[[1]])
acf(chain_px[[1]])
