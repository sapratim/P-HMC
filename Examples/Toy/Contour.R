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
# Testing

logpi_lam_dur <- function(x, lambda)
{
  prox <- sign(x)*pmax(abs(x) - lambda, 0)
  my <- abs(prox) + (x - prox)^2/(2*lambda)
  
  logpi_lam <- - my -sum((y - x)^2)/2
}

logpi_lam_px <- function(x, lambda)
{
  foo <- lambda/(1 + n*lambda)
  z <-  (x + lambda*n*mean(y))/(1 + n*lambda)
  prox <- sign(z)* pmax(abs(z) - foo, 0)
  
  my <- abs(prox) + (x - prox)^2/(2*lambda) + sum((y - prox)^2)/2
  
  logpi_lam <- - my
}

integrate_f <- function(x){
  sapply(x, function(t) exp(log_target(t)))
}

integrate_dur <- function(x, lambda){
  sapply(x, function(t) exp(logpi_lam_dur(t, lambda = lambda)))
}

integrate_px <- function(x, lambda){
  sapply(x, function(t) exp(logpi_lam_px(t, lambda = lambda)))
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

phmc_path <- function(q, p, lambda, eps, L)
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

ns_path <- function(q, p, lambda, eps, L)
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
  path_dur <- phmc_path(q, p, lambda = lambda, eps = eps, L = L)
  path_px <- ns_path(q, p, lambda = lambda, eps = eps, L = L)
  main_expr <- bquote(paste(lambda[g], " = ", lambda, " = ",
                            .(format(lambda, digits = 3))))
  contour(q_vec, p_vec, z, nlevels = 2, 
          xlim = c(.93, 1.15), ylim = c(-1.2, 1.2),
          main = main_expr, lwd = 2)
  points(path_dur, col = "purple", pch = 19)
  points(path_px, col = "orange", pch = 19)
}


eps <- .01
L <- 20
pdf("toy_contours.pdf", height = 5, width = 10)

# Layout: 8 plots + 1 row for top legend
layout(matrix(c(rep(9, 4), 1:8), nrow = 3, byrow = TRUE),
       heights = c(.12, 1, 1))  # first row for legend

# Tighter margins and minimal label-to-axis spacing
par(mar = c(2, 2, 1, 1), oma = c(2, 2, 2, 2), mgp = c(3, 0.3, 0))

q <- .947
p <- 0
exp(log_joint(q, p))
draw_contour(q, p, lambda = .001, eps = eps, L = L)
draw_contour(q, p, lambda = .01, eps = eps, L = L)
draw_contour(q, p, lambda = 1, eps = eps, L = L)
draw_contour(q, p, lambda = 10, eps = eps, L = L)
mtext(bquote("(" * x[0] * "," ~ p[0] * ")" == "(" * .(q) * "," ~ .(p) * ")"),
      side = 4, line = 0.5, outer = FALSE, cex = 0.8)

q <- 1.1
p <- -.8
exp(log_joint(q, p))

draw_contour(q, p, lambda = .001, eps = eps, L = L)
draw_contour(q, p, lambda = .01, eps = eps, L = L)
draw_contour(q, p, lambda = 1, eps = eps, L = L)
draw_contour(q, p, lambda = 10, eps = eps, L = L)
mtext(bquote("(" * x[0] * "," ~ p[0] * ")" == "(" * .(q) * "," ~ .(p) * ")"),
      side = 4, line = 0.5, outer = FALSE, cex = 0.8)

par(mar = c(0, 0, 0, 0))
plot.new()
legend("top", legend = c("ns-HMC", "p-HMC"),
       col = c("orange", "purple"), pch = 19, 
       horiz = TRUE, bty = "n", 
       inset = c(0, -.55),
       cex = 1.5)

# Add shared axis labels using outer margin text
mtext("Position: x", side = 1, line = .5, outer = TRUE, cex = 1)
mtext("Momentum: p", side = 2, line = .5, outer = TRUE, cex = 1)

dev.off()



############################
# Choosing lambda
q <- 1.04  # MAP
p <- -.5
eps <- 1e-7
lambda.seq <- seq(1e-3, 5, length = 1e3)

phmc_ham <- numeric(length = length(lambda.seq))
ns_ham <- numeric(length = length(lambda.seq))

ham <- function(q, p) - log_target(q) + p^2/2

for(i in 1:length(lambda.seq))
{
  phmc_state <- phmc_path(q, p, lambda.seq[i], eps, 1)[2, ]
  ns_state <- ns_path(q, p, lambda.seq[i], eps, 1)[2, ]

  phmc_ham[i] <- abs(ham(q,p) - ham(phmc_state[1], phmc_state[2]))/abs(ham(q,p) )
  ns_ham[i] <- abs(ham(q,p) - ham(ns_state[1], ns_state[2]))/abs(ham(q,p) )
}

pdf("lambda_toy.pdf", height = 3.5, width = 4.2)
plot(lambda.seq, phmc_ham, type = 'l', lwd = 2, 
  ylim = range(c(phmc_ham, 0)), xlab = expression(lambda[g]), ylab = expression(R[lambda[g]]))
abline(v = 1, lty = 2, lwd = 2)
legend("bottomright", legend = expression("Choice of " * lambda[g]),
       col = "black", lty = 2, lwd = 2, bty = "n")
dev.off()


# No need to run the chains
# ###########################
# # HMC with Durmus split
# durhmc <- function(lambda, iter, eps_hmc, L, start)
# {
#   samp.hmc <- numeric(length = iter)
#   
#   # starting value computations
#   samp <- start
#   samp.hmc[1] <- samp
#   
#   # For HMC
#   mom_mat <- rnorm(iter)
#   accept <- 0
#   
#   for (i in 2:iter) 
#   {
#     p_prop <- mom_mat[i]
#     U_samp <- -grad_log_pi_dur(samp, lambda)
#     p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
#     q_current <- samp
#     for (j in 1:L)
#     {
#       samp <- samp + eps_hmc*p_current   # full step for position
#       U_samp <- -grad_log_pi_dur(samp, lambda)
#       if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
#     }
#     p_current <- p_current - eps_hmc*U_samp/2
#     p_current <- - p_current  # negation to make proposal symmetric
#     
#     U_curr <- -log_target(q_current)
#     U_prop <- -log_target(samp)
#     K_curr <-  sum((p_prop^2)/2)
#     K_prop <-  sum((p_current^2)/2)
#     
#     log_acc_prob = U_curr - U_prop + K_curr - K_prop
#     
#     if(log(runif(1)) <= log_acc_prob )
#     {
#       samp.hmc[i] <- samp
#       accept <- accept + 1
#     }
#     else
#     {
#       samp.hmc[i] <- q_current
#       samp <- q_current
#     }
#     if(i %% (iter/10) == 0){
#       j <- accept/i
#       print(cat(i, j))}
#   } 
#   print(acc_rate <- accept/iter)
#   object <- list(samp.hmc, acc_rate)
#   return(object)
# }
# 
# 
# ###########################
# # HMC on whole proximal
# pxhmc <- function(lambda, iter, eps_hmc, L, start)
# {
#   samp.hmc <- numeric(length = iter)
#   
#   # starting value computations
#   samp <- start
#   samp.hmc[1] <- samp
#   
#   # For HMC
#   mom_mat <- rnorm(iter)
#   accept <- 0
#   
#   for (i in 2:iter) 
#   {
#     p_prop <- mom_mat[i]
#     U_samp <- -grad_log_pi_px(samp, lambda)
#     p_current <- p_prop - eps_hmc*U_samp /2  # half step for momentum
#     q_current <- samp
#     for (j in 1:L)
#     {
#       samp <- samp + eps_hmc*p_current   # full step for position
#       U_samp <- -grad_log_pi_px(samp, lambda)
#       if(j!=L) p_current <- p_current - eps_hmc*U_samp  # full step for momentum
#     }
#     p_current <- p_current - eps_hmc*U_samp/2
#     p_current <- - p_current  # negation to make proposal symmetric
#     
#     U_curr <- -log_target(q_current)
#     U_prop <- -log_target(samp)
#     K_curr <-  sum((p_prop^2)/2)
#     K_prop <-  sum((p_current^2)/2)
#     
#     log_acc_prob = U_curr - U_prop + K_curr - K_prop
#     
#     if(log(runif(1)) <= log_acc_prob )
#     {
#       samp.hmc[i] <- samp
#       accept <- accept + 1
#     }
#     else
#     {
#       samp.hmc[i] <- q_current
#       samp <- q_current
#     }
#     if(i %% (iter/10) == 0){
#       j <- accept/i
#       print(cat(i, j))}
#   } 
#   print(acc_rate <- accept/iter)
#   object <- list(samp.hmc, acc_rate)
#   return(object)
# }
# 
# 
# chain_dur <- durhmc(lambda = .01, iter = 1e4, eps_hmc = .1, L = 10, start = 2)
# chain_px <- pxhmc(lambda = .01, iter = 1e4, eps_hmc = .1, L = 10, start = 2)
# 
# plot(density(chain_dur[[1]]))
# plot(density(chain_px[[1]]))
# acf(chain_dur[[1]])
# acf(chain_px[[1]])
