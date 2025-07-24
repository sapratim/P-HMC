###############################
# File contains  Data generation
# Generating the CheckerBoard Image

# source this file. "y" stores the final image
##############################
library(ks)

# cameraman <- load.image("cameraman.png")
# plot(cameraman)
# dim(cameraman)
# 
# cam <- as.matrix(cameraman[, , , 1])
# plot(as.cimg(cam))

set.seed(8024248)
n <- 64
a <- 8
mat <- matrix(0, nrow = n, ncol = n)
vec.mat <- rep(c(1, 0), each = a, times = n/(2*a))
checker <- matrix(0, nrow = n, ncol = n)
for(j in seq(1, n/2, by = a))
{
  for(k in 1:a)
  {
    checker[j+k-1, ] <- vec.mat
  }
  vec.mat <- rev(vec.mat)
}
for(j in seq((n/2+1), n, by = a))
{
  for(k in 1:a)
  {
    checker[j+k-1, ] <- (vec.mat == 0)*(vec.mat)  + (vec.mat == 1)*(vec.mat - .50)
  }
  vec.mat <- rev(vec.mat)
}
n <- dim(checker)[2]
noise <- matrix(rnorm(n^2, 0, sqrt(0.01)), nrow = n, ncol = n)
image_mat <- checker + noise

x <- vec(checker)
y <- vec(image_mat)

# plot(as.cimg(y))


# sure_sv_threshold <- function(Y, sigma2, lambda) {
#   svd_Y <- svd(Y)
#   d <- svd_Y$d
#   r <- sum(d > lambda)
#   sure <- -length(Y) * sigma2 +
#     sum(pmin(d^2, lambda^2)) +
#     2 * sigma2 * r * (nrow(Y) + ncol(Y) - r)
#   return(sure)
# }
# 

# max_sig <- max(svd(image_mat)$d)
# ells <- seq(.1, 1.3, length = 1000)
# sures <- sapply(ells, function(l) sure_sv_threshold(Y = image_mat, sigma2 = .01, lambda = l))
# 
# ell.star <- ells[which.min(sures)]
# alpha_hat <- ell.star/sigma2_hat
sigma2_hat <- 0.01
alpha_hat <- 1.15/sigma2_hat

######### Display the checkerboard matrix as an image

# image(checker, col = gray.colors(4, start = 0, end = 1), axes = FALSE)
# image(image_mat, col = gray.colors(4, start = 0, end = 1), axes = FALSE)
