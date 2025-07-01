###############################
# File contains  Data generation
# Generating the CheckerBoard Image

# source this file. "y" stores the final image
##############################
library(ks)

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
noise <- matrix(rnorm(n^2, 0, sqrt(0.01)), nrow = n, ncol = n)
image_mat <- checker + noise

x <- vec(checker)
y <- vec(image_mat)

sigma2_hat <- 0.01
alpha_hat <- 1.15/sigma2_hat

######### Display the checkerboard matrix as an image

# image(checker, col = gray.colors(4, start = 0, end = 1), axes = FALSE)
# image(image_mat, col = gray.colors(4, start = 0, end = 1), axes = FALSE)
