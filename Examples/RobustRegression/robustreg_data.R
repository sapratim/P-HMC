####################
### For 1D synthetic signal

set.seed(123)
N <- 512
M <- 120
non_zero_t_count <- 20
t_index <- sample(c(1:512), size = non_zero_t_count)
w_truth <- numeric(length = N)
w_truth[t_index] <- rnorm(length(t_index))
Phi_mat <- t(apply(matrix(0, nrow = M, ncol = N), 1, function(x) {x <- rnorm(N)
                                                               x <- x/sqrt(sum(x^2))}))
add_noise <- rt(M, 5)
y <- Phi_mat %*% w_truth + add_noise
