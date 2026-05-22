####################
### For 1D synthetic signal

set.seed(100)
N <-  512
M <- 120

non_zero_t_count <- 20
t_index <- sort(sample(c(1:N), size = non_zero_t_count))
w_truth <- numeric(length = N)
w_truth[t_index] <- rt(length(t_index), 3)
Phi_mat <- t(apply(matrix(0, nrow = M, ncol = N), 1, function(x) {x <- rnorm(N)
                                                               x <- x/sqrt(sum(x^2))}))
nu <- 5
sigma <- .01                 # small but nonzero noise
add_noise <- sigma * rt(M, nu)
y <- Phi_mat %*% w_truth + add_noise

# Setting model parameters
alpha <- 10
B     <- Phi_mat