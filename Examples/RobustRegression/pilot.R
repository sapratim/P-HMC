################################################################################
## Pilot Run using RWM
################################################################################
set.seed(342)
library(Rcpp)
source("robustreg_data.R")
sourceCpp("pre_robustreg_functions.cpp")

MAP <- map_estimate(B, y, alpha, nu, sigma, w_truth)
w_start <- MAP + rnorm(length(MAP), 0, 0.01)

cat("--- Pilot RWM (estimating posterior diagonals) ---\n")
pilot_iter <- 2e6
pilot_burn <- 1e3

pilot_run <- rwm_cpp(B, y,
                     iter   = pilot_iter,
                     h      = 0.1,        # step size on the *preconditioned* scale
                     start  = w_start,
                     alpha  = alpha,
                     sigma  = sigma,
                     nu     = nu,
                     blather = TRUE)

pilot_samples <- pilot_run$samples[(pilot_burn + 1):pilot_iter, ]
post_var_diag <- apply(pilot_samples, 2, var)

print(summary(post_var_diag))

save("post_var_diag", file = "marginal_vars.RData")
