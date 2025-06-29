#install.packages(c("nloptr", "quadprog", "mvtnorm", "TruncatedNormal", "plot3D"))
library("lineqGPR")
library(MASS)
library(ggplot2)
source("cox_functions.R")

set.seed(1)

#Generating data
mu <- 46.47 #numerically obtained
n <- rpois(1, mu)
x <- numeric(n)

for(i in 1:n) 
{
  accept <- 0
  while(!accept)
  {
    y <- runif(1, min = 0, max = 50)
    if(runif(1, 0, 1) < lam1(y)/2.1)
    {
      x[i] <- y
      accept <- 1
    }
  }
}

#Interpolation
m <- 100
delta_m <- 50/(m-1)
t <- numeric(m)
c <- numeric(m)
phis <- numeric(m)
c[1] <- 0.5*delta_m
c[m] <- 0.5*delta_m


# Finding covariance parameters

sigfun <- function(x) return(lam1(50*x))
x <- seq(0, 1, 0.001); y <- sigfun(x)
DoE <- splitDoE(x, y, DoE.idx = seq(1, length(x), length = 10))
#### GP with nearly inactive boundedness constraints [-10,10] 
#### creating the "lineqGP" model
model <- lineqGPR::create(class = "lineqGP", x = DoE$xdesign, y = DoE$ydesign,
                          constrType = "boundedness")
model$localParam$m <- 100 # changing the (default) number of knots
model$bounds <- c(0, 2.1) # changing the (default) bounds
# sampling from the model
sim.model <- simulate(model, nsim = 1e3, seed = 1, xtest = DoE$xtest)
# ggplotLineqGPModel <- ggplot(sim.model)


set.seed(7)
#### GP with both boundedness constraints ####
model <- lineqGPR::create(class = "lineqGP", x = DoE$xdesign, y = DoE$ydesign,
                          constrType = c("boundedness"))
model$localParam$m <- 100 # changing the (default) number of knots
# modifying the bounds for first arg of "constrType" (boundedness)
model$bounds <- c(0,2.1)

sim.model <- simulate(model, nsim = 1e2, seed = 1, xtest = DoE$xtest)
message("Initial covariance parameters: ", model$kernParam$par[1],
        ", ", model$kernParam$par[2])
# ggplotLineqGPModel <- ggplot(sim.model, bounds = c(model$bounds))


# estimating the covariance parameter via MLE
model2 <- lineqGPOptim(model,
                       opts = list(algorithm = "NLOPT_LD_MMA",
                                   print_level = 0,
                                   ftol_abs = 1e-3, maxeval = 40,
                                   check_derivatives = TRUE,
                                   parfixed = c(FALSE, FALSE)),
                       lb = c(0.1, 0.01), ub = c(100, 0.8))
message("Estimated covariance parameters via MLE: ",
        model2$kernParam$par[1], ", ", model2$kernParam$par[2])
# evaluating the "optimal" model
sim.model2 <- simulate(model2, nsim = 1e2, seed = 1, xtest = DoE$xtest)
# ggplotLineqGPModel <- ggplot(sim.model2, bounds = c(model$bounds))

cov_params <- model2$kernParam$par
save(cov_params, file = "estimated-cov.RData")
