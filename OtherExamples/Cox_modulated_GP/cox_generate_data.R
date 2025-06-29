set.seed(1)
source("cox_functions.R")
load("estimated-cov.RData")

N0 <- 10 #Number of observations
mu <- 46.47 #obtained numerically

#Generating data
ns <- rpois(N0, mu) # list of values of n for each observation
x <- list() #list of data

for(j in 1:N0)
{
  temp <- numeric(ns[j])
  for(i in 1:ns[j]) 
  { 
    accept <- 0
    while(!accept)
    {
      y <- runif(1, min = 0, max = 50)
      if(runif(1, 0, 1) < lam1(y)/2.1)
      {
        temp[i] <- y
        accept <- 1
      }
    }
  }
  x[[j]] <- temp
}

#Linear interpolation
m <- 100
delta_m <- 50/(m-1) #the support is [0, 50]
t <- numeric(m)
c <- numeric(m)
phis <- numeric(m)
c[1] <- 0.5*delta_m
c[m] <- 0.5*delta_m
for(i in 1:m)
{
  t[i] <- (i-1)*delta_m
  if(i > 1 && i < m)
  {
    c[i] <- delta_m
  }
}


## Parameters of the covariance matrix obtained after tuning
cov <- matrix(0, m, m)
sigma <- sqrt(cov_params[1])
l <- cov_params[2]

for(i in 1:m)
{
  for(j in 1:m)
  {
    cov[i,j] <- (sigma^2)*(exp(-(t[i] - t[j])^2/(2*l^2) ) )
  }
}
cov[1,m]
cov <- cov

xn <- list(ns, x, c, t, cov, m, delta_m, N0, mu)
save(xn, file = "cox-data.RData")
