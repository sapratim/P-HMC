################### PIMA Indian diabetes dataset ###################
data <- MASS::Pima.tr
x <- as.matrix(data[,c(1:7)])
y <- as.matrix(ifelse(data$type == "Yes", 1, 0))
colnames(x) <- NULL
colnames(y) <- NULL
alpha <- 2


## Lipschitz coefficient for the smooth part
xtx <- t(x) %*% x

C_f <- eigen(xtx)$values[1]/4

1/C_f