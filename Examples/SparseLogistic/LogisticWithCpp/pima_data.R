################### PIMA Indian diabetes dataset ###################

data <- MASS::Pima.tr
x <- as.matrix(data[,c(1:7)])
y <- as.matrix(ifelse(data$type == "Yes", 1, 0))
colnames(x) <- NULL
colnames(y) <- NULL
alpha <- 2

colnames(data) = c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "y")