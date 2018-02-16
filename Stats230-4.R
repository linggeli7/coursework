library(ggplot2)

zeroPoisson <- function(y) {
  nIter <- 10
  theta <- 0.5
  lambda <- 1
  # total number of observations
  n <- length(y)
  # number of zero observations
  j <- sum(y == 0)
  for (i in 1:nIter) {
    # proportion of Poisson zero among zero observations
    prop <- theta / (theta + (1- theta) * exp(-lambda))
    # set lambda to mean of the rest
    lambda <- sum(y) / (n - j * prop)
    # set theta to proprotion of Poisson zero among all observations
    theta <- j * prop / n
  }
  return(c(theta, lambda))
}

zeroPoisson(roaches$count)

qplot(roaches$count, geom='histogram', binwidth=1)

fit <- glm(daysabs~., family='poisson', data=absent)
beta <- as.matrix(fit$coefficients)
X <- cbind(rep(1, 159), absent$male, absent$math, absent$langarts)
n <- 300
dist <- matrix(NA, nrow=n, ncol=4)
p <- exp(X %*% beta)
for (i in 1:n) {
  y <- rpois(159, lambda=p)
  dist[i, ] <- glm(y ~ 0 + X, family='poisson')$coefficients
}

temp <- data.frame(x=seq(2.2, 3, 0.01), p=dnorm(seq(2.2, 3, 0.01), 2.568596, 0.083403))
qplot(dist[, 1], geom='density', fill=I('#0072B2'), alpha=I(0.5), colour=I(NA)) +
  geom_line(data=temp, aes(x=x, y=p), colour='red')