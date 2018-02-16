library(ggplot2)
library(mvtnorm)

# rejection sampling
n <- 500
c <- 1
sample <- rep(NA, n)
accepted <- 0
for (i in 1:n) {
  # Gamma(1, 1) proposal distribution 
  proposed <- rgamma(1, shape=1, rate=1)
  u <- runif(1, min=0, max=proposed * c)
  if (u < dnorm(proposed, mean=-1, sd=2)) {
    accepted <- accepted + 1
    sample[accepted] <- proposed
  }
}
sample <- sample[1:accepted]
qplot(sample, geom='density', fill=I('red'), colour=I('red'), alpha=I(0.5))

# importance sampling
n <- 5000
proposed <- rgamma(n, shape=1, rate=1)
weights <- dnorm(proposed, mean=-1, sd=2) / dgamma(proposed, shape=1, rate=1)
sum((proposed^2) * weights) / sum(weights)

# metropolis hastings
poissonPostMH <- function(X, y, mu, sigma, initial) {
  nIterations <- 20000
  burnIn <- 2000
  draws <- matrix(0, nrow=nIterations - burnIn, ncol=4)
  beta <- initial
  accepted <- 0
  
  for(i in 1:nIterations) {
    betaProposed <- rep(0, 4)
    for (j in 1:4) {
      betaProposed[j] <- rnorm(1, mean=beta[j], sd=0.1)
    }
    
    postCurrent <- getPosterior(beta, X, y, mu, sigma)
    postProposed <- getPosterior(betaProposed, X, y, mu, sigma)
    
    u <- runif(1)  
    acceptanceRate <- min(1, exp(postProposed - postCurrent))    
    if (u < acceptanceRate) {
      beta <- betaProposed
      accepted <- accepted + 1
    }
    
    if (i > burnIn) {
      draws[i - burnIn, ] <- beta
    }
  }
  print('Acceptance probability: ')
  print(accepted / nIterations)
  return(draws)
}

getPosterior <- function(beta, X, y, mu, sigma) {
  n <- length(y)
  theta <- exp(X %*% as.matrix(beta))
  prior <- 1
  for (i in 1:4) {
    prior <- prior * dnorm(beta[i], mu, sigma)
  }
  loglik <- sum(y * log(theta) - theta)
  posterior <- log(prior) + loglik
  return(posterior)
}

y <- absent$daysabs
X <- cbind(absent$male, absent$math, absent$langarts)
X[, 2] <- (X[, 2] - mean(absent$math)) / sd(absent$math)
X[, 3] <- (X[, 3] - mean(absent$langarts)) / sd(absent$langarts)
n <- length(y)
X <- cbind(rep(1, n), X)

dens <- poissonPostMH(X, y, 0, 3, c(0, 0, 0, 0))
acf(dens[, 1])
summary(dens[ ,1])
hist(dens[, 1])



dens <- poissonSlice(X, y, 0, 3, c(0, 0, 0, 0))

poissonSlice <- function(X, y, mu, sigma, initial) {
  nIterations <- 2000
  draws <- matrix(0, nrow=nIterations, ncol=4)
  beta <- initial
  
  for (i in 1:nIterations) {
    for (j in 1:4) {
      z <- getPosterior(beta, X, y, mu, sigma) - rexp(1)
      slice <- stepout(beta, j, z, 0.5, 10)
      L <- slice[1]
      R <- slice[2]
      beta <- shrink(beta, j, z, L, R)
    }
    draws[i, ] <- beta
  }
  
  return(draws)
}


stepout <- function(beta, j, z, w, m) {
  currentParam <- beta[j]
  
  u <- runif(1)
  L <- currentParam - w * u
  R <- L + w
  v <- runif(1)
  J <- floor(m * v)
  K <- (m - 1) - J
  
  betaL <- beta
  betaL[j] <- L
  while (J > 0 && z < getPosterior(betaL, X, y, mu, sigma)) {
    L <- L - w
    J <- J - 1
  }
  
  betaR <- beta
  betaR[j] <- R
  while (K > 0 && z < getPosterior(betaR, X, y, mu, sigma)) {
    R <- R + w
    K <- K - 1
  }
  
  return(c(L, R))
}

shrink <- function(beta, j, z, L, R) {
  currentParam <- beta[j]
  betaProposed <- beta

  u <- runif(1)
  newParam <- L + u * (R - L)
  betaProposed[j] <- newParam
  
  while (z > getPosterior(betaProposed, X, y, mu, sigma)) {
    if (newParam < currentParam) {
      L <- newParam
    } else {
      R <- newParam
    }
    u <- runif(1)
    newParam <- L + u * (R - L)
    betaProposed[j] <- newParam
  }
  
  return(betaProposed)
}

postGradient <- function(beta, X, y, mu, sigma) {
  n <- length(y)
  theta <- exp(X %*% as.matrix(beta))
  grad <- rep(0, 4)
  for (j in 1:4) {
    grad[j] <- sum(y * X[, j] / theta - X[, j]) - (beta[j] - mu) / sigma^2
  }
  return(grad)
}



poissonHMC <- function(X, y, mu, sigma, initial) {
  eps <- 0.001
  L <- 10
  nIterations <- 2000
  draws <- matrix(NA, nrow=nIterations, ncol=4)
  beta <- initial
  accepted <- 0
  for (i in 1:nIterations) {
    momentum <- as.vector(rmvnorm(1, sigma=diag(4)))
    energy <- -getPosterior(beta, X, y, mu, sigma) - dmvnorm(momentum, sigma=diag(4), log=TRUE)
    
    betaProposed <- beta
    p <- momentum
    
    p <- p + 0.5 * eps * postGradient(betaProposed, X, y, mu, sigma)
    
    for (l in 1:L)
    {
      betaProposed <- betaProposed + eps * p
      if (l != L) {
        p <- p + eps * postGradient(betaProposed, X, y, mu, sigma)
      }
    }
    
    p <- p + 0.5 * eps * postGradient(betaProposed, X, y, mu, sigma)
    p <- -p
    
    energyProposed <- -getPosterior(betaProposed, X, y, mu, sigma) - dmvnorm(p, sigma=diag(4), log=TRUE)
    
    acceptProb <- min(1, exp(energy - energyProposed))
    
    if (is.finite(acceptProb) && runif(1) < acceptProb) {
      beta <- betaProposed
      accepted <- accepted + 1
    }
    draws[i, ] <- beta
  }
  print('Acceptance probability: ')
  print(accepted / nIterations)
  return(draws)
}

dens <- poissonHMC(X, y, 0, 3, c(0, 0, 0, 0))

normalGibbs <- function(y, mu0, tau0, v0, sigma0, initial) {
  n <- length(y)
  mu <- initial[1]
  sigma <- initial[2]
  nIter <- 1000
  burnIn <- 100
  draws <- matrix(NA, nrow=nIter - burnIn, ncol=2)
  for (i in 1:nIter) {
    m <- (mu0 / tau0 ^ 2 + sum(y) / sigma ^ 2) / (1 / tau0 ^ 2 + n / sigma ^ 2)
    s <- 1 / sqrt(1 / tau0 ^ 2 + n / sigma ^ 2)
    mu <- rnorm(1, mean=m, sd=s)
    v <- sum((y - mu) ^ 2) / n
    a <- v0 + n
    b <- (v0 * sigma0 ^ 2 + v * n) / (v0 + n)
    sigmaInv <- rgamma(1, 0.5 * a, 0.5 * a * b)
    sigma <- 1 / sigmaInv
    if (i > burnIn) {
      draws[i - burnIn, ] <- c(mu, sigma)
    }
  }
  return(draws)
}

dens <- normalGibbs(y, 40, 5, 1, 1, c(50, 10))
qplot(dens[, 1], geom='density', fill=I('orange'), colour=I('orange'), alpha=I(0.5))
qplot(dens[, 2], geom='density', fill=I('orange'), colour=I('orange'), alpha=I(0.5))