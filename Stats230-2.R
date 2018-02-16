# pima indians
X <- pima.indians.diabetes[, c(3, 6, 8)]
X <- as.matrix(X)
X <- scale(X)
X <- cbind(rep(1, 768), X)
y <- as.matrix(as.numeric(pima.indians.diabetes[, 9]))

glm(y ~ 0 + X, family=binomial)

# p for one observation
Pi <- function(Xi, beta) {
  temp <- exp(Xi %*% beta)
  return(temp / (1 + temp))
}

# p vector
Prob <- function(X, beta) {
  return(apply(X, 1, Pi, beta))
}

# log-likelihood
logL <- function(X, y, beta) {
  p <- Prob(X, beta)
  return(sum(y * log(p) + (1 - y) * log(1 - p)))
}

# gradient
Gradient <- function(X, y, beta) {
  p <- Prob(X, beta)
  k <- dim(X)[2]
  G <- rep(NA, k)
  G[1] <- sum(y - p)
  for (j in 2:k) {
    G[j] <- sum((y - p) * X[, j])
  }
  G <- as.matrix(G)
  return(G)
}

# drop AFTER CENTER AND SCALE
GradientDescent <- function(X, y) {
  # initialize
  d <- 0.5
  a <- 0.25
  k <- dim(X)[2]
  beta <- as.matrix(rep(0, k))
  G <- as.matrix(rep(1, k))
  # convergence criterion
  while(norm(G, type='2') > 0.001) {
    l <- -1.0 * logL(X, y, beta)
    G <- -1.0 * Gradient(X, y, beta)
    D <- -1.0 * G
    # new beta
    t <- 1.0
    betaNew <- beta + t * D
    lNew <- -1.0 * logL(X, y, betaNew)
    # backtrack
    while((lNew > as.numeric(l + a * t * t(G) %*% D)) | is.na(lNew)) {
      t <- d * t
      betaNew <- beta + t * D
      lNew <- -1.0 * logL(X, y, betaNew)
    }
    # move
    beta <- betaNew
    print(beta)
  }
  return(beta)
}

# hessian
Hessian <- function(X, y, beta) {
  p <- Prob(X, beta)
  k <- dim(X)[2]
  H <- matrix(NA, nrow=k, ncol=k)
  for (j in 1:k) {
    for (h in j:k) {
      H[j, h] <- -1.0 * sum(p * (1 - p) * X[, j] * X[, h])
      H[h, j] <- H[j, h]
    }
  }
  return(H)
}

# leibniz
Newton <- function(X, y) {
  # initialize
  k <- dim(X)[2]
  beta <- as.matrix(rep(0, k))
  a <- 0.25
  lambda2 <- 1
  # convergence criterion
  while(lambda2 > 0.001) {
    l <- -1.0 * logL(X, y, beta)
    G <- -1.0 * Gradient(X, y, beta)
    H <- -1.0 * Hessian(X, y, beta)
    D <- -1.0 * solve(H) %*% G
    lambda2 <- -1.0 * t(G) %*% D
    # new beta
    t <- 1.0
    betaNew <- beta + t * D
    lNew <- -1.0 * logL(X, y, betaNew)
    # backtrack
    while(lNew > as.numeric(l + a * t * t(G) %*% D)) {
      t <- 0.8 * t
      betaNew <- beta + t * D
      lNew <- -1.0 * logL(X, y, betaNew)
    }
    # move
    beta <- betaNew
    step <- t * D
  }
  return(beta)
}

# z response
zk <- function(X, y, beta, p) {
  z <- as.matrix(X %*% beta + (y - p) / (p * (1 - p)))
  return(z)
}

# weight vector
Wk <- function(p) {
  return(diag(p * (1 - p)))
}

# iterative
IterativeWeight <- function(X, y) {
  k <- dim(X)[2]
  beta <- as.matrix(rep(0, k))
  step <- 1
  while(step > 0.001) {
    p <- Prob(X, beta)
    W <- Wk(p)
    z <- zk(X, y, beta, p)
    betaOld <- beta
    beta <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% z
    step <- norm(beta - betaOld, type='2')
  }
  return(beta)
}