lm.QR <- function(X, y) {
  mn <- dim(X)
  Q <- matrix(0, nrow=mn[1], ncol=mn[2])
  R <- matrix(0, nrow=mn[2], ncol=mn[2])
  v <- as.matrix(X[, 1])
  R[1, 1] <- norm(v, type='2')
  Q[, 1] <- v / R[1, 1] 
  for (j in 2:mn[2]) {
    v <- as.matrix(X[, j])
    for (i in 1:(j - 1)) {
      R[i, j] <- t(Q[, i]) %*% X[, j]
      v <- v - R[i, j] * Q[, i]
    }
    R[j, j] <- norm(v, type='2')
    Q[, j] <- v / R[j, j]
  }
  return(solve(R) %*% t(Q) %*% y)
}

lm.Householder <- function(X, y) {
  np <- dim(X)
  U <- matrix(0, nrow=np[1], ncol=np[2])
  C <- cbind(X, y)
  for (k in 1:(np[2])) {
    w <- as.matrix(C[k:np[1], k])
    w[1] <- w[1] - norm(w, type='2')
    u <- w / norm(w, type='2')
    U[k:np[1], k] <- u
    C[k:np[1], k:(np[2] + 1)] <- C[k:np[1], k:(np[2] + 1)] - 2 * u %*% t(u) %*% C[k:np[1], k:(np[2] + 1)]
  }
  A <- C[1:np[2], 1:np[2]]
  b <- C[1:np[2], np[2] + 1]
  A[lower.tri(A)] <- 0
  return(solve(A) %*% b)
}

lm.Jacobi <- function(X, y, epsilon=0.01) {
  A <- t(X) %*% X
  b <- t(X) %*% y
  P <- diag(diag(A))
  p <- dim(X)[2]
  beta <- as.matrix(rep(0, p))
  step <- 1
  while(step > epsilon) {
    current <- beta
    beta <- (diag(p) - solve(P) %*% A) %*% beta + solve(P) %*% b
    step <- norm(current - beta, type='2')
  }
  return(beta)
}

lm.Iterative <- function(Xn, yn, Xk, yk) {
  k <- dim(Xk)[1]
  beta0 <- lm.QR(Xn, yn)
  Rn <- t(Xn) %*% Xn
  A <- solve(Rn)
  BCD <- A %*% t(Xk) %*% solve(diag(k) + Xk %*% A %*% t(Xk)) %*% Xk %*% A
  R <- A - BCD
  beta <- beta0 + R %*% t(Xk) %*% (yk - Xk %*% beta0)
  return(beta)
}
