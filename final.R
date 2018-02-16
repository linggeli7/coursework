# Emax model
emax <- function(C, E0, IC50, M) {
  return(E0/(1+(C/IC50)^M))
}

# Drug A similar to Rosner et al
E0 <- 100
IC50_A <- 40
M_A <- 1.5

C <- seq(0, 130, 0.1)
EY_A <- sapply(C, emax, E0, IC50_A, M_A)

X <- rep(c(0, 8, 16, 32, 64, 128), 4)
Y <- sapply(X, emax, E0, IC50_A, M_A)
Y <- Y + rnorm(24, 0, 4)
dat <- data.frame(dose=X, count=Y)

qplot(C, EY_A, geom='path', xlab='', ylab='') +
  geom_point(data=dat, aes(x=dose, y=count)) +
  scale_y_continuous(limits=c(0, 120)) +
  theme_bw()


E <- 100
C <- 30
M <- 0.8

Fisher <- function(x) {
  candy <- matrix(nrow = 2, ncol = 2)
  d1 <- E * (-M * x ^ M / C ^ (M + 1)) / (1 + (x / C) ^ M) ^ 2
  d2 <- E * (x / C) ^ M * log(x / C) / (1 + (x / C) ^ M) ^ 2
  candy[1, 1] <- d1 * d1
  candy[1, 2] <- d1 * d2
  candy[2, 1] <- d1 * d2
  candy[2, 2] <- d2 * d2
  return(candy)
}

D <- seq(21, 40, 1)

criterion <- rep(NA, 20)

for (i in 1:20) {
  Xd <- D[i]
  X <- c(0.25 * Xd, 0.5 * Xd, Xd, 2 * Xd, 4 * Xd)
  cov <- matrix(0, 2, 2)
  for (j in 1:5) {
    cov <- cov + Fisher(X[j])
  }
  criterion[i] <- det(solve(cov))
}

which.min(criterion)
which.max(criterion)
summary(criterion)
max(criterion) / min(criterion)

# create 
makeX <- function(conc, k) {
  X <- rep(conc, k)
  X <- sort(X)
  return(X)
}

generateY <- function(X, E0, IC50, M, l, s) {
  EY <- sapply(X, emax, E0, IC50, M)
  epsilon <- rnorm(length(X), 0, s)
  Y <- EY+epsilon*(EY^l)
  return(Y)
}

results <- rep(NA, 20)

## loop through design space
for (i in 1:20) {
  Xd <- D[i]
  X <- c(0.25 * Xd, 0.5 * Xd, Xd, 2 * Xd, 4 * Xd)
  conc <- makeX(X, 3)
  utility <- rep(NA, 300)
  ## sample data and fit model
  for (j in 1:300) {
    IC50 <- lognormal()
    M <- lognormal()
    y <- generateY(X, E0, IC50, M, 0, 4)
    fit <- sampling(sm, data=list(y=y, conc=conc), chains=1, show_messages=FALSE)
    distC <- extract(fit, c('C'))$C
    distM <- extract(fit, c('M'))$M
    utility[j] <- var(distC) * var(distM) - cov(distC, distM) * cov(distC, distM)
  }
  results[i] <- mean(utility)
}

qplot(21:40, criterion, xlab='', ylab='') + theme_bw()

IC50 <- rlnorm(2000, 3.4, 0.3)
qplot(IC50, geom="density", fill=I('red'), alpha=I(0.5), xlab='', ylab='') + theme_bw()

m <- rlnorm(1000, 0, 0.3)
qplot(m, geom="density", fill=I('blue'), alpha=I(0.5), xlab='', ylab='') + theme_bw()