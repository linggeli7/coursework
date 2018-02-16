# Data analysis
access.clean <- accessSurv[which(accessSurv$start == 0),]
summary(access.clean)
access.clean$female <- factor(access.clean$female)
access.clean$racegrp <- factor(access.clean$racegrp)
access.clean$smokegrp <- factor(access.clean$smokegrp)
access.clean$dx_diab <- factor(access.clean$dx_diab)
access.clean$fat <- ifelse(access.clean$bmi > 24.4, 1, 0)
access.clean$old <- ifelse(access.clean$age_ssd > 66, 1, 0)
access.clean$health <- ifelse(access.clean$ser_alb > 3.5, 1, 0)

library(stargazer)
stargazer(access.clean, median = TRUE,
          iqr = TRUE)
A <- access.clean[which(access.clean$acctype==3),]
summary(A)

fit=survfit(Surv(access.clean$stop, access.clean$failed) ~ access.clean$female)
plot(fit, mark.time=FALSE, lty=1:2, xlab="Time", ylab="Survival", lwd = 1.5)
legend( 700, .5, lty=1:3, legend=c("Male","Female"), bty="n" )

# COX
fit <- coxph(Surv(stop, failed) ~ factor(acctype) + age_ssd + female + factor(racegrp)
             + bmi + factor(smokegrp) + dx_diab + ser_alb, data = access.clean)
summary(fit)
# COX revision
access.clean$event <- ifelse(access.clean$revised < access.clean$stop, 1, access.clean$failed)
fit <- coxph(Surv(revised, event) ~ factor(acctype) + age_ssd + female + factor(racegrp)
             + bmi + factor(smokegrp) + dx_diab + ser_alb, data = access.clean)
summary(fit)

# problem 2
survdiff(Surv(Age, Death) ~ Male+strata(ID), data=CHD)

# problem 3
set.seed(12345)
larynx$t2death <- larynx$t2death + runif(90, -0.01, 0.01)
library(survival)
cox.fit <- coxph(Surv(t2death, death) ~ age + as.factor(stage), data=larynx)
X <- model.matrix(cox.fit)
Y <- list(y=larynx$t2death, delta=larynx$death)

beta <- c(0.0188,0.1293,0.6413,1.7169)
beta <- c(0,0,0,0)

partLklhd <- function(Y, X, beta) {
  # everybody at risk at time zero
  riskset.index <- rep(1,length(Y$y))
  # log likelihood
  l <- 0
  # score vector
  U <- rep(0, length(beta))
  # information matrix
  I <- matrix(0, length(beta), length(beta))
  for (n in 1:length(Y$y)) {
    # at each time
    t <- sort(Y$y)[n]
    # find jth individual
    j <- match(t, Y$y)
    # if death
    if (Y$delta[j] == 1) {
      # log likelihood
      l <- l + beta %*% X[j,] - log(sum(exp(X %*% beta)[riskset.index == 1]))
      # score
      for (k in 1:length(beta)) {
        # Xjk bar
        xjk.bar <- sum((X[,k]*exp(X %*% beta))[riskset.index == 1])
        xjk.bar <- xjk.bar/sum(exp(X %*% beta)[riskset.index == 1])
        U[k] <- U[k] + X[j,k] - xjk.bar
      }
      # information matrix
      for (h in 1:length(beta)) {
        for (k in 1:length(beta)) {
          xjk.bar <- sum((X[,k]*exp(X %*% beta))[riskset.index == 1])
          xjk.bar <- xjk.bar/sum(exp(X %*% beta)[riskset.index == 1])
          xjh.bar <- sum((X[,h]*exp(X %*% beta))[riskset.index == 1])
          xjh.bar <- xjh.bar/sum(exp(X %*% beta)[riskset.index == 1])
          term <- sum( ((X[,k]-xjk.bar)*(X[,h]-xjh.bar)*exp(X %*% beta))[riskset.index == 1] )
          term <- term/sum(exp(X %*% beta)[riskset.index == 1])
          I[h,k] <- I[h,k] + term
        }
      }
    }
    # remove from riskset
    riskset.index[j] <- 0
  }
  goodstuff <- list(lnlklhd=l, score=U, inform=I)
  return(goodstuff)
}

# problem 4
newtraph <- function(Y, X, beta, verbose=FALSE, eps=1e-8, iter.max=20) {
  beta0 <- beta
  lnlklhd0 <- partLklhd(Y,X,beta)$lnlklhd
  Delta <- rep(1, length(beta))
  iter <- 0
  prev.l <- lnlklhd0
  prev.Delta <- Delta
  while( sum( abs(Delta)>eps ) >= 1 && iter < iter.max) {
    Lik <- partLklhd(Y,X,beta)
    if (verbose) {
      print(beta)
      print(Lik) 
    }
    if (Lik$lnlklhd < prev.l) {
      Delta <- prev.Delta/2
      beta <- beta - Delta
    } else {
      Delta <- as.vector(solve(Lik$inform) %*% Lik$score)
      beta <- beta + Delta
    }
    iter <- iter + 1
    prev.l <- Lik$lnlklhd
    prev.Delta <- Delta
  }
  if (iter < iter.max) {
    converge <- TRUE
  } else {
    converge <- FALSE
  }
  goodstuff <- list(beta0=beta0, lnlklhd0=lnlklhd0, converge=converge,
                    nbriter=iter, beta=beta, lnlklhd=prev.l)
  return(goodstuff)
}

# Prostate cancer
psaData$event <- ifelse(psaData$inrem == "yes", 1, 0)
psaData$logpsa <- log10(psaData$nadirpsa)
fit <- coxph(Surv(obstime, event) ~ logpsa + ps + bss, data = psaData)
summary(fit)
X <- as.matrix(cbind(psaData$logpsa, psaData$ps, psaData$bss))
Y <- list(y=psaData$obstime, delta=psaData$event)
beta <- c(0,0,0)
newtraph(Y, X, beta)