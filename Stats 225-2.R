library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

accident <- c(24, 25, 31, 31, 22, 21, 26, 20, 16, 22)
death <- c(734, 516, 754, 877, 814, 362, 764, 809, 223, 1066)
rate <- c(0.19, 0.12, 0.15, 0.16, 0.14, 0.06, 0.13, 0.13, 0.03, 0.15)

Y <- accident
summary(Y)
X <- death / rate
X <- (X - mean(X)) / sd(X)

summary(lm(Y ~ X, weights = 1 / sqrt(Y)))

posterior <- function(alpha, beta, mu1, sd1, mu2, sd2) {
  prior <- dnorm(alpha, mu1, sd1) * dnorm(beta, mu2, sd2)
  likelihood <- 0
  lambda <- alpha + beta * X
  if (sum(lambda < 0) == 0) {
    likelihood <- prod(dpois(Y, lambda))
  } 
  return(likelihood * prior)
}

range.alpha <- seq(20, 30, length.out=200)
range.beta <- seq(-10, 0, length.out=200)

density <- matrix(NA, nrow=200, ncol=200)
for (i in 1:200) {
  for (j in 1:200) {
    density[i, j] <- posterior(range.alpha[i], range.beta[j], 25, 5, 0, 3)
  }
}

density <- c(density)
density <- density / sum(density)

pos <- sample(1:40000, 1000, prob=density)

draws <- matrix(NA, nrow=1000, ncol=2)
for (i in 1:1000) {
  pos.alpha <- ceiling(pos[i] / 200)
  draws[i, 1] <- range.alpha[pos.alpha]
  pos.beta <- pos[i] %% 200
  draws[i, 2] <- range.beta[pos.beta]
}

draws <- as.data.frame(draws)
colnames(draws) <- c('alpha', 'beta')

ggplot(data=draws, aes(x=alpha, y=beta)) + 
  geom_point(colour='blue', alpha=0.5) +
  geom_density_2d() 

x1986 <- (546 / 0.06 - mean(death / rate)) / sd(death / rate)

pred <- draws$alpha + x1986 * draws$beta
qplot(pred, geom='density', fill=I('red'), alpha=I(0.5))

chisq.observed <- rep(NA, 1000)
chisq.predictive <- rep(NA, 1000)
for (i in 1:1000){
  lambda <- draws$alpha[i] + X * draws$beta[i]
  Y.pred <- rpois(10, lambda)
  chisq.observed[i] <- sum((Y - lambda) ^ 2 / lambda)
  chisq.predictive[i] <- sum((Y.pred - lambda) ^ 2 / lambda)
}

qplot(chisq.observed, geom='density', fill=I('green'), alpha=I(0.5))
qplot(chisq.predictive, geom='density', fill=I('blue'), alpha=I(0.5))

discrepancy <- as.data.frame(cbind(chisq.observed, chisq.predictive))
ggplot(data=discrepancy, aes(x=chisq.observed, y=chisq.predictive)) + 
  geom_point(colour='blue', alpha=0.5) +
  geom_abline(intercept=0, slope=1, linetype=2) +
  scale_x_continuous(limits=c(0, 40)) +
  scale_y_continuous(limits=c(0, 40))

fit <- stan(file='/Users/linggeli/Documents/2017/Bayesian/passenger.stan', chains=1)
stan_dens(fit, c('alpha', 'beta'))
print(fit, c('alpha', 'beta'))