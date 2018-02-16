library(ggplot2)

gibbs <- function(obs, rats, n.draws, start) {
  groups <- unique(rats)
  J <- length(groups)
  I <- length(RatWeight$weight) / J
  draws <- matrix(NA, nrow=n.draws, ncol=34)
  draws[1, ] <- start
  for (t in 2:n.draws) {
    # sample tau2
    theta <- draws[t - 1, 33]
    draws[t, 34] <- 1 / rgamma(1, 0.5 * (J + 1), 0.5 * (sum((draws[t - 1, 1:J] - theta)^2) + 0.05))
    # sample theta
    tau2 <- draws[t, 34]
    variance <- 1 / (J / tau2 + 1 / 1000^2)
    draws[t, 33] <- rnorm(1, sum(draws[t - 1, 1:J]) * variance / tau2 , sqrt(variance))
    # sample sigma2.j
    for (j in 1:J) {
      mu.j <- draws[t - 1, j]
      draws[t, j + J] <- 1 / rgamma(1, 0.5 * (I + 1), 0.5 * (sum((obs[rats == groups[j]] - mu.j)^2) + 0.05))
    }
    # sample mu.j
    for (j in 1:J) {
      sigma2.j <- draws[t, j + J]
      theta <- draws[t, 33]
      tau2 <- draws[t, 34]
      variance <- 1 / (I / sigma2.j + 1 / tau2)
      mu <- sum(obs[rats == groups[j]]) / sigma2.j + theta / tau2
      mu <- mu * variance
      draws[t, j] <- rnorm(1, mu, sqrt(variance))
    }
  }
  return(draws)
}

pos <- gibbs(RatWeight$weight, RatWeight$Rat, 3000, rep(100, 34))

plot.ts(pos[, 1])
plot.ts(pos[, 17])
plot.ts(pos[, 33])
plot.ts(pos[, 34])

betas <- pos[2001:3000, 1:16]
temp <- cbind(c(betas), as.factor(sort(rep(1:16, 1000))))
temp <- as.data.frame(temp)
colnames(temp) <- c('beta', 'rat')
ggplot(temp) + geom_density(aes(x=beta), fill='red', alpha=0.5) +
  facet_wrap( ~ rat, ncol=4)

library(ggplot2)
prior <- c(rnorm(900, mean=1, sd=0.5), rnorm(100, mean=-1, sd=0.5))
qplot(prior, geom='density', fill=I('red'), alpha=0.5)
posterior <- c(rnorm(680, mean=(4 - 2.5) / 14, sd=sqrt(1/ 14)), 
               rnorm(320, mean=(-4 - 2.5) / 14, sd=sqrt(1 / 14)))
qplot(posterior, geom='density', fill=I('blue'), alpha=0.5)
