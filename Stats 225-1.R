library(ggplot2)

######################## BETA BINOMIAL ########################

# TREATMENT
# Beta(a, b) prior
# Beta(1, 1) is uniform
aT <- 1
bT <- 1
betaPriorT <- rbeta(1000, aT, bT)
qplot(betaPriorT, geom='density', fill=I('green'), alpha=I(0.5), 
      main='Treatment Prior')

# Beta(a + successes, b + failures) posterior
nT <- 680
successT <- 22
failT <- nT - successT
betaPosT <- rbeta(1000, aT + successT, bT + failT)
qplot(betaPosT, geom='density', fill=I('blue'), alpha=I(0.5), 
      main='Treatment Posterior')

# CONTROL
aC <- 1
bC <- 1
betaPriorC <- rbeta(1000, aC, bC)
qplot(betaPriorC, geom='density', fill=I('green'), alpha=I(0.5), 
      main='Control Prior')

# Beta(a + successes, b + failures) posterior
nC <- 674
successC <- 39
failC <- nC - successC
betaPosC <- rbeta(1000, aC + successC, bC + failC)
qplot(betaPosC, geom='density', fill=I('blue'), alpha=I(0.5), 
      main='Control Posterior')

# ODDS RATIO
oddsT <- betaPosT / (1 - betaPosT) 
oddsC <- betaPosC / (1 - betaPosC) 
oddsRatio <- oddsT / oddsC
summary(oddsRatio)
qplot(oddsRatio, geom='density', fill=I('red'), alpha=I(0.5), 
      main='Posterior Odds Ratio')


######################## CAUCHY ########################
likelihood <- function(Y, mu, sigma) {
# Cauchy likelihood
  return(prod(sigma / (pi * (sigma^2 + (Y - mu)^2))))
}

prior <- function(mu) {
# Uniform prior
  return(1)
}

posterior <- function(Y, mu, sigma) {
# Unnormalized posterior
  return(likelihood(Y, mu, sigma) * prior(mu))
}

summary(newcomb)

supp <- seq(20, 35, 0.01)

# standardized posterior density
muPos <- sapply(supp, posterior, Y=newcomb$V1, sigma=5)
muPos <- muPos / sum(muPos)

qplot(x=supp, y=muPos, geom='path', colour=I('blue'), 
      main='mu Posterior Density')

# posterior draws
muDraw <- sample(supp, 1000, replace=TRUE, prob=muPos)
qplot(muDraw, geom='histogram', bins=50, 
      main='mu Posterior Draw Histogram')

# posterior predictive density
yPred <- rcauchy(1000, location=muDraw, scale=5)
yPred <- yPred[yPred < 100 & yPred > -100]
qplot(yPred, geom='density', fill=I('red'), alpha=I(0.5), 
      main='y Posterior Predictive Density')