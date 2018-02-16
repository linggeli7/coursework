library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

wine <- read.delim("~/Documents/2017/Bayesian/wine.txt")
twins.dat <- read.csv("~/Documents/2017/Bayesian/twins.dat.txt", na.strings=".")

qplot(wine$WINE, wine$MORTALITY)

X <- wine$WINE
Y <- wine$MORTALITY
bayes.fit <- stan(file='/Users/linggeli/Documents/2017/Bayesian/wine.stan', chains=1)
stan_dens(bayes.fit, c('beta0', 'beta1'))
print(bayes.fit, c('beta0', 'beta1'))

qplot(twins.dat$DLHRWAGE, bins=20)
qplot(twins.dat$DEDUC1)
qplot(twins.dat$DEDUC1, twins.dat$DLHRWAGE)

qplot(as.factor(twins.dat$WHITEL), twins.dat$HRWAGEL, geom='boxplot') + coord_flip()
qplot(as.factor(twins.dat$MALEL), twins.dat$HRWAGEL, geom='boxplot') + coord_flip()
qplot(twins.dat$AGE, log(twins.dat$HRWAGEL))
qplot(twins.dat$EDUCL, log(twins.dat$HRWAGEL))

# twin 1
twin1 <- cbind(twins.dat$HRWAGEL, twins.dat$AGE, twins.dat$EDUCL,
           twins.dat$MALEL, twins.dat$WHITEL)
twin1 <- twin1[complete.cases(twin1), ]
twin1 <- as.data.frame(twin1)
colnames(twin1) <- c('HRWAGEL', 'AGE', 'EDUCL', 'MALEL', 'WHITEL')

logwage <- log(twin1$HRWAGEL)
age <- twin1$AGE / 10
education <- twin1$EDUCL
male <- twin1$MALEL
white <- twin1$WHITEL

twin1.fit <- stan(file='/Users/linggeli/Documents/2017/Bayesian/twin1.stan', chains=1)
stan_dens(twin1.fit, c('beta0', 'beta1', 'beta2', 'beta3', 'beta4'))
print(twin1.fit, c('beta0', 'beta1', 'beta2', 'beta3', 'beta4'))

# twindiff
twindiff <- cbind(twins.dat$DLHRWAGE, twins.dat$AGE, twins.dat$DEDUC1)
twindiff <- twindiff[complete.cases(twindiff), ]
twindiff <- as.data.frame(twindiff)
colnames(twindiff) <- c('DLHRWAGE', 'AGE', 'DEDUC1')

logwage <- twindiff$DLHRWAGE
age <- twindiff$AGE / 10
education <- twindiff$DEDUC1

twindiff.fit <- stan(file='/Users/linggeli/Documents/2017/Bayesian/twindiff.stan', chains=1)
stan_dens(twindiff.fit, c('beta0', 'beta1', 'beta2'))
print(twindiff.fit, c('beta0', 'beta1', 'beta2'))