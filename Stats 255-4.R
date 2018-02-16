(qnorm(0.975)-qnorm(0.1))^2/(log(0.75)^2*0.25)
508/(0.5*0.3843+0.5*0.3069)
6*4+7+9+10*2+11+13+16+17+19+20+22+23+25+32+32+34+35
17+15+5+1+22+4+1+8+2+11+2+8+4+5+8+12+12+3+8+11+23
c(qgamma(0.025, shape=12, rate=404), qgamma(0.975, shape=12, rate=404))
c(qgamma(0.025, shape=24, rate=227), qgamma(0.975, shape=24, rate=227))
library(ggplot2)
lambda1 <- rgamma(1000, shape=12, rate=404)
t1 <- rexp(1000, lambda1)
lambda2 <- rgamma(1000, shape=24, rate=227)
t2 <- rexp(1000, lambda2)
qplot(t1, geom="histogram") + scale_x_continuous(limits=c(0,200))
qplot(t2, geom="histogram") + scale_x_continuous(limits=c(0,200))