attach(splineExample)
plot(y~x1)
plot(y~x2)
plot(x1~x2)

# Smoothing splines
I <- sample(1:100,100)
fold <- rep(NA,100)
for (i in 1:5) {
  fold[I[(20*i-19):(20*i)]] <- i
}
N <- length(unique(x1))
RSS <- rep(0,51)
for (n in 2:51) {
  for (i in 1:5) {
    fit <- smooth.spline(x1[fold!=i],y[fold!=i],df=n)
    fit.y <- predict(fit,x1[fold==i])$y
    RSS[n] <- RSS[n]+sum((y[fold==i]-fit.y)^2)
  }
}
plot(RSS[2:51],type='l',xlab='Degrees of freedom',ylab='SSE',main='5-fold Cross Validation')

fit1 <- smooth.spline(x1,y,df=5)
fit2 <- smooth.spline(x1,y,df=10)
fit3 <- smooth.spline(x1,y,df=20)
plot(y~x1,main='Fitted Splines')
lines(fit1,lty=2)
lines(fit2,lty=1)
lines(fit3,lty=3)
legend(x=-3,y=2,c('df=5','df=10','df=20'),lty=c(2,1,3),bty='n')

# Natural cubic splines
library('splines')
fit <- lm(y~ns(x1,knots=c(-1.5,0,1.5))+ns(x2,knots=c(-4,0,4)))
summary(fit)
plot(y~fit$fitted.values,xlab='Fitted values',main='Natural Cubic Splines')
fit$coefficients
X1 <- as.matrix(ns(x1,knots=c(-1.5,0,1.5)))
b1 <- as.matrix(c(3.050015,-4.075312,-3.435617,-3.156741))
plot((X1%*%b1)[order(x1)]~x1[order(x1)],type='l',ylim=c(-4,2),xlab='x1',ylab='f(x1)')
points(x1,rep(-4,100),pch='|')
points(x1,y)
X2 <- as.matrix(ns(x2,knots=c(-4,0,4)))
b2 <- as.matrix(c(-12.608875,-4.344373,-1.717055,-7.262062))
plot((X2%*%b2)[order(x2)]~x2[order(x2)],type='l',ylim=c(-8,2),xlab='x2',ylab='f(x2)')
points(x2,rep(-8,100),pch='|')
points(x2,y)

# General additive model
library('gam')
fit <- gam(y~s(x1)+s(x2))
summary(fit)
plot(fit)