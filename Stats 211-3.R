# summary table
summary(salary95)
# salary difference between males and females
plot(salary ~ gender, data = salary95)
# salary associated with rank, field and admin
plot(salary ~ rank, data = salary95)
plot(salary ~ field, data = salary95)
plot(salary ~ as.factor(admin), data = salary95)
# rank associated with gender
plot(rank ~ gender, data = salary95)
# field and admin not so much (precision variables)
plot(as.factor(admin) ~ gender, data = salary95)
plot(field ~ gender, data = salary95)
# startyr and deg 
plot(salary ~ startyr, data = salary95)
plot(yrdeg ~ startyr, data = salary95)
plot(salary ~ deg, data = salary95)
# small difference with rank
model1 <- lm(salary ~ gender + field + startyr + deg + admin + rank, data=salary95)
summary(model1)
# big difference without rank
model2 <- lm(salary ~ gender + field + startyr + deg + admin, data=salary95)
summary(model2)
# not constant variance
plot(model2$residuals ~ model2$fitted.values)
# hella outliers
library(car)
influencePlot(model2,  id.method="noteworthy", main="Influence Plot")
# log transformation does not work very well
boxplot(salary95$salary)
model2 <- lm(log(salary) ~ gender + field + startyr + deg + admin, data=salary95)
summary(model2)
plot(model2$residuals ~ model2$fitted.values)
704.471+1.961458*101.595
qt(0.975,1589)

##
##### 	Compute robust (sandwich) variance-covariance estimate for a LM
##
robust.vcov.lm <- function( lm.obj ){
  X <- model.matrix( lm.obj )
  eps <- lm.obj$residuals
  robust.cov <- solve( t(X)%*%X ) %*%( t(X) %*% diag(eps^2) %*% X ) %*% solve( t(X)%*%X )
  dimnames( robust.cov ) <- dimnames( vcov(lm.obj) )
  return( robust.cov )
}

##
#####
#####	robust.se.lm() is a function to compute the Huber-White sandwich variance estimator
#####	for the linear regression model
#####	
##
robust.se.lm <- function( model) { 
  s <- summary( model) 
  X <- model.matrix( model )
  sandwich.cov <- robust.vcov.lm( model )
  sand.se <- sqrt( diag( sandwich.cov )) 
  t <- model$coefficients/sand.se
  p <- 2*pt( -abs( t ), dim(X)[1]-dim(X)[2] ) 
  ci95.lo <- model$coefficients - qt( .975, dim(X)[1]-dim(X)[2] ) * sand.se
  ci95.hi <- model$coefficients + qt( .975, dim(X)[1]-dim(X)[2] ) * sand.se
  rslt <- cbind( model$coefficients, sand.se, ci95.lo, ci95.hi, t, p ) 
  dimnames(rslt)[[2]] <- c( dimnames( s$coefficients )[[2]][1], "Robust SE", "ci95.lo", "ci95.hi", dimnames( s$coefficients )[[2]][3:4] ) 
  rslt 
} 

A <- matrix( nrow=1000, ncol=6) 
n <- 500
for( i in 1:1000 ){
  x <- rnorm( n, mean = 1, sd = 1 )
  mu <- 1 + 2*x
  # e <- rnorm( n, mean = 0, sd = 1 )
  e <- rnorm( n, mean = 0, sd = abs(mu) )
  y <- mu + e
  lm.fit <- lm(y ~ x)
  b0 <- lm.fit$coefficients[1]
  A[i,1] <- b0
  b1 <- lm.fit$coefficients[2]
  A[i,2] <- b1
  b0.se <- coef(summary(lm.fit))[, "Std. Error"][1]
  A[i,3] <- b0.se
  b1.se <- coef(summary(lm.fit))[, "Std. Error"][2]
  A[i,4] <- b1.se
  b0.robust_se <- robust.se.lm(lm.fit)[1,2]
  A[i,5] <- b0.robust_se
  b1.robust_se <- robust.se.lm(lm.fit)[2,2]
  A[i,6] <- b1.robust_se
}
c(mean(A[,1]),sd(A[,1]),mean(A[,3]),mean(A[,5]))
c(mean(A[,2]),sd(A[,2]),mean(A[,4]),mean(A[,6]))