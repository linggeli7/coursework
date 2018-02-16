##
#####
#####  Stat 211 -- Problem Set 2, Problem 4
#####
#####
##
#####
#####	getData() simulates nsims datasets with nobs observations per dataset
#####	with one of five error distributions
#####
##
getData <- function( nobs=10, nsims=1000, distn="normal", 
                     location=0, scale=1, minimum=-1, maximum=1 ){
  if ( !is.element( distn, c( "normal", "cauchy", "uniform", "dexp", "sexp" ) ) ) 
    stop( "Enter a valid distribution" )	
  beta0 <- 1
  beta1 <- 2
  x <- rnorm( n = nobs, mean = 10, sd = 3 )
  if ( distn == "normal" ) epsilon <- rnorm( nobs * nsims, location, scale )
  else if ( distn == "cauchy" ) epsilon <- rcauchy( nobs * nsims, location, scale )
  else if ( distn == "uniform" ) epsilon <- runif( nobs * nsims, minimum, maximum )
  else if ( distn == "dexp" ) epsilon <- rexp( nobs * nsims, 1/scale ) - rexp( nobs * nsims, 1/scale )
  else if ( distn == "sexp" ) epsilon <- rexp( nobs * nsims, 1/scale ) - scale 
  y <- matrix( beta0 + beta1 * x + epsilon, nrow = nobs, ncol=nsims )
  return( list( y=y, x=x ) )
}
##
##	Example call to simData : generate 1000 datasets with n=10 obs per data set
##	with errors from a Normal distribution with mean (location) 0 and scale 1
##
normData <- getData( nobs=10, nsims=1000, distn="normal", location=0, scale=1 )
normData$y[,1:5]
normData$x

##
#####
#####	OLSestimates() calculates the OLS estimates for each simulated dataset
#####
##
OLSestimates <- function( simData ){
  x <- simData$x
  y <- simData$y
  nobs <- length(x)
  x.bar <- mean(x)
  y.bar <- as.vector(rep(1, nobs) %*% y) / nobs
  s.xx <- sum(x^2) - nobs*(x.bar^2)
  s.yy <- as.vector(rep(1, nobs) %*% y^2) - nobs*(y.bar^2)
  s.xy <-  as.vector(rep(1, nobs) %*% (x*y)) - nobs*x.bar*y.bar
  b1.hat <- s.xy / s.xx
  b0.hat <- y.bar - b1.hat*x.bar
  sigma2.hat <- (1/(nobs-2)) * (s.yy - (s.xy^2/s.xx))
  sigma.hat <- sqrt(sigma2.hat)
  se.b1 <- sigma.hat * sqrt(1/s.xx)
  se.b0 <- sigma.hat * sqrt((s.xx + nobs*x.bar^2) / (nobs*s.xx))
  return( as.data.frame( cbind( b0.hat, b1.hat, se.b0, se.b1 ) ) )
}  
##
##	Example call to OLSestimates : Calculate the OLS estimates based upon 1000 
##	simulated datasets with n=10 obs per data set and errors from a Normal 
##  distribution with mean (location) 0 and scale 1
##
normOLS <- OLSestimates( normData )
normOLS[1:5,]
summary( lm( normData$y[,1] ~ normData$x ) )$coefficients

##
#####
#####	getCoverage() computes 95% CIs for each parameter and estimates
#####	the coverage probability of each CI across the simulations
#####
##
getCoverage <- function( estimates, nobs ){
  nsims <- dim( estimates )[1]
  b0Low <- estimates$b0.hat - ( qt( 0.975, nobs-2 ) * estimates$se.b0 )
  b0Up <- estimates$b0.hat + ( qt( 0.975, nobs-2 ) * estimates$se.b0 )
  b1Low <- estimates$b1.hat - ( qt( 0.975, nobs-2 ) * estimates$se.b1 )
  b1Up <- estimates$b1.hat + ( qt( 0.975, nobs-2 ) * estimates$se.b1 )
  b0Cover <- 1 - ( sum( 1 < b0Low ) + sum( 1 > b0Up ) ) / nsims
  b1Cover <- 1 - ( sum( 2 < b1Low ) + sum( 2 > b1Up ) ) / nsims
  return( c( b0Cover, b1Cover ) )
}
##
##	Example call to getCoverage : Calculate the coverage probability forthe previous 
##	simulated datasets with n=10 obs per data set and errors from a Normal 
##  distribution with mean (location) 0 and scale 1
##
getCoverage( estimates=normOLS, nobs=10 )

##
#####
##### 	getSumm() calculates the mean of the OLS estimates and variances and
#####	reports these along with the resulting coverage probabilities 
#####
##
getSumm <- function( estimates, nobs ){
  nsims <- dim( estimates )[1]
  mean.b0 <- mean( estimates$b0.hat )
  se.b0 <- mean( estimates$se.b0 )
  mean.b1 <- mean( estimates$b1.hat )
  se.b1 <- mean( estimates$se.b1 )
  coverage <- getCoverage( estimates, nobs = nobs )
  rslt <- round( c( mean.b0, se.b0, coverage[1], mean.b1, se.b1, coverage[2] ), 3 )
  names( rslt ) <- c( "mean_b0", "se_b0", "cov_b0", "mean_b1", "se_b1", "cov_b1" ) 
  return( rslt )
}
getSumm( normOLS, nobs=10 )

options(object.size = 9e6)
##
#####
#####  Case 1:  Normal errors with mean 0, n = 10, 100, 1000, scale = 1, 3, 5 #####
#####
##
nobs <- c( 10, 100, 1000 )
scale <- c( 1, 3, 5 )
summStats <- NULL

for ( i in scale){
  par( mfrow=c(3,4) )
  for ( j in nobs ){
    simData <- getData( nobs=j, nsims=1000, distn="sexp", location=0, scale=i)
    estimates <- OLSestimates( simData )
    summStats <- rbind( summStats, getSumm( estimates, nobs=j ) )
    #####  Graphical representations to investigate normality  #####
    qqnorm( estimates$b0.hat, ylab=expression( hat(beta)[0] ), 
            main=expression( paste( "Q-Q Plot : ", hat(beta)[0] ) ) )
    qqline( estimates$b0.hat )
    hist( 	estimates$b0.hat, xlab=expression( hat(beta)[0] ), 
           main=expression( paste( "Histogram of ", hat(beta)[0] ) ) )
    qqnorm( estimates$b1.hat, ylab=expression( hat(beta)[1] ), 
            main=expression( paste( "Q-Q Plot : ", hat(beta)[1] ) ) )
    qqline( estimates$b1.hat )
    hist( 	estimates$b1.hat, xlab=expression( hat(beta)[1] ), 
           main=expression( paste( "Histogram of ", hat(beta)[1] ) ) )			
  }
}
dimnames( summStats )[[1]] <- paste( "n=", rep(nobs, 3 ), ", scale=", rep(scale, each=3), " : ", sep="" )
summStats
