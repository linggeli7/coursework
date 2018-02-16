library(survival)
set.seed(12345)
true.time <- rexp(200, rate=.5)
cens.time <- runif(200, 1, 3)
obs.time <- ifelse( true.time <= cens.time, true.time, cens.time)
delta <- ifelse( obs.time==true.time, 1, 0)
group <- rep(0:1, each=100)


GstatFH <- function (time, event, group, rho = 1, gamma = 1) {
  fit <- survfit(Surv(time, event) ~ 1)
  fit.group <- survfit(Surv(time, event) ~ group)
  n <- length(time)
  n1 <- sum(group)
  n0 <- n - n1
  time0 <- fit.group$time[1:n0]
  time1 <- fit.group$time[(n0+1):n]
  ind0 <- rep(0,n)
  ind1 <- rep(0,n)
  n0.event <- rep(0,n)
  n1.event <- rep(0,n)
  n0.risk <- rep(0,n)
  n1.risk <- rep(0,n)
  n0.risk[1] <- n0
  n1.risk[1] <- n - n0
  for (i in 1:n) {
    pos <- match(fit$time[i], time0)
    if (!is.na(pos)) {
      ind0[i] <- 1
      n0.event[i] <- fit$n.event[i]
    } 
    pos <- match(fit$time[i], time1)
    if (!is.na(pos)) {
      ind1[i] <- 1
      n1.event[i] <- fit$n.event[i]
    }
    if (i > 1) {
      n0.risk[i] <- n0.risk[1] - sum(ind0[1:i-1])
      n1.risk[i] <- n1.risk[1] - sum(ind1[1:i-1])
    }
  }
  Gstat <- 0
  Var <- 0
  for (i in 1:n) {
    if (n0.risk[i] > 0 && n1.risk[i] > 0) {
      if (i == 1) {
        S <- 1
      } else {
        S <- fit$surv[i-1]
      }
      w <- (S^rho)*((1-S)^gamma)
      U <- n1.event[i] - n1.risk[i]*fit$n.event[i]/fit$n.risk[i]
      Gstat <- Gstat + w*U
      V <- n1.risk[i]*n0.risk[i]*fit$n.event[i]*(fit$n.risk[i]-fit$n.event[i])
      V <- V/(fit$n.risk[i]^2*(fit$n.risk[i]-1))
      Var <- Var + w^2*V 
    }
  }
  Zstat <- Gstat/(Var^0.5)
  Pval <- (1-pnorm(abs(Zstat)))*2
  rslt <- c(Gstat, Var, Zstat, Pval)
  names( rslt ) <- c("Gstat", "Var", "Zstat", "Pval")
  return(rslt)
}

group <- ifelse( tongue$dna==1, 0, 1)
GstatFH( time=tongue$timewks, event=tongue$idied, group=group, rho=1, gamma=1 )

ifelse1 <- function (test, yes, no){
  if (test) yes
  else no
}

qKM <- function( fit, p, xscale=1 ){
  fit$time <- fit$time / xscale
  if( is.null( fit$strata ) ) fit$strata <- length(fit$surv)
  nStrata <- length( fit$strata )
  rslt <- vector( "list", length=nStrata )
  names( rslt ) <- ifelse1( nStrata==1, "Survival Quantiles", names( fit$strata ) )
  for( i in 1:nStrata ){
    rslt[[i]] <- as.data.frame( cbind( p, matrix(NA, nrow=length(p), ncol=3) ) )
    names( rslt[[i]] ) <- c( "perc.surv", "km.est", "lower", "upper" ) 
    strata.fit <-  as.data.frame( cbind( fit$time, fit$surv, fit$lower, fit$upper )[ (c(0,cumsum(fit$strata))[i]+1):(c(0,cumsum(fit$strata))[i]+fit$strata[i]), ] )
    names( strata.fit ) <- c("time", "surv", "lower", "upper" )
    for(j in 1:length(p)){
      if( sum(strata.fit$surv<=p[j]) ){
        rslt[[i]][j,2] <- strata.fit$time[ min( which(strata.fit$surv <=p[j] ) ) ]
        rslt[[i]][j,3] <- ifelse( sum(strata.fit$lower>p[j], na.rm=TRUE) >= 1, strata.fit$time[ max(which(strata.fit$lower>p[j]))+1 ], 0 )
        rslt[[i]][j,4] <- ifelse( sum(strata.fit$upper<p[j], na.rm=TRUE) >= 1, strata.fit$time[ min(which(strata.fit$upper<p[j])) ], Inf )
      }
    }
  }
  return(rslt)
}	

obs.time <- bmt$tnodis
delta <- bmt$inodis
fit=survfit(Surv(obs.time,delta) ~ 1)
plot(fit, xlab = "Time in days", ylab="Disease free survival rate")
summary(fit)
qKM(fit = fit, p = 0.75)

obs.time1 <- bmt$tnodis[which(bmt$tnodis >= 90)]
delta1 <- bmt$inodis[which(bmt$tnodis >= 90)]
fit1=survfit(Surv(obs.time1, delta1) ~ 1)
plot(fit1)
summary(fit1)

fit$time
fit$n.risk[19]
bmt$g[which(bmt$tnodis >= 93)]

fit2=survfit(Surv(tongue$timewks, tongue$idied) ~ tongue$dna)
plot(fit2, lty=1:2, xlab = "Time", ylab = "Survival rate")
legend( 280, .9, lty=1:2, legend=c("DNA 1", "DNA 2") )
summary(fit2)
plot(fit2, fun="cloglog",lty=1:2, xlab = "Time", ylab = "Log cumulative hazard")
legend( 1.5, 0.2, lty=1:2, legend=c("DNA 1", "DNA 2") )
qKM(fit = fit2, p = 0.5)
survdiff(Surv(tongue$timewks, tongue$idied) ~ tongue$dna)

survdiff(Surv(obstime,death)~race+strata(gender), data=kidneytransplant)

survdiff(Surv(obstime[1:524],death[1:524])~race[1:524], data=kidneytransplant)
survdiff(Surv(obstime[525:863],death[525:863])~race[525:863], data=kidneytransplant)

