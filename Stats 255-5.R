library(survival)
# Bone marrow transplant
# Patients 127 and 37 
bmt$tcghd[127] <- 20
bmt[37,]$icghd <- 0
bmt[37,]$tcghd <- 110
# Duplicate each patient three times
bmt1 <- bmt
bmt2 <- bmt
bmt3 <- bmt
bmt1$start <- NA
bmt1$stop <- NA
bmt2$start <- NA
bmt2$stop <- NA
bmt3$start <- NA
bmt3$stop <- NA
# Set covariates
for (i in 1:137) {
  if (bmt[i,]$iaghd == 1) {
    if (bmt[i,]$icghd == 0) {
      # only A
      bmt1[i,]$start <- 0
      bmt1[i,]$stop <- bmt[i,]$taghd
      bmt1[i,]$irelpse <- 0
      bmt1[i,]$iaghd <- 0
      bmt2[i,]$start <- bmt[i,]$taghd
      bmt2[i,]$stop <- bmt[i,]$tnodis
    } else {
      # both A and C
      if (bmt[i,]$taghd <= bmt[i,]$tcghd) {
        # A before C
        bmt1[i,]$start <- 0
        bmt1[i,]$stop <- bmt[i,]$taghd
        bmt1[i,]$irelpse <- 0
        bmt1[i,]$iaghd <- 0
        bmt1[i,]$icghd <- 0
        bmt2[i,]$start <- bmt[i,]$taghd
        bmt2[i,]$stop <- bmt[i,]$tcghd
        bmt2[i,]$irelpse <- 0
        bmt2[i,]$icghd <- 0
        bmt3[i,]$start <- bmt[i,]$tcghd
        bmt3[i,]$stop <- bmt[i,]$tnodis
      } else {
        # C before A
        bmt1[i,]$start <- 0
        bmt1[i,]$stop <- bmt[i,]$tcghd
        bmt1[i,]$irelpse <- 0
        bmt1[i,]$iaghd <- 0
        bmt1[i,]$icghd <- 0
        bmt2[i,]$start <- bmt[i,]$tcghd
        bmt2[i,]$stop <- bmt[i,]$taghd
        bmt2[i,]$irelpse <- 0
        bmt2[i,]$iaghd <- 0
        bmt3[i,]$start <- bmt[i,]$taghd
        bmt3[i,]$stop <- bmt[i,]$tnodis
      }
    }
  } else if (bmt[i,]$icghd == 1) {
    # only C
    bmt1[i,]$start <- 0
    bmt1[i,]$stop <- bmt[i,]$tcghd
    bmt1[i,]$irelpse <- 0
    bmt1[i,]$icghd <- 0
    bmt2[i,]$start <- bmt[i,]$tcghd
    bmt2[i,]$stop <- bmt[i,]$tnodis
  } else {
    bmt1[i,]$start <- 0
    bmt1[i,]$stop <- bmt[i,]$tnodis
  }
}
# Combine data
bmt.ac <- rbind(bmt1, bmt2, bmt3)
bmt.ac <- bmt.ac[!is.na(bmt.ac$start),]
# Fit model
Surv(bmt.ac$start, bmt.ac$stop, bmt.ac$irelpse)
fit <- coxph(Surv(start, stop, irelpse)~iaghd+icghd, data=bmt.ac)
fit1 <- coxph(Surv(start, stop, irelpse)~iaghd*icghd, data=bmt.ac)

# Tongue
fit <- coxph(Surv(timewks,idied)~dna,data=tongue)
# Schoenfeld residuals
sres <- residuals(fit, type="scaledsch")
time <- sort(tongue$timewks[which(tongue$idied==1)])
plot(time, sres, xlab="Time", ylab="Scaled Schoenfeld Residual")
cox.zph(fit, transform="log")

# Create time-dependent covariate
tongue.t <- tongue[order(tongue$timewks),]
times <- unique(tongue.t$timewks[tongue.t$idied==1])
n <- length(times)
tongue.t <- tongue[rep(tongue$X,each=n),]
tongue.t$start <- rep(c(0,times[1:n-1]),80)
tongue.t$stop <- rep(times,80)
tongue.t <- tongue.t[tongue.t$timewks>tongue.t$start,]
tongue.t$idied <- ifelse(tongue.t$stop==tongue.t$timewks,1,0)

# Interaction with log(time)
fit2 <- coxph(Surv(start,stop,idied)~dna*log(stop)-log(stop), data=tongue.t)

# Stratified log-log
fit1 <- coxph(Surv(timewks,idied)~strata(dna),data=tongue)
plot(survfit(fit1), fun="cloglog", xlab="Time", 
     ylab="Log-Cumulative Hazard Function")

# Deviance residuals
dres <- residuals(fit, type="deviance")
plot(seq(1,80,1), dres, xlab="Observation Number",
     ylab="Deviance Residual")

# Delta betas
dfbeta <- residuals(fit, type="dfbeta")
plot(seq(1,80,1), dfbeta, xlab="Observation Number",
     ylab="Delta Beta")

# Burn
fit <- coxph(Surv(tinfect, iinfect)~ibodycln+as.factor(type)+area, data=burn)
mres <- residuals(fit, type="martingale")
lmfit <- lm(area~ibodycln+as.factor(type), data=burn)
res <- lmfit$residuals
ord <- order(burn$area)
mres <- mres[ord]
res <- res[ord]
plot(res, mres, xlab="Linear Regression Residual",
     ylab="Martingale Residual")
lines(res,fitted(lm(mres~res)),col="blue",lwd=2 )
lines(smooth.spline(res,mres,df=6),col="red",lwd=2 )
text(res[which(ord == 7)],mres[which(ord == 7)]+0.1,labels="7")
text(res[which(ord == 58)]-5,mres[which(ord == 58)],labels="58")

# Remove subjects 7 and 58 and take log(area)
burn1 <- burn[-c(7,58),]
fit <- coxph(Surv(tinfect, iinfect)~ibodycln+as.factor(type)+log(area), data=burn)
mres <- residuals(fit, type="martingale")
lmfit <- lm(log(area)~ibodycln+as.factor(type), data=burn1)
res <- lmfit$residuals
ord <- order(burn1$area)
mres <- mres[ord]
res <- res[ord]
plot(res, mres, xlab="Linear Regression Residual",
     ylab="Martingale Residual")
lines(res,fitted(lm(mres~res)),col="blue",lwd=2 )
lines(smooth.spline(res,mres,df=6),col="red",lwd=2 )

# Summary
fit <- coxph(Surv(tinfect, iinfect)~ibodycln+as.factor(type)+log(area), data=burn)

# Find cut-off value
sort(burn1$area)[sum(res<0)]
burn$la <- ifelse(burn$area>20,1,0)
fit <- coxph(Surv(tinfect, iinfect)~ibodycln+as.factor(type)+la, data=burn)

# Deviance residuals
dres <- residuals(fit, type="deviance")
lp <- predict(fit, type="lp")
plot(lp, dres, xlab="Linear Predictor",
     ylab="Deviance Residual")
abline(a=2.5,b=0)
burn[dres>2.2,]

# Delta betas
dfbeta <- residuals(fit, type="dfbeta")
plot(seq(1,154,1), dfbeta[,5], xlab="Observation Number",
     ylab="Delta Beta")