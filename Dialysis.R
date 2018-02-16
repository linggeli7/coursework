# Exploratory
summary(usrdsData)
library(Hmisc)
na.patterns <- naclus(usrdsData)
naplot(na.patterns)
plot(na.patterns)
A <- is.na(usrdsData$hist.cvd)+is.na(usrdsData$cholest)+is.na(usrdsData$trigly)
sum(A==3)
# albumin is heavy-tailed
boxplot(usrdsData$albumin.0)
boxplot(usrdsData$albumin.0~usrdsData$hist.cvd)
boxplot(usrdsData$albumin.0~usrdsData$undnour)
high <- usrdsData[which(usrdsData$albumin.0>=3.8),]
low <- usrdsData[which(usrdsData$albumin.0<3.8),]
summary(high)
mean(low$esrdtime,na.rm=TRUE)
sd(low$esrdtime,na.rm=TRUE)
sum(high$diabetes,na.rm=TRUE)
sum(high$diabetes,na.rm=TRUE)/1025

# Kaplan-Meier plot
stratified <- usrdsData
stratified$high <- ifelse(stratified$albumin.0>3.8,1,0)
fit <- coxph(Surv(tdeath,death)~strata(high),data=stratified)
plot(survfit(fit), xlab="Time", ylab="Cumulative Survival",
     mark.time=FALSE,lty=2:1)
legend(x=200,y=0.3,c("High serum albumin", "Low serum albumin"), lty = 1:2, 
       bty="n") 
# definitely not proportional
plot(survfit(fit), fun="cloglog", xlab="Time", 
     ylab="Log-Cumulative Hazard Function")
# difficult to tell interaction
card <- stratified[which(stratified$hist.cvd==1),]
fit <- coxph(Surv(tdeath,death)~strata(high),data=card)
plot(survfit(fit), xlab="Time", ylab="Cumulative Survival",
     mark.time=FALSE,lty=2:1)
legend(x=200,y=0.3,c("High serum albumin", "Low serum albumin"), lty = 1:2, 
       bty="n") 
nocard <- stratified[which(stratified$hist.cvd==0),]
fit <- coxph(Surv(tdeath,death)~strata(high),data=nocard)
plot(survfit(fit), xlab="Time", ylab="Cumulative Survival",
     mark.time=FALSE,lty=2:1)
legend(x=200,y=0.3,c("High serum albumin", "Low serum albumin"), lty = 1:2, 
       bty="n") 

# A priori model with complete data
library(survival)
complete <- usrdsData[,-c(14,15)]
complete <- complete[complete.cases(complete),]
fit1 <- coxph(Surv(tdeath,death)~albumin.0+age+female+as.factor(racegrp)
              +as.factor(smokegrp)+hist.cvd+diabetes+bmi+esrdtime
              +pst.sbp, data=complete)

fit2 <- coxph(Surv(tdeath,death)~albumin.0*hist.cvd+age+female+as.factor(racegrp)
              +as.factor(smokegrp)+diabetes+bmi+esrdtime
              +pst.sbp+undnour, data=complete)

anova(fit1,fit2)

# Naive imputation does not work very well
imputed <- usrdsData
impute1 <- glm(hist.cvd~cholest+trigly, family=binomial, data=usrdsData)
missing <- data.frame(cholest=imputed$cholest[which(is.na(imputed$hist.cvd))],
                      trigly=imputed$trigly[which(is.na(imputed$hist.cvd))])
fitted <- predict(impute1,missing)
imputed$hist.cvd[which(is.na(imputed$hist.cvd))] <- ifelse(fitted>0.5,1,0)
# undernourishment 
imputed$undnour[which(is.na(imputed$undnour))] <- 0
# smoking does not seem to matter
imputed$smokegrp[which(is.na(imputed$smokegrp))] <- 1
# BMI
imputed$bmi[which(is.na(imputed$bmi))] <- median(imputed$bmi)
# blood pressure
imputed$pst.sbp[which(is.na(imputed$pst.sbp))] <- median(imputed$pst.sbp)

fit1.imputed <- coxph(Surv(tdeath,death)~albumin.0+age+female+as.factor(racegrp)
              +as.factor(smokegrp)+hist.cvd+diabetes+bmi+esrdtime
              +pst.sbp, data=imputed)

fit2.imputed <- coxph(Surv(tdeath,death)~albumin.0*hist.cvd+age+female+as.factor(racegrp)
              +as.factor(smokegrp)+diabetes+bmi+esrdtime
              +pst.sbp, data=imputed)

anova(fit1.imputed,fit2.imputed)

# diagnostics
sres <- residuals(fit1, type="scaledsch")
time <- sort(complete$tdeath[which(complete$death==1)])
plot(time, sres[,1], xlab="Time", ylab="Scaled Schoenfeld Residual")
lines(smooth.spline(time, sres[,1]),col="red")
# non-proportional hazard
cox.zph(fit1, transform="log")

# martingales
mres <- residuals(fit1, type="martingale")
lmfit <- lm(log(albumin.0)~age+female+as.factor(racegrp)
            +as.factor(smokegrp)+hist.cvd+diabetes+bmi+esrdtime
            +pst.sbp, data=complete)
res <- lmfit$residuals
plot(res, mres, xlab="Linear Regression Residual",
     ylab="Martingale Residual")
lines(res,fitted(lm(mres~res)),col="blue",lwd=2 )
lines(smooth.spline(res,mres,df=6),col="red",lwd=2 )

# log-transformation does not seem to help
fit1.log <- coxph(Surv(tdeath,death)~log(albumin.0)+age+female+as.factor(racegrp)
              +as.factor(smokegrp)+hist.cvd+diabetes+bmi+esrdtime
              +pst.sbp, data=complete)
mres <- residuals(fit1.log, type="martingale")
lmfit <- lm(log(albumin.0)~age+female+as.factor(racegrp)
            +as.factor(smokegrp)+hist.cvd+diabetes+bmi+esrdtime
            +pst.sbp, data=complete)
res <- lmfit$residuals
plot(res, mres, xlab="Linear Regression Residual",
     ylab="Martingale Residual")
lines(res,fitted(lm(mres~res)),col="blue",lwd=2 )
lines(smooth.spline(res,mres,df=6),col="red",lwd=2 )

# influential points
dres <- residuals(fit1, type="deviance")
lp <- predict(fit1, type="lp")
# deviance residuals
plot(lp, dres, xlab="Linear Predictor",
     ylab="Deviance Residual")
abline(a=3,b=0)
# delta beta
dfbeta <- residuals(fit1, type="dfbeta")
plot(seq(1,1526,1), dfbeta[,1], xlab="Observation Number",
     ylab="Delta Beta")

# adding discrete time intervals is brilliant
stage1 <- complete
stage2 <- complete
stage3 <- complete
stage4 <- complete
stage1$stage <- 1
stage1$start <- 0
stage1$stop <- ifelse(complete$tdeath>100, 100, complete$tdeath)
stage1$died <- ifelse(complete$tdeath>100, 0, complete$death)
stage2$stage <- 2
stage2$start <- ifelse(complete$tdeath>100, 100, NA)
stage2$stop <- ifelse(complete$tdeath>200, 200, complete$tdeath)
stage2$died <- ifelse(complete$tdeath>200, 0, complete$death)
stage3$stage <- 3
stage3$start <- ifelse(complete$tdeath>200, 200, NA)
stage3$stop <- ifelse(complete$tdeath>300, 300, complete$tdeath)
stage3$died <- ifelse(complete$tdeath>300, 0, complete$death)
stage4$stage <- 4
stage4$start <- ifelse(complete$tdeath>300, 300, NA)
stage4$stop <- complete$tdeath
stage4$died <- complete$death
staged <- rbind(stage1,stage2,stage3,stage4)
staged <- staged[!is.na(staged$start),]

fit.staged1 <- coxph(Surv(start,stop,died)~albumin.0:as.factor(stage)+age+female+as.factor(racegrp)
                    +as.factor(smokegrp)+diabetes+bmi+esrdtime
                    +pst.sbp+albumin.0*hist.cvd, data=staged)

fit.staged2 <- coxph(Surv(start,stop,died)~age+female+as.factor(racegrp)
                     +as.factor(smokegrp)+diabetes+bmi+esrdtime
                     +pst.sbp+albumin.0+hist.cvd, data=staged)

# create a time-varying covariate for albumin
total <- merge(usrdsData,LongitudinalAlbumin,by="usrds.id")
total$start <- total$measday
total$stop <- NA
total$died <- NA
for (i in 1:11190) {
  if (total$start[i+1]==0) {
    total$stop[i] <- total$tdeath[i]
    total$died[i] <- total$death[i]
  } else {
    total$stop[i] <- total$start[i+1]
    total$died[i] <- 0
  }
}
total$stop[11191] <- 427
total$died[11191] <- 0

# Cox-proportional model
fit1 <- coxph(Surv(start,stop,died)~albumin+age+female+as.factor(racegrp)
              +as.factor(smokegrp)+hist.cvd+diabetes+bmi+esrdtime
              +pst.sbp, data=total)

fit2 <- coxph(Surv(start,stop,died)~albumin*hist.cvd+age+female+as.factor(racegrp)
              +as.factor(smokegrp)+diabetes+bmi+esrdtime
              +pst.sbp+undnour, data=total)

# martingales
mres <- residuals(fit1, type="martingale")
lmfit <- lm(albumin~age+female+as.factor(racegrp)
            +as.factor(smokegrp)+hist.cvd+diabetes+bmi+esrdtime
            +pst.sbp, data=total)
res <- lmfit$residuals
plot(res, mres, xlab="Linear Regression Residual",
     ylab="Martingale Residual")
lines(res,fitted(lm(mres~res)),col="blue",lwd=2 )
lines(smooth.spline(res,mres,df=6),col="red",lwd=2 )

# log-transformation again is not helpful
fit1.log <- coxph(Surv(start,stop,died)~log(albumin)+age+female+as.factor(racegrp)
                  +as.factor(smokegrp)+hist.cvd+diabetes+bmi+esrdtime
                  +pst.sbp, data=total)
mres <- residuals(fit1.log, type="martingale")
lmfit <- lm(log(albumin)~age+female+as.factor(racegrp)
            +as.factor(smokegrp)+hist.cvd+diabetes+bmi+esrdtime
            +pst.sbp, data=total)
res <- lmfit$residuals
plot(res, mres, xlab="Linear Regression Residual",
     ylab="Martingale Residual")
lines(res,fitted(lm(mres~res)),col="blue",lwd=2 )
lines(smooth.spline(res,mres,df=6),col="red",lwd=2 )

#######################################
ifelse1 <- function (test, yes, no){
  if (test) yes
  else no
}

linContr.coxph <- function( model, contr.names, contr.coef=rep(1,length(contr.names)) ){
  beta.hat <- model$coef 
  cov.beta <- vcov( model )
  
  contr.index <- match( contr.names, dimnames(cov.beta)[[1]] ) 
  beta.hat <- beta.hat[ contr.index ]
  cov.beta <- cov.beta[ contr.index,contr.index ]
  est <- contr.coef %*% beta.hat
  se.est <- sqrt( contr.coef %*% cov.beta %*% contr.coef )
  zStat <- est / se.est
  pVal <- 2*pnorm( abs(zStat), lower.tail=FALSE )
  ci95.lo <- exp( est - qnorm(.975)*se.est )
  ci95.hi <- exp( est + qnorm(.975)*se.est )
  est <- exp( est )
  
  cat( "\nTest of H_0: exp( " )
  for( i in 1:(length( contr.names )) ){
    if( i < length( contr.names ) ) cat( contr.coef[i], "*", contr.names[i], " + ", sep="" )
    else cat( contr.coef[i], "*", contr.names[i], " ) = 1 :\n\n", sep="" )	
  }
  
  rslt <- data.frame( est, se.est, zStat, pVal, ci95.lo, ci95.hi )
  colnames( rslt )[1] <- "exp( Est )"
  round( rslt, 3 )
}
