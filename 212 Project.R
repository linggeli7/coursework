################### Data cleaning ##########################
# Extract useful columns
mix.data.full <- mix.data.full[which(mix.data.full$Period=="ANN"),]
clean <- mix.data.full[,c('MFI.ID','MFI.name','Fiscal.Year','Diamonds',
                          'Gross.Loan.Portfolio','Country','Profit.status','Region',
                          'Scale','Operational.self.sufficiency')]
#clean$Profit.margin <- as.numeric(sub("%","",mix.data.full$Profit.margin))/100
#clean$Assets <- as.numeric(gsub(",","",mix.data.full$Assets))
#clean$Equity <- as.numeric(gsub(",","",mix.data.full$Equity))
#clean$Profit.margin <- as.numeric(sub("%","",mix.data.full$Profit.margin))/100
clean$Gross.Loan.Portfolio <- as.numeric(gsub(",","",mix.data.full$Gross.Loan.Portfolio))
clean$Yield.on.gross.portfolio..real. <- as.numeric(sub("%","",mix.data.full$Yield.on.gross.portfolio..real.))/100
clean$Operational.self.sufficiency <- as.numeric(sub("%","",clean$Operational.self.sufficiency))/100
clean$Portfolio.at.risk..gt..90.days <- as.numeric(sub("%","",clean$Portfolio.at.risk..gt..90.days))/100
clean <- clean[which(mix.data.full$Diamonds>3),]
clean <- clean[which(clean$Fiscal.Year>2008),]
clean <- clean[which(clean$Fiscal.Year<2015),]
clean$duration <- NA
for (i in 1:3687) {
  clean$duration[i] <- length(which(clean$MFI.ID==clean$MFI.ID[i]))
}
clean <- clean[which(clean$duration==6),]

clean <- clean[which(clean$Gross.Loan.Portfolio>0),]
clean$time <- clean$Fiscal.Year-2008
clean <- clean[which(clean$Profit.status=='Non-profit'|clean$Profit.status=='Profit'),]
clean <- clean[which(clean$Scale=='Small'|clean$Scale=='Medium'|clean$Scale=='Large'),]

length(unique(clean$MFI.ID))

####################### Exploratory ###########################
summary(Profit.status)
summary(Region)
summary(Scale)
summary(clean$Yield.on.gross.portfolio..real.)
summary(clean$Operational.self.sufficiency)
summary(clean[which(clean$Region=='South Asia'),])

# Make mad plots
library(ggplot2)

ggplot(data=clean, aes(x=as.factor(Fiscal.Year), y=Gross.Loan.Portfolio)) +
  geom_boxplot()

# Log is good
ggplot(data=clean, aes(x=as.factor(Fiscal.Year), y=Gross.Loan.Portfolio)) +
  geom_boxplot()

summary(clean$Region)

# Africa
a <- which(clean$Region=='Africa')

ggplot(data=clean[a,],aes(x=as.factor(Fiscal.Year),y=log10(Gross.Loan.Portfolio))) +
  geom_boxplot(fill="#E69F00")

ggplot(data=clean[a,],aes(x=time,y=log10(Gross.Loan.Portfolio),group=MFI.ID))+
  geom_line(colour="#E69F00")+geom_point(colour="#E69F00")

# East Asia and the Pacific
p <- which(clean$Region=='East Asia and the Pacific')

ggplot(data=clean[p,],aes(x=as.factor(Fiscal.Year),y=log10(Gross.Loan.Portfolio))) +
  geom_boxplot(fill="#CC79A7")

ggplot(data=clean[p,],aes(x=time,y=log10(Gross.Loan.Portfolio),group=MFI.ID,colour=Profit.status))+
  geom_line(colour="#CC79A7")+geom_point(colour="#CC79A7")

# Eastern Europe and Central Asia
e <- which(clean$Region=='Eastern Europe and Central Asia')

ggplot(data=clean[e,],aes(x=as.factor(Fiscal.Year),y=log10(Gross.Loan.Portfolio))) +
  geom_boxplot(fill="#009E73")

ggplot(data=clean[e,],aes(x=time,y=log10(Gross.Loan.Portfolio),group=MFI.ID,colour=Profit.status))+
  geom_line(colour="#009E73")+geom_point(colour="#009E73")

# Latin America and The Caribbean 
l <- which(clean$Region=='Latin America and The Caribbean')

ggplot(data=clean[l,],aes(x=as.factor(Fiscal.Year),y=log10(Gross.Loan.Portfolio))) +
  geom_boxplot(fill="#D55E00")

ggplot(data=clean[l,],aes(x=time,y=log10(Gross.Loan.Portfolio),group=MFI.ID,colour=Profit.status))+
  geom_line(colour="#D55E00")+geom_point(colour="#D55E00")

# Middle East and North Africa
m <- which(clean$Region=='Middle East and North Africa')

ggplot(data=clean[m,], aes(x=as.factor(Fiscal.Year), y=log10(Gross.Loan.Portfolio))) +
  geom_boxplot(fill='purple')

ggplot(data=clean[m,],aes(x=time,y=log10(Gross.Loan.Portfolio),group=MFI.ID,color=Profit.status))+
  geom_line(colour='purple')+geom_point(colour='purple')

# South Asia
s <- which(Region=='South Asia')

ggplot(data=clean[s,],aes(x=as.factor(Fiscal.Year),y=log10(Gross.Loan.Portfolio))) +
  geom_boxplot(fill="#0072B2")

ggplot(data=clean[s,],aes(x=time,y=log10(Gross.Loan.Portfolio),group=MFI.ID,colour=Profit.status))+
  geom_line(colour="#0072B2")+geom_point(colour="#0072B2")

# Covariates
ggplot(data=clean,aes(x=time,y=log10(Gross.Loan.Portfolio),group=MFI.ID,colour=sufficiency))+
  geom_line(alpha=0.5)+geom_point(alpha=0.5)

clean$sufficiency <- ifelse(clean$Operational.self.sufficiency>1.2,'High','Low')

ggplot(data=clean,aes(x=as.factor(Fiscal.Year),y=log10(Gross.Loan.Portfolio),fill=sufficiency)) +
  geom_boxplot()

ggplot(data=clean,aes(x=as.factor(Fiscal.Year),y=log10(Gross.Loan.Portfolio),fill=Scale)) +
  geom_boxplot()

ggplot(data=clean,aes(x=as.factor(Fiscal.Year),y=log10(Gross.Loan.Portfolio),fill=Profit.status)) +
  geom_boxplot()

# All together
ggplot(data=clean,aes(x=as.factor(Fiscal.Year),y=log10(Gross.Loan.Portfolio),fill=Region)) +
  geom_boxplot()

ggplot(data=clean,aes(x=time,y=log10(Gross.Loan.Portfolio),group=MFI.ID,colour=Region))+
  geom_line(alpha=0.5)+geom_point()

######################## lme4 ######################################
library(lme4)

# Linear regression model
fixed <- lm(log10(Gross.Loan.Portfolio)~Fiscal.Year*Region-Fiscal.Year-Region,data=clean)
fit1 <- lmer(log10(Gross.Loan.Portfolio)~time*Region-Region
            +Operational.self.sufficiency+Scale*time-time
            +Profit.status*time-time+(time|MFI.ID),data=clean)
fit2 <- lmer(log10(Gross.Loan.Portfolio)~time
             +Operational.self.sufficiency+Scale*time-time
             +Profit.status*time-time+(time|MFI.name)-1,data=clean)
plot(resid(fit1)~fitted(fit1))
qqnorm(resid(fit1),ylim=c(-3,3))
summary(fit1)
summary(log10(clean$Gross.Loan.Portfolio))
qqnorm(log10(clean$Gross.Loan.Portfolio)-7.2,ylim=c(-3,3))

V <- vcov(fit1)[c('time:RegionEast Asia and the Pacific','time:RegionSouth Asia'),
           c('time:RegionEast Asia and the Pacific','time:RegionSouth Asia')]
A <- t(as.matrix(c(1,-1)))
sd <- sqrt(A%*%V%*%t(A))

########################### nlme ######################################
library(nlme)

clean <- clean[which(!is.na(clean$Operational.self.sufficiency)),]

fit1 <- lme(fixed=log10(Gross.Loan.Portfolio)~time*Region-Region
            +Operational.self.sufficiency+Scale*time-time
            +Profit.status*time-time,
            random=~time|MFI.ID,data=clean)
summary(fit1)

fit2 <- lm(log10(Gross.Loan.Portfolio)~time*Region-Region
            +Operational.self.sufficiency+Scale*time-time
            +Profit.status*time-time,data=clean)
summary(fit2)

# AR1 is not appropriate
fit2 <- lme(fixed=log10(Gross.Loan.Portfolio)~time*(Region+Scale+Profit.status)-time-Region+
              Operational.self.sufficiency,random=~time|MFI.ID,
            data=clean,correlation=corAR1())
summary(fit2)
qqnorm(resid(fit2,resType='normalized'),ylim=c(-3,3))
plot(ACF(fit2,form=~time,resType='pearson'),alpha=0.05)
plot(ACF(fit3,form=~time,resType='normalized'),alpha=0.01)
plot(fit2)
qqnorm(fit1,grid=TRUE)
qqnorm(fit3,grid=TRUE)
plot(fit2,idResType='normalized')
plot(fit3,idResType='normalized')

# Plot residuals vs fitted values
clean$lme.fit <- fit1$fitted[,2]
clean$lme.resid <- resid(fit1,resType='normalized')
ggplot(data=clean,aes(x=lme.fit,y=lme.resid))+geom_point()

fit3 <- lme(fixed=log10(Gross.Loan.Portfolio)~time,random=~time|MFI.name,
            data=clean,na.action=na.exclude,correlation=NULL)

R <- fit3$residuals

### Calculate autocorrelation manually
MFI <- unique(clean$MFI.ID)
lag <- 5
ACF <- rep(NA,500)
for (i in 1:415) {
  index <- which(clean$MFI.ID==MFI[i])
  if (length(index)>lag) {
    numerator <- 0
    for (j in 1:(length(index)-lag)) {
      numerator <- numerator + R[index[j],2]*R[index[j+lag],2]
    }
    acf <- numerator/(length(index)-lag)/mean((R[index,2])^2)
    ACF[i] <- acf
  }
}
ACF <- ACF[which(!is.na(ACF))]
mean(ACF)

numerator/n1/mean((R[,2])^2)

######################### Two stage modeling ########################
Res.lme <- matrix(nrow=415,ncol=6)
for (i in 1:415) {
  index <- which(clean$MFI.ID==MFI[i])
  time <- clean$time[index]
  if (length(index)==6) {
    Res.lme[i,time] <- fit1$residuals[index,2]
  } 
}

Res2 <- matrix(nrow=415,ncol=6)
for (i in 1:10) {
  index <- which(clean$MFI.ID==MFI[i])
  time <- clean$time[index]
  loan <- log10(clean$Gross.Loan.Portfolio[index])
  fit <- lm(loan~time)
  plot(fit)
  Res2[i,time] <- fit$residuals 
}

Res.diff <- Res.lme-Res2
mean(Res.diff,na.rm=TRUE)
mean(Res.lme,na.rm=TRUE)
mean(Res2,na.rm=TRUE)
mean(abs(Res2),na.rm=TRUE)

# Stage I autocorrelation
Res <- na.omit(Res.lme)
ACF2 <- rep(0,6)
for (i in 1:347) {
  ACF2 <- ACF2+acf(Res[i,],plot=FALSE)$acf
}
auto <- data.frame(lag=0:5,acf=ACF2/347)
ggplot(data=auto, aes(x=lag, y=acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(aes(xend = lag, yend = 0))

# Fitted lines
loan.fitted <- rep(NA,2409)
for (i in 1:415) {
  index <- which(clean$MFI.ID==MFI[i])
  time <- clean$time[index]
  loan <- log10(clean$Gross.Loan.Portfolio[index])
  fit <- lm(loan~time)
  loan.fitted[index] <- fit$fitted.values
}
clean$fitted <- loan.fitted

# Coefficients 
MFI <- unique(clean$MFI.ID)
coef1.fitted <- rep(NA,400)
coef2.fitted <- rep(NA,400)
coef.region <- rep(NA,400)
for (i in 1:376) {
  index <- which(clean$MFI.ID==MFI[i])
  time <- clean$time[index]
  loan <- log10(clean$Gross.Loan.Portfolio[index])
  fit <- lm(loan~time)
  coef1.fitted[i] <- fit$coefficients[1]
  coef2.fitted[i] <- fit$coefficients[2]
  coef.region[i] <- clean$Region[index[1]]
}
mean(coef2.fitted[which(coef.region==1)])
mean(coef2.fitted[which(coef.region==2)])
mean(coef2.fitted[which(coef.region==3)])
mean(coef2.fitted[which(coef.region==4)])
mean(coef2.fitted[which(coef.region==5)])
mean(coef2.fitted[which(coef.region==6)])

lines <- data.frame(intercept=coef1.fitted,slope=coef2.fitted,region=as.factor(coef.region))
ggplot(data=lines,aes(x=intercept,y=slope,colour=region))+geom_point()

# Africa
ggplot(data=clean[a,],aes(x=time,y=fitted,group=MFI.ID))+
  geom_line(colour="#E69F00")+
  geom_point(aes(x=time,y=log10(Gross.Loan.Portfolio),group=MFI.ID),colour="#E69F00")

# East Asia and the Pacific
ggplot(data=clean[p,],aes(x=time,y=fitted,group=MFI.ID))+
  geom_line(colour="#CC79A7")+
  geom_point(aes(x=time,y=log10(Gross.Loan.Portfolio),group=MFI.ID),colour="#CC79A7")

# Eastern Europe and Central Asia
ggplot(data=clean[e,],aes(x=time,y=fitted,group=MFI.ID))+
  geom_line(colour="#009E73")+
  geom_point(aes(x=time,y=log10(Gross.Loan.Portfolio),group=MFI.ID),colour="#009E73")

# Latin America and the Caribbean
ggplot(data=clean[l[1:300],],aes(x=time,y=fitted,group=MFI.ID))+
  geom_line(colour="#D55E00",alpha=0.5)+
  geom_point(aes(x=time,y=log10(Gross.Loan.Portfolio),group=MFI.ID),colour="#D55E00")

# Middle East and North Africa
ggplot(data=clean[m,],aes(x=time,y=fitted,group=MFI.ID))+
  geom_line(colour='purple')+
  geom_point(aes(x=time,y=log10(Gross.Loan.Portfolio),group=MFI.ID),colour='purple')

# South Asia
ggplot(data=clean[s[1:205],],aes(x=time,y=fitted,group=MFI.ID))+
  geom_line(colour="#0072B2")+
  geom_point(aes(x=time,y=log10(Gross.Loan.Portfolio),group=MFI.ID),colour="#0072B2")


# Comparison
p1 <- which(clean$Region=='East Asia and the Pacific'|clean$Region=='South Asia')
ggplot(data=clean[p1,],aes(x=time,y=fitted,group=MFI.ID,colour=Region))+
  geom_line(alpha=0.8)+geom_point()

# Overall trends
summary(clean$Region)
region <- c(rep('Africa',6),rep('East Asia and the Pacific',6),
            rep('Eastern Europe and Central Asia',6),rep('Latin America and The Caribbean',6),
            rep('Middle East and North Africa',6),rep('South Asia',6))
time <- rep(1:6,6)
trend <- rep(NA,36)
for (i in 1:36) {
  trend[i] <- mean(clean$fitted[which(clean$time==time[i]&clean$Region==region[i])],na.rm=TRUE)
}
overall <- data.frame(time,trend,region)
ggplot(data=overall,aes(x=time,y=trend,group=region,colour=region))+
  geom_line()+geom_point()+theme(legend.position = "bottom")


ggplot()+geom_line(data=overall,aes(x=time,y=trend,group=region,colour=region))+
  geom_point(data=clean,aes(x=time,y=log10(Gross.Loan.Portfolio),colour=Region),alpha=0.1)

######################## Robust variance estimator ########################
clean <- clean[which(!is.na(clean$Operational.self.sufficiency)),]

fit <- lmer(log10(Gross.Loan.Portfolio)~time*Region-Region
             +Operational.self.sufficiency+Scale*time-time
             +Profit.status*time-time+(time|MFI.ID),data=clean)

fit1 <- lme(fixed=log10(Gross.Loan.Portfolio)~time*Region-Region
            +Operational.self.sufficiency+Scale*time-time
            +Profit.status*time-time,
            random=~time|MFI.ID,data=clean)

fit2 <- lme(fixed=log10(Gross.Loan.Portfolio)~time*Region-Region
            +Operational.self.sufficiency+Scale*time-time
            +Profit.status*time-time,
            random=~1|MFI.ID,data=clean)

W <- as.matrix(vcov(fit))
X <- model.matrix(fit)
R <- as.matrix(residuals(fit1,type='response'))
V <- as.matrix(bdiag(getVarCov(fit1,individuals=as.factor(unique(clean$MFI.ID)),type='marginal')))

robust <- W%*%t(X)%*%solve(V)%*%R%*%t(R)%*%solve(V)%*%X%*%W

# Simultaneous confidence intervals
cbind(summary(fit)$coef[6:11,1]-sqrt(diag(robust)[6:11])*qt(1-0.05/12,2000),
  summary(fit)$coef[6:11,1]+sqrt(diag(robust)[6:11])*qt(1-0.05/12,2000))
cbind(10^(summary(fit)$coef[6:11,1]-sqrt(diag(robust)[6:11])*qt(1-0.05/12,2000)),
      10^(summary(fit)$coef[6:11,1]+sqrt(diag(robust)[6:11])*qt(1-0.05/12,2000)))