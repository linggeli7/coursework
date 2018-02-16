# Mean imputation
missing1 <- which(is.na(communitiesCrimeRaw2$ViolentCrimesPerPop))
communitiesCrimeRaw2$ViolentCrimesPerPop[missing1] <- mean(communitiesCrimeRaw2$ViolentCrimesPerPop,na.rm=TRUE)
missing2 <- which(is.na(communitiesCrimeRaw2$nonViolPerPop))
communitiesCrimeRaw2$nonViolPerPop[missing2] <- mean(communitiesCrimeRaw2$nonViolPerPop,na.rm=TRUE)

# Collapse
library(doBy)
state <- summaryBy(.~state, data=communitiesCrimeRaw2, keep.names=TRUE)

# Principle component analysis
row.names(state) <- state[,1]
fit <- princomp(state[,2:20], cor=TRUE)
summary(fit) 
plot(fit,type="lines") 
biplot(fit,expand=1.2,cex=0.8,xlim=c(-0.6,0.7),ylim=c(-0.3,0.8))

# Factor analysis
fit <- factanal(state[,2:20], 2, scores="regression")
print(fit, digits=2, cutoff=.3, sort=TRUE)
x = fit$scores[,1:2]
y = fit$loadings[,1:2]
biplot(x,y,xlabs=state[,1],expand=1.2,cex=0.8,xlim=c(-2,4),ylim=c(-1,6))

# Independent component analysis
library(fastICA)
plot(SimulatedSignals$X.3~SimulatedSignals$Time,type="l")
a <- fastICA(SimulatedSignals[,2:4],2)
plot(a$S)
plot(a$S[,1]~SimulatedSignals$Time,type='l')
fit <- princomp(SimulatedSignals[,2:4], cor=TRUE)
plot(fit$scores[,1]~SimulatedSignals$Time,type='l')
plot(fit$scores[,2]~SimulatedSignals$Time,type='l')
summary(fit)