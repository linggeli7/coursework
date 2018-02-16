library(survival)
set.seed(12345)
true.time <- rexp(1000, rate = rep(c(1.2,1,.8,.6),each=250))
cens.time <- runif(1000,0,4)
obs.time <- ifelse( true.time <= cens.time, true.time, cens.time )
delta <- ifelse( obs.time==true.time, 1, 0 )
group <- rep( 1:4, each=250 )
kmPlot( fit=survfit(Surv(obs.time,delta) ~ group), xscale=1, lty=1:4,
          labelTimes=0:4, groupLabels=paste( "Group", 1:4) )
legend( 2, .9, lty=1:4, legend=paste("Group", 1:4), bty="n" )
fit=survfit(Surv(obs.time,delta) ~ group)
plot(fit)

kmPlot <- function(fit, xscale=1, lty=1, xlab="Time", ylab="Survival", labelTimes, groupLabels) {
  fit$time <- fit$time/xscale
  nStrata <- length( fit$strata )
  par(mar=c(nStrata+8,4,2,4))
  plot(fit, lty=lty, xlab=xlab, ylab=ylab, mark.time=FALSE)
  for (j in 1:length(labelTimes)) {
    start <- 1
    end <- 1
    riskTotal <- 0
    eventTotal <- 0
    for (i in 1:nStrata) {
      end <- start + fit$strata[i] - 1
      if (labelTimes[j] < fit$time[start]) {
        nrisk <- fit$n.risk[start]
        nevent <- 0
      } else if (labelTimes[j] > fit$time[end]) {
        nrisk <- 0
        nevent <- sum(fit$n.event[start:end])
      } else {
        index <- as.integer(start) + max( which(fit$time[start:end] <= labelTimes[j]))
        nrisk <- fit$n.risk[index]
        nevent <- sum(fit$n.event[start:index-1])
      }
      if (j == 1) {
        mtext(paste(groupLabels[i]," ",nrisk,"(",nevent,")"), side=1, line=i+4, at=labelTimes[j], cex=0.8)
      } else {
        mtext(paste(nrisk,"(",nevent,")"), side=1, line=i+4, at=labelTimes[j], cex=0.8)
      }
      riskTotal <- riskTotal + nrisk
      eventTotal <- eventTotal + nevent
      start <- start + fit$strata[i]
    }
    if (j == 1) {
      mtext(paste("Total ",riskTotal,"(",eventTotal,")"), side=1, line=nStrata+6, at=labelTimes[j], cex=0.8)
    } else {
      mtext(paste(riskTotal,"(",eventTotal,")"), side=1, line=nStrata+6, at=labelTimes[j], cex=0.8)
    }
  } 
}