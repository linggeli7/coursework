coindesk.bpi.USD.close_data.2016.02.13_2017.02.13 <- read.csv('~/Documents/2016/Time Series/coindesk-bpi-USD-close_data-2016-02-13_2017-02-13.csv')
bitcoin <- as.ts(coindesk.bpi.USD.close_data.2016.02.13_2017.02.13$Close.Price[4:367])
plot(bitcoin, ylab='Price')
acf(bitcoin, lag.max=50)
plot(diff(bitcoin))
plot(log(bitcoin))
plot(diff(log(bitcoin)))

acf(diff(log(bitcoin)), lag.max=50, ylim=c(-0.2, 0.2))
pacf(diff(log(bitcoin)), lag.max=50)
acf(diff(log(bitcoin)), lag.max=50)
pacf(diff(log(bitcoin)), lag.max=50)

periodogram(diff(log(bitcoin)))

library('astsa')
# sarima
sarima(log(bitcoin), 0, 1, 2, 1, 0, 0, 10, no.constant=TRUE, details=FALSE)
sarima.pred <- sarima.for(log(bitcoin), 30, 0, 1, 2, 1, 0, 0, 10, no.constant=TRUE)$pred
sarima.pred <- exp(sarima.pred)

library('quantmod')

getFX('USD/CNY', from='2016-02-13', to='2017-03-13')
getSymbols('^GSPC', src='yahoo', from='2016-02-13', to='2017-03-13')
getMetals(Metals='gold', base.currency='USD', from='2016-02-13', to='2017-03-13')

plot(as.ts(USDCNY[1:364]), ylab='USD/CNY')
acf(USDCNY)
plot(diff(USDCNY))
acf(diff(USDCNY))

plot(as.ts(XAUUSD[1:364]), ylab='Gold')
acf(XAUUSD)
plot(diff(XAUUSD))

SP500 <- GSPC$GSPC.Adjusted
plot(as.ts(X[4:367, 1]), ylab='SP500')


# match dates
X <- merge(SP500, USDCNY, XAUUSD)
X <- as.matrix(X)
# impute missing weekend/holiday values with previous values
temp1 <- c(0, X[1:394, 1])
temp2 <- c(0, 0, X[1:393, 1])
temp3 <- c(0, 0, 0, X[1:392, 1])
X[which(is.na(X[, 1])), 1] <- temp1[which(is.na(X[, 1]))]
X[which(is.na(X[, 1])), 1] <- temp2[which(is.na(X[, 1]))]
X[which(is.na(X[, 1])), 1] <- temp3[which(is.na(X[, 1]))]

X_past <- X[4:367, ]
X_future <- X[368:395, ]

plot(bitcoin ~ X_past[, 1])
plot(bitcoin ~ X_past[, 2])
plot(bitcoin ~ X_past[, 3])
plot(bitcoin ~ X_past[, 4])

past <- cbind(bitcoin, X_past)

pairs(past)

# linear regression plus sarima
fit <- lm(bitcoin ~ 0 + cbind(seq(1, 364, 1), X_past))
summary(fit)
res <- as.ts(fit$residuals)
plot(res)
acf(res, lag.max=50)
pacf(res, lag.max=50)
plot(diff(res))
acf(diff(res), lag.max=50)
pacf(diff(res), lag.max=50)

fit <- arima(res, c(2, 0, 0))
sarima.for(res, 30, 1, 1, 0)

pred <- cbind(365:392, X_future) %*% as.matrix(c(1.23666, 0.11610, 70.48203, -0.24736)) + 100
plot.ts(pred)

# vector autoregression
library('vars')

past <- cbind(bitcoin, X_past)
fitvar <- VAR(past, p=1, type='both')
summary(fitvar)
acf(residuals(fitvar)[, 1])

prd <- predict(fitvar, n.ahead=30, ci=0.95)
print(prd)
plot(prd, 'single')

# ARMA GARCH
library(fGarch)

garchres <- garchFit(formula=~arma(1, 1)+garch(1, 1), data=diff(log(bitcoin)), include.mean=FALSE)
plot(garchres)
predict(garchres, n.ahead=30, plot=TRUE)

garchdiff <- garchFit(formula=~arma(1, 1)+garch(1, 1), data=diff(bitcoin))
predict(garchdiff, n.ahead=30, plot=TRUE)

plot.ts(garchdiff@sigma.t^2)

# dynamic linear model
library(dlm)

quant <- dlmModReg(X_past)
dlmfit <- dlmSmooth(bitcoin, quant)

buildCapm <- function(u) {
  dlmModReg(as.ts(diff(SP500)), dV=exp(u[1]), dW=exp(u[2:3]))
}

outMLE <- dlmMLE(diff(bitcoin), parm=c(0,0,0), buildCapm)

# log return
plot(log(bitcoin))
acf(log(bitcoin), lag.max=50)
plot(diff(log(bitcoin)))
acf(diff(bitcoin), lag.max=50)
pacf(diff(bitcoin), lag.max=50)

######### FINAL PLOT ########
library(ggplot2)

bitcoin_future <- rep(0, 28)
for (i in 1:28) {
  bitcoin_future[i] <- coindesk.bpi.USD.close_data.2017.02.13_2017.03.12$Close.Price[1+24*(i-1)]
}

final <- cbind(bitcoin_future, pred[, 1], sarima.pred[1:28], var.pred[1:28])
final <- as.data.frame(final)

colnames(final) <- c('Bitcoin', 'Regression', 'SARIMA', 'VAR')
head(final)

final <- melt(final)
final$time <- rep(1:28, 4)
colnames(final)[1] <- 'model'

ggplot(data=final) +
  geom_line(aes(x=time, y=value, colour=model, linetype=model)) +
  theme_bw() +
  theme(legend.position='bottom')
