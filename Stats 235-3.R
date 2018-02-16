crime <- scale(communitiesCrimeRaw2[,4:21])
crime <- cbind(communitiesCrimeRaw2[,1:3],crime)
crime <- crime[which(!is.na(crime$ViolentCrimesPerPop)),]

# Linear regression with stepwise selection
fit1 = lm(ViolentCrimesPerPop~.-fold-communityname-state,data=crime)
summary(fit1)

RSS1 <- 0
for (i in 1:319) {
  test <- crime[i,4:21]
  train <- crime[-i,4:21]
  fit1 <- lm(ViolentCrimesPerPop~.,data = train)
  y1 <- predict(fit1,test)
  RSS1 <- RSS1+sum((test$ViolentCrimesPerPop-y1)^2)
}

RSS1 <- rep(NA,10)
for (i in 1:10) {
  test <- crime[which(crime$fold==i),4:21]
  train <- crime[which(crime$fold!=i),4:21]
  fit1 <- lm(ViolentCrimesPerPop~., data = train)
  y1 <- predict(fit1,test)
  RSS1[i] <- sum((test$ViolentCrimesPerPop-y1)^2)
}

# Principle component regression with spree plot
pca <- princomp(crime[,4:20], cor=TRUE)
summary(pca) 
plot(pca,type="lines") 
fit2 <- lm(crime$ViolentCrimesPerPop~pca$scores[,1:10])
as.matrix(pca$loadings[,1:10])%*%fit2$coefficients[2:11]

RSS2 <- 0
for (i in 1:319) {
  test <- as.data.frame(t(pca$scores[i,1:10]))
  train <- as.data.frame(pca$scores[-i,1:10])
  ViolentCrimes <- crime$ViolentCrimesPerPop[-i]
  train <- cbind(ViolentCrimes,train)
  fit2 <- lm(ViolentCrimes~.,data=train)
  y2 <- predict(fit2,test)
  RSS2 <- RSS2 +sum((crime$ViolentCrimesPerPop[i]-y2)^2)
}

# Partial least squares regression
library(pls)
fit3 <- plsr(ViolentCrimesPerPop~.-fold-communityname-state,ncomp=10,data=crime,model=TRUE)
sum((predict(fit3, ncomp=10)-crime$ViolentCrimesPerPop)^2)

RSS3 <- 0
for (i in 1:319) {
  test <- crime[i,]
  train <- crime[-i,]
  fit3 <- plsr(ViolentCrimesPerPop~.-fold-communityname-state,ncomp=10,data=train)
  y3 <- predict(fit3,ncomp=10,test)
  RSS3 <- RSS3+sum((test$ViolentCrimesPerPop-y3)^2)
}

# Ridge regression
library(MASS)
lambda <- seq(0,100,1)
fit4 <- lm.ridge(ViolentCrimesPerPop~.-fold-communityname-state,data=crime,lambda=lambda)
X <- as.matrix(crime[,4:20])
df <- NULL
for(i in 1:length(lambda)){
  df[i] <- sum(diag((X)%*%solve(t(X)%*%(X) + lambda[i]*diag(17))%*%t(X)))	
}
beta <- coef(fit4)
matplot(df, beta[, 1:18], ylim=c(-1,1), type='l', xlab=expression(df(lambda)), ylab=expression(beta),
        main="Ridge regressoin")
X <- cbind(rep(1,319),X)
RSS <- NULL
for (i in 1:100) {
  RSS[i] <- sum((X%*%as.matrix(beta[i,])-crime$ViolentCrimesPerPop)^2)
}
plot(RSS,type="l")

lambda <- seq(0,100,1)
RSS4 <- rep(NA,101)
for (l in lambda) {
  RSS4[l] <- 0
  for (i in 1:10) {
    test <- crime[which(crime$fold==i),]
    train <- crime[which(crime$fold!=i),]
    fit4 <- lm.ridge(ViolentCrimesPerPop~.-fold-communityname-state,data=train,lambda=l)
    beta <- coef(fit4)
    X <- as.matrix(test[,4:20])
    X <- cbind(1,X)
    RSS4[l] <- RSS4[l]+sum((X%*%as.matrix(beta)-test$ViolentCrimesPerPop)^2)
  }
}
plot(RSS4~lambda,type="l",ylab="10-fold Cross Validation",main="Ridge regression")

# LASSO
library(lars)
fit5 <- lars(as.matrix(crime[,4:20]),crime$ViolentCrimesPerPop,type="lasso")
plot(fit5)
summary(fit5)
fit5$lambda
fit5$RSS
predict(fit5,type="coefficients",s=0.1,mode="lambda")

cv.lars(as.matrix(crime[,4:20]),crime$ViolentCrimesPerPop,type="lasso")
RSS.lasso <- 0
for (i in 1:319) {
  test <- crime[i,4:20]
  train <- crime[-i,4:20]
  fit5 <- lars(as.matrix(train),crime$ViolentCrimesPerPop[-i],type="lasso")
  y5 <- predict(fit5,test,s=0,mode="lambda")
  RSS.lasso <- RSS.lasso+sum((crime$ViolentCrimesPerPop[i]-y5$fit)^2)
}
sum(RSS.lasso)