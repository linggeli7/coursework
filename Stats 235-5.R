diabetes <- pima.indians.diabetes
colnames(diabetes) <- c('Pregnancy','Glucose','BP','Tricep','Insulin','BMI','Pedigree','Age','Diabetes')
I <- sample(1:768,768)
train <- diabetes[I[1:400],]
test <- diabetes[I[401:768],]

# Conditional distribution
library(car)
scatterplotMatrix(~Glucose+Insulin+BP+BMI|Diabetes,data=train)

# Visualize decision boundaries
library(klaR)
partimat(as.factor(Diabetes)~Glucose+Insulin+BP+BMI,data=train,method="lda")
partimat(as.factor(Diabetes)~Glucose+Insulin+BP+BMI,data=train,method="qda")
partimat(as.factor(Diabetes)~Glucose+Insulin+BP+BMI,data=train,method="naiveBayes")

# Discriminant analysis
library(MASS)
fit.lda <- lda(Diabetes~.,data=train)
prob.lda <- predict(fit.lda,test[,-9])
fit.qda <- qda(Diabetes~.,data=train)
prob.qda <- predict(fit.qda,test[,-9])

# Logistic regression
fit.glm <- glm(Diabetes~.,family=binomial(),data=train)
prob.glm <- predict(fit.glm,test[,-9],type = "response")

# Naive bayes
library(e1071)
fit.nb <- naiveBayes(as.factor(Diabetes)~.,data=train)
prob.nb <- predict(fit.nb,test[,-9],type='raw')

# ROC curves
library(ROCR)
pred.lda <- prediction(prob.lda$posterior[,2],test[,9])
perf.lda <- performance(pred.lda,"tpr","fpr")
pred.qda <- prediction(prob.qda$posterior[,2],test[,9])
perf.qda <- performance(pred.qda,"tpr","fpr")
pred.glm <- prediction(prob.glm,test[,9])
perf.glm <- performance(pred.glm,"tpr","fpr")
pred.nb <- prediction(prob.nb[,2],test[,9])
perf.nb <- performance(pred.nb,"tpr","fpr")
plot(perf.lda,col=4)
plot(perf.qda,add=TRUE,col=6)
plot(perf.glm,add=TRUE,col=3)
plot(perf.nb,add=TRUE,col=2)
legend(x=0.6,y=0.4,c('LDA','QDA','Logistic','Naive Bayes'),col=c(4,6,3,2),lwd=1,bty='n')
performance(pred.lda,"auc")
performance(pred.qda,"auc")
performance(pred.glm,"auc")
performance(pred.nb,"auc")
