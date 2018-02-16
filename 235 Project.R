################ Wandering about ################
head(train[which(train$place_id==8523065625),])
length(unique(train$place_id))

# distributions seem to be uniform
hist(train$x)
hist(train$y)

# look at the first place
place1 <- train[which(train$place_id==8523065625),]
plot(place1$x,place1$y)
hist(place1$time)
hist(log10(place1$accuracy))
boxplot(place1$x)
boxplot(place1$y)
summary(place1$x)
summary(place1$y)

# Heatmaps
ggplot(data=train,aes(x=x,y=y))+geom_point(alpha=0.2)

# if time is measured in minute, we can get daily trends
place1$hour <- (place1$time/60)%%24
hist(place1$hour)
ggplot(data=block.train,aes(x=hour,group=as.factor(place_id),colour=as.factor(place_id)))+
  geom_density()+xlim(0,24)+theme(legend.position='none')

# what is accuracy
library(ggplot2)
place1$accuracy <- cut(place1$accuracy,c(0,56,64,83,1000))
xc <- median(place1$x)
yc <- median(place1$y)
place1$dist <- sqrt((place1$x-xc)^2+(place1$y-yc)^2)
ggplot(data=place1,aes(x=accuracy,y=dist,colour=precision))+
         geom_boxplot()+theme(legend.position='none')+ylim(0,0.1)
# look around then
surround <- train[which(train$x>0.7&train$x<0.9),]
surround <- surround[which(surround$y>9&surround$y<9.2),]
surround$place1 <- ifelse(surround$place_id==8523065625,1,0)

ggplot(data=surround,aes(x=x,y=y,colour=as.factor(place_id)))+
  geom_point()+theme(legend.position='none')

ggplot(data=surround,aes(x=x,y=y,colour=as.factor(place1)))+
  geom_point()+theme(legend.position='none')

surround$hour <- (surround$time/60)%%24
surround$ap <- cut(surround$hour,c(0,6,12,18,24))

ggplot(data=surround,aes(x=x,y=y,colour=as.factor(ap)))+
  geom_point()+theme(legend.position='none')

################ Battle Royale with Cheese ###############
block <- surround
block$hour <- (block$time/60)%%24

I <- sample(1:11000,11000)
block.train <- block[I[1:10000],]
block.test <- block[I[10001:11000],]
top <- as.numeric(names(summary(as.factor(block.train$place_id)))[1:50])
popping <- ifelse(block.train$place_id %in% top,1,0)
block.train <- block.train[which(popping==1),]

# Neural network performs poorly
library(nnet)
places <- class.ind(as.factor(block.train$place_id))
fit <- nnet(block.train[,c(2,3,7)],places,size=3,softmax=TRUE)
a <- predict(fit,block.test[,c(2,3,7)],type="class")
sum(as.numeric(a)==block.test$place_id)
block.test$correct <- ifelse(as.numeric(a)==block.test$place_id,1,0)
ggplot(data=block.test,aes(x=x,y=y,colour=as.factor(correct)))+
  geom_point()+theme(legend.position='none')

# Naive bayes is the current benchmark
library(e1071)
fit.nb <- naiveBayes(as.factor(place_id)~x+y+hour,data=block.train)
a <- predict(fit.nb,block.test[,-6],type='class')
result <- levels(a)[as.numeric(a)]
result <- as.numeric(result)
sum(result==block.test$place_id)
block.test$correct <- ifelse(result==block.test$place_id,1,0)
ggplot(data=block.test,aes(x=x,y=y,colour=as.factor(correct)))+
  geom_point()+theme(legend.position='none')

# SVM does better
fit.svm <- svm(as.factor(place_id)~x+y+hour+accuracy,data=block.train,kernel='radial')
a <- predict(fit.svm,block.test[,-6],type='class')
result <- levels(a)[as.numeric(a)]
result <- as.numeric(result)
sum(result==block.test$place_id)
block.test$correct <- ifelse(result==block.test$place_id,1,0)
ggplot(data=block.test,aes(x=x,y=y,colour=as.factor(correct)))+
  geom_point()+theme(legend.position='none')

# Random forest does even better
library(randomForest)
fit.rf <- randomForest(as.factor(place_id)~x+y+hour+accuracy,data=block.train)
a <- predict(fit.rf,block.test[,-6],type='prob')
result <- matrix(nrow=dim(a)[1],ncol=3)
for (k in 1:dim(a)[1]) {
  b <- sort(a[k,],decreasing=TRUE)
  result[k,] <- as.numeric(colnames(a)[match(b[1:3],a[k,])])
}
sum(block.test$place_id==result[,1])
sum(block.test$place_id==result[,2])
sum(block.test$place_id==result[,3])
block.test$correct <- ifelse(result[,1]==block.test$place_id,1,0)
ggplot(data=block.test,aes(x=x,y=y,colour=as.factor(correct)))+
  geom_point()+theme(legend.position='none')

train$hour <- (train$time/60)%%24
test$hour <- (test$time/60)%%24

################ Clustering #############
total <- rbind(train[,c(2,3)],test[,c(2,3)])
clus1 <- kmeans(train,centers=100,algorithm="Lloyd")
total$level1 <- clus1$cluster
total$level2 <- NA
for (i in 1:100) {
  index <- which(total$level1==i)
  clus2 <- kmeans(total[index,c(1,2)],centers=100,algorithm="Lloyd")
  total[index,4] <- clus2$cluster
}
area <- (total$level1-1)*100+total$level2
area.test <- area[29118022:37725251]
area <- area[1:29118021]

############### Divide into blocks ###################

east <- floor((train$x)/0.1)
east <- ifelse(east==0,1,east)
north <- floor((train$y)/0.1)
north <- ifelse(north==0,1,north)
area <- (east-1)*100+north

east.test <- floor((test$x)/0.1)
east.test <- ifelse(east.test==0,1,east.test)
north.test <- floor((test$y)/0.1)
north.test <- ifelse(north.test==0,1,north.test)
area.test <- (east.test-1)*100+north.test

############### For realz ###############

pred <- matrix(nrow=8607230,ncol=3)

# Naive bayes
for (i in 1:10000) {
  index <- which(area==i)
  index.test <- which(area.test==i)
  if (length(index)>2&length(index.test)>0) {
    hood <- train[index,]
    hood.test <- test[index.test,]
    top <- as.numeric(names(summary(as.factor(hood$place_id)))[1:20])
    popping <- ifelse(hood$place_id %in% top,1,0)
    hood.train <- hood[which(popping==1),]
    fit <- naiveBayes(as.factor(place_id)~x+y+hour,data=hood.train)
    a <- predict(fit,hood.test,type='raw')
    result <- matrix(nrow=dim(a)[1],ncol=3)
    for (k in 1:dim(a)[1]) {
      b <- sort(a[k,],decreasing=TRUE)
      result[k,] <- as.numeric(colnames(a)[match(b[1:3],a[k,])])
    }
    pred[index.test,] <- result
  }
  print(i)
}

# Random forest
for (i in 1:10000) {
  index <- which(area==i)
  index.test <- which(area.test==i)
  if (length(index)>2&length(index.test)>0) {
    hood <- train[index,]
    hood.test <- test[index.test,]
    top <- as.numeric(names(summary(as.factor(hood$place_id)))[1:20])
    popping <- ifelse(hood$place_id %in% top,1,0)
    hood.train <- hood[which(popping==1),]
    
    fit.rf <- randomForest(as.factor(place_id)~x+y+hour+accuracy,data=hood.train)
    a <- predict(fit.rf,hood.test,type='prob')
    result <- matrix(nrow=dim(a)[1],ncol=3)
    for (k in 1:dim(a)[1]) {
      b <- sort(a[k,],decreasing=TRUE)
      result[k,] <- as.numeric(colnames(a)[match(b[1:3],a[k,])])
    }
    
    pred[index.test,] <- result
  }
  print(i)
}

################## Results ##################
place <- apply(pred, 1, paste, collapse=" ") 
final <- data.frame(row_id=sample_submission$row_id,place_id=place)
head(final)
write.csv(final, file ="rando.csv",row.names=FALSE)
