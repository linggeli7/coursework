# Create data frame
data <- communities[,c(4,7,17,18,34,39,37,128)]
crime <- as.data.frame(data)
colnames(crime) <- c("name","householdsize","pctUrban","medIncome","PctPopUnderPov",
                     "PctEmploy","PctBSorMore","ViolentCrimesPerPop")

# K-means
K <- 4
clus <- NA
X <- crime[,2:8]
# Initializing the centroids 
Cent <- replicate(7, runif(K)) 
N <- 1994
D <- matrix(NA, N, K)

for(i in 1:500){
  # Calculating the distance of each data point from the K centroids.   
  for (k in 1:K){
    D[, k] <- rowSums((sweep(X, 2, Cent[k, ]))^2)
  }
  
  # Assigning the data points to the closest centroid  
  clus <- apply(D, 1, which.min)
  
  # Finding the new centroids   
  Cent <- by(X, INDICES=clus, FUN=colMeans)  
  Cent <- do.call(rbind, Cent)  
  
}  

crime[c(1939,967,1088,522),]
Cent

# Complete linkage seems to be the best
hc <- hclust(dist(X), "complete")
dend <- as.dendrogram(hc)

plot(dend, leaflab = "none")
abline(h=4, lty=2, lwd=2)
clus <- cutree(hc, k=20)
qplot(clus, geom="histogram") 


data <- data.frame(crime, clus = factor(clus))

p <- ggplot(data, aes(pctUrban, medIncome))
p + geom_point(aes(color=clus)) 

Cent <- by(X, INDICES=clus, FUN=colMeans)  
Cent <- do.call(rbind, Cent)  

# Soft K-means clustering with K=4
K = 4

# Initializing the parameters 
pi <- rep(1/K, K)
mu <- cbind(runif(K, min(X[, 1]), max(X[, 1])), runif(K, min(X[, 2]), max(X[, 2])) )
sigma <- array(rbind(c(1, 0), c(0, 1)), c(2, 2, K))

d <- matrix(NA, N, K)
p <- matrix(NA, N, K)

for(i in 1:20){
  # E-step   
  for (k in 1:K){
    d[, k] <- pi[k]*dmvnorm(X, mu[k, ], sigma[, , k])
  }
  
  total.d <- rowSums(d)  
  for (k in 1:K){
    p[, k] <- pi[k]*dmvnorm(X, mu[k, ], sigma[, , k])/total.d
  }
  
  
  # M-step
  pi <- colSums(p)/N
  
  sum.p <- colSums(p)  
  for (k in 1:K){
    mu[k, ] <- colSums(p[, k]*X)/sum.p[k]
  }
  
  for (k in 1:K){
    temp <- 0
    for(j in 1:N){
      temp <- temp + p[j, k]*(X[j, ] - mu[k, ])%*%t(X[j, ] - mu[k, ])
    }
    sigma[, , k] <- temp/sum.p[k]
  }
  
}  

# Gaussian mixture
mod1 = Mclust(crime[,2:8], G=4)
summary(mod1, parameters=TRUE)

