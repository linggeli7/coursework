chsData$sbp140 <- ifelse(chsData$sbp > 140, 1, 0)
chsData$sex <- ifelse(chsData$gender == 1, "M", "F")

# pkyrs, block0 and kcal0 heavily right skewed
summary(chsData)
# weight and sbp distributions
hist(chsData$weight)
hist(chsData$sbp)
abline(v=140, col = "red")
# Overall weight vs sbp140
plot(chsData$sbp ~ chsData$weight)
boxplot(weight ~ sbp140, xlab="High blood pressure", ylab="Weight", data=chsData[which(chsData$gender==1),])
boxplot(weight ~ sbp140, xlab="High blood pressure", ylab="Weight", data=chsData[which(chsData$gender==0),])

chsData$weightclass <-cut(chsData$weight, c(70,135,156,178,350))
barplot(table(chsData$sbp140, chsData$weightclass),col=c("lightblue","red"), 
        xlab="Weight", beside=TRUE)

# MAYBE EFFECT MODIFIER?
# Gender associated with weight but not so much SBP
boxplot(chsData$weight ~ chsData$sex, xlab="Gender", ylab="Weight")
table(chsData$sbp140, chsData$sex)
barplot(table(chsData$sbp140, chsData$sex),col=c("lightblue","red"),
        legend = c("Normal", "High"), xlab="Gender", beside=TRUE)

# PRECISION
# Age associated with SBP but not weight
chsData$agegroup <-cut(chsData$age, c(64.9,68,71,75,100))
table(chsData$sbp140, chsData$agegroup)
barplot(table(chsData$sbp140, chsData$agegroup),col=c("lightblue","red"),
        xlab="Age Group", legend = c("Normal", "High"), beside=TRUE)
boxplot(chsData$weight ~ chsData$agegroup, xlab="Age Group", ylab="Weight")

# THEORECTICAL CONFOUNDER
# Exercise probably measured best as kcal0
chsData$pe <-cut(chsData$kcal0, c(-1,214,735,1768,15000))
table(chsData$sbp140, chsData$pe)
barplot(table(chsData$sbp140, chsData$pe))
table(chsData$sbp140, chsData$exint0)
barplot(table(chsData$sbp140, chsData$exint0), col=c("lightblue","red"),
        legend = c("Normal", "High"), xlab="Exercise Intensity")
# Exercise does not seem to affect weight or blood pressure that much
boxplot(chsData$weight ~ chsData$pe)
boxplot(chsData$weight ~ chsData$exint0, xlab="Exercise Intensity", ylab="Weight")
boxplot(chsData$sbp ~ chsData$exint0)
plot(chsData$kcal0 ~ chsData$block0)
# exercise intensity not good
boxplot(chsData$block0 ~ chsData$exint0)
boxplot(chsData$kcal0 ~ chsData$exint0)

# CONFOUNDER
# Diabetics have higher blood pressure and weight
table(chsData$sbp140, chsData$diab)
barplot(table(chsData$sbp140, chsData$diab), col=c("lightblue","red"),
        legend = c("Normal", "High"), xlab="Diabetes")
boxplot(chsData$weight ~ chsData$diab, ylab="Weight", xlab="Diabetes")

# CONFOUNDER
# smoking increases SBP but lowers weight
chsData$smoke <-cut(chsData$pkyrs, c(-1,0,27,250))
table(chsData$sbp140, chsData$smoke)
barplot(table(chsData$sbp140, chsData$smoke), col=c("lightblue","red"),
        legend = c("Normal", "High"), xlab="Smoking status")
boxplot(chsData$weight ~ chsData$smoke, ylab="Weight", xlab="Smoking status")




