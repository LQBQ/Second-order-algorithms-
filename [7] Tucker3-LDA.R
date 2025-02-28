## Loading Required Packages
library(multiway)
library(ThreeWay)
library(R.matlab)
library(plot3D)
library(plotly)
library(prospectr)
library(MASS)
library(caret)

## Setting Working Directory
# In RStudio, go to Session > Set Working Directory > Choose Directory
dir <- "C:/Path/To/Your/Data"
setwd(dir)

## Loading Data
data <- readMat("SampleData.mat")

X <- data$x
A <- data.frame(matrix(1, nrow = 12, ncol = 1)) # nrow = number of samples class 1
colnames(A) <- "Class"
B <- data.frame(matrix(2, nrow = 15, ncol = 1)) # nrow = number of samples class 2
colnames(B) <- "Class"

Y <- rbind(A, B)

nmEX <- t(data.frame(seq(310, 430, by = 10)))
nmEM <- t(data.frame(seq(301, 700, by = 1)))

## Visualizing Data
# Matrix dimensions
mydim <- dim(X) # samples x emission x excitation
print(mydim)

# Average matrix of each class
data2 <- readMat("Class1.mat")
data3 <- readMat("Class2.mat")

X1 <- data2$x
X2 <- data3$x

X1m <- colMeans(X1)
X2m <- colMeans(X2)

# Visualization
filled.contour(X1m, color.palette = terrain.colors, main = "Class 1")
filled.contour(X2m, color.palette = terrain.colors, main = "Class 2")

matplot(t(nmEM), X1m, type="l", xlab="Emission Wavelength (nm)", ylab="Intensity", main="Class 1")
matplot(t(nmEX), t(X1m), type="l", xlab="Excitation Wavelength (nm)", ylab="Intensity", main="Class 1")

matplot(t(nmEM), X2m, type="l", xlab="Emission Wavelength (nm)", ylab="Intensity", main="Class 2")
matplot(t(nmEX), t(X2m), type="l", xlab="Excitation Wavelength (nm)", ylab="Intensity", main="Class 2")

## Tucker3 Model
set.seed(1)
nf <- 6 # Number of factors

Xr <- matrix(X, nrow = mydim[1])

model <- T3func(Xr, mydim[1], mydim[2], mydim[3], nf, nf, nf, 1, 1e-7)

# R2 Adjustment
print(model$fp)

## Plot Tucker3 Scores
matplot(model$A[Y == 1, 1], model$A[Y == 1, 2], pch="o", col="blue", xlab="Factor 1", ylab="Factor 2", main="TUCKER-3 Scores")
points(model$A[Y == 2, 1], model$A[Y == 2, 2], pch="x", col="red")
legend("topleft", legend = c("Class 1", "Class 2"), col = c("red", "blue"), lty = 1, lwd = 3)

## Plot Tucker3 Loadings
matplot(t(nmEM), model$B, type="l", xlab="Emission Wavelength (nm)", ylab="Intensity")
matplot(t(nmEX), model$C, type="l", xlab="Excitation Wavelength (nm)", ylab="Intensity")

## Training and Testing Sample Selection
perc <- 0.7 # 70% training

ntrain1 <- ceiling(perc * dim(X1)[1])
ntrain2 <- ceiling(perc * dim(X2)[1])

scores_X1 <- model$A[1:dim(X1)[1],]
scores_X2 <- model$A[(dim(X1)[1] + 1):dim(X)[1],]

sel1 <- kenStone(scores_X1, k = ntrain1)
sel2 <- kenStone(scores_X2, k = ntrain2)

train <- rbind(scores_X1[sel1$model,], scores_X2[sel2$model,])
test <- rbind(scores_X1[sel1$test,], scores_X2[sel2$test,])

group_train <- rbind(matrix(rep(1, ntrain1)), matrix(rep(2, ntrain2)))
group_test <- rbind(matrix(rep(1, dim(X1)[1] - ntrain1)), matrix(rep(2, dim(X2)[1] - ntrain2)))

## LDA Model
model_lda <- lda(train, group_train)
model_lda_cv <- lda(train, group_train, CV=TRUE)

## Predictions
pred_train <- predict(model_lda, train)
pred_test <- predict(model_lda, test)

## Performance Metrics
ac_train <- mean(pred_train$class == group_train)
spec_train <- mean(pred_train$class[1:ntrain1] == rep(1, ntrain1))
sens_train <- mean(pred_train$class[(ntrain1 + 1):(ntrain1 + ntrain2)] == rep(2, ntrain2))

ac_test <- mean(pred_test$class == group_test)
spec_test <- mean(pred_test$class[1:(dim(X1)[1] - ntrain1)] == rep(1, dim(X1)[1] - ntrain1))
sens_test <- mean(pred_test$class[(dim(X1)[1] - ntrain1 + 1):length(group_test)] == rep(2, dim(X2)[1] - ntrain2))

ac_cv <- mean(model_lda_cv$class == group_train)
spec_cv <- mean(model_lda_cv$class[1:ntrain1] == rep(1, ntrain1))
sens_cv <- mean(model_lda_cv$class[(ntrain1 + 1):(ntrain1 + ntrain2)] == rep(2, ntrain2))

cat("Accuracy:", ac_train, "\nSensitivity:", sens_train, "\nSpecificity:", spec_train, "\n")
cat("Accuracy:", ac_test, "\nSensitivity:", sens_test, "\nSpecificity:", spec_test, "\n")
cat("Accuracy:", ac_cv, "\nSensitivity:", sens_cv, "\nSpecificity:", spec_cv, "\n")

## Visualization
boxplot(pred_train$x[, 1] ~ group_train, col = c("blue", "red"), xlab = "Group", ylab = "LD1", main = "Boxplot (Training)")
boxplot(pred_test$x[, 1] ~ group_test, col = c("blue", "red"), xlab = "Group", ylab = "LD1", main = "Boxplot (Testing)")

# Confusion Matrices
confusion_train <- table(Predicted = pred_train$class, Actual = group_train)
print("Confusion Matrix - Training Set")
print(confusion_train)

confusion_test <- table(Predicted = pred_test$class, Actual = group_test)
print("Confusion Matrix - Test Set")
print(confusion_test)
