## Loading Required Packages
# Install packages if necessary
# install.packages("multiway")
# install.packages("ThreeWay")
# install.packages("R.matlab")
# install.packages("plot3D")
# install.packages("plotly")
# install.packages("prospectr")
# install.packages("MASS")
# install.packages("caret")
# install.packages("Stat2Data")
# install.packages("Metrics")
# install.packages("lintools")

library(multiway)
library(ThreeWay)
library(R.matlab)
library(plot3D)
library(plotly)
library(prospectr)
library(MASS)
library(caret)
library(Stat2Data)
library(Metrics)
library(lintools)

## Loading Data
setwd("C:/Path/To/Your/Data")
data <- readMat("Sample_calib.mat")

concentration <- c(0.00001, 0.75, 1.25, 1.5, 2.0, 2.75, 3.0, 3.49, 4.01, 4.48, 5.0, 6.01, 7.0, 7.94, 9.25, 9.97)
                  # Adjust according your data
x <- data$x
y <- concentration
nmEX <- t(data.frame(seq(310, 430, by = 10)))
nmEM <- t(data.frame(seq(301, 700, by = 1)))

## Data Visualization
# Matrix dimensions
mydim <- dim(x)
print(mydim)

# Average matrix
xm <- colMeans(x)

dev.new()
filled.contour(xm, color.palette = terrain.colors)

dev.new()
matplot(t(nmEM), xm, type = "l", xlab = "Emission Wavelength (nm)", ylab = "Intensity")

dev.new()
matplot(t(nmEX), t(xm), type = "l", xlab = "Excitation Wavelength (nm)", ylab = "Intensity")

## PARAFAC Model
set.seed(1)
nf <- 6 # Number of factors
model <- parafac(x, nfac = nf, nstart = 10, maxit = 500, ctol = 1e-5, parallel = FALSE, cl = NULL, output = c("best", "all"))

# R² Adjustment
model$Rsq

## PARAFAC Scores and Loadings Plots
dev.new()
matplot(model$A[, 1], model$A[, 2], pch = "o", col = "blue", xlab = "Factor 1", ylab = "Factor 2", main = "PARAFAC Scores")

dev.new()
matplot(t(nmEM), model$B, type = "l", xlab = "Emission Wavelength (nm)", ylab = "Intensity")

dev.new()
matplot(t(nmEX), model$C, type = "l", xlab = "Excitation Wavelength (nm)", ylab = "Intensity")

## Calibration and Prediction Split using Kennard-Stone
perc <- 0.7
A <- model$A
ntrain <- ceiling(perc * nrow(A))

combined_data <- cbind(A, y)
set.seed(2)
ks_result <- kenStone(combined_data, k = ntrain, metric = "euclid")

calibration_idx <- ks_result$model
prediction_idx <- ks_result$test

xcal <- A[calibration_idx, ]
ycal <- y[calibration_idx]

xpred <- A[prediction_idx, ]
ypred <- y[prediction_idx]

## Calibration and Prediction Split using Mutation-based Learning Method (MLM)
set.seed(1)
ks_result <- kenStone(combined_data, k = ntrain, metric = "euclid")
calibration_idx <- ks_result$model
prediction_idx <- ks_result$test

xcal <- A[calibration_idx, ]
ycal <- y[calibration_idx]

xpred <- A[prediction_idx, ]
ypred <- y[prediction_idx]

# Mutation
prob <- 0.2
p_cal <- ceiling(prob * nrow(xcal))
p_pred <- ceiling(prob * nrow(xpred))

cal_mut <- sample(1:nrow(xcal), size = p_cal)
pred_mut <- sample(1:nrow(xpred), size = p_pred)

xcal_new <- rbind(xcal[-cal_mut, ], xpred[pred_mut, ])
ycal_new <- c(ycal[-cal_mut], ypred[pred_mut])

xpred_new <- rbind(xpred[-pred_mut, ], xcal[cal_mut, ])
ypred_new <- c(ypred[-pred_mut], ycal[cal_mut])

# Outlier Removal
index_outlier <- which(ypred_new == 5.0)
xpred_clean <- xpred_new[-index_outlier, ]
ypred_clean <- ypred_new[-index_outlier]

## Regression Model (MLR)
xcal_df <- data.frame(xcal_new)
colnames(xcal_df) <- paste0("X", 1:ncol(xcal_new))

model_mlr <- lm(ycal_new ~ ., data = xcal_df)

ycal_calc <- predict(model_mlr, newdata = xcal_df)

xpred_df <- data.frame(xpred_clean)
colnames(xpred_df) <- paste0("X", 1:ncol(xpred_clean))
ypred_calc <- predict(model_mlr, newdata = xpred_df)

## Measured vs Predicted Concentration Plots
dev.new()
plot(ycal_new, ycal_calc, pch = 'o', col = "blue", xlab = "Measured Concentration (g/g)", ylab = "Predicted Concentration (g/g)", main = "Calibration")
lines(ycal_new, ycal_new, col = 'red')

dev.new()
plot(ypred_clean, ypred_calc, pch = 15, col = "red", xlab = "Measured Concentration (mg/L)", ylab = "Predicted Concentration (mg/L)", main = "Prediction")
lines(ypred_clean, ypred_clean, col = 'red')

## Model Evaluation
print(summary(model_mlr))

# Calibration Metrics
MAE <- mean(abs(ycal_new - ycal_calc))
R2cal <- cor(ycal_new, ycal_calc)^2
RMSEC <- rmse(ycal_new, ycal_calc)

# Prediction Metrics
MAE_pred <- mean(abs(ypred_clean - ypred_calc))
R2pred <- cor(ypred_clean, ypred_calc)^2
RMSEP <- rmse(ypred_clean, ypred_calc)