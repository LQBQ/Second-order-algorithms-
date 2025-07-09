#### NOTE: This script can only be used after the 3D matrix has been assembled (using "[5] Function to Create a 3D Array")
#### Note that the PARAFAC-MLR matrix considers just all of adulterated samples and one pure sample in this case.

### PARAFAC-MLR
### (Parallel factor analysis - multiple linear regression)

# Function to automatically install (if needed) and load required packages
load_required_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# Vector of required packages
required_packages <- c(
  "multiway", "ThreeWay", "R.matlab", "plot3D", "plotly", 
  "prospectr", "MASS", "caret", "Stat2Data", "Metrics", "lintools"
)

# Run the function to ensure all packages are installed and loaded
load_required_packages(required_packages)

# Define the input directory where the 3D matrix are attached
input_directory <- "path/to/your/directory"   # Adjust the path as needed
# Example: input_directory <- "C:/Users/Chemometric-PC/Documents/IC/Project - EVOO adulteration/Samples/Preprocessed samples"
# For compatibility, replace backslashes "\" with forward slashes "/" in your path.

concentration <- c(0.00001, 0.75, 1.25, 1.5, 2.0, 2.75, 3.0, 3.49, 4.01, 4.48, 5.0, 6.01, 7.0, 7.94, 9.25, 9.97)
# Here, enter the vector of concentrations according to the dataset you are working with.

x <- data$x
y <- concentration
nmEX <- t(data.frame(seq(310, 430, by = 10)))
nmEM <- t(data.frame(seq(301, 700, by = 1)))
# Here, the emission number has changed due to the clipping of the excitation-emission matrix performed in preprocessing. 
# Change the values ​​according to your dataset, if necessary.
# Make sure that the emission number is in accordance with the second dimension of the 3D matrix.

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

# Outlier Removal (opcional)
# Remove any sample if you identify that it is the only one that does not fit correctly to the regression line.
index_outlier <- which(ypred_new == 5.0) # In this case, the sample was removed.
xpred_clean <- xpred_new[-index_outlier, ]
ypred_clean <- ypred_new[-index_outlier]
# WARNING: This specific sample was removed considering the separation into training and testing 
# performed on the day of analysis by the authors of this work 
# and the results presented by the model for this specific set. 
# Consider that when separating the data into training and testing, the set may change. 
# Therefore, analyze the results before removing any outliers.

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
