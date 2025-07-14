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

## Load Data
data <- readMat("Sample_Data.mat") # Change the archive name to those that you chose 
mat <- data$x # Access the data matrix in the .mat file

### If you dataset have a mix of samples (regression and calibration samples), you must need perform a average of calibration samples, 
# and then you must concatenate the average samples matrix and the regression matrix (the other samples).

# Average the first 12 samples (pure samples) along the first dimension (row/samples)
# If you working with other dataset matrix, you must analyze your dataset and change the parameters 
average_m <- apply(Mat[1:12, , ], c(2, 3), mean)  # result: a 2D matrix (emission x excitation)

# Resize the 3D matrix with one sample (1 x emission x excitation)
average_3D <- array(average_m, dim = c(1, dim(X)[2], dim(X)[3]))

# Select the other samples (in this study, from 13 onwards)
reg_samples <- X[13:dim(X)[1], , ]

# Concatenate the avarege array into the remainder
x <- abind::abind(average_3D, reg_samples, along = 1)

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

################################################################################
## PARAFAC-MLR – Script Explanation
##
## This script implements a chemometric modeling pipeline combining
## Parallel Factor Analysis (PARAFAC) for feature extraction and
## Multiple Linear Regression (MLR) for quantitative prediction of analyte concentration.
##
## 1. Package Management:
##    - Loads and installs all necessary packages automatically for data handling,
##      modeling, and visualization.
##
## 2. Data Input:
##    - Loads a three-way data array (`x`) from a preprocessed fluorescence excitation-emission matrix (EEM).
##    - Defines the `concentration` vector representing known analyte concentrations.
##    - Sets emission (nmEM) and excitation (nmEX) axes.
##
## 3. Visualization:
##    - Visual inspection of the average sample spectrum using `filled.contour` and `matplot`
##      to assess general signal characteristics across emission and excitation modes.
##
## 4. PARAFAC Decomposition:
##    - Applies PARAFAC decomposition (`nfac = 6`) to reduce the three-way array into
##      three matrices (A: scores; B: emission loadings; C: excitation loadings).
##    - Visualizes sample distribution in the component space and spectral profiles.
##
## 5. Calibration/Test Split:
##    - Uses the Kennard-Stone algorithm to objectively divide the dataset into
##      calibration and prediction subsets, maintaining space-filling properties.
##
## 6. Mutation-Based Learning Method (MLM):
##    - Introduces a mutation step by swapping a fraction (20%) of calibration and
##      prediction samples to simulate variability and increase model robustness.
##
## 7.  Outlier Handling (Optional):
##    - Allows manual exclusion of samples (e.g., with known prediction deviation),
##      although such exclusions should be validated using objective residual/leverage analysis.
##
## 8.  Regression Model:
##    - Fits an MLR model using the calibration subset of PARAFAC scores (`A`) as predictors
##      and concentration (`ycal_new`) as response variable.
##
## 9.  Prediction & Validation:
##    - Predicts concentrations for both calibration and test samples.
##    - Produces diagnostic plots: measured vs. predicted concentrations.
##
## 10. Performance Metrics:
##    - Computes MAE, RMSE, and R² for both calibration and prediction datasets to assess
##      model performance quantitatively.
##
##  # Final Remarks:
##    - This script is structured to provide a complete end-to-end workflow for
##      spectroscopic regression using tensor decomposition and linear modeling.
##    - It assumes that the user understands the physical meaning of EEM matrices
##      and ensures consistent alignment between concentration vectors and sample dimensions.
##    - Adjust `nfac`, preprocessing, and mutation ratio according to your dataset characteristics.
##
################################################################################

