# Install and load required packages
install.packages("R.matlab")
library(R.matlab)

install.packages("pls")
library(pls)

install.packages("Stat2Data")
library(Stat2Data)

install.packages("Metrics")
library(Metrics)

install.packages("sNPLS")
library(sNPLS)

# Set working directory and load dataset
setwd("C:/Path/To/Your/Data")
data <- readMat("Sample_calib.mat")  # Adjust dataset name
concentration <- c(0.00001, 0.75 , 1.25, 1.5, 2.0, 2.75, 3.0, 3.49, 4.01, 4.48, 5.0, 6.01, 7.0, 7.94, 9.25, 9.97)
                  # Adjust according your data
# Extract data matrices
x <- data$x
y <- concentration
nmEX <- t(data.frame(seq(310, 430, by = 10)))  # Excitation wavelengths
nmEM <- t(data.frame(seq(300, 750, by = 1)))  # Emission wavelengths

# Check the structure of the data
print(x)

# Get matrix dimensions
mydim = dim(x)
print(mydim)

# Compute the average across columns
xm = colMeans(x)
print(xm)

# Plot contour of average matrix
dev.new()
filled.contour(xm, color.palette=terrain.colors())

# Plot intensity vs emission wavelengths
dev.new()
matplot(t(nmEM), xm, type="l", xlab="Emission Wavelength (nm)", ylab="Intensity")

# Plot intensity vs excitation wavelengths
dev.new()
matplot(t(nmEX), t(xm), type="l", xlab="Excitation Wavelength (nm)", ylab="Intensity")

# Set percentage for calibration and testing (70% calibration, 30% testing)
perc = 0.7
dim_x = dim(x)
size = 1:dim_x[1]
ntrain = ceiling(perc * dim_x[1])  # Number of calibration samples
sel = sample(size, ntrain)  # Random sample selection
xcal = x[sel, ,]  # Calibration data
ycal = matrix(y[sel])  # Calibration concentrations
xpred = x[-sel, ,]  # Testing data
ypred = matrix(y[-sel])  # Testing concentrations

# Cross-validation for nPLS model (can take a long time)
cv = cv_snpls(xcal, ycal, ncomp=1:3, keepJ=1:2, keepK=1:2, sample=10, parallel=FALSE)

# Set the number of components for the nPLS model
ncomp = 6

# Build the nPLS model
model_npls = sNPLS(xcal, ycal, ncomp=ncomp, keepJ=rep(2, ncomp), keepK=rep(1, ncomp))

# Predict values for calibration and testing sets
ycal_calc = predict(model_npls, xcal)
ypred_calc = predict(model_npls, xpred)

# Plot calibration results
dev.new()
plot(ycal, ycal_calc, pch='o', col="blue", xlab="Measured Concentration (g/g)", ylab="Predicted Concentration (g/g)", main="Calibration")
lines(ycal, ycal, col='red')

# Plot prediction results
dev.new()
plot(ypred, ypred_calc, pch=15, col="red", xlab="Measured Concentration (mg/L)", ylab="Predicted Concentration (mg/L)", main="Prediction")
lines(ypred, ypred, col='red')

# Combined calibration and prediction plot
dev.new()
plot(ycal, ycal_calc, pch='o', col="blue", xlab="Measured Concentration (g/g)", ylab="Predicted Concentration (g/g)", main="Calibration (blue) and Prediction (red)")
points(ypred, ypred_calc, pch=15, col="red")
lines(y, y, col="red")

# Calibration Metrics
MAE = mean(abs((ycal - ycal_calc))) * 100
R2cal = cor(ycal, ycal_calc)^2  # R-squared for calibration
RMSEC = rmse(ycal, ycal_calc)  # Root Mean Square Error of Calibration

# Prediction Metrics
MAE_pred = mean(abs((ypred - ypred_calc)))  # Mean Absolute Error for prediction
MAPEP = mean(abs((ypred - ypred_calc) / ypred)) * 100  # Mean Absolute Percentage Error
R2pred = cor(ypred, ypred_calc)^2  # R-squared for prediction
RMSEP = rmse(ypred, ypred_calc)  # Root Mean Square Error of Prediction

# Print calibration results
print("=========== Calibration =============")
print(paste("Mean Absolute Error Percentage (MAPE):", MAPEP))
print(paste("R-squared Calibration (R2cal):", R2cal))
print(paste("RMSEC (mg/L):", RMSEC))

# Print prediction results
print("=========== Test =============")
print(paste("Mean Absolute Error Percentage (MAPEP):", MAPEP))
print(paste("R-squared Prediction (R2pred):", R2pred))
print(paste("RMSEP (mg/L):", RMSEP))
