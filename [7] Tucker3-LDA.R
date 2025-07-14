#### NOTE: This script can only be used after the 3D matrix has been assembled (using "[5] Function to Create a 3D Array")
#### Note that the Tucker3-LDA matrix considers all samples in this case (pure and adulterated).

### Tucker3-LDA 
### (Tucker model - linear discriminant analysis)

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
  "prospectr", "MASS", "caret"
)

# Run the function to ensure all packages are installed and loaded
load_required_packages(required_packages)

# Define the input directory where the 3D matrix are attached
input_directory <- "path/to/your/directory"   # Adjust the path as needed
# Example: input_directory <- "C:/Users/Chemometric-PC/Documents/IC/Project - EVOO adulteration/Samples/Preprocessed samples"
# For compatibility, replace backslashes "\" with forward slashes "/" in your path.

## Load Data
data <- readMat("Sample_Data.mat") # Change the archive name to those that you chose 
X <- data$x # Access the data matrix in the .mat file 

nmEX <- t(data.frame(seq(310, 430, by = 10)))
nmEM <- t(data.frame(seq(301, 700, by = 1))) 
# Here, the emission number has changed due to the clipping of the excitation-emission matrix performed in preprocessing. 
# Change the values ​​according to your dataset, if necessary.
# Make sure that the emission number is in accordance with the second dimension of the 3D matrix.

# Total: 27 samples → 12 pure (class 1), 15 adulterated (class 2)
Y <- factor(c(rep(1, 12), rep(2, 15)))  # As a factor for modeling
# Again, if you working with other dataset, please, consider change the parameters 

# Accessing only class 1 samples (first 12)
X1 <- X[1:12,,]

# Accessing only class 2 samples (samples 13 to 27)
X2 <- X[13:27,,]

X1m <- colMeans(X1)
X2m <- colMeans(X2)

nmEX <- t(data.frame(seq(310, 430, by = 10)))
nmEM <- t(data.frame(seq(301, 700, by = 1)))
# Here, the emission number has changed due to the clipping of the excitation-emission matrix performed in preprocessing. 
# Change the values ​​according to your dataset, if necessary.
# Make sure that the emission number is in accordance with the second dimension of the 3D matrix.

## Visualizing Data
# Matrix dimensions
mydim <- dim(X) # samples x emission x excitation
print(mydim)

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

dev.new()
plot(pred_train$x[, 1], col = as.factor(group_train), pch = 19, 
     xlab = "Samples", ylab = "LD1", main = "Scatter Plot - Training")
legend("topleft", legend = c("Pure", "Adulterated"), col = c(1, 2), pch = 19)

dev.new()
plot(pred_test$x[, 1], col = as.factor(group_test), pch = 19, 
     xlab = "Samples", ylab = "LD1", main = "Scatter Plot - Test")
legend("topleft", legend = c("Pure", "Adulterated"), col = c(1, 2), pch = 19)

# Confusion Matrices
confusion_train <- table(Predicted = pred_train$class, Actual = group_train)
print("Confusion Matrix - Training Set")
print(confusion_train)

confusion_test <- table(Predicted = pred_test$class, Actual = group_test)
print("Confusion Matrix - Test Set")
print(confusion_test)

# This script performs a Tucker3 decomposition followed by Linear Discriminant Analysis (LDA) 
# for classification of fluorescence Excitation-Emission Matrices (EEMs).
#
# REQUIREMENTS:
# - The 3D matrix `X` (samples × emission × excitation) must be previously constructed and loaded.
# - In this case: 27 total samples → 12 pure (class 1), 15 adulterated (class 2).
#
# STEP-BY-STEP OVERVIEW:
# 1. PACKAGES:
#    - Automatically installs and loads the required libraries for tensor decomposition, visualization,
#      preprocessing, and classification.
#
# 2. DATA LOADING:
#    - Loads the fluorescence tensor `X` and defines excitation (nmEX) and emission (nmEM) axes.
#    - The emission range might be reduced due to preprocessing (e.g., Rayleigh scatter removal).
#    - Defines class labels `Y` as a factor (1 = pure, 2 = adulterated).
#
# 3. VISUALIZATION:
#    - Computes class-wise mean spectra (`X1m` and `X2m`) and visualizes them with `matplot` and `filled.contour`.
#
# 4. TUCKER3 MODEL:
#    - Tensor is unfolded into 2D for Tucker3 modeling using the `T3func()` function.
#    - Decomposition performed with equal number of components (`nf = 6`) across all three modes.
#    - Output matrices:
#       A = sample scores,
#       B = emission loadings,
#       C = excitation loadings,
#       G = core array (not used here).
#
# 5. SCORES AND CLASSIFICATION:
#    - Scores from `A` are split into class-wise matrices and subsets are selected for training/test 
#      using the KenStone algorithm (maximizing representativity).
#    - LDA is performed using the training set and evaluated on:
#       a) Training set
#       b) Test set
#       c) Cross-validation (leave-one-out)
#
# 6. METRICS:
#    - Accuracy, sensitivity (true positive rate), and specificity (true negative rate) are calculated 
#      for each evaluation strategy.
#    - Confusion matrices are printed for training and testing sets.
#
# 7. VISUALIZATION OF LDA RESULTS:
#    - Boxplots and scatter plots are used to explore the LD1 component across groups.
#
# NOTES:
# - This Tucker3-LDA pipeline assumes that the EEM dataset has been preprocessed, scaled, and 
#   shaped consistently.
# - The number of Tucker components (`nf`) should be optimized based on explained variance (fp), 
#   stability, or classification performance.
# - This implementation is suitable for detecting class differences in chemometric applications
#   such as adulteration detection, sample authentication, or quality control.
#
# TO ADAPT:
# - For other datasets, modify the class labels (`Y`), tensor dimensions, and wavelength axes 
#   accordingly.
