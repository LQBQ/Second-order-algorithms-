#### NOTE: This script can only be used after the 3D matrix has been assembled (using "[5] Function to Create a 3D Array")
#### Note that the PARAFAC-LDA matrix considers all samples in this case (pure and adulterated).

### PARAFAC-LDA 
### (Parallel factor analysis - linear discriminant analysis)

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
data <- readMat("Sample_Data.mat")
X <- data$x

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

# Define X and Y axis values
x_vals <- nmEM  
y_vals <- nmEX 

# Plot data
dev.new()
filled.contour(x = x_vals, y = y_vals, z = as.matrix(X1[1, , ]), 
               color.palette = topo.colors, main = "Class 1 Sample",
               xlab = "Emission (nm)", ylab = "Excitation (nm)")

dev.new()
filled.contour(X2m, color.palette = terrain.colors, main = "Class 2 Average")

dev.new()
matplot(t(nmEM), X1m, type="l", xlab = "Emission Wavelength (nm)", ylab = "Intensity", 
        main = "Class 1 - Normal")

dev.new()
matplot(t(nmEX), t(X1m), type="l", xlab = "Excitation Wavelength (nm)", ylab = "Intensity", 
        main = "Class 1 - Normal")

######## PARAFAC Model #########
set.seed(24)
num_factors <- 6  # Define number of factors

model <- parafac(X, nfac=num_factors, nstart=10, maxit=500, ctol=1e-5, parallel=FALSE)
# nstart, maxit and ctol need to be adjusted acoording your data. The values defined here are PARAFAC defaults

# R2 Adjustment
model$Rsq  

## Plot Parafac Scores 1x2
dev.new()
matplot(model$A[Y==1,1], model$A[Y==1,2], pch="o", col="blue", xlab="Factor 1", 
        ylab="Factor 2", main="Scores (o: Normal, x: Adulterated)")
points(model$A[Y==2,1], model$A[Y==2,2], pch="x", col="red")

## Plot Parafac Loadings
dev.new()
matplot(t(nmEM), model$B, type ="l", xlab = "Emission Wavelength (nm)", ylab = "Intensity")

dev.new()
matplot(t(nmEX), model$C, type ="l", xlab = "Excitation Wavelength (nm)", ylab = "Intensity")

######### Training and Testing Samples #########
train_ratio <- 0.7  # 70% for training
mutation_prob <- 0.2  # 20% mutation probability

dim_class1 <- dim(X1)
dim_class2 <- dim(X2)
dim_data <- dim(X)

# Calculate training samples per class
ntrain1 <- ceiling(train_ratio * dim_class1[1])
ntrain2 <- ceiling(train_ratio * dim_class2[1])

# Get scores
scores_X1 <- model$A[1:dim_class1[1], ]
scores_X2 <- model$A[(dim_class1[1] + 1):dim_data[1], ]

# Select training samples
sel1 <- kenStone(scores_X1, k = ntrain1)
sel2 <- kenStone(scores_X2, k = ntrain2)

train1 <- scores_X1[sel1$model, ]
train2 <- scores_X2[sel2$model, ]
train <- rbind(train1, train2)

group1 <- rep(1, dim_class1[1])
group2 <- rep(2, dim_class2[1])
group_train <- rbind(matrix(group1[sel1$model]), matrix(group2[sel2$model]))

test1 <- scores_X1[sel1$test, ]
test2 <- scores_X2[sel2$test, ]
test <- rbind(test1, test2)
group_test <- rbind(matrix(group1[sel1$test]), matrix(group2[sel2$test]))

# LDA Model
model_lda <- lda(train, group_train)
model_lda_cv <- lda(train, group_train, CV=TRUE)

# Predictions
pred_train <- predict(model_lda, train)
pred_test <- predict(model_lda, test)

# Performance Metrics
ac_train <- mean(pred_train$class == group_train) 
spec_train <- mean(pred_train$class[1:ntrain1] == group1[sel1$model]) 
sens_train <- mean(pred_train$class[(ntrain1 + 1):(ntrain1 + ntrain2)] == group2[sel2$model]) 

ac_cv <- mean(model_lda_cv$class == group_train)
spec_cv <- mean(model_lda_cv$class[1:ntrain1] == group1[sel1$model]) 
sens_cv <- mean(model_lda_cv$class[(ntrain1 + 1):(ntrain1 + ntrain2)] == group2[sel2$model]) 

ac_test <- mean(pred_test$class == group_test) 
spec_test <- mean(pred_test$class[1:ntrain1] == group1[sel1$test]) 
sens_test <- mean(pred_test$class[(ntrain1 + 1):(ntrain1 + ntrain2)] == group2[sel2$test]) 

cat("Training Accuracy:", ac_train, "\nSensitivity:", sens_train, "\nSpecificity:", spec_train, "\n")
cat("Test Accuracy:", ac_test, "\nSensitivity:", sens_test, "\nSpecificity:", spec_test, "\n")
cat("Cross-Validation Accuracy:", ac_cv, "\nSensitivity:", sens_cv, "\nSpecificity:", spec_cv, "\n")

# Confusion Matrices
confusion_train <- table(Predicted = pred_train$class, Actual = group_train)
confusion_test <- table(Predicted = pred_test$class, Actual = group_test)

print("Confusion Matrix - Training Set")
print(confusion_train)

print("Confusion Matrix - Test Set")
print(confusion_test)

# Visualization
dev.new()
boxplot(pred_train$x[, 1] ~ group_train, col = c("blue", "red"), 
        xlab = "Group", ylab = "LD1", main = "Boxplot (Training)")

dev.new()
plot(pred_train$x[, 1], col = as.factor(group_train), pch = 19, 
     xlab = "Samples", ylab = "LD1", main = "Scatter Plot - Training")
legend("topleft", legend = c("Pure", "Adulterated"), col = c(1, 2), pch = 19)

dev.new()
boxplot(pred_test$x[, 1] ~ group_test, col = c("blue", "red"), 
        xlab = "Group", ylab = "LD1", main = "Boxplot (Test)")

dev.new()
plot(pred_test$x[, 1], col = as.factor(group_test), pch = 19, 
     xlab = "Samples", ylab = "LD1", main = "Scatter Plot - Test")
legend("topleft", legend = c("Pure", "Adulterated"), col = c(1, 2), pch = 19)

# This script performs a classification analysis using the PARAFAC-LDA (Parallel Factor Analysis 
# combined with Linear Discriminant Analysis) approach, applied to Excitation-Emission Matrices (EEMs).
# 
# REQUIREMENTS:
# - The 3D data matrix 'X' must be pre-assembled and loaded in the workspace.
# - X must have dimensions [samples × emission × excitation], with all samples concatenated.
# - In this study: 12 pure samples (class 1) and 15 adulterated samples (class 2), total = 27.
#
# STEP-BY-STEP OVERVIEW:
# 1. Packages: Ensures all necessary packages are installed and loaded.
#
# 2. Data loading:
#    - Loads the 3D fluorescence matrix (`X`) and defines the excitation/emission wavelength axes.
#    - Defines class labels `Y` (1 for pure, 2 for adulterated).
#
# 3. Visualization:
#    - Visualizes raw and averaged fluorescence data for both classes.
#
# 4. PARAFAC modeling:
#    - Applies PARAFAC decomposition to the entire dataset (`X`).
#    - Parameters: `nfac = 6` (number of latent factors), with convergence and iteration control.
#    - Output matrices: 
#       A = scores (samples × components),
#       B = emission loadings (emission × components),
#       C = excitation loadings (excitation × components).
#
# 5. Classification preparation:
#    - Sample scores (from matrix A) are split by class.
#    - Training and test sets are selected using the KenStone algorithm for sample representativity.
#
# 6. LDA classification:
#    - Trains an LDA model on the selected training subset.
#    - Performs classification on training and test sets.
#    - Also performs cross-validation on the training set.
#    - Calculates performance metrics: accuracy, sensitivity, specificity.
#    - Constructs confusion matrices.
#
# 7. Visualization of classification results:
#    - Plots boxplots and scatter plots of LD1 scores for both training and test sets.
#
# NOTE:
# - The number of PARAFAC components (`nfac`), the number of parallel analysis starts (`nstart`), 
#   the number of maximum interactions (`maxit`), the convergence tolerance (`ctol`) and LDA parameters should be optimized based on 
#   explained variance (R²), core consistency, or validation performance.

