####
# SCRIPT [4]: PARAFAC-LDA
#
# OBJECTIVE:
# 1. Load processed data ('eem_cube_cleaned.mat') and metadata ('eem_data.mat').
# 2. Dynamically generate labels (Y) based on concentrations and names.
# 3. (Optional) Perform robust class balancing using sample names.
# 4. Build PARAFAC model and LDA classifier.
####

################################################################
# --- STEP 1: Load Required Packages ---
{
  load_required_packages <- function(packages) {
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
  }

  required_packages <- c(
    "multiway", "ThreeWay", "R.matlab", "plot3D", "plotly", 
    "prospectr", "MASS", "caret", "rstudioapi", "abind", "parallel"
  )
  
  load_required_packages(required_packages)
}

#####################################################################
# --- STEP 2: Load Data Files
{
  cat("Please select the 'eem_cube_cleaned.mat' file (contains pre-processed data)...\n")
  mat_cube_path <- rstudioapi::selectFile(caption = "Select eem_cube_cleaned.mat",
                                          filter = "MAT Files (*.mat)")
  if (!nzchar(mat_cube_path)) stop("Selection cancelled.")
  
  # Load matrix and data
  data_cube <- readMat(mat_cube_path)
  names(data_cube) <- gsub("\\.mat$", "", names(data_cube))
  
  X_original <- data_cube$cube.cleaned
  sample_names <- unlist(data_cube$sample.names)
  concentrations <- as.vector(data_cube$concentrations)
  nmEM <- as.vector(data_cube$nm.emission)
  nmEX <- as.vector(data_cube$nm.excitation)
  
  cat("Data loaded successfully!\n")
  cat(paste("Cube Dimensions (Samples x Emission x Excitation): [", 
            paste(dim(X_original), collapse = " x "), "]\n"))
}

######################################################################
# --- STEP 3: Label Generation and Data Balancing ---
{
  # 1. Generate labels based on names and concentrations
  
  # Binary label (Adulterated vs. Pure) for LDA
  Y_binario <- factor(ifelse(concentrations > 0, "Adulterated", "Pure"))
  
  # Multiclass label for visualization and selection
  Y_multiclasse <- factor(
    ifelse(grepl("^A", sample_names), "Brand A",
           ifelse(grepl("^B", sample_names), "Brand B",
                  ifelse(grepl("^C", sample_names), "Brand C",
                         ifelse(grepl("^D", sample_names), "Brand D", "Pure"))))
  )
  
  cat("Dynamically generated labels:\n")
  print(table(Y_binario, Y_multiclasse))
  
  # 2. Class Balancing
  # Original goal: Create a balanced cube with:
  # 10 samples from Brand A, 10 from B, 10 from C, 10 from D (Total 40 Adulterated)
  # 40 Pure samples (Total 40 Pure)
  
  set.seed(42) # For sampling reproducibility
  
  # Find INDICES for each group (regardless of order)
  indices_A <- which(Y_multiclasse == "Brand A")
  indices_B <- which(Y_multiclasse == "Brand B")
  indices_C <- which(Y_multiclasse == "Brand C")
  indices_D <- which(Y_multiclasse == "Brand D")
  indices_Puras <- which(Y_multiclasse == "Pure")
  
  # 2. Sample 10 samples from EACH adulterant stratum
  sel_indices_A <- sample(indices_A, size = 10, replace = FALSE)
  sel_indices_B <- sample(indices_B, size = 10, replace = FALSE)
  sel_indices_C <- sample(indices_C, size = 10, replace = FALSE)
  sel_indices_D <- sample(indices_D, size = 10, replace = FALSE)
  
  # 3. Combine selected indices
  indices_balanceados <- c(sel_indices_A, sel_indices_B, sel_indices_C, sel_indices_D, 
                           indices_Puras)
  
  # 4. Create the final balanced data cube
  X <- X_original[indices_balanceados, , ]
  
  # 5. Create final labels corresponding to the balanced cube
  Y <- Y_binario[indices_balanceados]
  Y_multi_balanceado <- Y_multiclasse[indices_balanceados]
  nomes_amostras <- sample_names[indices_balanceados]
  
  cat(paste("\nFinal balanced data cube created with dimensions:", 
            paste(dim(X), collapse = " x "), "\n"))
  print(table(Y)) # Should show 40 Adulterated, 40 Pure
  
  # 6. Access samples by class (for plots)
  # (Updated to use the new English labels Y)
  X1 <- X[Y == "Adulterated", , ]
  X2 <- X[Y == "Pure", , ]
  
  X1m <- colMeans(X1)
  X2m <- colMeans(X2)
}

################################################################
# --- STEP 4: Visualization

  dev.new()
  filled.contour(x = nmEM, y = nmEX, z = X1m, 
                 color.palette = topo.colors, main = "Class 1 Mean (Adulterated)",
                 xlab = "Emission (nm)", ylab = "Excitation (nm)")
  
  dev.new()
  filled.contour(x = nmEM, y = nmEX, z = X2m, 
                 color.palette = terrain.colors, main = "Class 2 Mean (Pure)",
                 xlab = "Emission (nm)", ylab = "Excitation (nm)")
  
  dev.new()
  matplot(nmEM, X1m, type="l", xlab = "Emission Wavelength (nm)", ylab = "Intensity", 
          main = "Class 1 - Normal (Emission Mode)")
  
  dev.new()
  matplot(nmEX, t(X1m), type="l", xlab = "Excitation Wavelength (nm)", ylab = "Intensity", 
          main = "Class 1 - Normal (Excitation Mode)")

#########################################################################
# --- STEP 5: PARAFAC Model (No logic changes) ---
{
  set.seed(24)
  num_factors <- 7  # Define number of factors
  
  
  cat(paste("Modeling PARAFAC with", num_factors, "factors..."))
 
  model <- parafac(X, nfac=num_factors, nstart=75, maxit=5000, 
                   ctol=1e-6, parallel=F, output = c("best", "all"))
  
  cat(paste("PARAFAC Model R²:", model$Rsq, "\n"))
}

###################################################################
# --- STEP 6: Outlier Diagnostics (RSS and Q vs T²) ---
{
  # --- Reconstruction error per sample (RSS) ---
  nfac <- ncol(model$A)
  n_samples <- dim(X)[1]
  rss <- numeric(n_samples)
  
  for(i in 1:n_samples) {
    recon_i <- matrix(0, nrow = dim(X)[2], ncol = dim(X)[3]) # emission x excitation
    for(k in 1:nfac) {
      recon_i <- recon_i + model$A[i,k] * (outer(model$B[,k], model$C[,k]))
    }
    resid_i <- X[i,,] - recon_i
    rss[i] <- sum(resid_i^2)
  }
  
  thr <- median(rss) + 3 * IQR(rss)
  outliers_rss <- which(rss > thr)
  
  dev.new()
  plot(rss, pch = 19, ylab = "RSS (Reconstruction Error)", xlab = "Sample Index",
       main = "Reconstruction Error per Sample (PARAFAC)")
  abline(h = thr, col = "red", lty = 2)
  text(outliers_rss, rss[outliers_rss], labels = nomes_amostras[outliers_rss], pos = 3, col = "red", cex=0.8)
  cat("Outliers by RSS (threshold median+3*IQR):", outliers_rss, "\n")
  
  # --- Hotelling T² vs Residual Q ---
  A <- model$A
  S_A <- cov(A)
  T2 <- mahalanobis(A, colMeans(A), S_A)
  
  nfac <- ncol(A)
  n_amostras <- nrow(A)
  alpha <- 0.05
  T2_lim <- nfac * (n_amostras - 1) / (n_amostras - nfac) * qf(1 - alpha, nfac, n_amostras - nfac)
  
  Q <- rss # Residual Q is the same as RSS calculated above
  Q_lim <- median(Q) + 3 * IQR(Q)
  
  dev.new()
  plot(T2, Q, pch = 19, col = ifelse(Y == "Adulterated", "red", "blue"),
       xlab = expression(Hotelling~T^2),
       ylab = expression(Q~"(Reconstruction Error)"),
       main = expression("Q vs T"^2~"Plot (Outlier Diagnosis)"))
  abline(v = T2_lim, col = "darkgreen", lty = 2)
  abline(h = Q_lim,  col = "darkgreen", lty = 2)
  text(T2, Q, labels = nomes_amostras, pos = 3, cex = 0.6)
  legend("topright",
         legend = c("Adulterated", "Pure", "95% Limit"),
         col = c("red", "blue", "darkgreen"), pch = c(19,19,NA), lty = c(NA,NA,2))
  
  out_QT <- which(T2 > T2_lim | Q > Q_lim)
  cat("Samples identified as outliers (Q–T²):", nomes_amostras[out_QT], "\n")
}

  # --- STEP 6.5 (Optional): Outlier Removal by NAME and Re-Modeling ---
  {
    # 1. DEFINE EXACT NAMES TO REMOVE HERE
    # Check the plot in Step 6 (red/blue labels).
    # Example: If sample "P C 0 5" and "Brand A 0 2" are outliers:
    
    nomes_para_remover <- c("B0900", "D0900", "PB02", "PD02", "PD03") 
    # Usage example: nomes_para_remover <- c("P C 0 5", "P B 1 0")
    
    # --- Search Logic (Name -> Index Decoding) ---
    if (length(nomes_para_remover) > 0) {
      
      # Identifies which indices (rows) correspond to the typed names
      # The %in% operator checks if the sample name is in your removal list
      outliers_indices <- which(nomes_amostras %in% nomes_para_remover)
      
      # Safety check: Warn if any name was typed incorrectly
      nomes_encontrados <- nomes_amostras[outliers_indices]
      nomes_errados <- setdiff(nomes_para_remover, nomes_encontrados)
      
      if (length(nomes_errados) > 0) {
        cat("\n[ATTENTION] The following names were not found in the dataset (check spaces or capitalization):\n")
        print(nomes_errados)
      }
      
      if (length(outliers_indices) > 0) {
        cat("\nRemoving", length(outliers_indices), "identified samples:\n")
        print(nomes_encontrados)
        
        # 2. REMOVE FROM DATA CUBE (X)
        # Use discovered indices to slice the cube
        X <- X[-outliers_indices, , ]
        
        # 3. UPDATE AUXILIARY VECTORS
        # It is crucial to also remove from Y and sample names vector
        Y <- Y[-outliers_indices]
        # If balancing was used, other class variables might need updating,
        # but Y and nomes_amostras are essential for the code to proceed.
        nomes_amostras <- nomes_amostras[-outliers_indices]
        
        cat("New cube dimensions:", paste(dim(X), collapse=" x "), "\n")
        
        # 4. RECALCULATE PARAFAC MODEL
        # The model needs to be redone as Loadings change without outliers
        cat("Recalculating PARAFAC model without outliers...\n")
        set.seed(24) 
        
        model <- parafac(X, nfac=num_factors, nstart=75, maxit=5000, 
                         ctol=1e-6, parallel=FALSE, output = c("best", "all"))
        
        cat(paste("New Clean Model R²:", model$Rsq, "\n"))
        
      } else {
        cat("No valid name was found for removal.\n")
      }
      
    } else {
      cat("Removal list is empty. Original model will be kept.\n")
    }
  }

############################################################################
# --- STEP 7: Scores and Loadings Visualization ---
# (Logic unchanged, uses 'nmEM', 'nmEX', 'Y', 'nomes_amostras' defined earlier)
{
  ## Plot Parafac Scores 1x2
  dev.new()
  fator_x <- 1 
  fator_y <- 2 
  
  plot(model$A[, fator_x], model$A[, fator_y], 
       type = "n",
       xlab = paste("Factor", fator_x), 
       ylab = paste("Factor", fator_y), 
       main = "PARAFAC Scores with Sample Names")
  
  text(x = model$A[Y == "Adulterated", fator_x], 
       y = model$A[Y == "Adulterated", fator_y], 
       labels = nomes_amostras[Y == "Adulterated"], 
       col = "red", cex = 0.7)
  
  text(x = model$A[Y == "Pure", fator_x], 
       y = model$A[Y == "Pure", fator_y], 
       labels = nomes_amostras[Y == "Pure"], 
       col = "blue", cex = 0.7)
  
  legend("topright", 
         legend = c("Adulterated", "Pure"), 
         col = c("red", "blue"), 
         pch = 19)
  grid()
  
  ## Plot Parafac Loadings
  dev.new()
  matplot(nmEM, model$B, type ="l", xlab = "Emission Wavelength (nm)", ylab = "Intensity")
  
  dev.new()
  matplot(nmEX, model$C, type ="l", xlab = "Excitation Wavelength (nm)", ylab = "Intensity")
}

###############################################################################
# --- STEP 8: Train/Test Split and LDA Model ---
# (Logic updated to use English labels "Adulterated"/"Pure")
{
  train_ratio <- 0.7 
  
  # Define numeric groups (1=Adulterated, 2=Pure) for LDA
  # (Matches original script logic but uses English check)
  group_numeric <- ifelse(Y == "Adulterated", 1, 2)
  
  # Scores by class
  scores_X1 <- model$A[group_numeric == 1, ] # Adulterated Scores
  scores_X2 <- model$A[group_numeric == 2, ] # Pure Scores
  
  # Calculate training samples per class
  ntrain1 <- ceiling(train_ratio * nrow(scores_X1))
  ntrain2 <- ceiling(train_ratio * nrow(scores_X2))
  
  # Select samples (KenStone)
  sel1 <- kenStone(scores_X1, k = ntrain1)
  sel2 <- kenStone(scores_X2, k = ntrain2)
  
  train1 <- scores_X1[sel1$model, ]
  train2 <- scores_X2[sel2$model, ]
  train <- rbind(train1, train2)
  
  group_train <- factor(c(rep(1, ntrain1), rep(2, ntrain2)))
  
  test1 <- scores_X1[sel1$test, ]
  test2 <- scores_X2[sel2$test, ]
  test <- rbind(test1, test2)
  group_test <- factor(c(rep(1, nrow(test1)), rep(2, nrow(test2))))
  
  # LDA Model
  model_lda <- lda(train, group_train)
  model_lda_cv <- lda(train, group_train, CV=TRUE)
}

############################################################################
# --- STEP 9: Performance Metrics and LDA Visualization ---
# (Logic unchanged, uses 'group_train' and 'group_test')
{
  # Predictions
  pred_train <- predict(model_lda, train)
  pred_test <- predict(model_lda, test)
  
  # Metrics
  # (Using robust form from your original script [PARAFAC-LDA.R])
  ac_train <- mean(pred_train$class == group_train) 
  spec_train <- mean(pred_train$class[group_train == 1] == 1) # Class 1 (Adulterated)
  sens_train <- mean(pred_train$class[group_train == 2] == 2) # Class 2 (Pure)
  
  ac_cv <- mean(model_lda_cv$class == group_train)
  spec_cv <- mean(model_lda_cv$class[group_train == 1] == 1)
  sens_cv <- mean(model_lda_cv$class[group_train == 2] == 2)
  
  ac_test <- mean(pred_test$class == group_test) 
  spec_test <- mean(pred_test$class[group_test == 1] == 1)
  sens_test <- mean(pred_test$class[group_test == 2] == 2)
  
  # Print results
  cat("\n--- LDA Model Metrics ---\n")
  cat(" (Class 1 = Adulterated, Class 2 = Pure)\n")
  cat("\n--- Training ---\n")
  cat("Accuracy:", ac_train, "\nSensitivity (Pure):", sens_train, "\nSpecificity (Adulterated):", spec_train, "\n\n")
  cat("--- Cross-Validation (CV) ---\n")
  cat("Accuracy:", ac_cv, "\nSensitivity (Pure):", sens_cv, "\nSpecificity (Adulterated):", spec_cv, "\n\n")
  cat("--- Test ---\n")
  cat("Accuracy:", ac_test, "\nSensitivity (Pure):", sens_test, "\nSpecificity (Adulterated):", spec_test, "\n\n")
  
  # Confusion Matrices
  confusion_train <- table(Predicted = pred_train$class, Actual = group_train)
  confusion_test <- table(Predicted = pred_test$class, Actual = group_test)
  
  print("Confusion Matrix - Training")
  print(confusion_train)
  
  print("Confusion Matrix - Test")
  print(confusion_test)
  
  # Visualization
  dev.new()
  boxplot(pred_train$x[, 1] ~ group_train, col = c("red", "blue"), 
          xlab = "Group (1=Adulterated, 2=Pure)", ylab = "LD1", main = "Boxplot (Training)")
  
  dev.new()
  plot(pred_train$x[, 1], col = as.factor(group_train), pch = 19, 
       xlab = "Samples", ylab = "LD1", main = "Scatter Plot - Training")
  legend("topleft", legend = c("Adulterated (1)", "Pure (2)"), col = c(2, 1), pch = 19) # Col 2=red, 1=blue
  
  dev.new()
  boxplot(pred_test$x[, 1] ~ group_test, col = c("red", "blue"), 
          xlab = "Group (1=Adulterated, 2=Pure)", ylab = "LD1", main = "Boxplot (Test)")
  
  dev.new()
  plot(pred_test$x[, 1], col = as.factor(group_test), pch = 19, 
       xlab = "Samples", ylab = "LD1", main = "Scatter Plot - Test")
  legend("topleft", legend = c("Adulterated (1)", "Pure (2)"), col = c(2, 1), pch = 19)
}
