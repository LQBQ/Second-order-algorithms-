####
# SCRIPT [5]: TUCKER3-LDA
#
# OBJECTIVE:
# 1. Load processed data ('eem_cube_cleaned.mat').
# 2. Dynamically generate labels (Y) based on concentrations and names.
# 3. Perform robust class balancing (Same as PARAFAC script).
# 4. Build TUCKER-3 model and LDA classifier.
####

###############################################################
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

########################################################################
# --- STEP 2: Load Data Files ---
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

########################################################################
# --- STEP 3: Label Generation and Data Balancing ---
{
  # 1. Generate labels based on names and concentrations
  
  # Binary label (Adulterated vs. Pure) for LDA
  # (Translated to English for consistency in plots)
  Y_binario <- factor(ifelse(concentrations > 0, "Adulterated", "Pure"))
  
  # Multiclass label for visualization and selection
  # (Translated to English)
  Y_multiclasse <- factor(
    ifelse(grepl("^A", sample_names), "Brand A",
           ifelse(grepl("^B", sample_names), "Brand B",
                  ifelse(grepl("^C", sample_names), "Brand C",
                         ifelse(grepl("^D", sample_names), "Brand D", "Pure"))))
  )
  
  cat("Dynamically generated labels:\n")
  print(table(Y_binario, Y_multiclasse))
  
  # 2. Class Balancing (Logic identical to PARAFAC-LDA)
  
  set.seed(42) # For reproducibility
  
  indices_A <- which(Y_multiclasse == "Brand A")
  indices_B <- which(Y_multiclasse == "Brand B")
  indices_C <- which(Y_multiclasse == "Brand C")
  indices_D <- which(Y_multiclasse == "Brand D")
  indices_Puras <- which(Y_multiclasse == "Pure")
  
  # Sample 10 samples from EACH adulterant stratum (or the maximum available)
  sel_indices_A <- sample(indices_A, size = min(10, length(indices_A)), replace = FALSE)
  sel_indices_B <- sample(indices_B, size = min(10, length(indices_B)), replace = FALSE)
  sel_indices_C <- sample(indices_C, size = min(10, length(indices_C)), replace = FALSE)
  sel_indices_D <- sample(indices_D, size = min(10, length(indices_D)), replace = FALSE)
  
  # Combine indices (Total ~80 samples)
  indices_balanceados <- c(sel_indices_A, sel_indices_B, sel_indices_C, sel_indices_D, 
                           indices_Puras)
  
  # 4. Create the final balanced data cube
  X <- X_original[indices_balanceados, , ]
  
  # 5. Create final labels
  Y <- Y_binario[indices_balanceados]
  Y_multi_balanceado <- Y_multiclasse[indices_balanceados]
  nomes_amostras <- sample_names[indices_balanceados]
  
  cat(paste("\nFinal balanced data cube created with dimensions:", 
            paste(dim(X), collapse = " x "), "\n"))
  print(table(Y))
  
  # 6. Means for visualization
  X1 <- X[Y == "Adulterated", , ]
  X2 <- X[Y == "Pure", , ]
  
  X1m <- colMeans(X1)
  X2m <- colMeans(X2)
}

#################################################################
# --- STEP 4: Visualization ---

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

###########################################################################
# --- STEP 5: TUCKER-3 Model ---
{
  set.seed(24)
  
  # DEFINE RANKS HERE (Samples, Emission, Excitation)
  # Unlike PARAFAC (which has 1 rank), Tucker has 3.
  # Suggestion: c(3, 3, 3) or c(4, 4, 4) depending on complexity.
  ranks_tucker <- c(7, 7, 7) 
  
  cat(paste("Modeling TUCKER-3 with ranks:", paste(ranks_tucker, collapse=","), "...\n"))
  
  # output="all" returns A, B, C and the Core Tensor (GA)
  model <- tucker(X, nfac=ranks_tucker, nstart=75, maxit=5000, ctol=1e-6, output = c("best", "all"))
  cat(paste("Tucker-3 Model R²:", model$Rsq, "\n"))
}

# --- STEP 6: Outlier Diagnostics (RSS and Q vs T²) ---
{
  # --- Reconstruction error per sample (RSS / Q) ---
  # In Tucker: X_hat_i = B * (Sum_p A[i,p] * Core[p,,]) * t(C)
  
  A <- model$A
  B <- model$B
  C <- model$C
  G <- model$G # Core Tensor
  
  # Ensure G is a 3D array [RankA x RankB x RankC]
  if(is.matrix(G)) G <- array(G, dim = ranks_tucker)
  
  n_samples <- dim(X)[1]
  rss <- numeric(n_samples)
  
  cat("Calculating Tucker residuals (this might take a while)...\n")
  
  for(i in 1:n_samples) {
    # 1. Calculate the Core slice weighted by sample i scores
    # Core_Slice = Sum (A[i, p] * G[p, , ])
    G_slice <- matrix(0, nrow = ranks_tucker[2], ncol = ranks_tucker[3])
    
    for(p in 1:ranks_tucker[1]) {
      G_slice <- G_slice + A[i, p] * G[p, , ]
    }
    
    # 2. Reconstruct sample i matrix
    recon_i <- B %*% G_slice %*% t(C)
    
    # 3. Calculate residual
    resid_i <- X[i,,] - recon_i
    rss[i] <- sum(resid_i^2)
  }
  
  thr <- median(rss) + 3 * IQR(rss)
  outliers_rss <- which(rss > thr)
  
  dev.new()
  plot(rss, pch = 19, ylab = "RSS (Reconstruction Error)", xlab = "Sample Index",
       main = "Reconstruction Error per Sample (Tucker-3)")
  abline(h = thr, col = "red", lty = 2)
  text(outliers_rss, rss[outliers_rss], labels = nomes_amostras[outliers_rss], pos = 3, col = "red", cex=0.8)
  
  # --- Hotelling T² vs Residual Q ---
  # We use matrix A (Sample Scores)
  S_A <- cov(A)
  T2 <- mahalanobis(A, colMeans(A), S_A)
  
  nfac_A <- ranks_tucker[1] # Degrees of freedom for Sample mode
  alpha <- 0.05
  T2_lim <- nfac_A * (n_samples - 1) / (n_samples - nfac_A) * qf(1 - alpha, nfac_A, n_samples - nfac_A)
  
  Q <- rss
  Q_lim <- median(Q) + 3 * IQR(Q)
  
  dev.new()
  plot(T2, Q, pch = 19, col = ifelse(Y == "Adulterated", "red", "blue"),
       xlab = expression(Hotelling~T^2),
       ylab = expression(Q~"(Reconstruction Error)"),
       main = expression("Tucker-3: Q vs T"^2~"Plot"))
  abline(v = T2_lim, col = "darkgreen", lty = 2)
  abline(h = Q_lim,  col = "darkgreen", lty = 2)
  text(T2, Q, labels = nomes_amostras, pos = 3, cex = 0.6)
  legend("topright",
         legend = c("Adulterated", "Pure", "95% Limit"),
         col = c("red", "blue", "darkgreen"), pch = c(19,19,NA), lty = c(NA,NA,2))
  
  out_QT <- which(T2 > T2_lim | Q > Q_lim)
  cat("Samples identified as outliers (Q–T²):", nomes_amostras[out_QT], "\n")
}

# --- STEP 6.5: Manual Outlier Removal by NAME and Re-Modeling (TUCKER) ---
{
  # 1. DEFINE EXACT NAMES TO REMOVE HERE
  # Example: nomes_para_remover <- c("B0900", "D0900", "PA02", "PB02", "PD02", "PD03")
  nomes_para_remover <- c("A0900", "C0900", "D0700", "D0010", "PB03", "PC02", "PC03") 
  
  if (length(nomes_para_remover) > 0) {
    
    outliers_indices <- which(nomes_amostras %in% nomes_para_remover)
    
    nomes_encontrados <- nomes_amostras[outliers_indices]
    nomes_errados <- setdiff(nomes_para_remover, nomes_encontrados)
    
    if (length(nomes_errados) > 0) {
      cat("\n[ATTENTION] Names not found:", nomes_errados, "\n")
    }
    
    if (length(outliers_indices) > 0) {
      cat("\nRemoving", length(outliers_indices), "samples:", nomes_encontrados, "\n")
      
      # 2. REMOVE FROM CUBE AND METADATA
      X <- X[-outliers_indices, , ]
      Y <- Y[-outliers_indices]
      nomes_amostras <- nomes_amostras[-outliers_indices]
      
      cat("New cube dimensions:", paste(dim(X), collapse=" x "), "\n")
      
      # 3. RECALCULATE TUCKER-3 MODEL
      cat("Recalculating TUCKER-3 model without outliers...\n")
      set.seed(24) 
      
      # Re-run Tucker with the same ranks
      model <- tucker(X, nfac=ranks_tucker, nstart=75, maxit=5000, ctol=1e-6, output = c("best", "all"))
      
      # Update R²
      cat(paste("Tucker-3 Model R²:", model$Rsq, "\n"))
      
    } else {
      cat("No valid name found.\n")
    }
    
  } else {
    cat("Removal list is empty. The original model will be kept.\n")
  }
}

#############################################################################
# --- STEP 7: Scores and Loadings Visualization (Adapted for Tucker) ---
{
  ## Plot Tucker Scores (Mode 1 - Samples)
  dev.new()
  fator_x <- 1 
  fator_y <- 2 
  
  plot(model$A[, fator_x], model$A[, fator_y], 
       type = "n",
       xlab = paste("Factor", fator_x), 
       ylab = paste("Factor", fator_y), 
       main = "TUCKER-3 Scores (Matrix A)")
  
  text(x = model$A[Y == "Adulterated", fator_x], 
       y = model$A[Y == "Adulterated", fator_y], 
       labels = nomes_amostras[Y == "Adulterated"], 
       col = "red", cex = 0.7)
  
  text(x = model$A[Y == "Pure", fator_x], 
       y = model$A[Y == "Pure", fator_y], 
       labels = nomes_amostras[Y == "Pure"], 
       col = "blue", cex = 0.7)
  
  legend("topright", legend = c("Adulterated", "Pure"), col = c("red", "blue"), pch = 19)
  grid()
  
  ## Plot Tucker Loadings (Mode 2 - Emission)
  dev.new()
  matplot(nmEM, model$B, type ="l", xlab = "Emission Wavelength (nm)", 
          ylab = "Loadings", main = "Tucker Loadings - Emission (Matrix B)")
  
  ## Plot Tucker Loadings (Mode 3 - Excitation)
  dev.new()
  matplot(nmEX, model$C, type ="l", xlab = "Excitation Wavelength (nm)", 
          ylab = "Loadings", main = "Tucker Loadings - Excitation (Matrix C)")
}

###############################################################################
# --- STEP 8: Train/Test Split and LDA Model ---
{
  train_ratio <- 0.7 
  
  # Define numeric groups (1=Adulterated, 2=Pure)
  group_numeric <- ifelse(Y == "Adulterated", 1, 2)
  
  # Class Scores (Using Tucker Matrix A)
  # Unlike PARAFAC, here we have 'ranks_tucker[1]' columns
  scores_X1 <- model$A[group_numeric == 1, ] # Adulterated Scores
  scores_X2 <- model$A[group_numeric == 2, ] # Pure Scores
  
  # Calculate training samples per class
  ntrain1 <- ceiling(train_ratio * nrow(scores_X1))
  ntrain2 <- ceiling(train_ratio * nrow(scores_X2))
  
  # Select samples (Kennard-Stone)
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

#########################################################################
# --- STEP 9: Performance Metrics and LDA Visualization ---
{
  # Predictions
  pred_train <- predict(model_lda, train)
  pred_test <- predict(model_lda, test)
  
  # Metrics
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
  cat("\n--- LDA Model Metrics (Tucker-3) ---\n")
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
          xlab = "Group (1=Adulterated, 2=Pure)", ylab = "LD1", main = "Boxplot (Training) - Tucker LDA")
  
  dev.new()
  plot(pred_train$x[, 1], col = as.factor(group_train), pch = 19, 
       xlab = "Samples", ylab = "LD1", main = "Scatter Plot - Training")
  legend("topleft", legend = c("Adulterated (1)", "Pure (2)"), col = c(2, 1), pch = 19)
  
  dev.new()
  boxplot(pred_test$x[, 1] ~ group_test, col = c("red", "blue"), 
          xlab = "Group (1=Adulterated, 2=Pure)", ylab = "LD1", main = "Boxplot (Test) - Tucker LDA")
  
  dev.new()
  plot(pred_test$x[, 1], col = as.factor(group_test), pch = 19, 
       xlab = "Samples", ylab = "LD1", main = "Scatter Plot - Test")
  legend("topleft", legend = c("Adulterated (1)", "Pure (2)"), col = c(2, 1), pch = 19)
}

