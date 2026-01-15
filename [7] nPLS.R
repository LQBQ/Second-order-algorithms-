####
# SCRIPT [9]: N-PLS (Regression)
#
# OBJECTIVE:
# 1. Load processed data ('eem_cube_cleaned.mat').
# 2. Balanced selection (1 Pure per brand + Adulterated).
# 3. N-PLS Modeling (sNPLS) and Outlier Diagnostics.
# 4. Train/Test Split via MLM (Kennard-Stone + Mutation).
# 5. Final N-PLS Regression with Cross-Validation.
####

##################################################################
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
    "sNPLS", "multiway", "ThreeWay", "R.matlab", "plot3D", "plotly", 
    "prospectr", "MASS", "caret", "rstudioapi", "abind", "parallel", "Metrics"
  )
  
  load_required_packages(required_packages)
}

#####################################################################
# --- STEP 2: Load Data Files ---
{
  cat("Please select the 'eem_cube_cleaned.mat' file...\n")
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
  cat(paste("Cube Dimensions: [", paste(dim(X_original), collapse = " x "), "]\n"))
}

#######################################################################
# --- STEP 3: Data Preparation (Balanced Selection) ---
{
  # Define brands (Translated to English for consistent logic)
  Y_multiclasse <- factor(
    ifelse(grepl("^A", sample_names), "Brand A",
           ifelse(grepl("^B", sample_names), "Brand B",
                  ifelse(grepl("^C", sample_names), "Brand C",
                         ifelse(grepl("^D", sample_names), "Brand D", "Pure")))
    )
  )
    
    set.seed(42) 
    
    # --- ADULTERATED SELECTION ---
    # Updated to match "Brand"
    indices_A <- which(Y_multiclasse == "Brand A")
    indices_B <- which(Y_multiclasse == "Brand B")
    indices_C <- which(Y_multiclasse == "Brand C")
    indices_D <- which(Y_multiclasse == "Brand D")
    
    sel_indices_A <- sample(indices_A, size = min(10, length(indices_A)), replace = FALSE)
    sel_indices_B <- sample(indices_B, size = min(10, length(indices_B)), replace = FALSE)
    sel_indices_C <- sample(indices_C, size = min(10, length(indices_C)), replace = FALSE)
    sel_indices_D <- sample(indices_D, size = min(10, length(indices_D)), replace = FALSE)
    
    # --- PURE SELECTION (1 from each brand) ---
    indices_todos_zeros <- which(concentrations == 0)
    nomes_zeros <- sample_names[indices_todos_zeros]
    
    # Try to identify brand in pure name (e.g., "P A 0 1")
    idx_pure_A <- indices_todos_zeros[grepl("A", nomes_zeros)]
    idx_pure_B <- indices_todos_zeros[grepl("B", nomes_zeros)]
    idx_pure_C <- indices_todos_zeros[grepl("C", nomes_zeros)]
    idx_pure_D <- indices_todos_zeros[grepl("D", nomes_zeros)]
    
    # Draw 1 from each if it exists
    sel_pure_A <- if(length(idx_pure_A) > 0) sample(idx_pure_A, 1) else NULL
    sel_pure_B <- if(length(idx_pure_B) > 0) sample(idx_pure_B, 1) else NULL
    sel_pure_C <- if(length(idx_pure_C) > 0) sample(idx_pure_C, 1) else NULL
    sel_pure_D <- if(length(idx_pure_D) > 0) sample(idx_pure_D, 1) else NULL
    
    indices_Puras_Final <- c(sel_pure_A, sel_pure_B, sel_pure_C, sel_pure_D)
    
    cat("Selected Pure Sample Indices:", indices_Puras_Final, "\n")
    
    # --- FINAL COMBINATION ---
    indices_balanceados <- c(sel_indices_A, sel_indices_B, sel_indices_C, sel_indices_D, 
                             indices_Puras_Final)
    
    X <- X_original[indices_balanceados, , ]
    Y_conc <- concentrations[indices_balanceados] 
    nomes_amostras <- sample_names[indices_balanceados]
    
    cat(paste("\nDataset prepared.\n"))
    cat("X Dimensions:", paste(dim(X), collapse=" x "), "\n")
}

####################################################################################
# --- STEP 4: Preliminary N-PLS Model ---
{
  set.seed(24)
  num_factors <- 9 
  
  cat(paste("Modeling N-PLS with", num_factors, "latent variables...\n"))
  
  # sNPLS is the package function. If threshold=0, acts as standard N-PLS (no variable selection)
  model <- sNPLS(X, as.matrix(Y_conc), ncomp = num_factors, threshold_j = 0, threshold_k = 0)
  
  # Extracting fitted values
  fitted_vals <- as.vector(fitted(model))
  r2_fit <- cor(Y_conc, fitted_vals)^2
  
  cat(paste("N-PLS Model R² (Prelim Fit):", round(r2_fit, 5), "\n"))
}

###################################################################################
# --- STEP 5: Complete Diagnostics and Outlier Removal ---
{
  # 1. HOTELLING T² CALCULATION (Score Space)
  # N-PLS generates scores (T) describing variation in X correlated with Y
  Scores <- model$T
  n_samples <- dim(X)[1]
  
  S_T <- cov(Scores)
  # Mahalanobis Distance on Scores
  T2 <- mahalanobis(Scores, colMeans(Scores), S_T)
  
  # Statistical Limit for T² (F-Distribution)
  alpha <- 0.05
  T2_lim <- num_factors * (n_samples - 1) / (n_samples - num_factors) * qf(1 - alpha, num_factors, n_samples - num_factors)
  
  # 2. PREDICTION ERROR CALCULATION (SPE-Y / Q-stat)
  # Unlike PARAFAC (which looks at X error), in N-PLS we look at how much it missed Y
  fitted_vals <- as.vector(fitted(model))
  residuals_y <- (Y_conc - fitted_vals)^2 # Squared prediction error
  Q <- residuals_y 
  
  # Limit for Q (Boxplot rule: Median + 3*IQR)
  Q_lim <- median(Q) + 3 * IQR(Q)
  outliers_Q <- which(Q > Q_lim)
  
  # Plot colors (Red = Adulterated, Blue = Pure)
  cores_plot <- ifelse(Y_conc > 0, "red", "blue")
  
  # PLOT 1: Prediction Error per Sample
  dev.new()
  plot(Q, pch = 19, col = cores_plot,
       ylab = expression(SPE[Y]~(Squared~Prediction~Error)), 
       xlab = "Sample Index",
       main = "Diagnostic 1: Fit Error (Y)")
  abline(h = Q_lim, col = "darkgreen", lty = 2)
  
  # Label only error outliers
  if(length(outliers_Q) > 0){
    text(outliers_Q, Q[outliers_Q], labels = nomes_amostras[outliers_Q], 
         pos = 3, cex=0.7, col="black")
  }
  cat("Potential outliers by Prediction Error (SPE-Y):", nomes_amostras[outliers_Q], "\n")
  
  # PLOT 2: Q vs T² (Influence Plot)
  dev.new()
  plot(T2, Q, pch = 19, col = cores_plot,
       xlab = expression(Hotelling~T^2~(Leverage)),
       ylab = expression(SPE[Y]~(Residuals)),
       main = expression("Diagnostic 2: N-PLS Outliers (Q vs T"^2*")"))
  
  abline(v = T2_lim, col = "darkgreen", lty = 2)
  abline(h = Q_lim,  col = "darkgreen", lty = 2)
  
  # Labels to identify samples on the plot
  text(T2, Q, labels = nomes_amostras, pos = 3, cex = 0.6, col = "black")
  
  legend("topright",
         legend = c("Adulterated (>0%)", "Pure (0%)", "95% Limits"),
         col = c("red", "blue", "darkgreen"), pch = c(19,19,NA), lty = c(NA,NA,2))
  
  out_QT <- which(T2 > T2_lim | Q > Q_lim)
  cat("Samples outside statistical limits (Q or T²):", nomes_amostras[out_QT], "\n")
  
}

{ 
 # 3. MANUAL REMOVAL AND RE-MODELING
  # ==========================================================
  # INSERT NAMES TO REMOVE HERE
  # Example: nomes_para_remover <- c("B0900", "PD02")
  
  nomes_para_remover <- c("B0800", "D0800") 
  
  # ==========================================================
  
  if (length(nomes_para_remover) > 0) {
    outliers_indices <- which(nomes_amostras %in% nomes_para_remover)
    
    # Check for not found names
    nomes_encontrados <- nomes_amostras[outliers_indices]
    nomes_errados <- setdiff(nomes_para_remover, nomes_encontrados)
    
    if(length(nomes_errados) > 0) {
      cat("\n[WARNING] Names not found:", nomes_errados, "\n")
    }
    
    if (length(outliers_indices) > 0) {
      cat("\nRemoving", length(outliers_indices), "samples and recalculating N-PLS...\n")
      
      # Remove from Cube and Vectors
      X <- X[-outliers_indices, , ]
      Y_conc <- Y_conc[-outliers_indices]
      nomes_amostras <- nomes_amostras[-outliers_indices]
      
      cat("New X Dimensions:", paste(dim(X), collapse=" x "), "\n")
      
      # Recalculate N-PLS model
      set.seed(24) 
      model <- sNPLS(X, as.matrix(Y_conc), ncomp = num_factors, threshold_j = 0, threshold_k = 0)
      
      # Check new fit
      fitted_vals <- as.vector(fitted(model))
      novo_r2 <- cor(Y_conc, fitted_vals)^2
      cat(paste("New R² (Clean Model):", round(novo_r2, 5), "\n"))
      
    } 
  } else {
    cat("\nNo samples listed for removal. Keeping original model.\n")
  }
}

###################################################################################
# --- STEP 6: Weight Visualization ---
{
  # N-PLS uses weights Wj (Emission) and Wk (Excitation)
  # The sNPLS object stores this in model$Wj and model$Wk
  
  dev.new()
  # Plotting weights for the first and second latent variable (LV1 and LV2)
  matplot(nmEM, model$Wj[, 1:2], type ="l", lty=1, col=c("black", "red"),
          xlab="Emission (nm)", ylab="Weights (Wj)", main="N-PLS Weights - Emission")
  legend("topright", legend=c("LV 1", "LV 2"), col=c("black", "red"), lty=1)
  
  dev.new()
  matplot(nmEX, model$Wk[, 1:2], type ="l", lty=1, col=c("black", "red"),
          xlab="Excitation (nm)", ylab="Weights (Wk)", main="N-PLS Weights - Excitation")
  legend("topright", legend=c("LV 1", "LV 2"), col=c("black", "red"), lty=1)
}

################################################################################
# --- STEP 7: Train/Test Split via MLM (Mutation) ---
{
  # 1. Initial Configuration and Kennard-Stone
  # For KS in 3D data, we unfold cube X into a matrix (Samples x Var)
  X_flat <- matrix(as.vector(X), nrow = dim(X)[1]) 
  
  perc <- 0.7
  ntrain <- ceiling(perc * nrow(X_flat))
  
  # Kennard-Stone (Euclidean on unfolded X)
  set.seed(2) 
  ks_result <- kenStone(X_flat, k = ntrain, metric = "euclid")
  
  calibration_idx <- ks_result$model
  prediction_idx <- ks_result$test
  
  cat("Initial KS Split: Train =", length(calibration_idx), "| Test =", length(prediction_idx), "\n")
  
  # 2. Mutation Logic (MLM)
  cat("Applying Mutation (Swapping 20% of samples)...\n")
  set.seed(1)
  
  prob <- 0.2
  p_cal <- ceiling(prob * length(calibration_idx))
  p_pred <- ceiling(prob * length(prediction_idx))
  
  # Relative indices for swap
  cal_mut_pos <- sample(1:length(calibration_idx), size = p_cal)
  pred_mut_pos <- sample(1:length(prediction_idx), size = p_pred)
  
  # Swapping real indices
  idx_to_pred <- calibration_idx[cal_mut_pos]
  idx_to_cal  <- prediction_idx[pred_mut_pos]
  
  final_cal_idx <- c(calibration_idx[-cal_mut_pos], idx_to_cal)
  final_pred_idx <- c(prediction_idx[-pred_mut_pos], idx_to_pred)
  
  # Creating Final Tensors
  # drop=FALSE maintains 3D structure even if there is only 1 sample
  X_cal <- X[final_cal_idx, , , drop=FALSE]
  Y_cal <- Y_conc[final_cal_idx]
  
  X_test <- X[final_pred_idx, , , drop=FALSE]
  Y_test <- Y_conc[final_pred_idx]
  
  cat("Final MLM Split: Train =", length(final_cal_idx), "| Test =", length(final_pred_idx), "\n")
}

#################################################################################
# --- STEP 8: Final N-PLS Regression with Cross-Validation ---
{
  set.seed(123)
  
  # 1. Final Calibration Model (sNPLS)
  final_model <- sNPLS(X_cal, as.matrix(Y_cal), ncomp = num_factors, threshold_j=0, threshold_k=0)
  
  # 2. Calibration Metrics (Fitted)
  pred_cal <- as.vector(fitted(final_model))
  
  rmse_cal <- rmse(Y_cal, pred_cal)
  mae_cal  <- mae(Y_cal, pred_cal)
  r2_cal   <- cor(Y_cal, pred_cal)^2
  
  # 3. Cross-Validation (10-Fold) - Manual Loop
  # (sNPLS does not have direct integration with 'caret' train, we do it manually for robustness)
  cat("Running Cross-Validation (10-Fold)...\n")
  k = 10
  folds <- createFolds(Y_cal, k = k, list = TRUE)
  
  cv_preds_vector <- numeric(length(Y_cal))
  # Store results in temporary dataframe to sort later
  temp_res <- data.frame(Index = integer(), Pred = numeric())
  
  for(i in 1:k) {
    idx_val <- folds[[i]]
    idx_train <- setdiff(1:length(Y_cal), idx_val)
    
    # Fold Subsets
    xc_fold <- X_cal[idx_train, , , drop=FALSE]
    yc_fold <- Y_cal[idx_train]
    xv_fold <- X_cal[idx_val, , , drop=FALSE]
    
    # CV Model
    mod_cv <- sNPLS(xc_fold, as.matrix(yc_fold), ncomp = num_factors, threshold_j=0, threshold_k=0)
    
    # CV Prediction
    # ATTENTION: predict in sNPLS requires newX with correct dimension
    pred_fold <- as.vector(predict(mod_cv, newX = xv_fold))
    
    temp_res <- rbind(temp_res, data.frame(Index = idx_val, Pred = pred_fold))
  }
  
  # Reorder predictions to align with Y_cal
  temp_res <- temp_res[order(temp_res$Index), ]
  cv_preds <- temp_res$Pred
  
  rmse_cv <- rmse(Y_cal, cv_preds)
  mae_cv  <- mae(Y_cal, cv_preds)
  r2_cv   <- cor(Y_cal, cv_preds)^2
  
  # 4. Prediction on External Test
  pred_test <- as.vector(predict(final_model, newX = X_test))
  
  rmse_test <- rmse(Y_test, pred_test)
  mae_test  <- mae(Y_test, pred_test)
  r2_test   <- cor(Y_test, pred_test)^2
  
  # --- Display Results ---
  cat("\n=============================================\n")
  cat(" FINAL RESULTS (N-PLS with MLM)\n")
  cat("=============================================\n")
  cat(sprintf("CALIBRATION (Mutated Train):\n   RMSE: %.4f | MAE: %.4f | R²: %.4f\n", rmse_cal, mae_cal, r2_cal))
  cat("---------------------------------------------\n")
  cat(sprintf("CROSS-VALIDATION (10-Fold):\n   RMSE: %.4f | MAE: %.4f | R²: %.4f\n", rmse_cv, mae_cv, r2_cv))
  cat("---------------------------------------------\n")
  cat(sprintf("PREDICTION (Mutated Test):\n   RMSE: %.4f | MAE: %.4f | R²: %.4f\n", rmse_test, mae_test, r2_test))
  cat("=============================================\n")
}

###########################################################################################
# --- STEP 9: Separate Plots ---
{
  lims <- range(c(0, 100, Y_conc))
  
  # Plot 1: Calibration
  dev.new()
  plot(Y_cal, pred_cal, pch=19, col="blue", xlim=lims, ylim=lims,
       xlab="Actual Concentration (%)", ylab="Predicted Concentration (%)",
       main="CALIBRATION N-PLS (MLM)")
  abline(0, 1, lty=2, lwd=2, col="gray")
  grid()
  legend("topleft", legend=paste0("R² Cal: ", round(r2_cal, 4)), bty="n")
  
  # Plot 2: Prediction
  dev.new()
  plot(Y_test, pred_test, pch=19, col="red", xlim=lims, ylim=lims,
       xlab="Actual Concentration (%)", ylab="Predicted Concentration (%)",
       main="PREDICTION N-PLS (MLM)")
  abline(0, 1, lty=2, lwd=2, col="gray")
  grid()
  legend("topleft", legend=paste0("R² Test: ", round(r2_test, 4)), bty="n", text.col="red")
  
  # Plot 3: Residuals
  dev.new()
  plot(Y_test, (Y_test - pred_test), pch=19, col="darkgreen",
       xlab="Actual Concentration (%)", ylab="Residual (Actual - Predicted)",
       main="Prediction Residuals (N-PLS)")
  abline(h=0, lty=2, col="black")
  grid()
}
