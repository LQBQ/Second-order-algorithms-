####
# SCRIPT [6]: PARAFAC-MLR (With MLM - Mutation Logic)
#
# OBJECTIVE:
# 1. Load processed data.
# 2. Balanced selection (1 Pure per brand + Adulterated).
# 3. PARAFAC Modeling and Outlier Removal.
# 4. Train/Test Split via MLM (Kennard-Stone + Mutation).
# 5. MLR Regression with Cross-Validation.
####

##############################################################
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
    "prospectr", "MASS", "caret", "rstudioapi", "abind", "parallel", "Metrics"
  )
  
  load_required_packages(required_packages)
}

########################################################################
# --- STEP 2: Load Data Files ---
{
  cat("Please select the 'eem_cube_cleaned.mat' file...\n")
  mat_cube_path <- rstudioapi::selectFile(caption = "Select eem_cube_cleaned.mat",
                                          filter = "MAT Files (*.mat)")
  if (!nzchar(mat_cube_path)) stop("Selection cancelled.")
  
  data_cube <- readMat(mat_cube_path)
  names(data_cube) <- gsub("\\.mat$", "", names(data_cube))
  
  X_original <- data_cube$cube.cleaned
  sample_names <- unlist(data_cube$sample.names)
  concentrations <- as.vector(data_cube$concentrations) 
  nmEM <- as.vector(data_cube$nm.emission)
  nmEX <- as.vector(data_cube$nm.excitation)
  
  cat("Data loaded successfully!\n")
}

#######################################################################
# --- STEP 3: Data Preparation (Balanced Selection) ---
{
  # Create brand labels
  Y_multiclasse <- factor(
    ifelse(grepl("^A", sample_names), "Brand A",
           ifelse(grepl("^B", sample_names), "Brand B",
                  ifelse(grepl("^C", sample_names), "Brand C",
                         ifelse(grepl("^D", sample_names), "Brand D", "Pure")))
    )
  )
    
    set.seed(42) 
    
    # Selection of Adulterated samples (Max available or 10)
    # Updated logic to match "Brand" instead of "Marca"
    sel_indices_A <- sample(which(Y_multiclasse == "Brand A"), size = min(10, sum(Y_multiclasse == "Brand A")))
    sel_indices_B <- sample(which(Y_multiclasse == "Brand B"), size = min(10, sum(Y_multiclasse == "Brand B")))
    sel_indices_C <- sample(which(Y_multiclasse == "Brand C"), size = min(10, sum(Y_multiclasse == "Brand C")))
    sel_indices_D <- sample(which(Y_multiclasse == "Brand D"), size = min(10, sum(Y_multiclasse == "Brand D")))
    
    # Selection of Pure samples (1 from each brand)
    indices_todos_zeros <- which(concentrations == 0)
    nomes_zeros <- sample_names[indices_todos_zeros]
    
    sel_pure_A <- if(any(grepl("A", nomes_zeros))) sample(indices_todos_zeros[grepl("A", nomes_zeros)], 1) else NULL
    sel_pure_B <- if(any(grepl("B", nomes_zeros))) sample(indices_todos_zeros[grepl("B", nomes_zeros)], 1) else NULL
    sel_pure_C <- if(any(grepl("C", nomes_zeros))) sample(indices_todos_zeros[grepl("C", nomes_zeros)], 1) else NULL
    sel_pure_D <- if(any(grepl("D", nomes_zeros))) sample(indices_todos_zeros[grepl("D", nomes_zeros)], 1) else NULL
    
    indices_balanceados <- c(sel_indices_A, sel_indices_B, sel_indices_C, sel_indices_D, 
                             sel_pure_A, sel_pure_B, sel_pure_C, sel_pure_D)
    
    X <- X_original[indices_balanceados, , ]
    Y_conc <- concentrations[indices_balanceados] 
    nomes_amostras <- sample_names[indices_balanceados]
    
    cat(paste("\nDataset prepared. X Dimensions:", paste(dim(X), collapse=" x "), "\n"))
}

##################################################################################
# --- STEP 4: PARAFAC Modeling ---
{
  set.seed(24)
  num_factors <- 8 
  
  cat(paste("Modeling PARAFAC with", num_factors, "factors...\n"))
  model <- parafac(X, nfac=num_factors, nstart=75, maxit=5000, 
                   ctol=1e-6, parallel=FALSE, output = c("best", "all"))
  cat(paste("PARAFAC Model R² (Fit):", round(model$Rsq, 5), "\n"))
}

###################################################################################
# --- STEP 5: Outlier Diagnostics (RSS and Q vs T²) ---
  # 1. RECONSTRUCTION ERROR CALCULATION (RSS / Q-statistic)
{ nfac <- ncol(model$A)
  n_samples <- dim(X)[1]
  rss <- numeric(n_samples)
  
  # Reconstruct each sample to calculate the residual
  for(i in 1:n_samples) {
    recon_i <- matrix(0, nrow = dim(X)[2], ncol = dim(X)[3]) 
    for(k in 1:nfac) {
      recon_i <- recon_i + model$A[i,k] * (outer(model$B[,k], model$C[,k]))
    }
    rss[i] <- sum((X[i,,] - recon_i)^2)
  }
  
  # Define threshold for RSS (Q)
  thr_rss <- median(rss) + 3 * IQR(rss)
  outliers_rss <- which(rss > thr_rss)
  
  # PLOT 1: RSS per Sample
  dev.new()
  plot(rss, pch = 19, ylab = "RSS (Reconstruction Error)", xlab = "Sample Index",
       main = "Diagnostic 1: Reconstruction Error (RSS)")
  abline(h = thr_rss, col = "red", lty = 2)
  if(length(outliers_rss) > 0){
    text(outliers_rss, rss[outliers_rss], labels = nomes_amostras[outliers_rss], 
         pos = 3, col = "red", cex=0.8)
  }
  cat("Potential outliers by RSS:", nomes_amostras[outliers_rss], "\n")
  
  # 2. HOTELLING T² CALCULATION (Score Space)
  A <- model$A
  S_A <- cov(A)
  T2 <- mahalanobis(A, colMeans(A), S_A)
  
  # T² Limit (F-Distribution)
  alpha <- 0.05
  T2_lim <- nfac * (n_samples - 1) / (n_samples - nfac) * qf(1 - alpha, nfac, n_samples - nfac)
  
  # Q Limit (Recalculated for visual consistency)
  Q <- rss
  Q_lim <- median(Q) + 3 * IQR(Q)
  
  # PLOT 2: Q vs T² (Influence Plot)
  # Define colors: Red if Concentration > 0 (Adulterated), Blue if 0 (Pure)
  cores_plot <- ifelse(Y_conc > 0, "red", "blue")
  
  dev.new()
  plot(T2, Q, pch = 19, col = cores_plot,
       xlab = expression(Hotelling~T^2),
       ylab = expression(Q~"(Residuals)"),
       main = expression("Diagnostic 2: Q vs T"^2))
  
  abline(v = T2_lim, col = "darkgreen", lty = 2)
  abline(h = Q_lim,  col = "darkgreen", lty = 2)
  
  # Add labels to all samples (small) or just outliers
  # Here we label everything for easier visual identification
  text(T2, Q, labels = nomes_amostras, pos = 3, cex = 0.6, col = "black")
  
  legend("topright",
         legend = c("Adulterated (>0%)", "Pure (0%)", "95% Limits"),
         col = c("red", "blue", "darkgreen"), pch = c(19,19,NA), lty = c(NA,NA,2))
  
  out_QT <- which(T2 > T2_lim | Q > Q_lim)
  cat("Samples outside limits (Q or T²):", nomes_amostras[out_QT], "\n")
}

{
  # 3. MANUAL REMOVAL AND RE-MODELING
  # ==========================================================
  # INSERT NAMES TO REMOVE HERE (Based on the plots above)
  # Example: nomes_para_remover <- c("B0900", "PD02")
  
  nomes_para_remover <- c("A0900", "A0800") 
  
  # ==========================================================
  
  if (length(nomes_para_remover) > 0) {
    # Search indices based on names
    outliers_indices <- which(nomes_amostras %in% nomes_para_remover)
    
    # Check for not found names
    nomes_encontrados <- nomes_amostras[outliers_indices]
    nomes_errados <- setdiff(nomes_para_remover, nomes_encontrados)
    
    if(length(nomes_errados) > 0) {
      cat("\n[WARNING] Names not found:", nomes_errados, "\n")
    }
    
    if (length(outliers_indices) > 0) {
      cat("\nRemoving", length(outliers_indices), "samples and recalculating the model...\n")
      
      # Remove from Cube and Auxiliary Vectors
      X <- X[-outliers_indices, , ]
      Y_conc <- Y_conc[-outliers_indices]
      nomes_amostras <- nomes_amostras[-outliers_indices]
      
      cat("New X Dimensions:", paste(dim(X), collapse=" x "), "\n")
      
      # Recalculate PARAFAC
      set.seed(24) 
      model <- parafac(X, nfac=num_factors, nstart=75, maxit=5000, 
                       ctol=1e-6, parallel=FALSE, output = c("best", "all"))
      
      cat(paste("New R² (Clean Model):", round(model$Rsq, 5), "\n"))
    } 
  } else {
    cat("\nNo samples listed for removal. Keeping original model.\n")
  }
}

############################################################################
# --- STEP 6: Loadings Visualization ---
{
  dev.new()
  par(mfrow=c(1,2))
  matplot(nmEM, model$B, type ="l", xlab="Emission (nm)", ylab="Loadings", main="Emission Mode")
  matplot(nmEX, model$C, type ="l", xlab="Excitation (nm)", ylab="Loadings", main="Excitation Mode")
  par(mfrow=c(1,1))
}

###############################################################################
# --- STEP 7: Train/Test Split via MLM (Mutation) ---
{
  # 1. Initial Configuration and Kennard-Stone
  perc <- 0.7
  A <- model$A
  y <- Y_conc
  ntrain <- ceiling(perc * nrow(A))
  
  combined_data <- cbind(A, y)
  
  # Initial Kennard-Stone (Euclidean)
  set.seed(2) 
  ks_result <- kenStone(combined_data, k = ntrain, metric = "euclid")
  
  calibration_idx <- ks_result$model
  prediction_idx <- ks_result$test
  
  xcal <- A[calibration_idx, ]
  ycal <- y[calibration_idx]
  
  xpred <- A[prediction_idx, ]
  ypred <- y[prediction_idx]
  
  cat("Initial KS Split: Train =", nrow(xcal), "| Test =", nrow(xpred), "\n")
  
  # 2. Mutation Logic (MLM)
  cat("Applying Mutation (Swapping 20% of samples)...\n")
  set.seed(2) # Seed to ensure mutation reproducibility
  
  prob <- 0.2
  p_cal <- ceiling(prob * nrow(xcal))
  p_pred <- ceiling(prob * nrow(xpred))
  
  cal_mut <- sample(1:nrow(xcal), size = p_cal)
  pred_mut <- sample(1:nrow(xpred), size = p_pred)
  
  # Creating new mutated sets
  xcal_new <- rbind(xcal[-cal_mut, ], xpred[pred_mut, ])
  ycal_new <- c(ycal[-cal_mut], ypred[pred_mut])
  
  xpred_new <- rbind(xpred[-pred_mut, ], xcal[cal_mut, ])
  ypred_new <- c(ypred[-pred_mut], ycal[cal_mut])
  
  cat("Final MLM Split: Train =", nrow(xcal_new), "| Test =", nrow(xpred_new), "\n")
  
  # Preparing Final DataFrames
  df_train_mlm <- data.frame(Y = ycal_new, xcal_new)
  df_test_mlm  <- data.frame(Y = ypred_new, xpred_new)
  
  # Adjusting column names
  colnames(df_train_mlm)[-1] <- paste0("Factor", 1:num_factors)
  colnames(df_test_mlm)[-1]  <- paste0("Factor", 1:num_factors)
}

##############################################################################
# --- STEP 8: MLR Regression and Cross-Validation ---
{
  # Caret Configuration (10-Fold CV)
  ctrl <- trainControl(method = "cv", number = 5)
  
  set.seed(123)
  # Train model using mutated dataset (df_train_mlm)
  model_caret <- train(Y ~ ., data = df_train_mlm, method = "lm", trControl = ctrl)
  
  # --- Calibration Metrics (Whole Training Set) ---
  pred_cal <- predict(model_caret, df_train_mlm)
  obs_cal  <- df_train_mlm$Y
  
  rmse_cal <- rmse(obs_cal, pred_cal)
  mae_cal  <- mae(obs_cal, pred_cal)
  r2_cal   <- cor(obs_cal, pred_cal)^2
  
  # --- Cross-Validation Metrics (CV) ---
  rmse_cv <- model_caret$results$RMSE
  mae_cv  <- model_caret$results$MAE
  r2_cv   <- model_caret$results$Rsquared
  
  # --- Prediction Metrics (External Test) ---
  pred_test <- predict(model_caret, df_test_mlm)
  obs_test  <- df_test_mlm$Y
  
  rmse_test <- rmse(obs_test, pred_test)
  mae_test  <- mae(obs_test, pred_test)
  r2_test   <- cor(obs_test, pred_test)^2
  
  # --- Display ---
  cat("\n=============================================\n")
  cat(" FINAL RESULTS (PARAFAC-MLR with MLM)\n")
  cat("=============================================\n")
  cat(sprintf("CALIBRATION (Mutated Train):\n   RMSE: %.4f | MAE: %.4f | R²: %.4f\n", rmse_cal, mae_cal, r2_cal))
  cat("---------------------------------------------\n")
  cat(sprintf("CROSS-VALIDATION (10-Fold):\n   RMSE: %.4f | MAE: %.4f | R²: %.4f\n", rmse_cv, mae_cv, r2_cv))
  cat("---------------------------------------------\n")
  cat(sprintf("PREDICTION (Mutated Test):\n   RMSE: %.4f | MAE: %.4f | R²: %.4f\n", rmse_test, mae_test, r2_test))
  cat("=============================================\n")
}

############################################################################
# --- STEP 9: Separate Plots ---
{
  lims <- range(c(0, 100, Y_conc))
  
  # Plot 1: Calibration
  dev.new()
  plot(obs_cal, pred_cal, pch=19, col="blue", xlim=lims, ylim=lims,
       xlab="Actual Concentration (%)", ylab="Predicted Concentration (%)",
       main="CALIBRATION (MLM)")
  abline(0, 1, lty=2, lwd=2, col="gray")
  grid()
  legend("topleft", legend=paste0("R² Cal: ", round(r2_cal, 4)), bty="n")
  
  # Plot 2: Prediction
  dev.new()
  plot(obs_test, pred_test, pch=19, col="red", xlim=lims, ylim=lims,
       xlab="Actual Concentration (%)", ylab="Predicted Concentration (%)",
       main="PREDICTION (MLM)")
  abline(0, 1, lty=2, lwd=2, col="gray")
  grid()
  legend("topleft", legend=paste0("R² Test: ", round(r2_test, 4)), bty="n", text.col="red")
  
  # Plot 3: Residuals
  dev.new()
  plot(obs_test, (obs_test - pred_test), pch=19, col="darkgreen",
       xlab="Actual Concentration (%)", ylab="Residual", main="Prediction Residuals")
  abline(h=0, lty=2, col="black")
  grid()
}
