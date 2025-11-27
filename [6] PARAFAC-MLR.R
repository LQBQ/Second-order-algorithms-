####
# SCRIPT [8] ATUALIZADO: PARAFAC-MLR (Com MLM - Mutation Logic)
#
# OBJETIVO:
# 1. Carregar dados processados.
# 2. Seleção balanceada (1 Pura por marca + Adulteradas).
# 3. Modelagem PARAFAC e Remoção de Outliers.
# 4. Divisão Treino/Teste via MLM (Kennard-Stone + Mutação).
# 5. Regressão MLR com Validação Cruzada.
####

# --- PASSO 1: Carregar Pacotes Necessários ---
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

# --- PASSO 2: Carregar Arquivos de Dados ---
{
  cat("Por favor, selecione o arquivo 'eem_cube_cleaned.mat'...\n")
  mat_cube_path <- rstudioapi::selectFile(caption = "Selecione eem_cube_cleaned.mat",
                                          filter = "MAT Files (*.mat)")
  if (!nzchar(mat_cube_path)) stop("Seleção cancelada.")
  
  data_cube <- readMat(mat_cube_path)
  names(data_cube) <- gsub("\\.mat$", "", names(data_cube))
  
  X_original <- data_cube$cube.cleaned
  sample_names <- unlist(data_cube$sample.names)
  concentrations <- as.vector(data_cube$concentrations) 
  nmEM <- as.vector(data_cube$nm.emission)
  nmEX <- as.vector(data_cube$nm.excitation)
  
  cat("Dados carregados com sucesso!\n")
}

# --- PASSO 3: Preparação dos Dados (Seleção Balanceada) ---
{
  Y_multiclasse <- factor(
    ifelse(grepl("^A", sample_names), "Marca A",
           ifelse(grepl("^B", sample_names), "Marca B",
                  ifelse(grepl("^C", sample_names), "Marca C",
                         ifelse(grepl("^D", sample_names), "Marca D", "Puro")))
    )
  )
    
    set.seed(42) 
    
    # Seleção de Adulterados (Máximo disponível ou 10)
    sel_indices_A <- sample(which(Y_multiclasse == "Marca A"), size = min(10, sum(Y_multiclasse == "Marca A")))
    sel_indices_B <- sample(which(Y_multiclasse == "Marca B"), size = min(10, sum(Y_multiclasse == "Marca B")))
    sel_indices_C <- sample(which(Y_multiclasse == "Marca C"), size = min(10, sum(Y_multiclasse == "Marca C")))
    sel_indices_D <- sample(which(Y_multiclasse == "Marca D"), size = min(10, sum(Y_multiclasse == "Marca D")))
    
    # Seleção de Puros (1 de cada marca)
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
    
    cat(paste("\nDataset preparado. Dimensões X:", paste(dim(X), collapse=" x "), "\n"))
}

# --- PASSO 4: Modelo PARAFAC ---
{
  set.seed(24)
  num_factors <- 7 
  
  cat(paste("Modelando PARAFAC com", num_factors, "fatores...\n"))
  model <- parafac(X, nfac=num_factors, nstart=75, maxit=3000, 
                   ctol=1e-6, parallel=FALSE, output = c("best", "all"))
  cat(paste("Modelo PARAFAC R² (Fit):", round(model$Rsq, 5), "\n"))
}

# --- PASSO 5: Diagnóstico e Remoção Manual ---
{
  # Diagnóstico Visual (Pule se quiser apenas rodar direto)
  rss <- numeric(dim(X)[1])
  for(i in 1:length(rss)) {
    recon_i <- matrix(0, nrow = dim(X)[2], ncol = dim(X)[3]) 
    for(k in 1:num_factors) recon_i <- recon_i + model$A[i,k] * (outer(model$B[,k], model$C[,k]))
    rss[i] <- sum((X[i,,] - recon_i)^2)
  }
  
  # LISTA DE REMOÇÃO MANUAL
  nomes_para_remover <- c() 
  
  if (length(nomes_para_remover) > 0) {
    outliers_indices <- which(nomes_amostras %in% nomes_para_remover)
    if (length(outliers_indices) > 0) {
      cat("\nRemovendo", length(outliers_indices), "amostras e recalculando...\n")
      X <- X[-outliers_indices, , ]
      Y_conc <- Y_conc[-outliers_indices]
      nomes_amostras <- nomes_amostras[-outliers_indices]
      
      set.seed(24) 
      model <- parafac(X, nfac=num_factors, nstart=50, maxit=3000, ctol=1e-6, parallel=F)
      cat(paste("Novo R²:", round(model$Rsq, 5), "\n"))
    } 
  }
}

# --- PASSO 6: Visualização Loadings ---
{
  dev.new()
  par(mfrow=c(1,2))
  matplot(nmEM, model$B, type ="l", xlab="Emissão (nm)", ylab="Loadings", main="Modo Emissão")
  matplot(nmEX, model$C, type ="l", xlab="Excitação (nm)", ylab="Loadings", main="Modo Excitação")
  par(mfrow=c(1,1))
}

# --- PASSO 7: Separação Treino/Teste via MLM (Mutation) ---
{
  # 1. Configuração Inicial e Kennard-Stone
  perc <- 0.7
  A <- model$A
  y <- Y_conc
  ntrain <- ceiling(perc * nrow(A))
  
  combined_data <- cbind(A, y)
  
  # Kennard-Stone Inicial (Euclidiana)
  set.seed(2) 
  ks_result <- kenStone(combined_data, k = ntrain, metric = "euclid")
  
  calibration_idx <- ks_result$model
  prediction_idx <- ks_result$test
  
  xcal <- A[calibration_idx, ]
  ycal <- y[calibration_idx]
  
  xpred <- A[prediction_idx, ]
  ypred <- y[prediction_idx]
  
  cat("Divisão KS Inicial: Treino =", nrow(xcal), "| Teste =", nrow(xpred), "\n")
  
  # 2. Lógica de Mutação (MLM)
  cat("Aplicando Mutação (Troca de 20% das amostras)...\n")
  set.seed(1) # Seed para garantir reprodutibilidade da mutação
  
  prob <- 0.2
  p_cal <- ceiling(prob * nrow(xcal))
  p_pred <- ceiling(prob * nrow(xpred))
  
  cal_mut <- sample(1:nrow(xcal), size = p_cal)
  pred_mut <- sample(1:nrow(xpred), size = p_pred)
  
  # Criando novos conjuntos mutados
  xcal_new <- rbind(xcal[-cal_mut, ], xpred[pred_mut, ])
  ycal_new <- c(ycal[-cal_mut], ypred[pred_mut])
  
  xpred_new <- rbind(xpred[-pred_mut, ], xcal[cal_mut, ])
  ypred_new <- c(ypred[-pred_mut], ycal[cal_mut])
  
  cat("Divisão MLM Final: Treino =", nrow(xcal_new), "| Teste =", nrow(xpred_new), "\n")
  
  # Preparando DataFrames Finais
  df_train_mlm <- data.frame(Y = ycal_new, xcal_new)
  df_test_mlm  <- data.frame(Y = ypred_new, xpred_new)
  
  # Ajustando nomes das colunas
  colnames(df_train_mlm)[-1] <- paste0("Factor", 1:num_factors)
  colnames(df_test_mlm)[-1]  <- paste0("Factor", 1:num_factors)
}

# --- PASSO 8: Regressão MLR e Validação Cruzada ---
{
  # Configuração do Caret (10-Fold CV)
  ctrl <- trainControl(method = "cv", number = 10)
  
  set.seed(123)
  # Treina o modelo usando o dataset mutado (df_train_mlm)
  model_caret <- train(Y ~ ., data = df_train_mlm, method = "lm", trControl = ctrl)
  
  # --- Métricas Calibração (Treino Inteiro) ---
  pred_cal <- predict(model_caret, df_train_mlm)
  obs_cal  <- df_train_mlm$Y
  
  rmse_cal <- rmse(obs_cal, pred_cal)
  mae_cal  <- mae(obs_cal, pred_cal)
  r2_cal   <- cor(obs_cal, pred_cal)^2
  
  # --- Métricas Validação Cruzada (CV) ---
  rmse_cv <- model_caret$results$RMSE
  mae_cv  <- model_caret$results$MAE
  r2_cv   <- model_caret$results$Rsquared
  
  # --- Métricas Predição (Teste Externo) ---
  pred_test <- predict(model_caret, df_test_mlm)
  obs_test  <- df_test_mlm$Y
  
  rmse_test <- rmse(obs_test, pred_test)
  mae_test  <- mae(obs_test, pred_test)
  r2_test   <- cor(obs_test, pred_test)^2
  
  # --- Exibição ---
  cat("\n=============================================\n")
  cat(" RESULTADOS FINAIS (PARAFAC-MLR com MLM)\n")
  cat("=============================================\n")
  cat(sprintf("CALIBRAÇÃO (Treino Mutado):\n   RMSE: %.4f | MAE: %.4f | R²: %.4f\n", rmse_cal, mae_cal, r2_cal))
  cat("---------------------------------------------\n")
  cat(sprintf("VALIDAÇÃO CRUZADA (10-Fold):\n   RMSE: %.4f | MAE: %.4f | R²: %.4f\n", rmse_cv, mae_cv, r2_cv))
  cat("---------------------------------------------\n")
  cat(sprintf("PREDIÇÃO (Teste Mutado):\n   RMSE: %.4f | MAE: %.4f | R²: %.4f\n", rmse_test, mae_test, r2_test))
  cat("=============================================\n")
}

# --- PASSO 9: Gráficos Separados ---
{
  lims <- range(c(0, 100, Y_conc))
  
  # Gráfico 1: Calibração
  dev.new()
  plot(obs_cal, pred_cal, pch=19, col="blue", xlim=lims, ylim=lims,
       xlab="Concentração Real (%)", ylab="Concentração Predita (%)",
       main="CALIBRAÇÃO (MLM)")
  abline(0, 1, lty=2, lwd=2, col="gray")
  grid()
  legend("topleft", legend=paste0("R² Cal: ", round(r2_cal, 4)), bty="n")
  
  # Gráfico 2: Predição
  dev.new()
  plot(obs_test, pred_test, pch=19, col="red", xlim=lims, ylim=lims,
       xlab="Concentração Real (%)", ylab="Concentração Predita (%)",
       main="PREDIÇÃO (MLM)")
  abline(0, 1, lty=2, lwd=2, col="gray")
  grid()
  legend("topleft", legend=paste0("R² Test: ", round(r2_test, 4)), bty="n", text.col="red")
  
  # Gráfico 3: Resíduos
  dev.new()
  plot(obs_test, (obs_test - pred_test), pch=19, col="darkgreen",
       xlab="Concentração Real (%)", ylab="Resíduo", main="Resíduos de Predição")
  abline(h=0, lty=2, col="black")
  grid()
}
