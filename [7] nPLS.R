####
# SCRIPT [9] ATUALIZADO: N-PLS (Regression)
# Baseado na estrutura do script [6] PARAFAC-MLR
#
# OBJETIVO:
# 1. Carregar dados processados ('eem_cube_cleaned.mat').
# 2. Seleção balanceada (1 Pura por marca + Adulteradas).
# 3. Modelagem N-PLS (sNPLS) e Diagnóstico de Outliers.
# 4. Divisão Treino/Teste via MLM (Kennard-Stone + Mutação).
# 5. Regressão N-PLS Final com Validação Cruzada.
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
    "sNPLS", "multiway", "ThreeWay", "R.matlab", "plot3D", "plotly", 
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
  
  # Carrega matriz e dados
  data_cube <- readMat(mat_cube_path)
  names(data_cube) <- gsub("\\.mat$", "", names(data_cube))
  
  X_original <- data_cube$cube.cleaned
  sample_names <- unlist(data_cube$sample.names)
  concentrations <- as.vector(data_cube$concentrations) 
  nmEM <- as.vector(data_cube$nm.emission)
  nmEX <- as.vector(data_cube$nm.excitation)
  
  cat("Dados carregados com sucesso!\n")
  cat(paste("Dimensões do Cubo: [", paste(dim(X_original), collapse = " x "), "]\n"))
}

# --- PASSO 3: Preparação dos Dados (Seleção Balanceada) ---
{
  # Definir marcas
  Y_multiclasse <- factor(
    ifelse(grepl("^A", sample_names), "Marca A",
           ifelse(grepl("^B", sample_names), "Marca B",
                  ifelse(grepl("^C", sample_names), "Marca C",
                         ifelse(grepl("^D", sample_names), "Marca D", "Puro")))
    )
  )
    
    set.seed(42) 
    
    # --- SELEÇÃO DE ADULTERADOS ---
    indices_A <- which(Y_multiclasse == "Marca A")
    indices_B <- which(Y_multiclasse == "Marca B")
    indices_C <- which(Y_multiclasse == "Marca C")
    indices_D <- which(Y_multiclasse == "Marca D")
    
    sel_indices_A <- sample(indices_A, size = min(10, length(indices_A)), replace = FALSE)
    sel_indices_B <- sample(indices_B, size = min(10, length(indices_B)), replace = FALSE)
    sel_indices_C <- sample(indices_C, size = min(10, length(indices_C)), replace = FALSE)
    sel_indices_D <- sample(indices_D, size = min(10, length(indices_D)), replace = FALSE)
    
    # --- SELEÇÃO DE PUROS (1 de cada marca) ---
    indices_todos_zeros <- which(concentrations == 0)
    nomes_zeros <- sample_names[indices_todos_zeros]
    
    # Tenta identificar marca no nome da pura (ex: "P A 0 1")
    idx_pure_A <- indices_todos_zeros[grepl("A", nomes_zeros)]
    idx_pure_B <- indices_todos_zeros[grepl("B", nomes_zeros)]
    idx_pure_C <- indices_todos_zeros[grepl("C", nomes_zeros)]
    idx_pure_D <- indices_todos_zeros[grepl("D", nomes_zeros)]
    
    # Sorteia 1 de cada se existir
    sel_pure_A <- if(length(idx_pure_A) > 0) sample(idx_pure_A, 1) else NULL
    sel_pure_B <- if(length(idx_pure_B) > 0) sample(idx_pure_B, 1) else NULL
    sel_pure_C <- if(length(idx_pure_C) > 0) sample(idx_pure_C, 1) else NULL
    sel_pure_D <- if(length(idx_pure_D) > 0) sample(idx_pure_D, 1) else NULL
    
    indices_Puras_Final <- c(sel_pure_A, sel_pure_B, sel_pure_C, sel_pure_D)
    
    cat("Índices das Puras Selecionadas:", indices_Puras_Final, "\n")
    
    # --- COMBINAÇÃO FINAL ---
    indices_balanceados <- c(sel_indices_A, sel_indices_B, sel_indices_C, sel_indices_D, 
                             indices_Puras_Final)
    
    X <- X_original[indices_balanceados, , ]
    Y_conc <- concentrations[indices_balanceados] 
    nomes_amostras <- sample_names[indices_balanceados]
    
    cat(paste("\nDataset preparado.\n"))
    cat("Dimensões X:", paste(dim(X), collapse=" x "), "\n")
}

# --- PASSO 4: Modelo N-PLS Preliminar ---
{
  set.seed(24)
  num_factors <- 6 
  
  cat(paste("Modelando N-PLS com", num_factors, "variáveis latentes...\n"))
  
  # sNPLS é a função do pacote. Se threshold=0, ele age como um N-PLS padrão (sem seleção de variável)
  model <- sNPLS(X, as.matrix(Y_conc), ncomp = num_factors, threshold_j = 0, threshold_k = 0)
  
  # Extraindo fitted values
  fitted_vals <- as.vector(fitted(model))
  r2_fit <- cor(Y_conc, fitted_vals)^2
  
  cat(paste("Modelo N-PLS R² (Fit Preliminar):", round(r2_fit, 5), "\n"))
}

# --- PASSO 5: Diagnóstico de Outliers (Q vs T²) ---
{
  # 1. Hotelling T² (Baseado nos Scores T do N-PLS)
  Scores <- model$T
  n_samples <- dim(X)[1]
  
  S_T <- cov(Scores)
  T2 <- mahalanobis(Scores, colMeans(Scores), S_T)
  
  # 2. Resíduos (Q_y)
  # No N-PLS foca-se muito no ajuste do Y. Usaremos o resíduo de Y quadrado.
  resid_y <- (Y_conc - fitted_vals)^2
  Q_y <- resid_y
  
  # Limites Estatísticos
  alpha <- 0.05
  T2_lim <- num_factors * (n_samples - 1) / (n_samples - num_factors) * qf(1 - alpha, num_factors, n_samples - num_factors)
  Q_lim <- median(Q_y) + 3 * IQR(Q_y)
  
  dev.new()
  colores <- ifelse(Y_conc > 0, "red", "blue")
  plot(T2, Q_y, pch = 19, col = colores,
       xlab = expression(Hotelling~T^2), ylab = expression(Q[Y]~(Erro~Quadrático~Y)),
       main = "Diagnóstico de Outliers (N-PLS)")
  abline(v = T2_lim, col = "darkgreen", lty = 2)
  abline(h = Q_lim,  col = "darkgreen", lty = 2)
  text(T2, Q_y, labels = nomes_amostras, pos = 3, cex = 0.6)
  legend("topright", legend = c("Adulterado", "Puro"), col = c("red", "blue"), pch = 19)
}

# --- PASSO 5.5: Remoção Manual de Outliers ---
{
  # 1. DEFINA AQUI OS NOMES EXATOS PARA REMOVER
  nomes_para_remover <- c() 
  
  if (length(nomes_para_remover) > 0) {
    outliers_indices <- which(nomes_amostras %in% nomes_para_remover)
    
    if (length(outliers_indices) > 0) {
      cat("\nRemovendo", length(outliers_indices), "amostras:", nomes_amostras[outliers_indices], "\n")
      
      X <- X[-outliers_indices, , ]
      Y_conc <- Y_conc[-outliers_indices]
      nomes_amostras <- nomes_amostras[-outliers_indices]
      
      cat("Recalculando N-PLS...\n")
      set.seed(24) 
      model <- sNPLS(X, as.matrix(Y_conc), ncomp = num_factors, threshold_j = 0, threshold_k = 0)
      
      fitted_vals <- as.vector(fitted(model))
      cat(paste("Novo R²:", round(cor(Y_conc, fitted_vals)^2, 5), "\n"))
    } 
  } else {
    cat("Nenhuma remoção solicitada.\n")
  }
}

# --- PASSO 6: Visualização dos Pesos (Weights) ---
{
  # N-PLS usa pesos (Weights) Wj (Emissão) e Wk (Excitação)
  # O objeto sNPLS armazena isso em model$Wj e model$Wk
  
  dev.new()
  # Plotando os pesos do primeiro e segundo fator latente (LV1 e LV2)
  matplot(nmEM, model$Wj[, 1:2], type ="l", lty=1, col=c("black", "red"),
          xlab="Emissão (nm)", ylab="Pesos (Weights Wj)", main="Pesos N-PLS - Emissão")
  legend("topright", legend=c("LV 1", "LV 2"), col=c("black", "red"), lty=1)
  
  dev.new()
  matplot(nmEX, model$Wk[, 1:2], type ="l", lty=1, col=c("black", "red"),
          xlab="Excitação (nm)", ylab="Pesos (Weights Wk)", main="Pesos N-PLS - Excitação")
  legend("topright", legend=c("LV 1", "LV 2"), col=c("black", "red"), lty=1)
}

# --- PASSO 7: Separação Treino/Teste via MLM (Mutation) ---
{
  # 1. Configuração Inicial e Kennard-Stone
  # Para KS em dados 3D, desdobramos o cubo X em uma matriz (Amostras x Var)
  X_flat <- matrix(as.vector(X), nrow = dim(X)[1]) 
  
  perc <- 0.7
  ntrain <- ceiling(perc * nrow(X_flat))
  
  # Kennard-Stone (Euclidiana no X desdobrado)
  set.seed(2) 
  ks_result <- kenStone(X_flat, k = ntrain, metric = "euclid")
  
  calibration_idx <- ks_result$model
  prediction_idx <- ks_result$test
  
  cat("Divisão KS Inicial: Treino =", length(calibration_idx), "| Teste =", length(prediction_idx), "\n")
  
  # 2. Lógica de Mutação (MLM)
  cat("Aplicando Mutação (Troca de 20% das amostras)...\n")
  set.seed(1)
  
  prob <- 0.2
  p_cal <- ceiling(prob * length(calibration_idx))
  p_pred <- ceiling(prob * length(prediction_idx))
  
  # Índices relativos para troca
  cal_mut_pos <- sample(1:length(calibration_idx), size = p_cal)
  pred_mut_pos <- sample(1:length(prediction_idx), size = p_pred)
  
  # Trocando os índices reais
  idx_to_pred <- calibration_idx[cal_mut_pos]
  idx_to_cal  <- prediction_idx[pred_mut_pos]
  
  final_cal_idx <- c(calibration_idx[-cal_mut_pos], idx_to_cal)
  final_pred_idx <- c(prediction_idx[-pred_mut_pos], idx_to_pred)
  
  # Criando os Tensores Finais
  # drop=FALSE mantém a estrutura 3D mesmo se houver 1 amostra
  X_cal <- X[final_cal_idx, , , drop=FALSE]
  Y_cal <- Y_conc[final_cal_idx]
  
  X_test <- X[final_pred_idx, , , drop=FALSE]
  Y_test <- Y_conc[final_pred_idx]
  
  cat("Divisão MLM Final: Treino =", length(final_cal_idx), "| Teste =", length(final_pred_idx), "\n")
}

# --- PASSO 8: Regressão N-PLS Final com Validação Cruzada ---
{
  set.seed(123)
  
  # 1. Modelo de Calibração Final (sNPLS)
  final_model <- sNPLS(X_cal, as.matrix(Y_cal), ncomp = num_factors, threshold_j=0, threshold_k=0)
  
  # 2. Métricas de Calibração (Fitted)
  pred_cal <- as.vector(fitted(final_model))
  
  rmse_cal <- rmse(Y_cal, pred_cal)
  mae_cal  <- mae(Y_cal, pred_cal)
  r2_cal   <- cor(Y_cal, pred_cal)^2
  
  # 3. Validação Cruzada (10-Fold) - Loop Manual
  # (sNPLS não tem integração direta com 'caret' train, fazemos manual para garantir robustez)
  cat("Rodando Validação Cruzada (10-Fold)...\n")
  folds <- createFolds(Y_cal, k = 10, list = TRUE)
  
  cv_preds_vector <- numeric(length(Y_cal))
  # Armazenar resultados em dataframe temporário para ordenar depois
  temp_res <- data.frame(Index = integer(), Pred = numeric())
  
  for(i in 1:10) {
    idx_val <- folds[[i]]
    idx_train <- setdiff(1:length(Y_cal), idx_val)
    
    # Subsets do Fold
    xc_fold <- X_cal[idx_train, , , drop=FALSE]
    yc_fold <- Y_cal[idx_train]
    xv_fold <- X_cal[idx_val, , , drop=FALSE]
    
    # Modelo CV
    mod_cv <- sNPLS(xc_fold, as.matrix(yc_fold), ncomp = num_factors, threshold_j=0, threshold_k=0)
    
    # Predição CV
    # ATENÇÃO: predict no sNPLS exige newX com dimensão correta
    pred_fold <- as.vector(predict(mod_cv, newX = xv_fold))
    
    temp_res <- rbind(temp_res, data.frame(Index = idx_val, Pred = pred_fold))
  }
  
  # Reordenar predições para alinhar com Y_cal
  temp_res <- temp_res[order(temp_res$Index), ]
  cv_preds <- temp_res$Pred
  
  rmse_cv <- rmse(Y_cal, cv_preds)
  mae_cv  <- mae(Y_cal, cv_preds)
  r2_cv   <- cor(Y_cal, cv_preds)^2
  
  # 4. Predição no Teste Externo
  pred_test <- as.vector(predict(final_model, newX = X_test))
  
  rmse_test <- rmse(Y_test, pred_test)
  mae_test  <- mae(Y_test, pred_test)
  r2_test   <- cor(Y_test, pred_test)^2
  
  # --- Exibição dos Resultados ---
  cat("\n=============================================\n")
  cat(" RESULTADOS FINAIS (N-PLS com MLM)\n")
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
  plot(Y_cal, pred_cal, pch=19, col="blue", xlim=lims, ylim=lims,
       xlab="Concentração Real (%)", ylab="Concentração Predita (%)",
       main="CALIBRAÇÃO N-PLS (MLM)")
  abline(0, 1, lty=2, lwd=2, col="gray")
  grid()
  legend("topleft", legend=paste0("R² Cal: ", round(r2_cal, 4)), bty="n")
  
  # Gráfico 2: Predição
  dev.new()
  plot(Y_test, pred_test, pch=19, col="red", xlim=lims, ylim=lims,
       xlab="Concentração Real (%)", ylab="Concentração Predita (%)",
       main="PREDIÇÃO N-PLS (MLM)")
  abline(0, 1, lty=2, lwd=2, col="gray")
  grid()
  legend("topleft", legend=paste0("R² Test: ", round(r2_test, 4)), bty="n", text.col="red")
  
  # Gráfico 3: Resíduos
  dev.new()
  plot(Y_test, (Y_test - pred_test), pch=19, col="darkgreen",
       xlab="Concentração Real (%)", ylab="Resíduo (Real - Predito)",
       main="Resíduos de Predição (N-PLS)")
  abline(h=0, lty=2, col="black")
  grid()
}
