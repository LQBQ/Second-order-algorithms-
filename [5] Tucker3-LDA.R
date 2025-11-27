####
# SCRIPT [7] ATUALIZADO: TUCKER3-LDA
#
# OBJETIVO:
# 1. Carregar dados processados ('eem_cube_cleaned.mat').
# 2. Gerar rótulos (Y) dinamicamente a partir das concentrações e nomes.
# 3. Realizar balanceamento de classe robusto (Igual ao script PARAFAC).
# 4. Construir modelo TUCKER-3 e classificador LDA.
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
    "prospectr", "MASS", "caret", "rstudioapi", "abind", "parallel"
  )
  
  load_required_packages(required_packages)
}

# --- PASSO 2: Carregar Arquivos de Dados ---
{
  cat("Por favor, selecione o arquivo 'eem_cube_cleaned.mat' (contém dados pré-processados)...\n")
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
  cat(paste("Dimensões do Cubo (Amostras x Emissão x Excitação): [", 
            paste(dim(X_original), collapse = " x "), "]\n"))
}

# --- PASSO 3: Geração de Rótulos e Balanceamento de Dados ---
{
  # 1. Gerar rótulos baseados nos nomes e concentrações
  
  # Rótulo binário (Adulterado vs. Puro) para o LDA
  Y_binario <- factor(ifelse(concentrations > 0, "Adulterado", "Puro"))
  
  # Rótulo multiclasse para visualização e seleção
  Y_multiclasse <- factor(
    ifelse(grepl("^A", sample_names), "Marca A",
           ifelse(grepl("^B", sample_names), "Marca B",
                  ifelse(grepl("^C", sample_names), "Marca C",
                         ifelse(grepl("^D", sample_names), "Marca D", "Puro"))))
  )
  
  cat("Rótulos gerados dinamicamente:\n")
  print(table(Y_binario, Y_multiclasse))
  
  # 2. Balanceamento de Classe (Lógica idêntica ao PARAFAC-LDA)
  
  set.seed(42) # Para reprodutibilidade
  
  indices_A <- which(Y_multiclasse == "Marca A")
  indices_B <- which(Y_multiclasse == "Marca B")
  indices_C <- which(Y_multiclasse == "Marca C")
  indices_D <- which(Y_multiclasse == "Marca D")
  indices_Puras <- which(Y_multiclasse == "Puro")
  
  # Sortear 10 amostras de CADA estrato adulterante (ou o máximo disponível)
  sel_indices_A <- sample(indices_A, size = min(10, length(indices_A)), replace = FALSE)
  sel_indices_B <- sample(indices_B, size = min(10, length(indices_B)), replace = FALSE)
  sel_indices_C <- sample(indices_C, size = min(10, length(indices_C)), replace = FALSE)
  sel_indices_D <- sample(indices_D, size = min(10, length(indices_D)), replace = FALSE)
  
  # Combinar índices (Total ~80 amostras)
  indices_balanceados <- c(sel_indices_A, sel_indices_B, sel_indices_C, sel_indices_D, 
                           indices_Puras)
  
  # 4. Criar o cubo de dados final e balanceado
  X <- X_original[indices_balanceados, , ]
  
  # 5. Criar os rótulos finais
  Y <- Y_binario[indices_balanceados]
  Y_multi_balanceado <- Y_multiclasse[indices_balanceados]
  nomes_amostras <- sample_names[indices_balanceados]
  
  cat(paste("\nCubo de dados final balanceado criado com dimensões:", 
            paste(dim(X), collapse = " x "), "\n"))
  print(table(Y))
  
  # 6. Médias para visualização
  X1 <- X[Y == "Adulterado", , ]
  X2 <- X[Y == "Puro", , ]
  
  X1m <- colMeans(X1)
  X2m <- colMeans(X2)
}

# --- PASSO 4: Visualização ---

  dev.new() +
  filled.contour(x = nmEM, y = nmEX, z = X1m, 
                 color.palette = topo.colors, main = "Média Classe 1 (Adulterado)",
                 xlab = "Emissão (nm)", ylab = "Excitação (nm)")
  
  dev.new() +
  filled.contour(x = nmEM, y = nmEX, z = X2m, 
                 color.palette = terrain.colors, main = "Média Classe 2 (Puro)",
                 xlab = "Emissão (nm)", ylab = "Excitação (nm)")
  
  dev.new() +
  matplot(nmEM, X1m, type="l", xlab = "Comprimento de Onda Emissão (nm)", ylab = "Intensidade", 
          main = "Classe 1 - Normal (Modo Emissão)")
  
  dev.new() +
  matplot(nmEX, t(X1m), type="l", xlab = "Comprimento de Onda Excitação (nm)", ylab = "Intensidade", 
          main = "Classe 1 - Normal (Modo Excitação)")


# --- PASSO 5: Modelo TUCKER-3 ---
{
  set.seed(24)
  
  # DEFINA OS RANKS AQUI (Amostras, Emissão, Excitação)
  # Diferente do PARAFAC (que tem 1 rank), o Tucker tem 3.
  # Sugestão: c(3, 3, 3) ou c(4, 4, 4) dependendo da complexidade.
  ranks_tucker <- c(7, 7, 7) 
  
  cat(paste("Modelando TUCKER-3 com ranks:", paste(ranks_tucker, collapse=","), "...\n"))
  
  # output="all" retorna A, B, C e o Core Tensor (GA)
  model <- tucker(X, nfac=ranks_tucker, nstart=75, maxit=5000, ctol=1e-6, output = c("best", "all"))
  cat(paste("Modelo Tucker-3 R²:", model$Rsq, "\n"))
}

# --- PASSO 6: Diagnóstico de Outliers (RSS e Q vs T²) ---
{
  # --- Erro de reconstrução por amostra (RSS / Q) ---
  # No Tucker: X_hat_i = B * (Somatorio_p A[i,p] * Core[p,,]) * t(C)
  
  A <- model$A
  B <- model$B
  C <- model$C
  G <- model$G # Core Tensor
  
  # Garantir que G seja um array 3D [RankA x RankB x RankC]
  if(is.matrix(G)) G <- array(G, dim = ranks_tucker)
  
  n_samples <- dim(X)[1]
  rss <- numeric(n_samples)
  
  cat("Calculando resíduos do Tucker (isso pode demorar um pouco)...\n")
  
  for(i in 1:n_samples) {
    # 1. Calcular a fatia do Core ponderada pelos scores da amostra i
    # Core_Slice = Somatorio (A[i, p] * G[p, , ])
    G_slice <- matrix(0, nrow = ranks_tucker[2], ncol = ranks_tucker[3])
    
    for(p in 1:ranks_tucker[1]) {
      G_slice <- G_slice + A[i, p] * G[p, , ]
    }
    
    # 2. Reconstruir a matriz da amostra i
    recon_i <- B %*% G_slice %*% t(C)
    
    # 3. Calcular resíduo
    resid_i <- X[i,,] - recon_i
    rss[i] <- sum(resid_i^2)
  }
  
  thr <- median(rss) + 3 * IQR(rss)
  outliers_rss <- which(rss > thr)
  
  dev.new()
  plot(rss, pch = 19, ylab = "RSS (erro de reconstrução)", xlab = "Índice da Amostra",
       main = "Erro de reconstrução por amostra (Tucker-3)")
  abline(h = thr, col = "red", lty = 2)
  text(outliers_rss, rss[outliers_rss], labels = nomes_amostras[outliers_rss], pos = 3, col = "red", cex=0.8)
  
  # --- T² Hotelling vs Q residual ---
  # Usamos a matriz A (Scores das Amostras)
  S_A <- cov(A)
  T2 <- mahalanobis(A, colMeans(A), S_A)
  
  nfac_A <- ranks_tucker[1] # Graus de liberdade do modo Amostra
  alpha <- 0.05
  T2_lim <- nfac_A * (n_samples - 1) / (n_samples - nfac_A) * qf(1 - alpha, nfac_A, n_samples - nfac_A)
  
  Q <- rss
  Q_lim <- median(Q) + 3 * IQR(Q)
  
  dev.new()
  plot(T2, Q, pch = 19, col = ifelse(Y == "Adulterado", "red", "blue"),
       xlab = expression(Hotelling~T^2),
       ylab = expression(Q~"(Erro de reconstrução)"),
       main = expression("Tucker-3: Gráfico Q vs T"^2))
  abline(v = T2_lim, col = "darkgreen", lty = 2)
  abline(h = Q_lim,  col = "darkgreen", lty = 2)
  text(T2, Q, labels = nomes_amostras, pos = 3, cex = 0.6)
  legend("topright",
         legend = c("Adulterado", "Puro", "Limite 95%"),
         col = c("red", "blue", "darkgreen"), pch = c(19,19,NA), lty = c(NA,NA,2))
  
  out_QT <- which(T2 > T2_lim | Q > Q_lim)
  cat("Amostras identificadas como outliers (Q–T²):", nomes_amostras[out_QT], "\n")
}

# --- PASSO 6.5: Remoção Manual de Outliers por NOME e Re-Modelagem (TUCKER) ---
{
  # 1. DEFINA AQUI OS NOMES EXATOS PARA REMOVER
  # Exemplo: nomes_para_remover <- c("B0900", "D0900", "PA02", "PB02", "PD02", "PD03")
  nomes_para_remover <- c("A0900", "C0900", "D0700", "D0010", "PB03", "PC02", "PC03") 
  
  if (length(nomes_para_remover) > 0) {
    
    outliers_indices <- which(nomes_amostras %in% nomes_para_remover)
    
    nomes_encontrados <- nomes_amostras[outliers_indices]
    nomes_errados <- setdiff(nomes_para_remover, nomes_encontrados)
    
    if (length(nomes_errados) > 0) {
      cat("\n[ATENÇÃO] Nomes não encontrados:", nomes_errados, "\n")
    }
    
    if (length(outliers_indices) > 0) {
      cat("\nRemovendo", length(outliers_indices), "amostras:", nomes_encontrados, "\n")
      
      # 2. REMOVER DO CUBO E METADADOS
      X <- X[-outliers_indices, , ]
      Y <- Y[-outliers_indices]
      nomes_amostras <- nomes_amostras[-outliers_indices]
      
      cat("Novas dimensões do cubo:", paste(dim(X), collapse=" x "), "\n")
      
      # 3. RECALCULAR O MODELO TUCKER-3
      cat("Recalculando modelo TUCKER-3 sem outliers...\n")
      set.seed(24) 
      
      # Re-roda o Tucker com os mesmos ranks
      model <- tucker(X, nfac=ranks_tucker, nstart=75, maxit=5000, ctol=1e-6, output = c("best", "all"))
      
      # Atualiza R²
      cat(paste("Modelo Tucker-3 R²:", model$Rsq, "\n"))
      
    } else {
      cat("Nenhum nome válido encontrado.\n")
    }
    
  } else {
    cat("Lista de remoção vazia. O modelo original será mantido.\n")
  }
}

# --- PASSO 7: Visualização de Scores e Loadings (Adaptado para Tucker) ---
{
  ## Plot Tucker Scores (Modo 1 - Amostras)
  dev.new()
  fator_x <- 1 
  fator_y <- 2 
  
  plot(model$A[, fator_x], model$A[, fator_y], 
       type = "n",
       xlab = paste("Fator", fator_x), 
       ylab = paste("Fator", fator_y), 
       main = "Scores TUCKER-3 (Matriz A)")
  
  text(x = model$A[Y == "Adulterado", fator_x], 
       y = model$A[Y == "Adulterado", fator_y], 
       labels = nomes_amostras[Y == "Adulterado"], 
       col = "red", cex = 0.7)
  
  text(x = model$A[Y == "Puro", fator_x], 
       y = model$A[Y == "Puro", fator_y], 
       labels = nomes_amostras[Y == "Puro"], 
       col = "blue", cex = 0.7)
  
  legend("topright", legend = c("Adulterado", "Puro"), col = c("red", "blue"), pch = 19)
  grid()
  
  ## Plot Tucker Loadings (Modo 2 - Emissão)
  dev.new()
  matplot(nmEM, model$B, type ="l", xlab = "Comprimento de Onda Emissão (nm)", 
          ylab = "Loadings", main = "Tucker Loadings - Emissão (Matriz B)")
  
  ## Plot Tucker Loadings (Modo 3 - Excitação)
  dev.new()
  matplot(nmEX, model$C, type ="l", xlab = "Comprimento de Onda Excitação (nm)", 
          ylab = "Loadings", main = "Tucker Loadings - Excitação (Matriz C)")
}

# --- PASSO 8: Separação Treino/Teste e Modelo LDA ---
{
  train_ratio <- 0.7 
  
  # Definir grupos numéricos (1=Adulterado, 2=Puro)
  group_numeric <- ifelse(Y == "Adulterado", 1, 2)
  
  # Scores das classes (Usando Matriz A do Tucker)
  # Diferente do PARAFAC, aqui temos 'ranks_tucker[1]' colunas
  scores_X1 <- model$A[group_numeric == 1, ] # Scores Adulterados
  scores_X2 <- model$A[group_numeric == 2, ] # Scores Puros
  
  # Calcular amostras de treino por classe
  ntrain1 <- ceiling(train_ratio * nrow(scores_X1))
  ntrain2 <- ceiling(train_ratio * nrow(scores_X2))
  
  # Selecionar amostras (Kennard-Stone)
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
  
  # Modelo LDA
  model_lda <- lda(train, group_train)
  model_lda_cv <- lda(train, group_train, CV=TRUE)
}

# --- PASSO 9: Métricas de Performance e Visualização LDA ---
{
  # Predições
  pred_train <- predict(model_lda, train)
  pred_test <- predict(model_lda, test)
  
  # Métricas
  ac_train <- mean(pred_train$class == group_train) 
  spec_train <- mean(pred_train$class[group_train == 1] == 1) # Classe 1 (Adulterado)
  sens_train <- mean(pred_train$class[group_train == 2] == 2) # Classe 2 (Puro)
  
  ac_cv <- mean(model_lda_cv$class == group_train)
  spec_cv <- mean(model_lda_cv$class[group_train == 1] == 1)
  sens_cv <- mean(model_lda_cv$class[group_train == 2] == 2)
  
  ac_test <- mean(pred_test$class == group_test) 
  spec_test <- mean(pred_test$class[group_test == 1] == 1)
  sens_test <- mean(pred_test$class[group_test == 2] == 2)
  
  # Imprimir os resultados
  cat("\n--- Métricas do Modelo LDA (Tucker-3) ---\n")
  cat(" (Classe 1 = Adulterado, Classe 2 = Puro)\n")
  cat("\n--- Treinamento ---\n")
  cat("Acurácia:", ac_train, "\nSensibilidade (Puro):", sens_train, "\nEspecificidade (Adulterado):", spec_train, "\n\n")
  cat("--- Validação Cruzada (CV) ---\n")
  cat("Acurácia:", ac_cv, "\nSensibilidade (Puro):", sens_cv, "\nEspecificidade (Adulterado):", spec_cv, "\n\n")
  cat("--- Teste ---\n")
  cat("Acurácia:", ac_test, "\nSensibilidade (Puro):", sens_test, "\nEspecificidade (Adulterado):", spec_test, "\n\n")
  
  # Matrizes de Confusão
  confusion_train <- table(Previsto = pred_train$class, Real = group_train)
  confusion_test <- table(Previsto = pred_test$class, Real = group_test)
  
  print("Matriz de Confusão - Treino")
  print(confusion_train)
  
  print("Matriz de Confusão - Teste")
  print(confusion_test)
  
  # Visualização
  dev.new()
  boxplot(pred_train$x[, 1] ~ group_train, col = c("red", "blue"), 
          xlab = "Grupo (1=Adulterado, 2=Puro)", ylab = "LD1", main = "Boxplot (Treino) - Tucker LDA")
  
  dev.new()
  plot(pred_train$x[, 1], col = as.factor(group_train), pch = 19, 
       xlab = "Amostras", ylab = "LD1", main = "Scatter Plot - Treino")
  legend("topleft", legend = c("Adulterado (1)", "Puro (2)"), col = c(2, 1), pch = 19)
  
  dev.new()
  boxplot(pred_test$x[, 1] ~ group_test, col = c("red", "blue"), 
          xlab = "Grupo (1=Adulterado, 2=Puro)", ylab = "LD1", main = "Boxplot (Teste) - Tucker LDA")
  
  dev.new()
  plot(pred_test$x[, 1], col = as.factor(group_test), pch = 19, 
       xlab = "Amostras", ylab = "LD1", main = "Scatter Plot - Teste")
  legend("topleft", legend = c("Adulterado (1)", "Puro (2)"), col = c(2, 1), pch = 19)
}