####
# SCRIPT [6] ATUALIZADO: PARAFAC-LDA
#
# OBJETIVO:
# 1. Carregar dados processados ('eem_cube_cleaned.mat') e metadados ('eem_data.mat').
# 2. Gerar rótulos (Y) dinamicamente a partir das concentrações e nomes.
# 3. (Opcional) Realizar balanceamento de classe robusto usando nomes.
# 4. Construir modelo PARAFAC e classificador LDA.
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
  
  # Pacotes do script original PARAFAC-LDA.R]
  required_packages <- c(
    "multiway", "ThreeWay", "R.matlab", "plot3D", "plotly", 
    "prospectr", "MASS", "caret", "rstudioapi", "abind", "parallel"
  )
  
  load_required_packages(required_packages)
}

# --- PASSO 2: Carregar Arquivos de Dados (Novo Fluxo) ---
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

# --- PASSO 3: Geração de Rótulos e Balanceamento de Dados (Lógica ATUALIZADA) ---
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
  
  # 2. Balanceamento de Classe (Recriando a lógica do script [6] PARAFAC-LDA.R] de forma robusta)
  # O objetivo original era criar um cubo balanceado com:
  # 10 amostras da Marca A, 10 da B, 10 da C, 10 da D (Total 40 Adulteradas)
  # 40 amostras Puras (Total 40 Puras)
  
  set.seed(42) # Para reprodutibilidade da amostragem
  
  # Encontra os ÍNDICES de cada grupo (independente da ordem)
  indices_A <- which(Y_multiclasse == "Marca A")
  indices_B <- which(Y_multiclasse == "Marca B")
  indices_C <- which(Y_multiclasse == "Marca C")
  indices_D <- which(Y_multiclasse == "Marca D")
  indices_Puras <- which(Y_multiclasse == "Puro")
  
  # 2. Sortear 10 amostras de CADA estrato adulterante
  sel_indices_A <- sample(indices_A, size = 10, replace = FALSE)
  sel_indices_B <- sample(indices_B, size = 10, replace = FALSE)
  sel_indices_C <- sample(indices_C, size = 10, replace = FALSE)
  sel_indices_D <- sample(indices_D, size = 10, replace = FALSE)
  
  # 3. Combinar os índices selecionados
  indices_balanceados <- c(sel_indices_A, sel_indices_B, sel_indices_C, sel_indices_D, 
                           indices_Puras)
  
  # 4. Criar o cubo de dados final e balanceado
  X <- X_original[indices_balanceados, , ]
  
  # 5. Criar os rótulos finais que correspondem ao cubo balanceado
  Y <- Y_binario[indices_balanceados]
  Y_multi_balanceado <- Y_multiclasse[indices_balanceados]
  nomes_amostras <- sample_names[indices_balanceados]
  
  cat(paste("\nCubo de dados final balanceado criado com dimensões:", 
            paste(dim(X), collapse = " x "), "\n"))
  print(table(Y)) # Deve mostrar 40 Adulterado, 40 Puro
  
  # 6. Acessar amostras por classe (para plots)
  # (Atualizado para usar os novos rótulos Y)
  X1 <- X[Y == "Adulterado", , ]
  X2 <- X[Y == "Puro", , ]
  
  X1m <- colMeans(X1)
  X2m <- colMeans(X2)
}

# --- PASSO 4: Visualização (Sem alteração de lógica) ---

  dev.new()
  filled.contour(x = nmEM, y = nmEX, z = X1m, 
                 color.palette = topo.colors, main = "Média Classe 1 (Adulterado)",
                 xlab = "Emissão (nm)", ylab = "Excitação (nm)")
  
  dev.new()
  filled.contour(x = nmEM, y = nmEX, z = X2m, 
                 color.palette = terrain.colors, main = "Média Classe 2 (Puro)",
                 xlab = "Emissão (nm)", ylab = "Excitação (nm)")
  
  dev.new()
  matplot(nmEM, X1m, type="l", xlab = "Comprimento de Onda Emissão (nm)", ylab = "Intensidade", 
          main = "Classe 1 - Normal (Modo Emissão)")
  
  dev.new()
  matplot(nmEX, t(X1m), type="l", xlab = "Comprimento de Onda Excitação (nm)", ylab = "Intensidade", 
          main = "Classe 1 - Normal (Modo Excitação)")

# --- PASSO 5: Modelo PARAFAC (Sem alteração de lógica) ---
{
  set.seed(24)
  num_factors <- 7  # Define number of factors
  
  
  cat(paste("Modelando PARAFAC com", num_factors))
 
  model <- parafac(X, nfac=num_factors, nstart=75, maxit=5000, 
                   ctol=1e-6, parallel=F, output = c("best", "all"))
  
  cat(paste("Modelo PARAFAC R²:", model$Rsq, "\n"))
}

# --- PASSO 6: Diagnóstico de Outliers (RSS e Q vs T²) ---
# (Sem alteração de lógica, já usa 'nomes_amostras' e 'Y' que definimos)
{
  # --- Erro de reconstrução por amostra (RSS) ---
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
  plot(rss, pch = 19, ylab = "RSS (erro de reconstrução)", xlab = "Índice da Amostra",
       main = "Erro de reconstrução por amostra (PARAFAC)")
  abline(h = thr, col = "red", lty = 2)
  text(outliers_rss, rss[outliers_rss], labels = nomes_amostras[outliers_rss], pos = 3, col = "red", cex=0.8)
  cat("Outliers por RSS (threshold median+3*IQR):", outliers_rss, "\n")
  
  # --- T² Hotelling vs Q residual ---
  A <- model$A
  S_A <- cov(A)
  T2 <- mahalanobis(A, colMeans(A), S_A)
  
  nfac <- ncol(A)
  n_amostras <- nrow(A)
  alpha <- 0.05
  T2_lim <- nfac * (n_amostras - 1) / (n_amostras - nfac) * qf(1 - alpha, nfac, n_amostras - nfac)
  
  Q <- rss # Q residual é o mesmo que o RSS calculado acima
  Q_lim <- median(Q) + 3 * IQR(Q)
  
  dev.new()
  plot(T2, Q, pch = 19, col = ifelse(Y == "Adulterado", "red", "blue"),
       xlab = expression(Hotelling~T^2),
       ylab = expression(Q~"(Erro de reconstrução)"),
       main = expression("Gráfico Q vs T"^2~"(Diagnóstico de outliers)"))
  abline(v = T2_lim, col = "darkgreen", lty = 2)
  abline(h = Q_lim,  col = "darkgreen", lty = 2)
  text(T2, Q, labels = nomes_amostras, pos = 3, cex = 0.6)
  legend("topright",
         legend = c("Adulterado", "Puro", "Limite 95%"),
         col = c("red", "blue", "darkgreen"), pch = c(19,19,NA), lty = c(NA,NA,2))
  
  out_QT <- which(T2 > T2_lim | Q > Q_lim)
  cat("Amostras identificadas como outliers (Q–T²):", nomes_amostras[out_QT], "\n")
}

  # --- PASSO 6.5 (opcional): Remoção de Outliers por NOME e Re-Modelagem ---
  {
    # 1. DEFINA AQUI OS NOMES EXATOS PARA REMOVER
    # Olhe no gráfico do Passo 6 (os rótulos vermelhos/azuis).
    # Exemplo: Se a amostra "P C 0 5" e "Marca A 0 2" forem outliers:
    
    nomes_para_remover <- c("B0900", "D0900", "PB02", "PD02", "PD03") 
    # Exemplo de uso: nomes_para_remover <- c("P C 0 5", "P B 1 0")
    
    # --- Lógica de Busca (Decodificação Nome -> Índice) ---
    if (length(nomes_para_remover) > 0) {
      
      # Identifica quais índices (linhas) correspondem aos nomes digitados
      # O operador %in% verifica se o nome da amostra está na sua lista de remoção
      outliers_indices <- which(nomes_amostras %in% nomes_para_remover)
      
      # Verificação de segurança: Avise se algum nome foi digitado errado
      nomes_encontrados <- nomes_amostras[outliers_indices]
      nomes_errados <- setdiff(nomes_para_remover, nomes_encontrados)
      
      if (length(nomes_errados) > 0) {
        cat("\n[ATENÇÃO] Os seguintes nomes não foram encontrados no conjunto de dados (verifique espaços ou maiúsculas):\n")
        print(nomes_errados)
      }
      
      if (length(outliers_indices) > 0) {
        cat("\nRemovendo", length(outliers_indices), "amostras identificadas:\n")
        print(nomes_encontrados)
        
        # 2. REMOVER DO CUBO DE DADOS (X)
        # Usamos os índices descobertos para cortar o cubo
        X <- X[-outliers_indices, , ]
        
        # 3. ATUALIZAR VETORES AUXILIARES
        # É crucial remover também do vetor Y e do vetor de nomes
        Y <- Y[-outliers_indices]
        # Se você usou balanceamento, pode ter outras variáveis de classe para atualizar, 
        # mas Y e nomes_amostras são os essenciais para o código seguir.
        nomes_amostras <- nomes_amostras[-outliers_indices]
        
        cat("Novas dimensões do cubo:", paste(dim(X), collapse=" x "), "\n")
        
        # 4. RECALCULAR O MODELO PARAFAC
        # O modelo precisa ser refeito pois os Loadings mudam sem os outliers
        cat("Recalculando modelo PARAFAC sem os outliers...\n")
        set.seed(24) 
        
        model <- parafac(X, nfac=num_factors, nstart=75, maxit=5000, 
                         ctol=1e-6, parallel=FALSE, output = c("best", "all"))
        
        cat(paste("Novo R² do Modelo Limpo:", model$Rsq, "\n"))
        
      } else {
        cat("Nenhum nome válido foi encontrado para remoção.\n")
      }
      
    } else {
      cat("Lista de remoção vazia. O modelo original será mantido.\n")
    }
  }
  
# --- PASSO 7: Visualização de Scores e Loadings ---
# (Lógica sem alterações, usa 'nmEM', 'nmEX', 'Y', 'nomes_amostras' que definimos)
{
  ## Plot Parafac Scores 1x2
  dev.new()
  fator_x <- 1 
  fator_y <- 2 
  
  plot(model$A[, fator_x], model$A[, fator_y], 
       type = "n",
       xlab = paste("Fator", fator_x), 
       ylab = paste("Fator", fator_y), 
       main = "Scores PARAFAC com Nomes das Amostras")
  
  text(x = model$A[Y == "Adulterado", fator_x], 
       y = model$A[Y == "Adulterado", fator_y], 
       labels = nomes_amostras[Y == "Adulterado"], 
       col = "red", cex = 0.7)
  
  text(x = model$A[Y == "Puro", fator_x], 
       y = model$A[Y == "Puro", fator_y], 
       labels = nomes_amostras[Y == "Puro"], 
       col = "blue", cex = 0.7)
  
  legend("topright", 
         legend = c("Adulterado", "Puro"), 
         col = c("red", "blue"), 
         pch = 19)
  grid()
  
  ## Plot Parafac Loadings
  dev.new()
  matplot(nmEM, model$B, type ="l", xlab = "Comprimento de Onda Emissão (nm)", ylab = "Intensidade")
  
  dev.new()
  matplot(nmEX, model$C, type ="l", xlab = "Comprimento de Onda Excitação (nm)", ylab = "Intensidade")
}

# --- PASSO 8: Separação Treino/Teste e Modelo LDA ---
# (Lógica atualizada para usar Y e classes "1" e "2" como no original)
{
  train_ratio <- 0.7 
  
  # Definir grupos numéricos (1=Adulterado, 2=Puro) para o LDA
  # (O script original usava 1 e 2 PARAFAC-LDA.R])
  group_numeric <- ifelse(Y == "Adulterado", 1, 2)
  
  # Scores das classes
  scores_X1 <- model$A[group_numeric == 1, ] # Scores Adulterados
  scores_X2 <- model$A[group_numeric == 2, ] # Scores Puros
  
  # Calcular amostras de treino por classe
  ntrain1 <- ceiling(train_ratio * nrow(scores_X1))
  ntrain2 <- ceiling(train_ratio * nrow(scores_X2))
  
  # Selecionar amostras
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
# (Lógica sem alterações, usa 'group_train' e 'group_test')
{
  # Predições
  pred_train <- predict(model_lda, train)
  pred_test <- predict(model_lda, test)
  
  # Métricas
  # (Usando a forma robusta do seu script original PARAFAC-LDA.R])
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
  cat("\n--- Métricas do Modelo LDA ---\n")
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
          xlab = "Grupo (1=Adulterado, 2=Puro)", ylab = "LD1", main = "Boxplot (Treino)")
  
  dev.new()
  plot(pred_train$x[, 1], col = as.factor(group_train), pch = 19, 
       xlab = "Amostras", ylab = "LD1", main = "Scatter Plot - Treino")
  legend("topleft", legend = c("Adulterado (1)", "Puro (2)"), col = c(2, 1), pch = 19) # Col 2=red, 1=blue
  
  dev.new()
  boxplot(pred_test$x[, 1] ~ group_test, col = c("red", "blue"), 
          xlab = "Grupo (1=Adulterado, 2=Puro)", ylab = "LD1", main = "Boxplot (Teste)")
  
  dev.new()
  plot(pred_test$x[, 1], col = as.factor(group_test), pch = 19, 
       xlab = "Amostras", ylab = "LD1", main = "Scatter Plot - Teste")
  legend("topleft", legend = c("Adulterado (1)", "Puro (2)"), col = c(2, 1), pch = 19)
}
  