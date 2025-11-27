####
# SCRIPT [4] ATUALIZADO (v2): Pré-processamento EEM Completo
#
# NOME: [4] preprocessing_EEM.R
#
# OBJETIVO:
# 1. Carregar o cubo 3D bruto ('eem_data.mat').
# 2. (ETAPA A) Aplicar a remoção de espalhamento (Rayleigh, Raman).
# 3. (ETAPA B) Aplicar correção de baseline.
# 4. (ETAPA C) Aplicar suavização Savitzky-Golay.
# 5. (NOVO) Oferecer visualização imediata do resultado.
# 6. Salvar o cubo 3D final e limpo ('eem_cube_cleaned.mat').
####

# --- PASSO 1: Instalar e Carregar Pacotes Necessários ---
{
  load_required_packages <- function(packages) {
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
  }
  
  # Adicionamos 'plotly' (para a função de plot 3D) e 'baseline'
  required_packages <- c("zoo", "signal", "R.matlab", "rstudioapi", 
                         "baseline", "plotly")
  
  load_required_packages(required_packages)
}

# --- PASSO 2: Função de Remoção de Espalhamento (ETAPA A) ---
# (Função 'removeScatteringR' original preprocessing.R])
removeScatteringR <- function(eem_matrix, exc_wavelengths, em_wavelengths, 
                              Rayleigh_1st, Rayleigh_2nd, Rayleigh_3rd, 
                              Raman_1st, Raman_2nd, raman_offset) {
  
  if (length(dim(eem_matrix)) != 3) {
    stop('eem_matrix must be a 3D array: [samples x emissions x excitations].')
  }
  dimensions <- dim(eem_matrix)
  nSamples <- dimensions[1]
  nEm <- dimensions[2]
  nEx <- dimensions[3]
  if (nEx != length(exc_wavelengths) || nEm != length(em_wavelengths)) {
    stop("Matrix dimensions must match wavelength vector lengths (Emissions: dim2, Excitations: dim3).")
  }
  
  eem_cleaned <- eem_matrix
  
  # ... (O código de remoção de espalhamento e interpolação PCHIP preprocessing.R] continua aqui) ...
  
  # --- 2. Scattering Removal (Masking with NaN) ---
  for (sampleIdx in 1:nSamples) {
    for (i in 1:length(exc_wavelengths)) {
      exc_wl <- exc_wavelengths[i]
      
      rayleigh1_lower = exc_wl - Rayleigh_1st
      rayleigh1_upper = exc_wl + Rayleigh_1st
      rayleigh2_lower = 2 * exc_wl - Rayleigh_2nd
      rayleigh2_upper = 2 * exc_wl + Rayleigh_2nd
      rayleigh3_lower = 3 * exc_wl - Rayleigh_3rd
      rayleigh3_upper = 3 * exc_wl + Rayleigh_3rd
      
      raman1_lower = exc_wl + raman_offset - Raman_1st
      raman1_upper = exc_wl + raman_offset + Raman_1st
      raman2_lower = exc_wl + 2 * raman_offset - Raman_2nd
      raman2_upper = 2 * raman_offset + Raman_2nd
      
      idx_rayleigh1 = em_wavelengths >= rayleigh1_lower & em_wavelengths <= rayleigh1_upper
      idx_rayleigh2 = em_wavelengths >= rayleigh2_lower & em_wavelengths <= rayleigh2_upper
      idx_rayleigh3 = em_wavelengths >= rayleigh3_lower & em_wavelengths <= rayleigh3_upper
      idx_raman1 = em_wavelengths >= raman1_lower & em_wavelengths <= raman1_upper
      idx_raman2 = em_wavelengths >= raman2_lower & em_wavelengths <= raman2_upper
      
      eem_cleaned[sampleIdx, idx_rayleigh1, i] = NaN
      eem_cleaned[sampleIdx, idx_rayleigh2, i] = NaN
      eem_cleaned[sampleIdx, idx_rayleigh3, i] = NaN
      eem_cleaned[sampleIdx, idx_raman1, i] = NaN
      eem_cleaned[sampleIdx, idx_raman2, i] = NaN
    }
  }
  
  # --- 3. Interpolation (PCHIP) by Sample and Excitation (Column-wise) ---
  cat("ETAPA A: Iniciando remoção de espalhamento e interpolação PCHIP...\n")
  pb <- txtProgressBar(min = 0, max = nSamples, style = 3) 
  
  for (sampleIdx in 1:nSamples) {
    for (ex in 1:nEx) {
      y <- eem_cleaned[sampleIdx, , ex]
      nan_idx <- is.na(y)
      if (!any(nan_idx)) next 
      
      W <- !nan_idx 
      valid_idx <- which(W)
      
      if (length(valid_idx) < 2) { 
        next
      }
      
      x <- em_wavelengths
      y_interp <- y
      x_known <- x[valid_idx]
      y_known <- y[valid_idx]
      x_fill  <- x[nan_idx]
      
      x_fill <- x_fill[x_fill >= min(x_known) & x_fill <= max(x_known)]
      
      if (length(x_fill) > 0) {
        y_fill <- signal::interp1(x = x_known, y = y_known, xi = x_fill, method = "pchip")
        y_interp[match(x_fill, x)] <- y_fill
      }
      
      first_valid_idx <- valid_idx[1]
      last_valid_idx <- valid_idx[length(valid_idx)]
      
      if (first_valid_idx > 1) {
        y_interp[1:(first_valid_idx - 1)] <- y_interp[first_valid_idx] 
      }
      if (last_valid_idx < nEm) {
        y_interp[(last_valid_idx + 1):nEm] <- y_interp[last_valid_idx] 
      }
      
      eem_cleaned[sampleIdx, , ex] <- y_interp
    }
    setTxtProgressBar(pb, sampleIdx) 
  }
  close(pb)
  cat("\nETAPA A: Interpolação PCHIP concluída.\n")
  return(eem_cleaned)
}

# --- PASSO 3: Função para Baseline e Suavização (ETAPAS B e C) ---
# (Função 'apply_baseline_and_smoothing' que criamos anteriormente)
apply_baseline_and_smoothing <- function(eem_cube, 
                                         baseline_method = 'modpoly', 
                                         baseline_degree = 4,
                                         sgolay_p = 3, 
                                         sgolay_n = 11) {
  
  if (sgolay_n %% 2 == 0) {
    stop(paste("O tamanho da janela S-G (sgolay_n =", sgolay_n, ") deve ser ímpar."))
  }
  
  eem_final <- eem_cube
  
  dimensions <- dim(eem_cube)
  nSamples <- dimensions[1]
  nEm <- dimensions[2]
  nEx <- dimensions[3]
  
  cat("ETAPA B/C: Iniciando correção de baseline e suavização S-G...\n")
  pb <- txtProgressBar(min = 0, max = nSamples * nEx, style = 3) 
  
  for (i in 1:nSamples) {
    for (j in 1:nEx) {
      
      spectrum <- eem_cube[i, , j]
      spec_matrix <- matrix(spectrum, ncol = 1)
      
      # ETAPA B: Correção de Baseline
      bl_correction <- try(baseline::baseline(spec_matrix, 
                                              method = baseline_method, 
                                              degree = baseline_degree), 
                           silent = TRUE)
      
      if (inherits(bl_correction, "try-error")) {
        spectrum_corrected <- spectrum 
      } else {
        spectrum_corrected <- baseline::getCorrected(bl_correction)
      }
      
      # ETAPA C: Suavização Savitzky-Golay
      spectrum_smoothed <- signal::sgolayfilt(spectrum_corrected, 
                                              p = sgolay_p, 
                                              n = sgolay_n)
      
      spectrum_smoothed[spectrum_smoothed < 0] <- 0
      
      eem_final[i, , j] <- spectrum_smoothed
      setTxtProgressBar(pb, (i - 1) * nEx + j)
    }
  }
  
  close(pb)
  cat("\nETAPA B/C: Baseline e suavização concluídos.\n")
  return(eem_final)
}

# --- PASSO 4: Funções de Visualização (COPIADAS DO SCRIPT [3]) ---

#' Função 1: Visualiza a superfície 3D (plotly)
visualize_3D_surface <- function(eem_matrix, em_wavelengths, exc_wavelengths, sample_name = NULL) {
  if (is.null(eem_matrix) || nrow(eem_matrix) == 0 || ncol(eem_matrix) == 0) {
    stop("A matriz de dados deve ser não-nula e ter dimensões maiores que zero.")
  }
  title_text <- ifelse(is.null(sample_name),
                       "Superfície 3D da Amostra EEM",
                       paste("Superfície 3D - Amostra:", sample_name))
  
  plot_ly() %>%
    add_surface(z = ~eem_matrix, 
                x = ~exc_wavelengths, 
                y = ~em_wavelengths, 
                colorscale = 'Viridis') %>%
    layout(
      title = title_text,
      scene = list(
        xaxis = list(title = "Excitação (nm)"),
        yaxis = list(title = "Emissão (nm)"),
        zaxis = list(title = "Intensidade")
      )
    )
}

#' Função 2: Visualiza o contorno 2D (filled.contour)
visualize_2D_contour <- function(eem_matrix, em_wavelengths, exc_wavelengths, sample_name = NULL) {
  if (is.null(eem_matrix) || nrow(eem_matrix) == 0 || ncol(eem_matrix) == 0) {
    stop("A matriz de dados deve ser não-nula e ter dimensões maiores que zero.")
  }
  title_text <- ifelse(is.null(sample_name),
                       "Contorno 2D (Mapa de Calor) da Amostra EEM",
                       paste("Contorno 2D - Amostra:", sample_name))
  
  filled.contour(
    x = em_wavelengths, 
    y = exc_wavelengths, 
    z = eem_matrix, # Transpõe a matriz para Z[x, y]
    color.palette = topo.colors,
    main = title_text,
    xlab = "Emission (nm)",
    ylab = "Excitation (nm)"
  )
}

# --- PASSO 5: Carregamento e Execução do Pré-processamento ---

# 1. Defina seus parâmetros de pré-processamento (AJUSTE CONFORME NECESSÁRIO)
params_scatter <- list(
  Rayleigh_1st = 15,  
  Rayleigh_2nd = 15, 
  Rayleigh_3rd = 15, 
  Raman_1st = 0,     
  Raman_2nd = 0,    
  raman_offset = 0  
)
params_smooth <- list(
  baseline_method = 'modpoly', 
  baseline_degree = 3,       
  sgolay_p = 3,              
  sgolay_n = 11              
)

# 2. Selecione e carregue o arquivo .mat bruto
cat("Por favor, selecione o arquivo 'eem_data.mat' bruto...\n")
mat_file_path <- rstudioapi::selectFile(caption = "Selecione o arquivo eem_data.mat",
                                        filter = "MAT Files (*.mat)")

if (!nzchar(mat_file_path)) {
  stop("Nenhum arquivo selecionado. Script terminado.")
}

eem_data <- readMat(mat_file_path)
names(eem_data) <- gsub("\\.mat$", "", names(eem_data))

eem_cube_raw <- eem_data$cube
em_wavelengths <- as.vector(eem_data$nm.emission)
exc_wavelengths <- as.vector(eem_data$nm.excitation)
concentrations <- as.vector(eem_data$concentrations)

# Carrega os nomes para usar nos títulos dos gráficos
sample_names <- unlist(eem_data$sample.names) 

cat(paste("Dados brutos carregados. Dimensões:", 
          paste(dim(eem_cube_raw), collapse = " x "), "\n"))

# 3. ETAPA A: Execute a remoção de espalhamento
eem_cube_scatter_removed <- removeScatteringR(
  eem_matrix = eem_cube_raw, 
  exc_wavelengths = exc_wavelengths, 
  em_wavelengths = em_wavelengths, 
  Rayleigh_1st = params_scatter$Rayleigh_1st, 
  Rayleigh_2nd = params_scatter$Rayleigh_2nd, 
  Rayleigh_3rd = params_scatter$Rayleigh_3rd, 
  Raman_1st = params_scatter$Raman_1st, 
  Raman_2nd = params_scatter$Raman_2nd, 
  raman_offset = params_scatter$raman_offset
)

# 4. ETAPA B/C: Execute a correção de baseline e suavização
eem_cube_cleaned <- apply_baseline_and_smoothing(
  eem_cube = eem_cube_scatter_removed,
  baseline_method = params_smooth$baseline_method,
  baseline_degree = params_smooth$baseline_degree,
  sgolay_p = params_smooth$sgolay_p,
  sgolay_n = params_smooth$sgolay_n
)

cat("\n--- Pré-processamento completo --- \n")

# --- PASSO 6: Verificação Visual (Opcional) ---
{
  prompt <- "Deseja visualizar uma amostra pré-processada antes de salvar? (s/n): "
  resposta <- tolower(trimws(readline(prompt)))
  
  if (resposta == "s" || resposta == "sim") {
    
    # Pergunta pelo índice
    n_samples_total <- dim(eem_cube_cleaned)[1]
    prompt_idx <- paste0("Digite o NÚMERO (índice) da amostra que deseja ver (1 a ", n_samples_total, "): ")
    idx_ver <- as.integer(readline(prompt_idx))
    
    # Valida o índice
    if (!is.na(idx_ver) && idx_ver >= 1 && idx_ver <= n_samples_total) {
      
      # Extrai os dados para plotar
      amostra_para_plotar <- eem_cube_cleaned[idx_ver, , ]
      nome_da_amostra <- sample_names[idx_ver]
      titulo_plot <- paste(nome_da_amostra, "(Pré-processado)")
      
      cat(paste("Gerando gráficos para a amostra:", nome_da_amostra, "\n"))
      
      # Gera Gráfico 3D (interativo)
      print(visualize_3D_surface(
        eem_matrix = amostra_para_plotar,
        em_wavelengths = em_wavelengths,
        exc_wavelengths = exc_wavelengths,
        sample_name = titulo_plot
      ))
      
      # Gera Gráfico 2D (estático)
      cat("Pressione [Enter] no console para ver o gráfico 2D...\n")
      readline() # Pausa para o usuário ver o gráfico 3D
      
      dev.new() # Abre uma nova janela gráfica para o contorno
      visualize_2D_contour(
        eem_matrix = amostra_para_plotar,
        em_wavelengths = em_wavelengths,
        exc_wavelengths = exc_wavelengths,
        sample_name = titulo_plot
      )
      
    } else {
      cat("Índice inválido ou não numérico. Pulando visualização.\n")
    }
  } else {
    cat("Pulando visualização.\n")
  }
}

# --- PASSO 7: Salvar o Cubo Limpo ---
cat("\nPronto para salvar. Por favor, selecione a pasta de SAÍDA para salvar o arquivo limpo.\n")
output_directory <- selectDirectory(caption = "SELECIONE A PASTA DE SAÍDA")

if (nzchar(output_directory)) {
  output_filename <- "eem_cube_cleaned.mat"
  output_filepath <- file.path(output_directory, output_filename)
  
  # Salva o cubo limpo E os vetores de comprimento de onda
  writeMat(output_filepath, 
           cube_cleaned = eem_cube_cleaned,
           nm_emission = em_wavelengths,
           nm_excitation = exc_wavelengths,
           concentrations = concentrations,
           sample_names = sample_names
           
  )
  cat(paste("Cubo de dados limpo salvo com sucesso em:", output_filepath, "\n"))
  
} else {
  cat("Seleção de diretório de saída cancelada. O arquivo limpo NÃO foi salvo.\n")
}
