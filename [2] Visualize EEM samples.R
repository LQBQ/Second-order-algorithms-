####
# SCRIPT [3] ATUALIZADO: Visualização de Amostras EEM (3D e Contorno 2D)
#
# OBJETIVO:
# 1. Carregar o cubo de dados .mat processado.
# 2. Fornecer funções para visualizar amostras EEM individuais.
#    - visualize_3D_surface: Gráfico de superfície 3D interativo.
#    - visualize_2D_contour: Gráfico de contorno 2D (mapa de calor) interativo.
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
  
  # plotly para gráficos interativos, R.matlab para carregar dados
  required_packages <- c("plotly", "rstudioapi", "R.matlab")
  
  load_required_packages(required_packages)
}

# --- PASSO 2: Funções de Visualização ---

#' Função 1: Visualiza a superfície 3D (como no seu script original)
#' @param eem_matrix Uma única matriz 2D [emissão x excitação].
#' @param em_wavelengths Vetor de comprimentos de onda de emissão (eixo Y).
#' @param exc_wavelengths Vetor de comprimentos de onda de excitação (eixo X).
#' @param sample_name (Opcional) Nome/ID da amostra para o título.
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

#' Função 2: Visualiza o contorno 2D (mapa de calor) (NOVA FUNÇÃO)
#' @param eem_matrix Uma única matriz 2D [emissão x excitação].
#' @param em_wavelengths Vetor de comprimentos de onda de emissão (eixo Y).
#' @param exc_wavelengths Vetor de comprimentos de onda de excitação (eixo X).
#' @param sample_name (Opcional) Nome/ID da amostra para o título.
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
    z = eem_matrix,
    color.palette = topo.colors, 
    main = title_text,
    xlab = "Emissão (nm)",
    ylab = "Excitação (nm)"
  )
}

# --- PASSO 3: Exemplo de Uso ---

# 1. Selecione o arquivo .mat que criamos no script consolidado
{
mat_file_path <- rstudioapi::selectFile(caption = "Selecione o arquivo eem_data.mat",
                                        filter = "MAT Files (*.mat)")

if (!nzchar(mat_file_path)) {
  stop("Nenhum arquivo selecionado. Script terminado.")
}
  }
# 2. Carregue os dados do arquivo .mat
eem_data <- readMat(mat_file_path)

# 3. Extraia os componentes (usando os nomes que definimos)
#    Nota: readMat pode adicionar '.mat' ao nome do objeto, removemos se existir
{
names(eem_data) <- gsub("\\.mat$", "", names(eem_data))
eem_cube_raw <- eem_data$cube
sample_names <- unlist(eem_data$sample.names) # Garante que seja um vetor
em_wavelengths <- as.vector(eem_data$nm.emission) # Garante que seja um vetor
exc_wavelengths <- as.vector(eem_data$nm.excitation) # Garante que seja um vetor
}

cat("Dados carregados com sucesso!\n")
cat(paste(length(sample_names), "amostras encontradas.\n"))


# 4. Escolha qual amostra visualizar
#    Pode ser por índice (ex: 1) ou por nome (ex: "D0005")

# --- Exemplo por Índice ---
{
Sample_ID_Index <- 10 # Mude para o índice da amostra que deseja ver
Sample_Name <- sample_names[Sample_ID_Index]
eem_to_plot <- eem_cube_raw[Sample_ID_Index, , ]
}

# --- Exemplo por Nome (descomente para usar) ---
{
 Sample_Name <- "PA05" # Mude para o nome da amostra que deseja ver
 Sample_ID_Index <- which(sample_names == Sample_Name)
 if(length(Sample_ID_Index) == 0) stop("Nome da amostra não encontrado.")
 eem_to_plot <- eem_cube_raw[Sample_ID_Index, , ]
}


# 5. Gere os gráficos para a amostra selecionada
cat(paste("Gerando gráficos para a amostra:", Sample_Name, "(Índice:", Sample_ID_Index, ")\n"))

# Gráfico 3D
grafico_3d <- visualize_3D_surface(
  eem_matrix = eem_to_plot,
  em_wavelengths = em_wavelengths,
  exc_wavelengths = exc_wavelengths,
  sample_name = Sample_Name
)
print(grafico_3d)

# Gráfico 2D (Mapa de Calor)
dev.new()
grafico_2d <- visualize_2D_contour(
  eem_matrix = eem_to_plot,
  em_wavelengths = em_wavelengths,
  exc_wavelengths = exc_wavelengths,
  sample_name = Sample_Name
)
