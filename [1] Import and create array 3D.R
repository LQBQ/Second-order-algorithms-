####
# SCRIPT [1] IMPORTING, COMBINING AND CREATING 3D CUBE EEM
#
# OBJETIVO:
# 1. Ler arquivos .ASC brutos de um espectrofluorímetro.
  # 2. Combinar arquivos de excitação para formar matrizes 2D (EEMs) por amostra.
    # 3. Empilhar todas as amostras em um único array 3D (cubo).
      # 4. Mapear e associar um vetor de concentração.
        # 5. Salvar o cubo 3D, nomes das amostras e concentrações em um arquivo .mat
###########################################################################################
# --- STEP 1: Load necessary packages
{
  pacotes_necessarios <- c("dplyr", "rstudioapi", "R.matlab")
  
  for (pkg in pacotes_necessarios) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}
###########################################################################################
# --- STEP 2: Auxiliary Functions

#' process_sample: Reads a group of .ASC files and combines them into a 2D array.
process_sample <- function(sample_files, input_directory) {
  data_frames <- list()
  
  for (file in sample_files) {
    file_path <- file.path(input_directory, file) 
    lines <- readLines(file_path)
    lines <- gsub("([0-9]),([0-9])", "\\1.\\2", lines)
    
    temp_file <- tempfile(fileext = ".txt")
    writeLines(lines, temp_file)
    
    data <- read.delim(temp_file, header = FALSE, sep = ",", dec = ".")
    data <- data[, -1] 
    
    data_frames[[file]] <- data
    unlink(temp_file)
  }
  
 
  combined_data <- bind_cols(data_frames)
  return(as.matrix(combined_data))
}

#' rename_sample: Standardize the sample name.
rename_sample <- function(name) {
  if (grepl("^[ABCD][0-9%]+$", name, ignore.case = TRUE)) {
    
    prefix <- substr(name, 1, 1) 
    name_uppercase <- toupper(name)
    
    if (grepl("[ABCD]%[0-9]+$", name_uppercase)) {
      base_value <- 5 
    } else {
      num_str <- substr(name, 2, nchar(name)) 
      base_perc <- as.numeric(num_str) 
      base_value <- base_perc * 10
    }
    formatted_value <- sprintf("%04d", base_value) 
    formatted <- paste0(toupper(prefix), formatted_value)
    
  } else if (grepl("^P[ABCD][0-9]+$", name, ignore.case = TRUE)) {
    formatted <- toupper(name)
  } else {
    formatted <- name 
  }
  return(formatted)
}

############################################################################################
# --- STEP 3: Directory Selection 
cat("Por favor, selecione a pasta de ENTRADA contendo os arquivos .ASC.\n")
input_directory <- selectDirectory(caption = "SELECIONE A PASTA DE ENTRADA (.ASC Files)")
if (!nzchar(input_directory)) {
  stop("Seleção de diretório de entrada cancelada. Script terminado.")
}

cat("Por favor, selecione a pasta de SAÍDA onde o arquivo .mat será salvo.\n")
output_directory <- selectDirectory(caption = "SELECIONE A PASTA DE SAÍDA (para .mat)")
if (!nzchar(output_directory)) {
  stop("Seleção de diretório de saída cancelada. Script terminado.")
}
setwd(output_directory) 

# --- STEP 4: Define Import Parameters ---
files_per_sample <- 13  # Adjust as needed (e.g., 13 .ASC files per EEM)
cat(paste("Usando", files_per_sample, "Files per samples.\n"))

# Defining the axes 
# Adjust the seq() values ​​if your scans have changed
nmEM_vector <- seq(300, 800, by = 1)
nmEX_vector <- seq(310, 430, by = 10) 

# --- STEP 5: Processing .ASC Files
files <- list.files(path = input_directory, pattern = "\\.ASC$", ignore.case = TRUE)

if (length(files) == 0) {
  stop("Nenhum arquivo .ASC encontrado no diretório de entrada.")
}
if (length(files) %% files_per_sample != 0) {
  stop(paste("Erro: O número total de arquivos .ASC (", length(files), 
             ") não é divisível pelo número de arquivos por amostra (", files_per_sample, ")."))
}

cat(paste("Encontrados", length(files), "arquivos .ASC, correspondendo a", 
          length(files) / files_per_sample, "amostras.\n"))

eem_matrix_list <- list()

for (i in seq(1, length(files), by = files_per_sample)) {
  sample_files <- files[i:(i + files_per_sample - 1)]
  
  raw_name <- sub("^((?:[ABCD][0-9%]+)|(?:P[ABCD][0-9]+)).*", "\\1", sample_files[1], ignore.case = TRUE)
  sample_name <- rename_sample(raw_name)
  
  cat(paste("Processando:", sample_name, "...\n"))
  
  eem_matrix_list[[sample_name]] <- process_sample(sample_files, input_directory)
}

# --- PASSO 6: Creating a 3D Cube

# Ensures the alphabetical order of the samples.
sample_names_ordered <- sort(names(eem_matrix_list))
eem_matrix_list_sorted <- eem_matrix_list[sample_names_ordered]

dimensions <- sapply(eem_matrix_list_sorted, dim)
if (!all(apply(dimensions, 1, function(x) length(unique(x)) == 1))) {
  dims_df <- data.frame(Sample = names(eem_matrix_list_sorted), 
                        Rows_Em = dimensions[1,], 
                        Cols_Ex = dimensions[2,])
  print(dims_df)
  stop("As matrizes possuem dimensões inconsistentes (ver tabela acima).")
}

# Creates the 3D array.
dim_em <- dim(eem_matrix_list_sorted[[1]])[1]
dim_ex <- dim(eem_matrix_list_sorted[[1]])[2]
n_samples <- length(eem_matrix_list_sorted)

# Validate the dimensions using the wavelength vectors.
if (dim_em != length(nmEM_vector) || dim_ex != length(nmEX_vector)) {
  stop(paste("Dimensões da matriz (", dim_em, "x", dim_ex, 
             ") não correspondem aos vetores de comprimento de onda (", 
             length(nmEM_vector), "x", length(nmEX_vector), ")."))
}
               
array_3d_temp <- array(
  data = unlist(eem_matrix_list_sorted),
  dim = c(dim_em, dim_ex, n_samples)
)
               
eem_cube <- aperm(array_3d_temp, perm = c(3, 1, 2))

# Assigns names to the dimensions.
dimnames(eem_cube) <- list(
  Amostras = sample_names_ordered,
  Emissao = nmEM_vector,
  Excitacao = nmEX_vector
)

cat("Cubo 3D 'eem_cube' criado com sucesso.\n")
cat(paste("Dimensões Finais: [", paste(dim(eem_cube), collapse = " x "), "]\n"))


# --- STEP 7: Defining the Concentration Vector
####!!!--- Define concentrations here.---!!!####

  # - STEP A: Define the concentration vectors.
{
Y_A <- c(0.56, 1.1, 3.1, 5.1, 7.1, 9.9, 19.7, 30.0, 42.2, 50.0, 60.3, 72.0, 80.1, 90.9)
Y_B <- c(0.73, 1.0, 3.3, 5.0, 6.7, 11.5, 20.0, 30.4, 40.0, 50.2, 60.0, 70.2, 80.2, 89.3)
Y_C <- c(0.47, 1.0, 2.9, 5.0, 7.1, 9.8, 20.3, 33.3, 40.8, 50.1, 59.9, 70.0, 80.3, 90.6)
Y_D <- c(1.2, 1.2, 3.0, 5.2, 7.5, 9.8, 20.0, 29.9, 41.5, 50.0, 61.0, 69.8, 79.9, 90.3)

  # - STEP B: Define the corresponding SAMPLE NAMES (after renaming)
nomes_A <- c("A0005", "A0010", "A0030", "A0050", "A0070", "A0100", "A0200", 
             "A0300", "A0400", "A0500", "A0600", "A0700", "A0800", "A0900")

nomes_B <- c("B0005", "B0010", "B0030", "B0050", "B0070", "B0100", "B0200", 
             "B0300", "B0400", "B0500", "B0600", "B0700", "B0800", "B0900")

nomes_C <- c("C0005", "C0010", "C0030", "C0050", "C0070", "C0100", "C0200", 
             "C0300", "C0400", "C0500", "C0600", "C0700", "C0800", "C0900")

nomes_D <- c("D0005", "D0010", "D0030", "D0050", "D0070", "D0100", "D0200", 
             "D0300", "D0400", "D0500", "D0600", "D0700", "D0800", "D0900")

  # - STEP C: Define the names and values ​​of the PURE samples
    # (Assuming 40 pure samples, 10 of each brand PA, PB, PC, PD)
nomes_Puras <- c(
  paste0("PA", sprintf("%02d", 1:10)), # PA01, PA02... PA10
  paste0("PB", sprintf("%02d", 1:10)), # PB01, PB02... PB10
  paste0("PC", sprintf("%02d", 1:10)), # PC01, PC02... PC10
  paste0("PD", sprintf("%02d", 1:10))  # PD01, PD02... PD10
)
valores_Puras <- as.list(rep(0.0, 40))
mapa_Puras <- setNames(valores_Puras, nomes_Puras)

  # - STEP D: Create the maps for the adulterated ones.
mapa_A <- setNames(as.list(Y_A), nomes_A)
mapa_B <- setNames(as.list(Y_B), nomes_B)
mapa_C <- setNames(as.list(Y_C), nomes_C)
mapa_D <- setNames(as.list(Y_D), nomes_D)

  # - STEP E: Combine EVERYTHING into a single map
concentration_map <- c(mapa_A, mapa_B, mapa_C, mapa_D, mapa_Puras)
}

# ------------------------------------------------
# WARNING: If 'concentration_map' is not defined, the script will create a
# vector of NAs (missing values) as a placeholder.
if (!exists("concentration_map")) {
  cat("AVISO: 'concentration_map' não foi definido. Criando vetor de concentração com NAs.\n")
  cat("Por favor, edite este script na seção 'DEFINIR CONCENTRAÇÕES AQUI' para resultados corretos.\n")
  concentration_map <- setNames(as.list(rep(NA, length(sample_names_ordered))), sample_names_ordered)
}

# Creates the final concentration vector, ensuring the same order as the 3D cube.
concentration_vector <- unlist(concentration_map[sample_names_ordered])
               
if (any(is.na(concentration_vector))) {
  cat("AVISO: O vetor de concentração final contém valores NA (ausentes).\n")
  cat("Amostras com NA:\n")
  print(sample_names_ordered[is.na(concentration_vector)])
}

if (length(concentration_vector) != n_samples) {
  stop("Erro: O vetor de concentração final não tem o mesmo número de amostras que o cubo EEM.")
}


# --- STEP 8: Save file .mat final ---
{
output_filename <- "eem_data.mat"
output_filepath <- file.path(output_directory, output_filename)

writeMat(output_filepath, 
         cube = eem_cube, 
         concentrations = concentration_vector,
         sample_names = sample_names_ordered,
         nm_emission = nmEM_vector,
         nm_excitation = nmEX_vector
)


cat(paste("\n--- SUCESSO! ---\n"))
cat(paste("Arquivo .mat salvo em:", output_filepath, "\n"))
cat("O arquivo contém:\n")
cat(" 1. 'cube': O array 3D [", paste(dim(eem_cube), collapse = " x "), "]\n")
cat(" 2. 'concentrations': Vetor numérico de concentrações (", length(concentration_vector), " valores)\n")
cat(" 3. 'sample_names': Vetor de nomes das amostras (", length(sample_names_ordered), " nomes)\n")
cat(" 4. 'nm_emission': Vetor de comprimentos de onda de emissão\n")
cat(" 5. 'nm_excitation': Vetor de comprimentos de onda de excitação\n")
}

