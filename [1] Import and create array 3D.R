#### NOTE: Before using this script, please download the sample .ASC files.
#### These files are available at the following repository:
#### https://github.com/LQBQ/Second-order-algorithms-.git

####
#  SCRIPT: [1] IMPORT, COMBINE, AND CREATE 3D EEM CUBE
#
# OBJECTIVE:
# 1. Read raw .ASC files from a spectrofluorometer.
# 2. Combine excitation files to form 2D matrices (EEMs) per sample.
# 3. Stack all samples into a single 3D array (cube).
# 4. Map and associate a concentration vector.
# 5. Save the 3D cube, sample names, and concentrations into a .mat file.
####

###########################################################################
# --- STEP 1: Load Required Packages ---
# (Installs if not present)
{
  pacotes_necessarios <- c("dplyr", "rstudioapi", "R.matlab")
  
  for (pkg in pacotes_necessarios) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

#############################################################################
# --- STEP 2: Helper Functions ---

#' process_sample: Reads a group of .ASC files and combines them into a 2D matrix.
process_sample <- function(sample_files, input_directory) {
  data_frames <- list()
  
  for (file in sample_files) {
    file_path <- file.path(input_directory, file) # Ensures correct path
    lines <- readLines(file_path)
    lines <- gsub("([0-9]),([0-9])", "\\1.\\2", lines)
    
    temp_file <- tempfile(fileext = ".txt")
    writeLines(lines, temp_file)
    
    data <- read.delim(temp_file, header = FALSE, sep = ",", dec = ".")
    data <- data[, -1] # Remove the first column (wavelength)
    
    data_frames[[file]] <- data
    unlink(temp_file)
  }
  
  # Combine all columns (excitations)
  combined_data <- bind_cols(data_frames)
  return(as.matrix(combined_data)) # Return as matrix
}

#' rename_sample: Standardizes the sample name.
rename_sample <- function(name) {
  # 1. Adulterated Samples: Start with A/B/C/D followed by numbers, may contain '%'
  if (grepl("^[ABCD][0-9%]+$", name, ignore.case = TRUE)) {
    
    prefix <- substr(name, 1, 1) # Type letter (e.g., D)
    name_uppercase <- toupper(name)
    
    if (grepl("[ABCD]%[0-9]+$", name_uppercase)) {
      base_value <- 5 # Special Case: D%05 (0.5%)
    } else {
      num_str <- substr(name, 2, nchar(name)) 
      base_perc <- as.numeric(num_str) 
      base_value <- base_perc * 10
    }
    formatted_value <- sprintf("%04d", base_value) # 4 digits, zero-padded integer
    formatted <- paste0(toupper(prefix), formatted_value)
    
  } else if (grepl("^P[ABCD][0-9]+$", name, ignore.case = TRUE)) {
    # 2. Pure Samples: PA01, PB02, etc. (Standardize case)
    formatted <- toupper(name)
  } else {
    formatted <- name # Fallback
  }
  return(formatted)
}

########################################################################
# --- STEP 3: Directory Selection ---
cat("Please select the INPUT folder containing the .ASC files.\n")
input_directory <- selectDirectory(caption = "SELECT INPUT FOLDER (.ASC Files)")
if (!nzchar(input_directory)) {
  stop("Input directory selection cancelled. Script terminated.")
}

cat("Please select the OUTPUT folder where the .mat file will be saved.\n")
output_directory <- selectDirectory(caption = "SELECT OUTPUT FOLDER (for .mat)")
if (!nzchar(output_directory)) {
  stop("Output directory selection cancelled. Script terminated.")
}
setwd(output_directory) # Sets the working directory for saving

################################################################################
# --- STEP 4: Define Import Parameters ---
files_per_sample <- 13  # Adjust as necessary (e.g., 13 .ASC files per EEM)
cat(paste("Using", files_per_sample, "files per sample.\n"))

# Definition of axes
# Adjust seq() values if your scan ranges changed
nmEM_vector <- seq(300, 800, by = 1) # 501 emission points
nmEX_vector <- seq(310, 430, by = 10) # 13 excitation points

################################################################################
# --- STEP 5: Process .ASC Files ---
files <- list.files(path = input_directory, pattern = "\\.ASC$", ignore.case = TRUE)

if (length(files) == 0) {
  stop("No .ASC files found in the input directory.")
}
if (length(files) %% files_per_sample != 0) {
  stop(paste("Error: The total number of .ASC files (", length(files), 
             ") is not divisible by the number of files per sample (", files_per_sample, ")."))
}

cat(paste("Found", length(files), ".ASC files, corresponding to", 
          length(files) / files_per_sample, "samples.\n"))

eem_matrix_list <- list() # List to store the 2D matrices

for (i in seq(1, length(files), by = files_per_sample)) {
  sample_files <- files[i:(i + files_per_sample - 1)]
  
  # Extract the base name of the first file in the group
  raw_name <- sub("^((?:[ABCD][0-9%]+)|(?:P[ABCD][0-9]+)).*", "\\1", sample_files[1], ignore.case = TRUE)
  sample_name <- rename_sample(raw_name)
  
  cat(paste("Processing:", sample_name, "...\n"))
  
  # Process the file group and store the 2D matrix in the list
  eem_matrix_list[[sample_name]] <- process_sample(sample_files, input_directory)
}

#############################################################################
# --- STEP 6: Create 3D Cube ---

# Ensure alphabetical order of samples
sample_names_ordered <- sort(names(eem_matrix_list))
eem_matrix_list_sorted <- eem_matrix_list[sample_names_ordered]

# Check dimensions
dimensions <- sapply(eem_matrix_list_sorted, dim)
if (!all(apply(dimensions, 1, function(x) length(unique(x)) == 1))) {
  dims_df <- data.frame(Sample = names(eem_matrix_list_sorted), 
                        Rows_Em = dimensions[1,], 
                        Cols_Ex = dimensions[2,])
  print(dims_df)
  stop("Matrices have inconsistent dimensions (see table above).")
}

# Create the 3D array
dim_em <- dim(eem_matrix_list_sorted[[1]])[1]
dim_ex <- dim(eem_matrix_list_sorted[[1]])[2]
n_samples <- length(eem_matrix_list_sorted)

# Validate dimensions with wavelength vectors
if (dim_em != length(nmEM_vector) || dim_ex != length(nmEX_vector)) {
  stop(paste("Matrix dimensions (", dim_em, "x", dim_ex, 
             ") do not match the wavelength vectors (", 
             length(nmEM_vector), "x", length(nmEX_vector), ")."))
}

# Unlist and fill the array
array_3d_temp <- array(
  data = unlist(eem_matrix_list_sorted),
  dim = c(dim_em, dim_ex, n_samples)
)

# Permute dimensions to [Samples x Emission x Excitation]
eem_cube <- aperm(array_3d_temp, perm = c(3, 1, 2))

# Assign names to dimensions
dimnames(eem_cube) <- list(
  Amostras = sample_names_ordered,
  Emissao = nmEM_vector,
  Excitacao = nmEX_vector
)

cat("3D Cube 'eem_cube' created successfully.\n")
cat(paste("Final Dimensions: [", paste(dim(eem_cube), collapse = " x "), "]\n"))

##############################################################################
# --- STEP 7: Define Concentration Vector ---
####!!!--- DEFINE CONCENTRATIONS HERE ---!!!####
#
  # STEP A: Define concentration vectors
{
Y_A <- c(0.60, 1.0, 3.1, 5.0, 6.9, 10.2, 20.0, 30.0, 42.2, 50.6, 60.3, 72.0, 80.1, 90.0)
Y_B <- c(0.50, 1.1, 3.2, 5.0, 7.1, 10.3, 20.1, 30.2, 40.0, 50.0, 60.0, 70.2, 80.2, 90.0)
Y_C <- c(0.47, 1.0, 2.9, 5.0, 7.1, 9.8, 20.3, 33.3, 40.8, 50.1, 59.9, 70.0, 80.3, 90.6)
Y_D <- c(1.2, 1.2, 3.0, 5.2, 7.5, 9.8, 20.0, 29.9, 41.5, 50.0, 61.0, 69.8, 79.9, 90.3)

  # STEP B: Define corresponding SAMPLE NAMES (after renaming)
#          (!!! THIS IS AN EXAMPLE - REPLACE WITH YOUR REAL NAMES !!!)
#          The order here must correspond to the order of values in Y_A, Y_B...
nomes_A <- c("A0005", "A0010", "A0030", "A0050", "A0070", "A0100", "A0200", 
             "A0300", "A0400", "A0500", "A0600", "A0700", "A0800", "A0900")

nomes_B <- c("B0005", "B0010", "B0030", "B0050", "B0070", "B0100", "B0200", 
             "B0300", "B0400", "B0500", "B0600", "B0700", "B0800", "B0900")

nomes_C <- c("C0005", "C0010", "C0030", "C0050", "C0070", "C0100", "C0200", 
             "C0300", "C0400", "C0500", "C0600", "C0700", "C0800", "C0900")

nomes_D <- c("D0005", "D0010", "D0030", "D0050", "D0070", "D0100", "D0200", 
             "D0300", "D0400", "D0500", "D0600", "D0700", "D0800", "D0900") # Watch out for duplicate names

  # STEP C: Define names and values for PURE samples
#          (Assuming 40 pure samples, 10 for each brand PA, PB, PC, PD)
nomes_Puras <- c(
  paste0("PA", sprintf("%02d", 1:10)), # PA01, PA02... PA10
  paste0("PB", sprintf("%02d", 1:10)), # PB01, PB02... PB10
  paste0("PC", sprintf("%02d", 1:10)), # PC01, PC02... PC10
  paste0("PD", sprintf("%02d", 1:10))  # PD01, PD02... PD10
)
valores_Puras <- as.list(rep(0.0, 40))
mapa_Puras <- setNames(valores_Puras, nomes_Puras)

  # STEP D: Create maps for adulterated samples
mapa_A <- setNames(as.list(Y_A), nomes_A)
mapa_B <- setNames(as.list(Y_B), nomes_B)
mapa_C <- setNames(as.list(Y_C), nomes_C)
mapa_D <- setNames(as.list(Y_D), nomes_D)

  # STEP E: Combine EVERYTHING into a single map
concentration_map <- c(mapa_A, mapa_B, mapa_C, mapa_D, mapa_Puras)
}

# ------------------------------------------------
# ATTENTION: If 'concentration_map' is not defined, the script will create a
# vector of NAs (missing values) as a placeholder.
if (!exists("concentration_map")) {
  cat("WARNING: 'concentration_map' was not defined. Creating concentration vector with NAs.\n")
  cat("Please edit this script in the 'DEFINE CONCENTRATIONS HERE' section for correct results.\n")
  concentration_map <- setNames(as.list(rep(NA, length(sample_names_ordered))), sample_names_ordered)
}

# Create the final concentration vector, ensuring the same order as the 3D cube
concentration_vector <- unlist(concentration_map[sample_names_ordered])

# Final verification
if (any(is.na(concentration_vector))) {
  cat("WARNING: The final concentration vector contains NA (missing) values.\n")
  cat("Samples with NA:\n")
  print(sample_names_ordered[is.na(concentration_vector)])
}

if (length(concentration_vector) != n_samples) {
  stop("Error: The final concentration vector does not have the same number of samples as the EEM cube.")
}


# --- STEP 8: Save Final .mat File ---
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


cat(paste("\n--- SUCCESS! ---\n"))
cat(paste(".mat file saved at:", output_filepath, "\n"))
cat("The file contains:\n")
cat(" 1. 'cube': The 3D array [", paste(dim(eem_cube), collapse = " x "), "]\n")
cat(" 2. 'concentrations': Numeric concentration vector (", length(concentration_vector), " values)\n")
cat(" 3. 'sample_names': Vector of sample names (", length(sample_names_ordered), " names)\n")
cat(" 4. 'nm_emission': Emission wavelength vector\n")
cat(" 5. 'nm_excitation': Excitation wavelength vector\n")
}
