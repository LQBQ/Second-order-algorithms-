#### NOTE: This script is intended solely for importing samples in .ASC format.

### Step 1: Install and load the required package
# If the package is not yet installed in your R environment, use the command below:
install.packages("dplyr")

# After installation, load the package:
library(dplyr)

### Step 2: Import the sample files into the R environment
# First, download the sample files available in the repository:
# https://github.com/LQBQ/Second-order-algorithms-.git
# The samples are attached in a folder called "ACS Files.rar"

# Set the directory where the .ASC files are located on your computer:
input_directory <- "path/to/your/directory"  # Replace with the correct path
# Example: input_directory <- "C:/Users/Chemometric-PC/Documents/IC/Project - EVOO adulteration/Samples"
# For compatibility, replace backslashes "\" with forward slashes "/" in your path.

setwd(input_directory)

# List all .ASC files in the directory
files <- list.files(pattern = "\\.ASC$", ignore.case = TRUE)

# Define the number of .ASC files associated with each sample.
# In this study, each sample was excited at 13 different wavelengths 
# (from 310 to 430 nm, in 10 nm increments), resulting in 13 .ASC files per sample.
# This value is used to group and process files correctly for multivariate analysis.
files_per_sample <- 13  # If using this script with other datasets, adjust this number accordingly.

# Ensure that the total number of files is divisible by the number of files per sample
if (length(files) %% files_per_sample != 0) {
  stop("Error: The total number of files is not evenly divisible by the number of files per sample. Please check the directory.")
}

# Create a list to store the data organized by sample
data_list <- list()

# Function to process and combine data from a sample
process_sample <- function(sample_files) {
  data_frames <- list()
  
  for (file in sample_files) {
    # Read file content as text lines
    lines <- readLines(file)
    
    # Replace decimal commas with dots
    lines <- gsub("([0-9]),([0-9])", "\\1.\\2", lines)
    
    # Create a temporary file with corrected decimal format
    temp_file <- tempfile(fileext = ".txt")
    writeLines(lines, temp_file)
    
    # Read the corrected data
    data <- read.delim(temp_file, header = FALSE, sep = ",", dec = ".")
    
    # Remove the excitation column (first column)
    data <- data[, -1]
    
    data_frames[[file]] <- data
    unlink(temp_file)
  }
  
  # Combine all data frames column-wise
  combined_data <- bind_cols(data_frames)
  return(combined_data)
}

# Process files in groups according to the defined sample size
for (i in seq(1, length(files), by = files_per_sample)) {
  sample_files <- files[i:(i + files_per_sample - 1)]
  
  # Extract the sample name (e.g., T1, C2) from the file name
  sample_name <- sub("^((T[0-9]+)|(C[0-9]+)).*", "\\1", sample_files[1], ignore.case = TRUE)
  
  # Process and store the data
  data_list[[sample_name]] <- process_sample(sample_files)
}

# Set the output directory to save the processed CSV files
output_directory <- "path/to/your/directory"
# Example: output_directory <- "C:/Users/Chemometric-PC/Documents/IC/Project - EVOO adulteration/Samples/CSV samples"
# Again, use forward slashes "/" for compatibility.

# Create the output directory if it does not exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# Save each sample's data as a CSV file
for (sample in names(data_list)) {
  file_path <- file.path(output_directory, paste0(sample, ".csv"))
  write.csv(data_list[[sample]], file = file_path, row.names = FALSE)
  cat("File saved:", file_path, "\n")
}
