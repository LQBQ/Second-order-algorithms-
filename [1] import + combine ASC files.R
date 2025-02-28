# Load necessary packages
library(dplyr)

# Set the directory where the .ASC files are located
input_directory <- "path/to/your/directory"  # Change this to the correct directory path
setwd(input_directory)

# List all .ASC files in the directory
files <- list.files(pattern = "\\.ASC$", ignore.case = TRUE)

# Define the number of files per sample (adjust based on the experiment setup)
files_per_sample <- as.integer(readline(prompt = "Enter the number of files per sample: "))

# Ensure that the total number of files is divisible by the number of files per sample
if (length(files) %% files_per_sample != 0) {
  stop("Error: The total number of files is not evenly divisible by the number of files per sample. Please check the directory.")
}

# Create a list to store the data organized by sample
data_list <- list()

# Function to process and combine data from a sample
process_sample <- function(sample_files) {
  # List to store the dataframes from the files
  data_frames <- list()
  
  for (file in sample_files) {
    # Read the file content as text lines
    lines <- readLines(file)
    
    # Replace decimal commas with dots
    lines <- gsub("([0-9]),([0-9])", "\\1.\\2", lines)
    
    # Create a temporary file to store the corrected data
    temp_file <- tempfile(fileext = ".txt")
    writeLines(lines, temp_file)
    
    # Read the temporary file with corrected data
    data <- read.delim(temp_file, header = FALSE, sep = ",", dec = ".")
    
    # Remove the excitation column (first column)
    data <- data[,-1]
    
    # Add the dataframe to the list
    data_frames[[file]] <- data
    
    # Remove the temporary file
    unlink(temp_file)
  }
  
  # Combine the dataframes side by side (by columns)
  combined_data <- bind_cols(data_frames)
  
  return(combined_data)
}

# Process the files in groups based on the defined sample size
for (i in seq(1, length(files), by = files_per_sample)) {
  # Determine the files corresponding to this sample
  sample_files <- files[i:(i + files_per_sample - 1)]
  
  # Name the sample based on the first file in the group
  sample_name <- sub("^((T[0-9]+)|(C[0-9]+)).*", "\\1", sample_files[1], ignore.case = TRUE)
  
  # Process the sample data
  data_list[[sample_name]] <- process_sample(sample_files)
}

# Set the directory to save the CSV files
output_directory <- readline(prompt = "Enter the output directory path: ")

# Create the directory if it does not exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)  # "recursive = TRUE" creates subdirectories if necessary
}

# Loop to save each sample as a CSV file
for (sample in names(data_list)) {
  # Define the output file path
  file_path <- file.path(output_directory, paste0(sample, ".csv"))
  
  # Save the data in CSV format
  write.csv(data_list[[sample]], file = file_path, row.names = FALSE)
  
  cat("File saved:", file_path, "\n")  # Confirmation message
}
