#### NOTE: This script is intended solely for importing samples in .ASC format.
#### It requires the samples to have the correct naming convention.

### Step 1: Install and load the required package
# If the package is not yet installed in your R environment, use the command below:
install.packages("dplyr")
install.packages("rstudioapi")

  # After installation, load the package:
library(dplyr)
library(rstudioapi) 

### Step 2: Import the sample files into the R environment
# First, download the sample files available in the repository:
# https://github.com/LQBQ/Second-order-algorithms-.git
# The samples are attached in a folder called "ACS Files.rar"

### Step 3: Set the directory where the .ASC files are located on your computer:
{
  input_directory <- selectDirectory(caption = "SELECT INPUT FOLDER (.ASC Files)")

  if (input_directory != "") 
  {
    setwd(input_directory) } 
  else {
    stop("Directory selection cancelled. Script terminated.")
  }
}
 
### Step 4: List all .ASC files in the directory
files <- list.files(pattern = "\\.ASC$", ignore.case = TRUE)

### Step 5: Define the number of .ASC files associated with each sample.
# In this study, each sample was excited at 13 different wavelengths 
# (from 310 to 430 nm, in 10 nm increments), resulting in 13 .ASC files per sample.
# This value is used to group and process files correctly for multivariate analysis.
{
  files_per_sample <- 13  # If using this script with other datasets, adjust this number accordingly.

# Ensure that the total number of files is divisible by the number of files per sample
if (length(files) %% files_per_sample != 0) {
  stop("Error: The total number of files is not evenly divisible by the number of files per sample. Please check the directory.")
}
 }

### Step 6: Create a list to store the data organized by sample
data_list <- list()

### Step 7: Function to process and combine data from a sample
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

### Step 8: Process files in groups according to the defined sample size
for (i in seq(1, length(files), by = files_per_sample)) {
  sample_files <- files[i:(i + files_per_sample - 1)]
  
  # Extract the sample name (e.g., D%05, D050, PA01) from the file name
  # It now includes the '%' character for special cases (0.5%)
  raw_name <- sub("^((?:[ABCD][0-9%]+)|(?:P[ABCD][0-9]+)).*", "\\1", sample_files[1], ignore.case = TRUE)
  
  rename_sample <- function(name) {
    # 1. Adulterated Samples: Start with A/B/C/D followed by numbers, may contain '%'
    if (grepl("^[ABCD][0-9%]+$", name, ignore.case = TRUE)) {
      
      prefix <- substr(name, 1, 1) # Type letter (e.g., D)
      name_uppercase <- toupper(name)
      
      # --- Determine the concentration value (numeric) ---
      if (grepl("[ABCD]%[0-9]+$", name_uppercase)) {
        # Special Case: D%05 (0.5%)
        # Extract only the last digit (5)
        num_str <- sub("^[ABCD]%", "", name_uppercase) 
        last_digit <- sub(".*([0-9])$", "\\1", num_str)
        perc <- as.numeric(last_digit) / 10  # 5 / 10 = 0.5
        
      } else {
        # General Cases: D005 (5%), D050 (50%), D100 (100%)
        
        # Extract the numeric string (e.g., "050")
        num_str <- substr(name, 2, nchar(name)) 
        
        # Convert to number (005 -> 5, 050 -> 50, 100 -> 100)
        perc <- as.numeric(num_str) 
      }
      
      # --- Format the name using the D000.0% pattern ---
      # We use sprintf() to format the number as a string of length 5 (e.g., "100.0", "005.0", "000.5")
      # This ensures alphabetical order matches numerical order.
      # The format "%05.1f" means: fill with zeros to a length of 5, with 1 decimal place (e.g., 000.5, 005.0)
      
      formatted_perc <- sprintf("%05.1f", perc)
      
      # Assemble final name
      formatted <- paste0(toupper(prefix), formatted_perc, "%")
      
    } else if (grepl("^P[ABCD][0-9]+$", name, ignore.case = TRUE)) {
      # 2. Pure Samples: PA01, PB02, etc. (Standardize case)
      formatted <- toupper(name)
      
    } else {
      # Fallback
      formatted <- name 
    }
    
    return(formatted)
  }
  
  sample_name <- rename_sample(raw_name)
  
  # Process and store the data
  data_list[[sample_name]] <- process_sample(sample_files)
}

### Step 9: Set the output directory where you are going to save the spectra EEM
{
  output_directory <- selectDirectory(caption = "SELECT OUTPUT FOLDER")
  
  if (input_directory != "") 
  {
    setwd(input_directory) } 
  else {
    stop("Directory selection cancelled. Script terminated.")
  }
}

### Step 10: Save each sample's data as a CSV file
for (sample in names(data_list)) {
  file_path <- file.path(output_directory, paste0(sample, ".csv"))
  write.csv(data_list[[sample]], file = file_path, row.names = FALSE)
  cat("File saved:", file_path, "\n")
}
