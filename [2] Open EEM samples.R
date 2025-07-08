#### NOTE: This script is intended solely for samples already combined using the script [1] ("import + combine ASC files")

### Open EEM Samples

# Set the directory where the .ASC files are located on your computer:
input_directory <- "path/to/your/directory"  # Replace with the correct path
# Example: input_directory <- "C:/Users/Chemometric-PC/Documents/IC/Project - EVOO adulteration/Samples/CSV samples"
# For compatibility, replace backslashes "\" with forward slashes "/" in your path.

# Get a list of all .csv files in the directory
csv_files <- list.files(pattern = "\\.csv$")

# Loop to read each file individually
for (file in csv_files) {
  # Create an object name from the file name (without extension)
  object_name <- tools::file_path_sans_ext(file)
  
  # Read the file, ignoring the first row
  data <- read.csv(file, header = TRUE, sep = ",", dec = ".", stringsAsFactors = FALSE)
  
  # Convert all data to numeric type
  data <- as.data.frame(lapply(data, as.numeric))
  
  # Store the dataset in the environment with the created object name
  assign(object_name, data)
}

# Check the loaded data
ls()

# After running this script, each .csv file will be loaded into the R environment
# as a separate data frame named according to its file name (e.g., T1, C1, etc.).
#
# Each data frame represents a single sample and corresponds to an Excitation-Emission Matrix (EEM).
#
# Based on this study's setup:
# - Excitation wavelengths: from 310 to 430 nm in 10 nm steps → 13 columns (variables)
# - Emission wavelengths: from 300 to 800 nm in 1 nm steps → 501 rows
#
# Thus, each data frame will have dimensions 501 × 13 (emission × excitation),
# where each row is a fluorescence emission spectrum for a given excitation wavelength.
