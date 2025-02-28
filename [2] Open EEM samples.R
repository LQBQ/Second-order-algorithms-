### Open EEM Samples

# Prompt the user to define the directory where the .csv files are located
input_directory <- readline(prompt = "Enter the directory path containing the .csv files: ")
setwd(input_directory)

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
