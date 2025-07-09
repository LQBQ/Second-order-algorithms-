#### NOTE: This script is intended for use only after pre-processing of EEM samples (using "[4] preprocessing").
#### It requires the samples to be open in the R environment and to have the correct naming convention.

### Function to Create a 3D Array
create_3d_array <- function() {
  # Step 1: List objects in the global environment matching the pattern "PXXT" or "AXXT"
  objects <- ls(envir = .GlobalEnv)
  pattern <- "^(PU|A)\\d{1,}T$"  # Adjusted pattern (e.g., A10T, PU1T, P10T, ...)
  matrix_names <- grep(pattern, objects, value = TRUE)
  
  # Check if there are any matrices matching the pattern
  if (length(matrix_names) == 0) {
    stop("No matrices matching the pattern '(P|A)XXT' were found.")
  }
  
  # Step 2: Retrieve objects and convert data.frames to matrices
  matrix_list <- lapply(matrix_names, function(name) {
    obj <- get(name, envir = .GlobalEnv)
    
    # Convert to matrix if it is a data.frame
    if (is.data.frame(obj)) {
      return(as.matrix(obj))
    } else if (is.matrix(obj)) {
      return(obj)
    } else {
      stop(paste("The object", name, "is neither a data.frame nor a matrix."))
    }
  })
  
  # Step 3: Ensure all matrices have the same dimensions
  dimensions <- sapply(matrix_list, dim)
  if (!all(apply(dimensions, 1, function(x) length(unique(x)) == 1))) {
    stop("The matrices have inconsistent dimensions.")
  }
  
  # Step 4: Create the 3D array
  array_3d <- array(
    data = unlist(matrix_list),
    dim = c(dim(matrix_list[[1]]), length(matrix_list))
  )
  
  # Adjust the order of dimensions (layers, rows, columns)
  array_3d <- aperm(array_3d, perm = c(3, 1, 2))
  
  return(array_3d)
}

# Execute the function
Samples <- create_3d_array()

### Save and Load the 3D Array in .mat Format
# If the package is not yet installed in your R environment, use the command below:
install.packages("R.matlab")

# After installation, load the package:
library(R.matlab)

# Define the output directory dynamically
output_directory <- "path/to/your/directory"   # Adjust the path as needed
# Example: output_directory <- "C:/Users/Chemometric-PC/Documents/IC/Project - EVOO adulteration/Samples/Preprocessed samples"
# For compatibility, replace backslashes "\" with forward slashes "/" in your path. 

# Save as .mat format to preserve the 3D array dimensions
writeMat("Samples.mat", x = Samples)

# Confirmation message
cat("3D array saved successfully in:", output_directory, "\n")

# Load the .mat file back into R
loaded_data <- readMat("Samples.mat")

# Access the array
X <- loaded_data$x  # The variable name inside the .mat file should match the saved object name

# This script defines a function to create a 3D array from multiple EEM samples 
# previously loaded into the R environment. It is designed for matrices named 
# according to the pattern "AXXT" or "PXXT" (e.g., A1T, P2T), where each object 
# represents an Excitation-Emission Matrix (EEM) for one sample.
#
# The function performs the following steps:
# 1. Searches the global environment for objects matching the pattern (e.g., A1T, P1T).
# 2. Converts all matched data.frames or matrices into a list of numeric matrices.
# 3. Verifies that all matrices have identical dimensions.
# 4. Constructs a 3D array with dimensions: [samples × emission × excitation].
#
# ⚠️ NOTE ON SAMPLE ORDER:
# The sample objects are retrieved using `ls()` followed by `grep()`, 
# which means the ordering in the 3D array follows **alphabetical order of the object names**.
# Therefore, samples starting with "A" (e.g., A1T, A2T, ...) will appear **before** samples 
# starting with "P" (e.g., P1T, P2T, ...), regardless of their experimental grouping.
#
# If a specific sample order is required (e.g., group-wise or numerically sorted),
# the user must reorder the `matrix_names` vector **before** constructing the array.
#
# Example:
#   matrix_names <- c("P1T", "P2T", "A1T", "A2T")  # Custom order
#
# OUTPUT:
# - The resulting array has shape [n_samples × n_excitation × n_emission],
#   ready for use in multilinear modeling (e.g., PARAFAC, Tucker3) or machine learning workflows.
