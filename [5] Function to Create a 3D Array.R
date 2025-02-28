### Function to Create a 3D Array
create_3d_array <- function() {
  # Step 1: List objects in the global environment matching the pattern "PXXT" or "AXXT"
  objects <- ls(envir = .GlobalEnv)
  pattern <- "^(P|A)\\d{1,}T$"  # Adjusted pattern (e.g., A1T, A10T)
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
install.packages("R.matlab")
library(R.matlab)

# Define the output directory dynamically
output_directory <- "C:/Path/To/Save/Directory"  # Change this to your desired path
setwd(output_directory)

# Save as .mat format to preserve the 3D array dimensions
writeMat("Samples.mat", x = Samples)

# Confirmation message
cat("3D array saved successfully in:", output_directory, "\n")

# Load the .mat file back into R
loaded_data <- readMat("Samples.mat")

# Access the array
X <- loaded_data$x  # The variable name inside the .mat file should match the saved object name
