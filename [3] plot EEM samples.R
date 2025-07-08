#### NOTE: This script is intended for visualizing EEM samples that have already been loaded
#### into the R environment using the script "Open EEM Samples".

### Step 1: Install and load the required package
# If the package is not yet installed in your R environment, use the command below:
install.packages("plotly")

# After installation, load the package:
library(plotly)

# Function to visualize the fluoresence surface 
visualize_3D_surface <- function(data_matrix) {
  # Check if the matrix is valid
  if (is.null(data_matrix) || nrow(data_matrix) == 0 || ncol(data_matrix) == 0) {
    stop("The data matrix must be non-null and have dimensions greater than zero.")
  }
  
  # Create a grid for the axes
  x <- 1:ncol(data_matrix)   # X-axis: Excitations
  y <- 1:nrow(data_matrix)   # Y-axis: Wavelengths
  z <- data_matrix           # Z-axis: Intensity
  
  # Create the 3D surface plot
  plot_ly() %>%
    add_surface(z = ~z, x = ~x, y = ~y, colorscale = 'Viridis') %>%
    layout(
      title = "3D Surface of the Data Matrix",
      scene = list(
        xaxis = list(title = "Excitations"),
        yaxis = list(title = "Wavelengths"),
        zaxis = list(title = "Intensity")
      )
    )
}

# Check if the sample exists in the dataset list
if (!exists("data_list") || !sample_name %in% names(data_list)) {
  stop("Error: The specified sample does not exist in the dataset.")
}

# Use the function to visualize the data
visualize_3D_surface(as.matrix(sample_name)) #Example: visualize_3D_surface(as.matrix(A05)), change the "A05" accordingly the sample that you want to see

# The function `visualize_3D_surface()` generates a 3D surface plot of a given EEM matrix 
# using the `plotly` package.
#
# INPUT:
# - data_matrix: a numeric matrix with dimensions [emission × excitation], 
#   where columns correspond to excitation wavelengths (e.g., 310–430 nm) 
#   and rows to emission wavelengths (e.g., 300–800 nm).
#
# OUTPUT:
# - An interactive 3D surface plot, with:
#   - X-axis: excitation channels (column indices)
#   - Y-axis: emission channels (row indices)
#   - Z-axis: fluorescence intensity values
#
# USAGE:
# - Ensure that the object corresponding to a sample (e.g., `A05`) is present in the environment.
# - Convert the sample object to a matrix and pass it to the function:
#
#     visualize_3D_surface(as.matrix(A05))
#
# NOTE:
# - Axis labels in this example use index positions. For physical wavelength labeling, 
#   you may replace x/y with actual wavelength vectors.
#
# - Example of how to generate wavelength vectors (optional enhancement):
#     x <- seq(300, 800, by = 1)   # Excitation wavelengths (13 columns)
#     y <- seq(310, 430, by = 10) # Emission wavelengths (501 rows)
