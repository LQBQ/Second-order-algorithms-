### Function to visualize the matrix in 3D
# Load necessary packages
install.packages("plotly")
library(plotly)

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
visualize_3D_surface(as.matrix(sample_name))
