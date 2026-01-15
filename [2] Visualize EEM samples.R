####
# SCRIPT [2]: EEM Sample Visualization (3D and 2D Contour)
#
# OBJECTIVE:
# 1. Load the processed .mat data cube.
# 2. Provide functions to visualize individual EEM samples.
#    - visualize_3D_surface: Interactive 3D surface plot.
#    - visualize_2D_contour: Interactive 2D contour plot (heatmap).
####

# --- STEP 1: Install and Load Required Packages ---
{
  load_required_packages <- function(packages) {
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
  }
  
  # plotly for interactive plots, R.matlab for loading data
  required_packages <- c("plotly", "rstudioapi", "R.matlab")
  
  load_required_packages(required_packages)
}

# --- STEP 2: Visualization Functions ---

#' Function 1: Visualizes the 3D surface (as in your original script)
#' @param eem_matrix A single 2D matrix [emission x excitation].
#' @param em_wavelengths Vector of emission wavelengths (Y-axis).
#' @param exc_wavelengths Vector of excitation wavelengths (X-axis).
#' @param sample_name (Optional) Name/ID of the sample for the title.
visualize_3D_surface <- function(eem_matrix, em_wavelengths, exc_wavelengths, sample_name = NULL) {
  if (is.null(eem_matrix) || nrow(eem_matrix) == 0 || ncol(eem_matrix) == 0) {
    stop("Data matrix must be non-null and have dimensions greater than zero.")
  }
  
  title_text <- ifelse(is.null(sample_name),
                       "3D Surface of EEM Sample",
                       paste("3D Surface - Sample:", sample_name))
  
  plot_ly() %>%
    add_surface(z = ~eem_matrix, 
                x = ~exc_wavelengths, 
                y = ~em_wavelengths, 
                colorscale = 'Viridis') %>%
    layout(
      title = title_text,
      scene = list(
        xaxis = list(title = "Excitation (nm)"),
        yaxis = list(title = "Emission (nm)"),
        zaxis = list(title = "Intensity")
      )
    )
}

#' Function 2: Visualizes the 2D contour (heatmap) (NEW FUNCTION)
#' @param eem_matrix A single 2D matrix [emission x excitation].
#' @param em_wavelengths Vector of emission wavelengths (Y-axis).
#' @param exc_wavelengths Vector of excitation wavelengths (X-axis).
#' @param sample_name (Optional) Name/ID of the sample for the title.
visualize_2D_contour <- function(eem_matrix, em_wavelengths, exc_wavelengths, sample_name = NULL) {
  if (is.null(eem_matrix) || nrow(eem_matrix) == 0 || ncol(eem_matrix) == 0) {
    stop("Data matrix must be non-null and have dimensions greater than zero.")
  }
  
  title_text <- ifelse(is.null(sample_name),
                       "2D Contour (Heatmap) of EEM Sample",
                       paste("2D Contour - Sample:", sample_name))
  
  filled.contour(
    x = em_wavelengths, 
    y = exc_wavelengths, 
    z = eem_matrix,
    color.palette = topo.colors, 
    main = title_text,
    xlab = "Emission (nm)",
    ylab = "Excitation (nm)"
  )
}

# --- STEP 3: Usage Example ---

# 1. Select the .mat file created in the consolidated script
{
mat_file_path <- rstudioapi::selectFile(caption = "Select the eem_data.mat file",
                                        filter = "MAT Files (*.mat)")

if (!nzchar(mat_file_path)) {
  stop("No file selected. Script terminated.")
}
  }
# 2. Load data from the .mat file
eem_data <- readMat(mat_file_path)

# 3. Extract components (using the names we defined)
#    Note: readMat might add '.mat' to the object name; we remove it if it exists
{
names(eem_data) <- gsub("\\.mat$", "", names(eem_data))
eem_cube_raw <- eem_data$cube
sample_names <- unlist(eem_data$sample.names) # Ensure it is a vector
em_wavelengths <- as.vector(eem_data$nm.emission) # Ensure it is a vector
exc_wavelengths <- as.vector(eem_data$nm.excitation) # Ensure it is a vector
}

cat("Data loaded successfully!\n")
cat(paste(length(sample_names), "samples found.\n"))


# 4. Choose which sample to visualize
#    Can be by index (e.g., 1) or by name (e.g., "D0005")

# --- Example by Index ---
{
Sample_ID_Index <- 10 # Change to the index of the sample you want to view
Sample_Name <- sample_names[Sample_ID_Index]
eem_to_plot <- eem_cube_raw[Sample_ID_Index, , ]
}

# --- Example by Name (uncomment to use) ---
{
 # Sample_Name <- "PA05" # Change to the name of the sample you want to view
 # Sample_ID_Index <- which(sample_names == Sample_Name)
 # if(length(Sample_ID_Index) == 0) stop("Sample name not found.")
 # eem_to_plot <- eem_cube_raw[Sample_ID_Index, , ]
}


# 5. Generate plots for the selected sample
cat(paste("Generating plots for sample:", Sample_Name, "(Index:", Sample_ID_Index, ")\n"))

# 3D Plot
grafico_3d <- visualize_3D_surface(
  eem_matrix = eem_to_plot,
  em_wavelengths = em_wavelengths,
  exc_wavelengths = exc_wavelengths,
  sample_name = Sample_Name
)
print(grafico_3d)

# 2D Plot (Heatmap)
dev.new()
grafico_2d <- visualize_2D_contour(
  eem_matrix = eem_to_plot,
  em_wavelengths = em_wavelengths,
  exc_wavelengths = exc_wavelengths,
  sample_name = Sample_Name
)

