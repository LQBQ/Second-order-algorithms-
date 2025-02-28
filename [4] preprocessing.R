### Function to Remove Harmonic Scattering
# Install and load necessary package
install.packages("zoo")
library(zoo)

# Function to remove and interpolate scattered regions
remove_and_fill_scattering <- function(eem_matrix, 
                                       exc_wavelengths, 
                                       em_wavelengths, 
                                       exc_limit, 
                                       em_limit) {
  # Copy matrix to avoid modifying the original
  eem_cleaned <- eem_matrix
  
  # Loop to adjust the removal of scattering
  for (i in 1:length(exc_wavelengths)) {
    exc_wl <- exc_wavelengths[i]
    
    # Define the scattering region (1st and 2nd order Rayleigh scattering)
    rayleigh1 <- exc_wl
    rayleigh2 <- 2 * exc_wl
    
    # Remove only data near the scattering regions
    eem_cleaned[em_wavelengths <= rayleigh1 + em_limit & em_wavelengths >= rayleigh1 - em_limit, i] <- NA
    eem_cleaned[em_wavelengths <= rayleigh2 + em_limit & em_wavelengths >= rayleigh2 - em_limit, i] <- NA
  }
  
  # Fill removed values using interpolation
  for (i in 1:ncol(eem_cleaned)) {
    eem_cleaned[, i] <- na.approx(eem_cleaned[, i], na.rm = FALSE)  # Interpolates NAs
  }
  
  return(eem_cleaned)
}

# Example usage:
exc_wavelengths <- seq(310, 430, length.out = 13)  # Define excitation wavelengths
em_wavelengths <- seq(300, 800, length.out = 501)  # Define emission wavelengths
exc_limit <- 5  # Excitation scattering limit
em_limit <- 5   # Emission scattering limit

# Apply the function
sampleT <- remove_and_fill_scattering(raw_sample, exc_wavelengths, em_wavelengths, exc_limit, em_limit)

visualize_3D_surface(as.matrix(sampleT)) #remenber to apply the visualize_3d_surface to plot the cleaned_sample

# Define the output directory where files will be saved
output_dir <- "~/Your/Path/Here"  # Adjust the path as needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Get a list of objects in the environment that match the pattern "T" (processed samples)
processed_objects <- ls(pattern = "T$")  # Selects objects ending with "T"

# Loop to save each object as a .csv file
for (obj in processed_objects) {
  # Retrieve the data from the object
  data <- get(obj)
  
  # Ensure the object is a matrix or dataframe before saving
  if (is.matrix(data) || is.data.frame(data)) {
    # Build the output file path
    file_path <- file.path(output_dir, paste0(obj, ".csv"))
    
    # Save the file as .csv
    write.csv(data, file = file_path, row.names = FALSE)
    
    # Confirmation message
    cat("File saved:", file_path, "\n")
  } else {
    cat("The object", obj, "is not a matrix or dataframe and was ignored.\n")
  }
}