### Function to Remove Harmonic Scattering
# If the package is not yet installed in your R environment, use the command below:
install.packages("zoo")

# After installation, load the package:
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
#Examploe to use: A05T <- remove_and_fill_scattering(A05, exc_wavelengths, em_wavelengths, exc_limit, em_limit)

visualize_3D_surface(as.matrix(sampleT)) #remenber to apply the visualize_3d_surface to plot the cleaned_sample
# Example to use: visualize_3D_surface(as.matrix(A05T))

# Define the output directory where files will be saved
# Set the directory where the .ASC files are located on your computer:
input_directory <- "path/to/your/directory"   # Adjust the path as needed
# Example: input_directory <- "C:/Users/Chemometric-PC/Documents/IC/Project - EVOO adulteration/Samples/Preprocessed samples"
# For compatibility, replace backslashes "\" with forward slashes "/" in your path. 

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

# This script defines and applies a function to remove first- and second-order Rayleigh scattering 
# from EEM data (excitation-emission matrices), commonly present in fluorescence spectroscopy.
#
# The function `remove_and_fill_scattering()`:
# - Identifies spectral regions affected by scattering based on the excitation wavelengths.
# - Masks these regions (sets them to NA) within a defined ±limit around each scattering center.
# - Reconstructs the missing values using linear interpolation (`na.approx`, from the `zoo` package).
#
# INPUTS:
# - eem_matrix: the raw EEM data matrix (emission × excitation).
# - exc_wavelengths: vector of excitation wavelengths (e.g., 310–430 nm).
# - em_wavelengths: vector of emission wavelengths (e.g., 300–800 nm).
# - exc_limit / em_limit: numeric range defining how close to the Rayleigh lines data should be removed.
#
# USAGE EXAMPLE:
# A05T <- remove_and_fill_scattering(A05, exc_wavelengths, em_wavelengths, exc_limit, em_limit)
# visualize_3D_surface(as.matrix(A05T))  # Optional: visualize the cleaned matrix
#
# After cleaning, all processed objects ending with "T" are saved as individual .csv files in the specified directory.
# The script filters objects in the environment, validates them as matrices/data.frames, and writes them to disk.
