####
# SCRIPT [4]: Complete EEM Pre-processing
#
# OBJECTIVE:
# 1. Load the raw 3D cube ('eem_data.mat').
# 2. (STEP A) Apply scattering removal (Rayleigh, Raman).
# 3. (STEP B) Apply baseline correction.
# 4. (STEP C) Apply Savitzky-Golay smoothing.
# 5. (NEW) Offer immediate visualization of the result.
# 6. Save the final clean 3D cube ('eem_cube_cleaned.mat').
####

###########################################################################
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
  
  required_packages <- c("zoo", "signal", "R.matlab", "rstudioapi", 
                         "baseline", "plotly")
  
  load_required_packages(required_packages)
}

###########################################################################
# --- STEP 2: Scattering Removal Function (STEP A) ---
# (Original 'removeScatteringR' function from preprocessing.R)
removeScatteringR <- function(eem_matrix, exc_wavelengths, em_wavelengths, 
                              Rayleigh_1st, Rayleigh_2nd, Rayleigh_3rd, 
                              Raman_1st, Raman_2nd, raman_offset) {
  
  if (length(dim(eem_matrix)) != 3) {
    stop('eem_matrix must be a 3D array: [samples x emissions x excitations].')
  }
  dimensions <- dim(eem_matrix)
  nSamples <- dimensions[1]
  nEm <- dimensions[2]
  nEx <- dimensions[3]
  if (nEx != length(exc_wavelengths) || nEm != length(em_wavelengths)) {
    stop("Matrix dimensions must match wavelength vector lengths (Emissions: dim2, Excitations: dim3).")
  }
  
  eem_cleaned <- eem_matrix
  
  # --- 2. Scattering Removal (Masking with NaN) ---
  for (sampleIdx in 1:nSamples) {
    for (i in 1:length(exc_wavelengths)) {
      exc_wl <- exc_wavelengths[i]
      
      rayleigh1_lower = exc_wl - Rayleigh_1st
      rayleigh1_upper = exc_wl + Rayleigh_1st
      rayleigh2_lower = 2 * exc_wl - Rayleigh_2nd
      rayleigh2_upper = 2 * exc_wl + Rayleigh_2nd
      rayleigh3_lower = 3 * exc_wl - Rayleigh_3rd
      rayleigh3_upper = 3 * exc_wl + Rayleigh_3rd
      
      raman1_lower = exc_wl + raman_offset - Raman_1st
      raman1_upper = exc_wl + raman_offset + Raman_1st
      raman2_lower = exc_wl + 2 * raman_offset - Raman_2nd
      raman2_upper = 2 * raman_offset + Raman_2nd
      
      idx_rayleigh1 = em_wavelengths >= rayleigh1_lower & em_wavelengths <= rayleigh1_upper
      idx_rayleigh2 = em_wavelengths >= rayleigh2_lower & em_wavelengths <= rayleigh2_upper
      idx_rayleigh3 = em_wavelengths >= rayleigh3_lower & em_wavelengths <= rayleigh3_upper
      idx_raman1 = em_wavelengths >= raman1_lower & em_wavelengths <= raman1_upper
      idx_raman2 = em_wavelengths >= raman2_lower & em_wavelengths <= raman2_upper
      
      eem_cleaned[sampleIdx, idx_rayleigh1, i] = NaN
      eem_cleaned[sampleIdx, idx_rayleigh2, i] = NaN
      eem_cleaned[sampleIdx, idx_rayleigh3, i] = NaN
      eem_cleaned[sampleIdx, idx_raman1, i] = NaN
      eem_cleaned[sampleIdx, idx_raman2, i] = NaN
    }
  }
  
  # --- 3. Interpolation (PCHIP) by Sample and Excitation (Column-wise) ---
  cat("STEP A: Starting scattering removal and PCHIP interpolation...\n")
  pb <- txtProgressBar(min = 0, max = nSamples, style = 3) 
  
  for (sampleIdx in 1:nSamples) {
    for (ex in 1:nEx) {
      y <- eem_cleaned[sampleIdx, , ex]
      nan_idx <- is.na(y)
      if (!any(nan_idx)) next 
      
      W <- !nan_idx 
      valid_idx <- which(W)
      
      if (length(valid_idx) < 2) { 
        next
      }
      
      x <- em_wavelengths
      y_interp <- y
      x_known <- x[valid_idx]
      y_known <- y[valid_idx]
      x_fill  <- x[nan_idx]
      
      x_fill <- x_fill[x_fill >= min(x_known) & x_fill <= max(x_known)]
      
      if (length(x_fill) > 0) {
        y_fill <- signal::interp1(x = x_known, y = y_known, xi = x_fill, method = "pchip")
        y_interp[match(x_fill, x)] <- y_fill
      }
      
      first_valid_idx <- valid_idx[1]
      last_valid_idx <- valid_idx[length(valid_idx)]
      
      if (first_valid_idx > 1) {
        y_interp[1:(first_valid_idx - 1)] <- y_interp[first_valid_idx] 
      }
      if (last_valid_idx < nEm) {
        y_interp[(last_valid_idx + 1):nEm] <- y_interp[last_valid_idx] 
      }
      
      eem_cleaned[sampleIdx, , ex] <- y_interp
    }
    setTxtProgressBar(pb, sampleIdx) 
  }
  close(pb)
  cat("\nSTEP A: PCHIP interpolation completed.\n")
  return(eem_cleaned)
}

#########################################################################
# --- STEP 3: Baseline and Smoothing Function (STEPS B and C) ---
# (Function 'apply_baseline_and_smoothing' previously created)
apply_baseline_and_smoothing <- function(eem_cube, 
                                         baseline_method = 'modpoly', 
                                         baseline_degree = 4,
                                         sgolay_p = 3, 
                                         sgolay_n = 11) {
  
  if (sgolay_n %% 2 == 0) {
    stop(paste("The S-G window size (sgolay_n =", sgolay_n, ") must be odd."))
  }
  
  eem_final <- eem_cube
  
  dimensions <- dim(eem_cube)
  nSamples <- dimensions[1]
  nEm <- dimensions[2]
  nEx <- dimensions[3]
  
  cat("STEP B/C: Starting baseline correction and S-G smoothing...\n")
  pb <- txtProgressBar(min = 0, max = nSamples * nEx, style = 3) 
  
  for (i in 1:nSamples) {
    for (j in 1:nEx) {
      
      spectrum <- eem_cube[i, , j]
      spec_matrix <- matrix(spectrum, ncol = 1)
      
      # STEP B: Baseline Correction
      bl_correction <- try(baseline::baseline(spec_matrix, 
                                              method = baseline_method, 
                                              degree = baseline_degree), 
                           silent = TRUE)
      
      if (inherits(bl_correction, "try-error")) {
        spectrum_corrected <- spectrum 
      } else {
        spectrum_corrected <- baseline::getCorrected(bl_correction)
      }
      
      # STEP C: Savitzky-Golay Smoothing
      spectrum_smoothed <- signal::sgolayfilt(spectrum_corrected, 
                                              p = sgolay_p, 
                                              n = sgolay_n)
      
      spectrum_smoothed[spectrum_smoothed < 0] <- 0
      
      eem_final[i, , j] <- spectrum_smoothed
      setTxtProgressBar(pb, (i - 1) * nEx + j)
    }
  }
  
  close(pb)
  cat("\nSTEP B/C: Baseline and smoothing completed.\n")
  return(eem_final)
}

#########################################################################
# --- STEP 4: Visualization Functions (COPIED FROM SCRIPT [3]) ---

#' Function 1: Visualizes the 3D surface (plotly)
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

#' Function 2: Visualizes the 2D contour (filled.contour)
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
    z = eem_matrix, # Transposes the matrix to Z[x, y]
    color.palette = topo.colors,
    main = title_text,
    xlab = "Emission (nm)",
    ylab = "Excitation (nm)"
  )
}

######################################################################
# --- STEP 5: Loading and Execution of Pre-processing ---

# 1. Define your pre-processing parameters (ADJUST AS NEEDED)
params_scatter <- list(
  Rayleigh_1st = 15,  
  Rayleigh_2nd = 15, 
  Rayleigh_3rd = 15, 
  Raman_1st = 0,     
  Raman_2nd = 0,    
  raman_offset = 0  
)
params_smooth <- list(
  baseline_method = 'modpoly', 
  baseline_degree = 3,       
  sgolay_p = 3,              
  sgolay_n = 11              
)

# 2. Select and load the raw .mat file
cat("Please select the raw 'eem_data.mat' file...\n")
mat_file_path <- rstudioapi::selectFile(caption = "Select the eem_data.mat file",
                                        filter = "MAT Files (*.mat)")

if (!nzchar(mat_file_path)) {
  stop("No file selected. Script terminated.")
}

eem_data <- readMat(mat_file_path)
names(eem_data) <- gsub("\\.mat$", "", names(eem_data))

eem_cube_raw <- eem_data$cube
em_wavelengths <- as.vector(eem_data$nm.emission)
exc_wavelengths <- as.vector(eem_data$nm.excitation)
concentrations <- as.vector(eem_data$concentrations)

# Loads names to use in plot titles
sample_names <- unlist(eem_data$sample.names) 

cat(paste("Raw data loaded. Dimensions:", 
          paste(dim(eem_cube_raw), collapse = " x "), "\n"))

# 3. STEP A: Execute scattering removal
eem_cube_scatter_removed <- removeScatteringR(
  eem_matrix = eem_cube_raw, 
  exc_wavelengths = exc_wavelengths, 
  em_wavelengths = em_wavelengths, 
  Rayleigh_1st = params_scatter$Rayleigh_1st, 
  Rayleigh_2nd = params_scatter$Rayleigh_2nd, 
  Rayleigh_3rd = params_scatter$Rayleigh_3rd, 
  Raman_1st = params_scatter$Raman_1st, 
  Raman_2nd = params_scatter$Raman_2nd, 
  raman_offset = params_scatter$raman_offset
)

# 4. STEP B/C: Execute baseline correction and smoothing
eem_cube_cleaned <- apply_baseline_and_smoothing(
  eem_cube = eem_cube_scatter_removed,
  baseline_method = params_smooth$baseline_method,
  baseline_degree = params_smooth$baseline_degree,
  sgolay_p = params_smooth$sgolay_p,
  sgolay_n = params_smooth$sgolay_n
)

cat("\n--- Pre-processing complete --- \n")

######################################################################
# --- STEP 6: Visual Verification (Optional) ---
{
  prompt <- "Do you want to visualize a pre-processed sample before saving? (y/n): "
  resposta <- tolower(trimws(readline(prompt)))
  
  if (resposta == "y" || resposta == "yes") {
    
    # Ask for the index
    n_samples_total <- dim(eem_cube_cleaned)[1]
    prompt_idx <- paste0("Enter the NUMBER (index) of the sample you want to view (1 to ", n_samples_total, "): ")
    idx_ver <- as.integer(readline(prompt_idx))
    
    # Validate the index
    if (!is.na(idx_ver) && idx_ver >= 1 && idx_ver <= n_samples_total) {
      
      # Extract data to plot
      amostra_para_plotar <- eem_cube_cleaned[idx_ver, , ]
      nome_da_amostra <- sample_names[idx_ver]
      titulo_plot <- paste(nome_da_amostra, "(Pre-processed)")
      
      cat(paste("Generating plots for sample:", nome_da_amostra, "\n"))
      
      # Generate 3D Plot (interactive)
      print(visualize_3D_surface(
        eem_matrix = amostra_para_plotar,
        em_wavelengths = em_wavelengths,
        exc_wavelengths = exc_wavelengths,
        sample_name = titulo_plot
      ))
      
      # Generate 2D Plot (static)
      cat("Press [Enter] in the console to view the 2D plot...\n")
      readline() # Pause for user to see the 3D plot
      
      dev.new() # Open a new graphics window for the contour
      visualize_2D_contour(
        eem_matrix = amostra_para_plotar,
        em_wavelengths = em_wavelengths,
        exc_wavelengths = exc_wavelengths,
        sample_name = titulo_plot
      )
      
    } else {
      cat("Invalid or non-numeric index. Skipping visualization.\n")
    }
  } else {
    cat("Skipping visualization.\n")
  }
}

####################################################################
# --- STEP 7: Save Clean Cube ---
cat("\nReady to save. Please select the OUTPUT folder to save the clean file.\n")
output_directory <- selectDirectory(caption = "SELECT OUTPUT FOLDER")

if (nzchar(output_directory)) {
  output_filename <- "eem_cube_cleaned.mat"
  output_filepath <- file.path(output_directory, output_filename)
  
  # Save the clean cube AND the wavelength vectors
  writeMat(output_filepath, 
           cube_cleaned = eem_cube_cleaned,
           nm_emission = em_wavelengths,
           nm_excitation = exc_wavelengths,
           concentrations = concentrations,
           sample_names = sample_names
           
  )
  cat(paste("Clean data cube saved successfully at:", output_filepath, "\n"))
  
} else {
  cat("Output directory selection cancelled. The clean file was NOT saved.\n")
}
