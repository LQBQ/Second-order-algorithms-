# 🧪 Second-Order Calibration & Classification Workflow

This section provides a quick tutorial on how to use the updated scripts for high-order chemometric analysis (PARAFAC, Tucker3, and N-PLS) using the provided EEM fluorescence data.

**⚠️ Note on Reproducibility:** To ensure stability and avoid compatibility issues previously reported with Jupyter Notebooks, these scripts are provided as standard **R scripts (.R)**. We strongly recommend running them in **RStudio**.

## 📂 Dataset Structure

The repository includes the necessary data in two formats:

* **`eem_cube_cleaned.mat`**: The pre-processed 3D tensor (Samples × Emission × Excitation). **This is the main input file required for all scripts.**
* **`Raw_Samples.rar`**: A compressed archive containing the original `.ASC` instrument files (raw EEMs), provided for users who wish to reproduce the pre-processing steps from scratch.

## 🚀 How to Run the Analysis

The scripts are numbered and organized by algorithm and objective. They include automatic package installation checks.

### Available Scripts

| Script | Algorithm | Objective | Key Features |
| :--- | :--- | :--- | :--- |
| **`[4] PARAFAC-LDA.R`** | PARAFAC | Classification | Adulterated vs. Pure discrimination using Linear Discriminant Analysis. |
| **`[6] PARAFAC-MLR.R`** | PARAFAC | Regression | Concentration prediction using Multiple Linear Regression |
| **`[7] Tucker3-LDA.R`** | Tucker3 | Classification | Similar to [4] but using Tucker decomposition (flexible ranks). |
| **`[9] nPLS.R`** | N-PLS | Regression | N-way Partial Least Squares (via `sNPLS`) |

### Step-by-Step Tutorial

#### 1. Setup & Loading Data
1. Open the desired script (e.g., `[4] PARAFAC-LDA.R`) in **RStudio**.
2. **Important:** Ensure the file `eem_cube_cleaned.mat` is downloaded to your local machine.
3. Run the code. An interactive dialog box will appear asking you to select the `.mat` file.
   * *Troubleshooting:* If the dialog box does not appear or you are running in a non-interactive environment, locate the line `file.choose()` in the script and manually replace it with the full path to the file:
     `data <- readMat("C:/path/to/your/eem_cube_cleaned.mat")`

#### 2. Outlier Diagnosis & Manual Removal
Each script includes an interactive step for outlier management.
1. The script generates a **Hotelling's T² vs. Residuals (Q)** plot.
2. Identify outliers visually (e.g., a sample named `"P C 0 5"` with high residuals).
3. Locate the **"Manual Outlier Removal"** block in the code (usually Step X.5) and add the sample names to the vector:

```r
# Example: Removing specific samples based on the diagnosis plot
nomes_para_remover <- c("PC05", "PB10")
