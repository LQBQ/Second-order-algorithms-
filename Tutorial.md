# 🧪 Second-Order Calibration & Classification Workflow

This section provides a quick tutorial on how to use the updated scripts for high-order chemometric analysis (PARAFAC, Tucker3, and N-PLS) using the provided EEM fluorescence data.

## 📂 Dataset Structure

The repository includes the necessary data in two formats:

* **`eem_cube_cleaned.mat`**: The pre-processed 3D tensor (Samples × Emission × Excitation). **This is the main input file required for all scripts.**
* **`Raw_Samples.rar`**: A compressed archive containing the original `.ASC` instrument files (raw EEMs), provided for users who wish to reproduce the pre-processing steps from scratch.

## 🚀 How to Run the Analysis

The scripts are numbered and organized by algorithm and objective (Classification or Regression). They are designed to be interactive and self-contained (automatic package installation included).

### Available Scripts

| Script | Algorithm | Objective | Key Features |
| :--- | :--- | :--- | :--- |
| **`[4] PARAFAC-LDA.R`** | PARAFAC | Classification | Adulterated vs. Pure discrimination using Linear Discriminant Analysis. |
| **`[6] PARAFAC-MLR.R`** | PARAFAC | Regression | Concentration prediction using Multiple Linear Regression |
| **`[7] Tucker3-LDA.R`** | Tucker3 | Classification | Similar to [4] but using Tucker decomposition (flexible ranks). |
| **`[9] nPLS.R`** | N-PLS | Regression | N-way Partial Least Squares (via `sNPLS`) |

### Step-by-Step Tutorial

#### 1. Loading the Data
Open any of the scripts in **RStudio**. When you run the code, a dialog box will appear. Select the file **`eem_cube_cleaned.mat`**. The script will automatically parse the sample names and concentrations.

#### 2. Outlier Diagnosis & Manual Removal
Each script includes an interactive step for outlier management.
1.  The script generates a **Hotelling's T² vs. Residuals (Q)** plot.
2.  Identify outliers visually (e.g., a sample named `"P C 0 5"` with high residuals).
3.  Locate the **"Manual Outlier Removal"** block in the code (usually Step X.5) and add the sample names to the list:

```r
# Example: Removing specific samples based on the diagnosis plot
nomes_para_remover <- c("PC05", "PB10")
