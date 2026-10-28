# Tutorial: Second-Order Calibration and Classification Workflow

This tutorial explains how to run the R scripts in this repository for second-order chemometric analysis (PARAFAC, Tucker3, and N-PLS) of the provided EEM fluorescence data.

**Note on reproducibility.** To avoid compatibility problems previously reported with Jupyter Notebooks, the code is provided as standard R scripts (`.R`). Run them in **RStudio**: the scripts use `rstudioapi` dialogs to select files and folders.

## 1. Dataset structure

All data files are in the `Samples/` folder:

* **`Raw files.rar`**: compressed archive with the original `.ASC` instrument files (raw EEMs). Use it only if you want to reproduce the pre-processing from scratch.
* **`eem_data.mat`**: raw 3D array (Samples x Emission x Excitation), produced by script [1] and used as input of script [3].
* **`eem_cube_cleaned.mat`**: pre-processed 3D array, produced by script [3]. **This is the main input file for scripts [4] to [7].**

Sample names follow a fixed pattern, which you will need for outlier removal:

* Adulterated: brand letter (A-D) + soybean oil concentration x 10, zero-padded to four digits (`B0900` = brand B, 90.0%; `D0010` = brand D, 1.0%; `D0005` = brand D, 0.5%).
* Pure: `P` + brand letter + two-digit replicate (`PB02` = pure brand B, replicate 2).

## 2. Scripts and execution order

| Script | Algorithm / task | Objective | Input | Output |
| :--- | :--- | :--- | :--- | :--- |
| `[1] Import and create array 3D.R` | Data import | Reads the `.ASC` files and stacks them into a 3D array | Folder with `.ASC` files | `eem_data.mat` |
| `[2] Visualize EEM samples.R` | Visualization | Interactive 3D surface and 2D contour plots of one sample | `eem_data.mat` | Plots |
| `[3] Preprocessing samples.R` | Pre-processing | Rayleigh/Raman scattering removal, baseline correction, Savitzky-Golay smoothing | `eem_data.mat` | `eem_cube_cleaned.mat` |
| `[4] PARAFAC-LDA.R` | PARAFAC + LDA | Classification of pure vs. adulterated samples | `eem_cube_cleaned.mat` | Plots and metrics |
| `[5] Tucker3-LDA.R` | Tucker3 + LDA | Same as [4], using Tucker3 (flexible ranks per mode) | `eem_cube_cleaned.mat` | Plots and metrics |
| `[6] PARAFAC-MLR.R` | PARAFAC + MLR | Regression of soybean oil concentration | `eem_cube_cleaned.mat` | Plots and metrics |
| `[7] nPLS.R` | N-PLS (`sNPLS`) | Regression of soybean oil concentration | `eem_cube_cleaned.mat` | Plots and metrics |

**Two ways to use the scripts**

* **Modeling only (recommended for the classroom):** run any of [4] to [7] directly with `eem_cube_cleaned.mat`.
* **Full pipeline from raw files:** extract `Raw files.rar`, then run [1], optionally [2], then [3], and finally [4] to [7].

All scripts check for the required packages and install the missing ones automatically (internet connection needed on the first run).

## 3. Step-by-step

### 3.1 Loading the data

1. Open the desired script (for example, `[4] PARAFAC-LDA.R`) in RStudio.
2. Make sure `eem_cube_cleaned.mat` is on your computer.
3. Run the script. A dialog box titled "Select eem_cube_cleaned.mat" appears; select the file.

*Troubleshooting.* If the dialog does not appear (for example, in a non-interactive session), find the line that assigns `mat_cube_path` using `rstudioapi::selectFile(...)` and replace it with the full path to the file:

```r
mat_cube_path <- "C:/path/to/your/eem_cube_cleaned.mat"
```

Scripts [1] and [3] also ask for input and output folders through dialogs; the same approach applies if they do not appear.

### 3.2 Running script [1] (only for the full pipeline)

1. Extract `Raw files.rar` to a folder.
2. Run the script and select that folder as the input folder, then choose an output folder.
3. The script assumes `files_per_sample <- 13` (one `.ASC` file per excitation wavelength). Change this value if your acquisition used a different number of excitation wavelengths.

### 3.3 Outlier diagnosis and manual removal

Scripts [4] to [7] include an interactive outlier step.

1. The script plots **Hotelling's T² vs. Q residuals** for the model.
2. Identify samples that fall outside the limits (for example, a sample named `PC05` with high residuals).
3. Find the vector `nomes_para_remover` in the outlier block of the script and list the sample names to remove:

```r
# Example: removing samples identified in the diagnosis plot
nomes_para_remover <- c("PC05", "PB10")
```

4. Re-run from that block. The script rebuilds the model without those samples and reports any names it could not find.

**Important:** each script ships with an example list of names already filled in. These lists were defined for the original dataset and class activity. Review them against your own diagnosis plots before accepting the results.

### 3.4 Choosing model parameters

The scripts use baseline values as a starting point. They are not optimized, and choosing them is part of the learning activity:

| Script | Parameter | Baseline value |
| :--- | :--- | :--- |
| [4] PARAFAC-LDA | `num_factors` | 7 |
| [5] Tucker3-LDA | `ranks_tucker` | `c(7, 7, 7)` |
| [6] PARAFAC-MLR | `num_factors` | 8 |
| [7] nPLS | `num_factors` (latent variables) | 9 |

Other settings you can change: the training proportion (`train_ratio` in [4] and [5], `perc` in [6] and [7]; baseline 0.7) and the random seeds (`set.seed`). Changing the number of factors or the data subset changes the figures of merit, which is intended: compare calibration, cross-validation, and external prediction results to discuss overfitting and the limit of detection.

### 3.5 Train/test partitioning

* Scripts [4] and [5] (classification) use Kennard-Stone on the scores.
* Scripts [6] and [7] (regression) use the MLM approach: Kennard-Stone followed by a mutation step, with a fixed seed so that the split is reproducible.

### 3.6 Reading the results

* **Classification ([4], [5]):** accuracy, sensitivity, and specificity for the training, cross-validation, and test sets, plus LD1 plots for calibration and test samples.
* **Regression ([6], [7]):** R², RMSE, and MAE for calibration, cross-validation, and prediction, plus predicted vs. actual concentration plots.

These are the figures of merit discussed in the manuscript.

## 4. Reproducibility notes

* Random seeds are fixed in the scripts, but results can differ slightly across R and package versions.
* The scripts were developed with R 4.3.1.
* The models are teaching examples. Use them to learn the workflow, not as a validated method for authenticating olive oil.

## 5. Contact

For questions or suggestions, contact kassio.lima@ufrn.br.
