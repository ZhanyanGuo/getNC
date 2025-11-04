
# GetNC

## Description

**GetNC** (*Gaussian Experiment Testing and Knockout Calculator*) is an R package for modeling and validating **gene regulatory relationships** using multivariate Gaussian inference. It fits sparse Gaussian graphical models from single-cell expression data, predicts **combinational knockout outcomes** via closed-form conditional Gaussian formulas, and validates model reliability against observed perturbation experiments.

**Key innovation:**
GetNC bridges the gap between network inference and biological validation — providing a **mathematically interpretable**, simulation-free method to predict and assess perturbation effects directly from fitted graphical models.

Developed under **R version 4.5.2 (2023-10-31)** on **macOS 12.5**.

---

## Installation

To install the latest version of the package:

```r
install.packages("devtools")
library("devtools")
devtools::install_github("ZhanyanGuo/getNC", build_vignettes = TRUE)
library("getNC")
```

To run the Shiny app: *Under construction.*

---

## Overview

### Quick start

```r
ls("package:getNC")
data(package = "getNC")
browseVignettes("getNC")
```

### Core user functions

| Function                                 | Purpose                                                                                |
| ---------------------------------------- | -------------------------------------------------------------------------------------- |
| `fit_glasso()`                           | Fit sparse Gaussian graphical model from gene expression data (using graphical lasso). |
| `predict_conditional_knockout()`         | Compute conditional mean and variance for a target gene given knockout set.            |
| `predict_knockout_from_fit()`            | Wrapper for conditional prediction using a fitted `fit_glasso()` object.               |
| `plot_partner_knockout_densities_dual()` | Generate 3D Gaussian density plots for top-K mean/variance partners.                   |
| `validate_knockout_group_from_fit()`     | Statistical validation of model predictions using observed knockout replicates.        |

### Workflow diagram

![getNC Workflow](man/figures/getNC_overview.png)

**Input:** Normalized single-cell RNA-seq matrix (e.g., Seurat object).
**Output:** Conditional Gaussian prediction and validation summary.

For a step-by-step walkthrough:

```r
browseVignettes("getNC")
```

---

## Contributions

**Author:** Zhanyan Guo

* Designed core statistical framework for conditional Gaussian modeling of combinational knockouts.
* Implemented covariance recovery, conditional prediction, and validation workflow.
* Authored vignette and documentation.

**External packages used**

| Package      | Role                                                          | Reference                                                                                                                                       |
| ------------ | ------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| **Seurat**   | Normalization, variable feature selection, scaling            | Hao, Y. et al. (2021). *Integrated analysis of multimodal single-cell data*. Cell 184(13):3573–3587.                                            |
| **glasso**   | Graphical lasso estimation for sparse inverse covariance      | Friedman, J., Hastie, T., & Tibshirani, R. (2008). *Sparse inverse covariance estimation with the graphical lasso*. Biostatistics 9(3):432–441. |
| **plotly**   | 3D visualization of conditional densities                     | Sievert, C. (2020). *Interactive Web-Based Data Visualization with R, plotly, and shiny.* CRC Press.                                            |
| **stats**    | Multivariate normal conditioning and base statistical methods | R Core Team (2024). *R: A Language and Environment for Statistical Computing.*                                                                  |
| **Matrix**   | Efficient linear algebra and sparse matrix support            | Bates, D., Maechler, M., & Jagan, M. (2023). *Matrix: Sparse and Dense Matrix Classes and Methods*. CRAN.                                       |
| **testthat** | Automated testing framework                                   | Wickham, H. (2011). *testthat: Get Started with Testing.* The R Journal 3(1):5–10.                                                              |
| **usethis**  | Package development utilities                                 | Wickham, H., & Bryan, J. (2023). *usethis: Automate Package and Project Setup.* CRAN.                                                           |
| **devtools** | Installation and build tools for R packages                   | Wickham, H., & Chang, W. (2023). *devtools: Tools to Make Developing R Packages Easier.* CRAN.                                                  |

**Generative AI assistance**

* ChatGPT (OpenAI GPT-5) was used for:

  * drafting documentation templates,
  * optimizing code clarity and roxygen comments,
  * generating part of the test and vignette narrative.

---

## References

* Friedman, J., Hastie, T., & Tibshirani, R. (2008). *Sparse inverse covariance estimation with the graphical lasso.* Biostatistics, 9(3):432–441.
* Hao, Y. et al. (2021). *Integrated analysis of multimodal single-cell data.* *Cell*, 184(13), 3573–3587.
* C. Sievert (2020). *Interactive Web-Based Data Visualization with R, plotly, and shiny.* CRC Press.
* Bates, D., Maechler, M., & Jagan, M. (2023). *Matrix: Sparse and Dense Matrix Classes and Methods.* CRAN.
* Wickham, H. (2011). *testthat: Get Started with Testing.* *The R Journal*, 3(1), 5–10.
* Wickham, H., & Bryan, J. (2023). *usethis: Automate Package and Project Setup.* CRAN.
* Wickham, H., & Chang, W. (2023). *devtools: Tools to Make Developing R Packages Easier.* CRAN.
* R Core Team (2024). *R: A Language and Environment for Statistical Computing.* Vienna, Austria.

---

## Acknowledgements

This package was developed as part of an assessment for the **2025 BCB410H: Applied Bioinformatics** course
at the **University of Toronto, Toronto, CANADA**.

**GetNC** welcomes issues, enhancement requests, and contributions.
To submit an issue, please use the **GitHub Issues** page.

---

```r
sessionInfo()
```
