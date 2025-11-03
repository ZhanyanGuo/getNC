#' Quality control and normalization
#' for a given seurat single cell expression matrix
#' (Seurat Object)
#'
#' @param obj Seurat object or matrix. If Seurat, QC is applied.
#' @param min_genes minimum genes per cell (default 200)
#' @param max_genes maximum genes per cell (default 2500)
#' @param max_mt    maximum percent mitochondrial genes (default 5)
#' @param normalize boolean for optional normalization (default True)
#' @param sf scale.factor for normalization (defult 10000)
#' @return subset of the seurat object
#' @import Seurat
#' @export

preprocess_seurat_data <- function(obj,
                                   min_genes = 200,
                                   max_genes = 2500,
                                   max_mt = 5,
                                   normalize = TRUE,
                                   sf = 10000) {
    if (!"percent.mt" %in% colnames(obj[[]])) {
        obj[["percent.mt"]] <-
            Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
    }
    obj <- subset(obj,
        subset = obj$nFeature_RNA > min_genes &
            obj$nFeature_RNA < max_genes &
            obj$percent.mt < max_mt
    )
    if (normalize) {
        obj <- Seurat::NormalizeData(obj,
            normalization.method = "LogNormalize",
            scale.factor = sf
        )
    }
    return(obj)
}

#' Z-score Seurat data (internal)
#'
#' Selects variable features, scales them, and returns the scaled matrix.
#' Assumes the object is already normalized.
#'
#' @param obj Seurat object.
#' @param nfeatures number of variable features to select.
#'
#' @return numeric matrix (features x cells)
#' @keywords internal
zscore_seurat <- function(obj, nfeatures = 2000) {
    obj <- Seurat::FindVariableFeatures(
        obj,
        nfeatures = nfeatures
    )

    feats <- Seurat::VariableFeatures(obj)
    obj <- Seurat::ScaleData(obj, features = feats)
    mat <- Seurat::GetAssayData(obj, slot = "scale.data")
    return(mat[feats, , drop = FALSE])
}

#' Run graphical lasso on a Seurat object
#'
#' This is a convenience wrapper that (1) selects variable features,
#' (2) z-scores them, and (3) fits a sparse Gaussian graphical model
#' using glasso. It is meant for single-cell data stored in a Seurat object.
#'
#' @param obj Seurat object that has already been normalized
#'   (i.e. \code{Seurat::NormalizeData()} has been run).
#' @param nfeatures Number of variable features to use (default 2000).
#' @param rho Regularization parameter for \code{glasso::glasso()}.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{omega}: precision matrix (inverse covariance)
#'     \item \code{sigma}: covariance matrix
#'     \item \code{features}: character vector of genes used
#'     \item \code{glasso}: the raw glasso fit object
#'   }
#'
#' @import Seurat
#' @import glasso
#' @export
run_glasso_seurat <- function(obj,
                              nfeatures = 2000,
                              rho = 0.1) {
    if (!inherits(obj, "Seurat")) {
        stop("run_glasso_seurat() expects a Seurat object.")
    }


    zmat <- zscore_seurat(obj, nfeatures = nfeatures)
    S <- stats::cov(t(zmat))
    fit <- glasso::glasso(S, rho = rho)

    list(
        omega    = fit$wi,
        sigma    = fit$w,
        features = rownames(zmat),
        glasso   = fit
    )
}


#' Fit graphical lasso on a Seurat object (with QC + default data)
#'
#' This is a high-level wrapper that
#' 1) preprocesses a Seurat object (QC + optional normalize)
#' 2) runs graphical lasso on variable, z-scored features.
#'
#' If \code{obj} is NULL, it will load the example dataset shipped with the
#' package (e.g. \code{pbmc_small}) and run on that.
#'
#' @param obj Seurat object or NULL.
#' @param min_genes minimum genes per cell for QC.
#' @param max_genes maximum genes per cell for QC.
#' @param max_mt maximum mitochondrial percent for QC.
#' @param normalize logical, whether to run Seurat::NormalizeData().
#' @param sf scale factor for normalization.
#' @param nfeatures number of variable features to use for glasso.
#' @param rho regularization parameter for glasso.
#'
#' @return list with precision, covariance, features, and glasso fit.
#'
#' @export
fit_glasso <- function(obj = NULL,
                       min_genes = 200,
                       max_genes = 2500,
                       max_mt    = 5,
                       normalize = TRUE,
                       sf        = 10000,
                       nfeatures = 20,
                       rho       = 0.1) {

  if (is.null(obj)) {
    data("pbmc_small", package = "getNC", envir = environment())
    obj <- get("pbmc_small", envir = environment())
  }

  if (!inherits(obj, "Seurat")) {
    stop("fit_glasso() expects a Seurat object or NULL (to use example data).")
  }

  # preprocess (QC + optional normalize)
  obj <- preprocess_seurat_data(
    obj,
    min_genes = min_genes,
    max_genes = max_genes,
    max_mt    = max_mt,
    normalize = normalize,
    sf        = sf
  )

  # run glasso
  res <- run_glasso_seurat(
    obj,
    nfeatures = nfeatures,
    rho       = rho
  )

  return(res)
}