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

# helper for Seurat v4/v5 compatibility
.get_data <- function(obj, which = c("data", "scale.data")) {
  which <- match.arg(which)
  if ("layer" %in% names(formals(Seurat::GetAssayData))) {
    Seurat::GetAssayData(obj, layer = which)
  } else {
    Seurat::GetAssayData(obj, slot  = which)
  }
}

#' Variable-feature selection + record mean/sd + z-score
#' (returns z, mu, sd, and features)
#' @param obj Seurat object (already normalized)
#' @param nfeatures number of variable features to select
#' @return list(z, mu, sd, features)
#' @keywords internal
zscore_seurat_with_params <- function(obj, nfeatures = 2000) {
  stopifnot(inherits(obj, "Seurat"))

  obj   <- Seurat::FindVariableFeatures(obj, nfeatures = nfeatures)
  feats <- Seurat::VariableFeatures(obj)

  # pull normalized data (compat for Seurat v4/v5)
  mat <- if ("layer" %in% names(formals(Seurat::GetAssayData))) {
    Seurat::GetAssayData(obj, layer = "data")
  } else {
    Seurat::GetAssayData(obj, slot  = "data")
  }

  # subset to features and coerce to plain matrix
  mat <- mat[feats, , drop = FALSE]
  mat <- as.matrix(mat)

  # ensure at least 2 features for covariance/glasso
  if (nrow(mat) < 2L) {
    stop("After variable-feature selection, only ", nrow(mat),
         " feature remained; need at least 2 for covariance/glasso. ",
         "Try increasing `nfeatures` or relaxing QC.")
  }

  mu <- rowMeans(mat)
  sd <- apply(mat, 1L, stats::sd)
  sd[sd == 0 | is.na(sd)] <- 1

  z <- (mat - mu) / sd

  list(
    z        = z,
    mu       = mu,
    sd       = sd,
    features = feats
  )
}
#' Run graphical lasso on a Seurat object
#'
#' This wrapper selects variable features, records per-gene mean and standard
#' deviation from the normalized expression, produces a z-scored matrix for
#' those features, and fits a sparse Gaussian graphical model via
#' \code{glasso::glasso()}.
#'
#' @param obj A Seurat object that has already been normalized
#'   (e.g., via \code{Seurat::NormalizeData()}).
#' @param nfeatures Integer. Number of variable features to use (default 2000).
#' @param rho Numeric. \eqn{\ell_1} regularization parameter for
#'   \code{glasso::glasso()}.
#'
#' @details Variable features are selected with \code{Seurat::FindVariableFeatures()}.
#' Per-gene mean (\code{mu}) and standard deviation (\code{sd}) are computed from the
#' normalized (unscaled) expression for the selected features. The covariance
#' matrix used by glasso is computed from the z-scored matrix of those features.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{omega}: estimated precision matrix (inverse covariance) for the selected features.
#'   \item \code{sigma}: estimated covariance matrix for the selected features.
#'   \item \code{mu}: numeric vector of per-gene means (order matches \code{features}).
#'   \item \code{sd}: numeric vector of per-gene standard deviations (order matches \code{features}).
#'   \item \code{features}: character vector of feature (gene) names used.
#'   \item \code{glasso}: the raw \code{glasso} fit object.
#' }
#'
#' @import Seurat
#' @import glasso
#' @importFrom stats cov
#' @export
run_glasso_seurat <- function(obj, nfeatures = 2000, rho = 0.1) {
  stopifnot(inherits(obj, "Seurat"))

  prep <- zscore_seurat_with_params(obj, nfeatures = nfeatures)
  zmat <- prep$z
  S    <- stats::cov(t(zmat))
  fit  <- glasso::glasso(S, rho = rho)

  list(
    omega    = fit$wi,
    sigma    = fit$w,
    mu       = prep$mu,
    sd       = prep$sd,
    features = prep$features,
    glasso   = fit
  )
}


#' Fit graphical lasso on a Seurat object (with QC + default data)
#'
#' This high-level wrapper (1) preprocesses a Seurat object (QC + optional
#' normalization) and (2) runs graphical lasso on variable, z-scored features.
#' If \code{obj} is NULL, it loads a small example dataset shipped with the
#' package (e.g., \code{pbmc_small}) and runs on that.
#'
#' @param obj Seurat object or NULL.
#' @param min_genes Minimum genes per cell for QC.
#' @param max_genes Maximum genes per cell for QC.
#' @param max_mt Maximum mitochondrial percent for QC.
#' @param normalize Logical; whether to run \code{Seurat::NormalizeData()}.
#' @param sf Scale factor for normalization.
#' @param nfeatures Number of variable features to use for glasso.
#' @param rho Regularization parameter for glasso.
#'
#' @return A list (from \code{run_glasso_seurat()}) with:
#' \itemize{
#'   \item \code{omega}: precision matrix (inverse covariance)
#'   \item \code{sigma}: covariance matrix
#'   \item \code{mu}: per-gene means (order matches \code{features})
#'   \item \code{sd}: per-gene standard deviations (order matches \code{features})
#'   \item \code{features}: character vector of genes used
#'   \item \code{glasso}: the raw \code{glasso} fit object
#' }
#'
#' @importFrom utils data
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
    utils::data("pbmc_small", package = "getNC", envir = environment())
    obj <- get("pbmc_small", envir = environment())
  }

  if (!inherits(obj, "Seurat")) {
    stop("fit_glasso() expects a Seurat object or NULL (to use example data).")
  }

  # QC + normalization
  obj <- preprocess_seurat_data(
    obj,
    min_genes = min_genes,
    max_genes = max_genes,
    max_mt    = max_mt,
    normalize = normalize,
    sf        = sf
  )

  # Glasso (returns omega, sigma, mu, sd, features, glasso)
  res <- run_glasso_seurat(
    obj,
    nfeatures = nfeatures,
    rho       = rho
  )

  res
}