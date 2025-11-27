#' Quality control and normalization
#' for a given Seurat single-cell expression matrix
#' (Seurat Object)
#'
#' @param obj Seurat object or matrix. If Seurat, QC is applied.
#' @param min_genes minimum genes per cell (default 200)
#' @param max_genes maximum genes per cell (default 2500)
#' @param max_mt    maximum percent mitochondrial genes (default 5)
#' @param normalize boolean for optional normalization (default TRUE)
#' @param sf scale.factor for normalization (default 10000)
#' @return subset of the Seurat object (post-QC; optionally normalized)
#'
#' @examples
#' \donttest{
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   set.seed(1)
#'   m <- matrix(rpois(2000, 5), nrow = 200, ncol = 10,
#'               dimnames = list(paste0("G",1:200), paste0("C",1:10)))
#'   obj <- Seurat::CreateSeuratObject(m)
#'   obj_qc <- preprocess_seurat_data(obj,
#'                                    min_genes = 0,
#'                                    max_genes = 5000,
#'                                    max_mt    = 100,
#'                                    normalize = TRUE,
#'                                    sf        = 10000)
#'   obj_qc
#' }
#' }
#'
#' @references
#' Hao, Y. et al. (2021). Integrated analysis of multimodal single-cell data.
#'   \emph{Cell} 184(13):3573–3587.
#'
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
                  obj$percent.mt < max_mt)
  if (normalize) {
    obj <- Seurat::NormalizeData(obj,
                                 normalization.method = "LogNormalize",
                                 scale.factor = sf)
  }
  return(obj)
}

# helper for Seurat v4/v5 compatibility
#' @keywords internal
.get_data <- function(obj, which = c("data", "scale.data")) {
  which <- match.arg(which)
  if ("layer" %in% names(formals(Seurat::GetAssayData))) {
    Seurat::GetAssayData(obj, layer = which)
  } else {
    Seurat::GetAssayData(obj, slot  = which)
  }
}

#' Variable-feature selection + record mean/sd + z-score (internal)
#'
#' Returns a list with z-scored matrix, per-gene mean/sd, and the feature set.
#'
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

#' @keywords internal
run_glasso_seurat_helper <- function(obj, nfeatures, rho) {
  stopifnot(inherits(obj, "Seurat"))

  prep <- zscore_seurat_with_params(obj, nfeatures = nfeatures)

  zmat <- prep$z                     # genes × cells
  S    <- stats::cov(t(zmat))        # covariance across genes
  fit  <- glasso::glasso(S, rho = rho)

  list(
    omega    = fit$wi,
    sigma    = fit$w,
    mu       = prep$mu,
    sd       = prep$sd,
    features = prep$features,
    glasso   = fit,
    type     = "seurat"
  )
}

#' @keywords internal
run_glasso_matrix_helper <- function(mat, nfeatures, rho) {

  stopifnot(is.matrix(mat), is.numeric(mat))

  # HVG + z-scoring (cells × genes)
  prep <- zscore_matrix_with_params(mat, nfeatures = nfeatures)

  # covariance across genes
  S <- stats::cov(prep$z)
  fit <- glasso::glasso(S, rho = rho)

  list(
    omega    = fit$wi,
    sigma    = fit$w,
    mu       = prep$mu,
    sd       = prep$sd,
    features = prep$features,
    glasso   = fit,
    type     = "matrix"
  )
}

#' Run Graphical Lasso on a Seurat Object or Raw Count Matrix
#'
#' This unified wrapper automatically detects whether the input is a
#' \code{Seurat} object or a raw count matrix with rows representing cells
#' and columns representing genes. It applies the appropriate preprocessing
#' pipeline (QC + normalization for raw matrices, variable feature
#' selection + z-scoring for Seurat) before estimating a sparse precision
#' matrix using the graphical lasso.
#'
#' @param x A Seurat object (already normalized) or a numeric matrix
#'          with \strong{rows = cells} and \strong{columns = genes}.
#' @param nfeatures Integer. Number of highly variable genes to select.
#'        Default: 2000.
#' @param rho Numeric. L1-regularization parameter passed to
#'        \code{glasso::glasso()}. Default: 0.1.
#'
#' @details
#' \strong{If the input is a Seurat object}  
#' The function:
#' \enumerate{
#'   \item selects variable features via \code{Seurat::FindVariableFeatures()},
#'   \item extracts normalized data (layer/slot = "data"),
#'   \item computes per-gene mean and standard deviation,
#'   \item z-scores genes,
#'   \item computes covariance across selected genes,
#'   \item runs \code{glasso::glasso()}.
#' }
#'
#' \strong{If the input is a raw matrix}  
#' The function:
#' \enumerate{
#'   \item performs QC (min/max genes, mitochondrial%, etc.),
#'   \item normalizes by library size and applies log1p,
#'   \item selects variable genes by variance,
#'   \item z-scores genes,
#'   \item computes covariance,
#'   \item runs \code{glasso::glasso()}.
#' }
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{omega} — estimated precision matrix (inverse covariance)
#'   \item \code{sigma} — estimated covariance matrix
#'   \item \code{mu} — per-gene means (of normalized values)
#'   \item \code{sd} — per-gene standard deviations
#'   \item \code{features} — vector of gene names used
#'   \item \code{glasso} — raw \code{glasso} fit object
#'   \item \code{type} — either \code{"seurat"} or \code{"matrix"}
#' }
#'
#' @examples
#' \donttest{
#' if (requireNamespace("Seurat", quietly = TRUE) &&
#'     requireNamespace("glasso", quietly = TRUE)) {
#'   set.seed(1)
#'   m <- matrix(rpois(3000, 5), nrow = 100, ncol = 30,
#'               dimnames = list(paste0("G",1:100), paste0("C",1:30)))
#'   obj <- Seurat::CreateSeuratObject(m)
#'   obj <- Seurat::NormalizeData(obj)
#'
#'   fit1 <- run_glasso(obj, nfeatures = 200, rho = 0.2)
#'   fit2 <- run_glasso(m,   nfeatures = 200, rho = 0.2)
#' }
#' }
#'
#' @references
#' Friedman, J., Hastie, T., & Tibshirani, R. (2008).
#'   Sparse inverse covariance estimation with the graphical lasso.
#'   \emph{Biostatistics} 9(3), 432–441.
#'
#' Hao, Y. et al. (2021). Integrated analysis of multimodal single-cell data.
#'   \emph{Cell} 184(13):3573–3587.
#' 
#' @import Seurat
#' @import glasso
#' @importFrom stats cov
#' @export
run_glasso <- function(x = NULL, nfeatures = 2000, rho = 0.1) {
  auto_fit_glasso(x, nfeatures = 2000, rho = 0.1)
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
#' @return A list (from \code{run_glasso()}) with:
#' \itemize{
#'   \item \code{omega}: precision matrix (inverse covariance)
#'   \item \code{sigma}: covariance matrix
#'   \item \code{mu}: per-gene means (order matches \code{features})
#'   \item \code{sd}: per-gene standard deviations (order matches \code{features})
#'   \item \code{features}: character vector of genes used
#'   \item \code{glasso}: the raw \code{glasso} fit object
#' }
#'
#' @examples
#' \donttest{
#' if (requireNamespace("Seurat", quietly = TRUE) &&
#'     requireNamespace("glasso", quietly = TRUE)) {
#'   # Use internal example when obj = NULL (expects pbmc_small bundled in package)
#'   fit <- fit_glasso(obj = NULL, nfeatures = 50, rho = 0.15)
#'   names(fit)
#' }
#' }
#'
#' @references
#' Hao, Y. et al. (2021). Integrated analysis of multimodal single-cell data.
#'   \emph{Cell} 184(13):3573–3587.
#' Friedman, J., Hastie, T., & Tibshirani, R. (2008).
#'   Sparse inverse covariance estimation with the graphical lasso.
#'   \emph{Biostatistics}, 9(3), 432–441.
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
  res <- run_glasso(
    obj,
    nfeatures = nfeatures,
    rho       = rho
  )

  res
}

#' Fit graphical lasso on a raw count matrix (cells × genes)
#'
#' This high-level wrapper (1) preprocesses a raw single-cell count matrix
#' (QC + optional normalization) and (2) runs graphical lasso on variable,
#' z-scored features. The interface mirrors \code{fit_glasso()} so that both
#' Seurat and raw-matrix workflows use identical parameters.
#'
#' @param mat Numeric matrix with \strong{rows = cells} and \strong{columns = genes}.
#' @param min_genes Minimum genes per cell for QC (default 200).
#' @param max_genes Maximum genes per cell for QC (default 2500).
#' @param max_mt Maximum mitochondrial percent for QC (default 5).
#' @param normalize Logical; whether to perform library-size normalization
#'        and log1p transform. Default TRUE.
#' @param sf Scale factor for normalization (default 10000).
#' @param nfeatures Number of variable features (genes) to use for glasso.
#' @param rho Regularization parameter for glasso.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{omega}: precision matrix (inverse covariance)
#'   \item \code{sigma}: covariance matrix
#'   \item \code{mu}: per-gene means
#'   \item \code{sd}: per-gene standard deviations
#'   \item \code{features}: selected genes
#'   \item \code{glasso}: raw \code{glasso} fit object
#'   \item \code{type}: "matrix"
#' }
#'
#' @examples
#' \donttest{
#' if (requireNamespace("glasso", quietly = TRUE)) {
#'   m <- matrix(rpois(3000, 4), nrow = 100, ncol = 30,
#'               dimnames = list(paste0("C",1:100), paste0("G",1:30)))
#'   fit <- fit_glasso_raw(m,
#'  min_genes = 0,
#'  max_genes = 1e6,
#'  max_mt    = 100,  
#'  nfeatures = 5,
#' rho       = 0.05)
#' }
#' }
#'
#' @export
fit_glasso_raw <- function(mat,
                           min_genes = 200,
                           max_genes = 2500,
                           max_mt    = 5,
                           normalize = TRUE,
                           sf        = 10000,
                           nfeatures = 20,
                           rho       = 0.1) {

  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("fit_glasso_raw() expects a numeric matrix with rows=cells and cols=genes.")
  }

  # Step 1 — QC + normalization
  mat_qc <- preprocess_matrix_raw(
    mat,
    min_genes = min_genes,
    max_genes = max_genes,
    max_mt    = max_mt,
    normalize = normalize,
    sf        = sf
  )

  # Step 2 — Call the unified dispatcher (matrix case)
  res <- run_glasso(
    mat_qc,
    nfeatures = nfeatures,
    rho       = rho
  )

  res
}

#' Automatically dispatch to Seurat or raw-matrix graphical lasso
#'
#' This function detects whether the input is NULL, a Seurat object, or
#' a raw count matrix (rows = cells, columns = genes), and calls the
#' appropriate wrapper:
#' \itemize{
#'   \item \code{fit_glasso()} for NULL or Seurat input
#'   \item \code{fit_glasso_raw()} for raw matrices
#' }
#'
#' @param x A Seurat object, a numeric matrix, or NULL.
#' @param min_genes Minimum genes per cell for QC.
#' @param max_genes Maximum genes per cell for QC.
#' @param max_mt Maximum mitochondrial percent for QC.
#' @param normalize Whether to normalize data.
#' @param sf Scale factor for normalization.
#' @param nfeatures Number of variable genes for glasso.
#' @param rho Regularization parameter for glasso.
#'
#' @return Output of either \code{fit_glasso()} or \code{fit_glasso_raw()}.
#'
#' @export
auto_fit_glasso <- function(x = NULL,
                            min_genes = 200,
                            max_genes = 2500,
                            max_mt    = 5,
                            normalize = TRUE,
                            sf        = 10000,
                            nfeatures = 20,
                            rho       = 0.1) {
  # Case 1: NULL or Seurat → use fit_glasso()
  if (is.null(x) || inherits(x, "Seurat")) {
    return(
      fit_glasso(
        obj        = x,
        min_genes  = min_genes,
        max_genes  = max_genes,
        max_mt     = max_mt,
        normalize  = normalize,
        sf         = sf,
        nfeatures  = nfeatures,
        rho        = rho
      )
    )
  }

  # Case 2: raw matrix → use fit_glasso_raw()
  if (is.matrix(x) && is.numeric(x)) {
    return(
      fit_glasso_raw(
        mat        = x,
        min_genes  = min_genes,
        max_genes  = max_genes,
        max_mt     = max_mt,
        normalize  = normalize,
        sf         = sf,
        nfeatures  = nfeatures,
        rho        = rho
      )
    )
  }
  stop("auto_fit_glasso(): input must be NULL, a Seurat object, or a numeric matrix.")
}




#' Quality control + normalization for a raw count matrix (cells × genes)
#'
#' @param mat numeric matrix, cells × genes (raw counts)
#' @param min_genes minimum detected genes per cell
#' @param max_genes maximum detected genes per cell
#' @param max_mt maximum mitochondrial percent (default 5)
#' @param normalize whether to normalize (default TRUE)
#' @param sf scale factor for normalization (default 10000)
#'
#' @return a normalized log1p matrix (cells × genes)
#' @export
preprocess_matrix_raw <- function(mat,
                                  min_genes = 200,
                                  max_genes = 2500,
                                  max_mt    = 5,
                                  normalize = TRUE,
                                  sf        = 10000) {

  stopifnot(is.matrix(mat), is.numeric(mat))

  # detect mitochondrial genes
  mt_genes <- grepl("^MT-", colnames(mat), ignore.case = TRUE)
  if (all(!mt_genes)) {
    warning("No mitochondrial genes detected using '^MT-' pattern.")
  }

  # cell-level QC
  ngenes <- rowSums(mat > 0)
  libsize <- rowSums(mat)

  percent_mt <- if (any(mt_genes)) {
    rowSums(mat[, mt_genes, drop = FALSE]) / libsize * 100
  } else {
    rep(0, nrow(mat))
  }

  keep_cells <- ngenes > min_genes &
                ngenes < max_genes &
                percent_mt < max_mt

  mat <- mat[keep_cells, , drop = FALSE]

  if (!normalize) return(mat)

  # library size normalization
  lib <- rowSums(mat)
  mat_norm <- sweep(mat, 1L, lib, "/") * sf

  # log1p to log the variance for numeric stability
  mat_log <- log1p(mat_norm)

  mat_log
}

#' Variable gene selection & Z-scoring for raw matrix
#'
#' @param mat normalized log1p matrix, cells × genes
#' @param nfeatures number of highly variable genes to select
#'
#' @return list(z, mu, sd, features)
#' @export
zscore_matrix_with_params <- function(mat, nfeatures = 2000) {

  stopifnot(is.matrix(mat), is.numeric(mat))

  #compute per-gene variance
  gene_var <- apply(mat, 2L, stats::var)

  # rank genes by variance
  nfeatures <- min(nfeatures, ncol(mat))
  feats <- names(sort(gene_var, decreasing = TRUE))[1:nfeatures]

  submat <- mat[, feats, drop = FALSE]

  #compute mean and sd
  mu <- colMeans(submat)
  sd <- apply(submat, 2L, stats::sd)

  # avoid division by zero
  sd[sd == 0 | is.na(sd)] <- 1

  # z-score -
  z <- sweep(submat, 2L, mu, "-")
  z <- sweep(z,      2L, sd, "/")

  list(
    z        = z,
    mu       = mu,
    sd       = sd,
    features = feats
  )
}

# Gen AI used for documentation and input verify
