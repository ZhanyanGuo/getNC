#' Conditional prediction for combinational knockouts
#'
#' Compute \emph{Expectation target | knocked = 0} and \emph{Var target | knocked = 0} given
#' a covariance matrix and gene names under a multivariate Normal model.
#'
#' @param Sigma Covariance matrix of the Gaussian model (features x features).
#' @param genes Character vector of feature (gene) names aligned with Sigma.
#' @param target Character or integer vector: genes to predict.
#' @param knocked Character or integer vector: genes fixed to 0.
#' @param mu Optional mean vector (same order as `genes`). If NULL, zeros are used.
#'
#' @return list with `mean_cond` and `cov_cond`.
#'
#' @examples
#' # Small synthetic example (4 genes). Predict G1 given {G3, G4} knocked to 0.
#' set.seed(123)
#' genes <- c("G1","G2","G3","G4")
#' R <- matrix(c(1, 0.4, 0.2, 0,
#'               0.4, 1, 0.1, 0.3,
#'               0.2, 0.1, 1, -0.25,
#'               0,   0.3, -0.25, 1), 4, 4, dimnames=list(genes, genes))
#' sdv   <- c(1.0, 0.8, 1.2, 0.9)
#' Sigma <- diag(sdv) %*% R %*% diag(sdv)
#' mu    <- rep(0, length(genes))
#' ans   <- predict_conditional_knockout(
#'   Sigma, genes, target = "G1", knocked = c("G3","G4"), mu = mu
#' )
#' ans$mean_cond   # conditional mean for G1
#' ans$cov_cond    # conditional variance matrix for G1 (1x1 here)
#'
#' @seealso \code{\link{recover_covariance}}, \code{\link{predict_knockout_from_fit}}
#'
#' @references
#' Anderson, T.W. (2003) *An Introduction to Multivariate Statistical Analysis* (3rd ed.). Wiley.
#' Mardia, K.V., Kent, J.T., & Bibby, J.M. (1979) *Multivariate Analysis*. Academic Press.
#'
#' @export
predict_conditional_knockout <- function(Sigma,
                                         genes,
                                         target,
                                         knocked,
                                         mu = NULL) {
  if (is.null(Sigma) || is.null(genes)) {
    stop("Provide both `Sigma` (covariance) and `genes` (feature names).")
  }
  if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma)) {
    stop("`Sigma` must be a square matrix.")
  }
  if (length(genes) != nrow(Sigma)) {
    stop("`genes` length must match nrow(Sigma).")
  }

  # name/index resolution
  to_index <- function(sel) if (is.numeric(sel)) sel else match(sel, genes)

  a_idx <- to_index(target)
  b_idx <- to_index(knocked)

  if (any(is.na(a_idx))) stop("Some `target` genes not found in `genes`.")
  if (any(is.na(b_idx))) stop("Some `knocked` genes not found in `genes`.")

  # partitions
  Sigma_aa <- Sigma[a_idx, a_idx, drop = FALSE]
  Sigma_ab <- Sigma[a_idx, b_idx, drop = FALSE]
  Sigma_bb <- Sigma[b_idx, b_idx, drop = FALSE]

  # mean handling 
  if (is.null(mu)) {
    mu_full <- rep(0, length(genes))
  } else {
    if (length(mu) != length(genes)) {
      stop("`mu` length must equal length of `genes`.")
    }
    mu_full <- mu
  }
  mu_a <- mu_full[a_idx]
  mu_b <- mu_full[b_idx]

  # condition on X_b = 0
  x_b <- rep(0, length(b_idx))

  # conditional moments
  Sigma_bb_inv <- solve(Sigma_bb)

  mean_cond <- as.vector(mu_a + Sigma_ab %*% Sigma_bb_inv %*% (x_b - mu_b))
  cov_cond  <- Sigma_aa - Sigma_ab %*% Sigma_bb_inv %*% t(Sigma_ab)

  names(mean_cond)   <- genes[a_idx]
  dimnames(cov_cond) <- list(genes[a_idx], genes[a_idx])

  list(mean_cond = mean_cond, cov_cond = cov_cond)
}

#' Recover covariance in original units from a glasso fit on correlation scale
#'
#' Given the output list from \code{run_glasso_seurat()} (which contains
#' \code{sigma} estimated on the correlation scale and per-gene \code{sd}),
#' reconstruct the covariance matrix in the original units via
#' \eqn{\Sigma = D \, R \, D}, where \eqn{D = \mathrm{diag}(\text{sd})}
#' and \eqn{R = \mathrm{cov2cor}(\text{sigma})}.
#'
#' @param fit A list with components \code{sigma}, \code{sd}, and \code{features}
#'   as returned by \code{run_glasso_seurat()}.
#' @param renormalize Logical; if TRUE (default), first convert \code{fit$sigma}
#'   to an exact correlation matrix with \code{stats::cov2cor()} before scaling.
#'
#' @return A covariance matrix in original units (dimnames set to \code{features}).
#'
#' @examples
#' # Suppose a fit list stores correlation (sigma) and per-gene sds:
#' \dontrun{
#' set.seed(1)
#' genes <- paste0("G", 1:5)
#' R <- diag(5); R[1,2] <- R[2,1] <- 0.5
#' sdv <- runif(5, 0.7, 1.3)
#' fit_demo <- list(sigma = R, sd = sdv, features = genes)
#' S_cov <- recover_covariance(fit_demo, renormalize = TRUE)
#' dimnames(S_cov)
#' }
#'
#' @seealso \code{\link{predict_knockout_from_fit}}, \code{\link{run_glasso_seurat}}
#'
#' @references
#' Anderson, T.W. (2003) *An Introduction to Multivariate Statistical Analysis* (3rd ed.). Wiley.
#'
#' @importFrom stats cov2cor
#' @keywords internal
recover_covariance <- function(fit, renormalize = TRUE) {
  req <- c("sigma", "sd", "features")
  if (!all(req %in% names(fit))) {
    stop("fit must contain: ", paste(req, collapse = ", "), ".")
  }
  if (length(fit$sd) != nrow(fit$sigma)) {
    stop("Length of sd must match nrow(sigma).")
  }

  R <- if (isTRUE(renormalize)) stats::cov2cor(fit$sigma) else fit$sigma
  D <- diag(fit$sd, nrow = length(fit$sd), ncol = length(fit$sd))

  S_cov <- D %*% R %*% D
  dimnames(S_cov) <- list(fit$features, fit$features)
  S_cov
}

#' Predict conditional mean/covariance from a glasso-style fit
#'
#' Uses a fit list produced by \code{run_glasso_seurat()} / \code{fit_glasso()}
#' to compute \emph{Expectation target | knocked = 0]} and \emph{Var target | knocked = 0}.
#' You can choose to work on original units (recover covariance via sd)
#' or on the correlation scale.
#'
#' @param fit List with at least `sigma`, `features`, optionally `sd`, `mu`.
#' @param target Character or integer vector: genes to predict.
#' @param knocked Character or integer vector: genes set to 0.
#' @param use_original_units Logical; if TRUE (default) reconstruct covariance
#'   in original units via `recover_covariance(fit)`. If FALSE, use correlation scale.
#' @param renormalize Logical; if TRUE and `use_original_units=TRUE`,
#'   first renormalize `fit$sigma` to exact correlation via `cov2cor()` before scaling.
#'
#' @return list with `mean_cond` and `cov_cond`.
#'
#' @examples
#' # Build a tiny synthetic "fit" object and compute conditional for a target:
#' \dontrun{
#' set.seed(10)
#' genes <- paste0("G", 1:6)
#' R <- diag(6); R[1,2] <- R[2,1] <- 0.35; R[1,3] <- R[3,1] <- -0.25
#' sdv <- runif(6, 0.8, 1.2)
#' fit_small <- list(
#'   sigma    = R,           # here treated as correlation for the demo
#'   sd       = sdv,         # needed to recover covariance in original units
#'   mu       = rep(0, 6),   # optional; default 0 if absent
#'   features = genes
#' )
#' # Predict target G1 given G2 knocked to 0:
#' out <- predict_knockout_from_fit(
#'   fit_small, target = "G1", knocked = "G2",
#'   use_original_units = TRUE
#' )
#' out$mean_cond
#' out$cov_cond
#'}
#' 
#' @seealso \code{\link{predict_conditional_knockout}}, \code{\link{recover_covariance}},
#'   \code{\link{run_glasso_seurat}}, \code{\link{fit_glasso}}
#'
#' @references
#' Mardia, K.V., Kent, J.T., & Bibby, J.M. (1979) *Multivariate Analysis*. Academic Press.
#'
#' @export
predict_knockout_from_fit <- function(fit,
                                      target,
                                      knocked,
                                      use_original_units = TRUE,
                                      renormalize = TRUE) {
  if (is.null(fit$sigma) || is.null(fit$features)) {
    stop("`fit` must contain `sigma` and `features`.")
  }

  # choose covariance matrix to use
  if (isTRUE(use_original_units)) {
    if (is.null(fit$sd)) {
      stop("`fit$sd` is required to recover covariance in original units. ",
           "Set `use_original_units = FALSE` to operate on correlation scale.")
    }
    Sigma <- recover_covariance(fit, renormalize = renormalize)
  } else {
    # operate on correlation scale directly
    Sigma <- fit$sigma
    # if you want exact unit diagonals on corr scale:
    # Sigma <- stats::cov2cor(Sigma)
  }

  genes <- fit$features
  mu    <- if (!is.null(fit$mu)) fit$mu else rep(0, length(genes))

  predict_conditional_knockout(
    Sigma  = Sigma,
    genes  = genes,
    target = target,
    knocked = knocked,
    mu     = mu
  )
}

# Gen Ai used for documentation and input verify