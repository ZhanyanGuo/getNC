#' Sweep partner knockouts for a single target
#'
#' For each candidate partner gene g (excluding the target and already-knocked genes),
#' compute E[target | knocked ∪ {g} = 0] and Var[target | knocked ∪ {g} = 0].
#' Inputs mirror `predict_conditional_knockout()`.
#'
#' @param Sigma Covariance matrix (features x features).
#' @param genes Character vector of feature names aligned with `Sigma`.
#' @param target Character or integer (must specify exactly one target gene).
#' @param knocked Character or integer vector of genes already set to 0.
#' @param mu Optional mean vector aligned with `genes`; if NULL, zeros used.
#'
#' @return list with two named numeric vectors (same names/order):
#'   \itemize{
#'     \item \code{mean}: conditional means of the target under each partner knockout
#'     \item \code{var}:  conditional variances of the target under each partner knockout
#'   }
#' @export
sweep_partner_knockouts <- function(Sigma,
                                    genes,
                                    target,
                                    knocked = character(0),
                                    mu = NULL) {
  # basic checks
  if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma)) {
    stop("`Sigma` must be a square matrix.")
  }
  if (length(genes) != nrow(Sigma)) {
    stop("`genes` length must match nrow(Sigma).")
  }

  # resolve indices
  to_index <- function(sel) if (is.numeric(sel)) sel else match(sel, genes)

  t_idx <- to_index(target)
  if (length(t_idx) != 1L || is.na(t_idx)) {
    stop("`target` must refer to exactly one valid gene.")
  }

  b_idx <- integer(0)
  if (length(knocked)) {
    b_idx <- to_index(knocked)
    if (any(is.na(b_idx))) stop("Some `knocked` genes not found in `genes`.")
  }

  # candidate partners: all genes except target and already-knocked
  all_idx <- seq_along(genes)
  cand_idx <- setdiff(all_idx, c(t_idx, b_idx))
  if (length(cand_idx) == 0L) {
    warning("No candidate partner genes found (everything excluded).")
    return(list(mean = numeric(0), var = numeric(0)))
  }

  # pre-allocate
  partner_names <- genes[cand_idx]
  means <- numeric(length(cand_idx))
  vars  <- numeric(length(cand_idx))
  names(means) <- partner_names
  names(vars)  <- partner_names

  # loop over partners, add each to the knockout set, call conditioning
  base_knocked <- if (length(b_idx)) genes[b_idx] else character(0)
  for (k in seq_along(cand_idx)) {
    gname <- genes[cand_idx[k]]
    res <- predict_conditional_knockout(
      Sigma  = Sigma,
      genes  = genes,
      target = genes[t_idx],
      knocked = c(base_knocked, gname),
      mu     = mu
    )
    # extract scalar mean/var for the (single) target
    means[k] <- unname(res$mean_cond[1])
    vars[k]  <- unname(res$cov_cond[1, 1])
  }

  list(mean = means, var = vars)
}

#' 3D density plots from a glasso fit, ranked by mean and by variance
#'
#' Takes a glasso-style fit list and produces two 3D plots:
#' (1) top-K partners by conditional mean and (2) top-K partners by conditional variance.
#' Each plot shows the Normal pdf of the target under knockout set {knocked ∪ g}
#' for each selected partner g.
#'
#' @param fit List from run_glasso_seurat()/fit_glasso() (needs $sigma, $features; for
#'   original units also needs $sd, and optionally $mu).
#' @param target Single gene (name or index) to predict.
#' @param knocked Vector of already-knocked genes (names or indices).
#' @param k Integer, how many partners to plot for each ranking (default 15).
#' @param use_original_units TRUE to recover covariance via sd (Σ = D R D). If FALSE, use corr scale.
#' @param renormalize TRUE to cov2cor(fit$sigma) before scaling when recovering Σ.
#' @param xlim Optional length-2 numeric x-range for the pdf curves; NULL = auto.
#' @param n Number of x points per curve (default 200).
#'
#' @return list(mean_plot = plotly_object, var_plot = plotly_object)
#' @export
plot_partner_knockout_densities_dual <- function(
  fit, target, knocked = character(0),
  k = 15,
  use_original_units = TRUE,
  renormalize = TRUE,
  xlim = NULL,
  n = 200
) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. install.packages('plotly')")
  }

  # Build Sigma, genes, mu from the fit
  genes <- fit$features
  if (isTRUE(use_original_units)) {
    if (is.null(fit$sd)) stop("fit$sd missing; set use_original_units=FALSE or refit storing sd.")
    Sigma <- recover_covariance(fit, renormalize = renormalize)
    mu    <- if (!is.null(fit$mu)) fit$mu else rep(0, length(genes))
  } else {
    Sigma <- stats::cov2cor(fit$sigma)
    mu    <- rep(0, length(genes))   # z-scored/corr scale => zero mean
  }

  # Get sweep results (means/vars per partner)
  sv <- sweep_partner_knockouts(Sigma, genes, target, knocked, mu)
  if (length(sv$mean) == 0L) stop("No candidate partners to plot.")

  partners_all <- names(sv$mean)
  m_all <- as.numeric(sv$mean)
  v_all <- as.numeric(sv$var)
  sd_all <- sqrt(pmax(v_all, 0))

  # Helper to build a plot given chosen partner indices
  .build_plot <- function(partners, m, sd, title_suffix) {
    # x-range
    if (is.null(xlim)) {
      x_min <- min(m - 4 * sd); x_max <- max(m + 4 * sd)
      if (!is.finite(x_min) || !is.finite(x_max) || x_min == x_max) {
        x_min <- -3; x_max <- 3
      }
      xrng <- c(x_min, x_max)
    } else {
      xrng <- xlim
    }
    xs <- seq(xrng[1], xrng[2], length.out = n)

    p <- plotly::plot_ly()
    for (i in seq_along(partners)) {
      dens <- stats::dnorm(xs, mean = m[i], sd = ifelse(sd[i] > 0, sd[i], 1e-8))
      p <- plotly::add_trace(
        p,
        x = xs,
        y = rep(i, length(xs)),   # stack by index; label with partner names
        z = dens,
        type = "scatter3d",
        mode = "lines",
        name = sprintf("%s (sd=%.3f)", partners[i], sd[i]),
        line = list(width = 4)
      )
    }
    targ_name <- if (is.numeric(target)) genes[target] else target
    plotly::layout(
      p,
      title = paste("Conditional density of", targ_name, title_suffix),
      scene = list(
        xaxis = list(title = paste0("Target ", targ_name)),
        yaxis = list(title = "Partner",
                     tickmode = "array",
                     tickvals = seq_along(partners),
                     ticktext = partners),
        zaxis = list(title = "Density f(x)")
      ),
      legend = list(title = list(text = "Partner (sd)"))
    )
  }

  # Select top-K by mean (largest means) and by variance (largest variances)
  ord_mean <- order(m_all, decreasing = TRUE)
  idx_mean <- ord_mean[seq_len(min(k, length(ord_mean)))]
  mean_plot <- .build_plot(
    partners = partners_all[idx_mean],
    m = m_all[idx_mean],
    sd = sd_all[idx_mean],
    title_suffix = "(top-K by conditional mean)"
  )

  ord_var <- order(v_all, decreasing = TRUE)
  idx_var <- ord_var[seq_len(min(k, length(ord_var)))]
  var_plot <- .build_plot(
    partners = partners_all[idx_var],
    m = m_all[idx_var],
    sd = sd_all[idx_var],
    title_suffix = "(top-K by conditional variance)"
  )

  list(mean_plot = mean_plot, var_plot = var_plot)
}