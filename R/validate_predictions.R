#' Two p-values for a 1D Gaussian check (t-mean + chi-square)
#'
#' Given observations y and a model mean/variance (mu, var):
#' 1) t-test on the mean using sample sd (decoupled from `var`) — two-sided.
#' 2) Chi-square on sum((y - mu)^2 / var) with df = n — two-sided.
#'
#' @param y   numeric vector of observations.
#' @param mu  scalar mean under H0.
#' @param var scalar variance under H0 (>= 0) used only for the chi-square part.
#'
#' @return named numeric vector: c(p_tmean = ..., p_chisq = ...)
#' @keywords internal
gaussian_1d_group_pvals <- function(y, mu, var) {
  y <- as.numeric(y)
  if (!length(y)) stop("`y` must have at least one value.")
  if (!is.finite(var) || var < 0) stop("`var` must be finite and >= 0.")

  n  <- length(y)
  m  <- mean(y)
  s2 <- if (n > 1) stats::var(y) else NA_real_
  s  <- sqrt(s2)

  # t-test for mean using sample sd (decoupled from model var)
  if (n < 2 || !is.finite(s) || s == 0) {
    # t-mean undefined when n<2 or zero sample variance; fall back to exact match
    p_mean <- as.numeric(m == mu)
  } else {
    tstat   <- sqrt(n) * (m - mu) / s
    p_mean <- 2 * stats::pt(-abs(tstat), df = n - 1)
  }

  # Chi-square on sum of squared deviations scaled by model var
  if (var == 0) {
    p_chisq <- as.numeric(all(y == mu))
  } else {
    z2_sum  <- sum((y - mu)^2 / var)
    p_one   <- stats::pchisq(z2_sum, df = n)
    p_chisq <- 2 * min(p_one, 1 - p_one)
    p_chisq <- min(1, max(0, p_chisq))
  }

  c(p_mean = p_mean, p_chisq = p_chisq)
}

#' Group tests for a knockout using a glasso-style fit (1D target)
#'
#' For a single target gene and a set of knocked genes, compute the model's
#' conditional Normal N(mu_cond, var_cond) from `fit`, then evaluate a vector
#' of observations `y` with two omnibus p-values:
#'  - p_mean  : t test on mean(z), two-sided
#'  - p_chisq : Chi-square test on sum(z^2), two-sided
#'
#' @param fit   List from run_glasso_seurat()/fit_glasso() (needs $sigma, $features;
#'              for original units also needs $sd and optionally $mu).
#' @param target  Single gene (name or index) whose expression was measured.
#' @param knocked Character or integer vector of genes set to 0 in the experiment.
#' @param y       Numeric vector of observed target expression under this knockout set.
#' @param use_original_units Logical; if TRUE (default) recover covariance in original units;
#'                           if FALSE, operate on correlation scale.
#' @param renormalize Logical; passed to recover_covariance() when use_original_units=TRUE.
#'
#' @return A list with:
#'   \describe{
#'     \item{target}{target gene name}
#'     \item{knocked}{character vector of knocked genes}
#'     \item{mean_cond}{scalar conditional mean from the model}
#'     \item{var_cond}{scalar conditional variance from the model}
#'     \item{p_mean}{two-sided t test p-value on mean(z)}
#'     \item{p_chisq}{two-sided Chi-square test p-value on sum(z^2)}
#'     \item{n}{number of observations}
#'   }
#' @export
validate_knockout_group_from_fit <- function(
  fit,
  target,
  knocked,
  y,
  use_original_units = TRUE,
  renormalize = TRUE
) {
  # Get conditional 1D Gaussian for the target under the knockout set
  cond <- predict_knockout_from_fit(
    fit,
    target    = target,
    knocked   = knocked,
    use_original_units = use_original_units,
    renormalize = renormalize
  )

  # enforce scalar target
  if (length(cond$mean_cond) != 1L || nrow(cond$cov_cond) != 1L || ncol(cond$cov_cond) != 1L) {
    stop("`target` must correspond to exactly one gene for 1D validation.")
  }

  mu_c  <- as.numeric(cond$mean_cond[1])
  var_c <- as.numeric(cond$cov_cond[1, 1])

  # Run the two omnibus tests on the data vector
  pv <- gaussian_1d_group_pvals(y, mu = mu_c, var = var_c)

  # names
  genes <- fit$features
  tgt_name <- if (is.numeric(target)) genes[target] else as.character(target)
  knocked_names <- if (is.numeric(knocked)) genes[knocked] else as.character(knocked)

  list(
    target    = tgt_name,
    knocked   = knocked_names,
    mean_cond = mu_c,
    var_cond  = var_c,
    p_mean    = unname(pv["p_mean"]),
    p_chisq   = unname(pv["p_chisq"]),
    n         = length(as.numeric(y))
  )
}