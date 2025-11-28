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
#' @return named numeric vector: c(p_mean = ..., p_chisq = ...)
#'
#' @examples
#' \dontrun{
#'   # Synthetic data under N(0, 1)
#'   set.seed(1)
#'   y <- rnorm(30, mean = 0, sd = 1)
#'   pv <- gaussian_1d_group_pvals(y, mu = 0, var = 1)
#'   print(pv)
#' }
#'
#' @references
#'   Student (1908). "The probable error of a mean." \emph{Biometrika}, 6(1), 1–25.
#'
#'   Fisher, R.A. (1924). "On a distribution yielding the error functions of
#'   several well known statistics." \emph{Proceedings of the International
#'   Congress of Mathematics}, 2, 805–813. (Chi-square usage)
#'
#'   Casella, G. & Berger, R.L. (2002). \emph{Statistical Inference} (2nd ed.).
#'   Duxbury. (Classical t and chi-square results)
#'
#' @keywords internal
#' @noRd
#' @importFrom stats var pchisq pt rnorm
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
#'  - p_mean  : two-sided t test on the sample mean (using the sample sd)
#'  - p_chisq : two-sided chi-square test on sum((y - mu_cond)^2 / var_cond)
#'
#' @param fit   List from auto_fit_glasso()/fit_glasso() (needs $sigma, $features;
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
#'     \item{p_mean}{two-sided t test p-value on mean}
#'     \item{p_chisq}{two-sided Chi-square test p-value on sum of squares}
#'     \item{n}{number of observations}
#'   }
#'
#' @examples
#' \dontrun{
#'   suppressPackageStartupMessages(library(getNC))
#'
#'   # 1) Fit a model using the packaged example data
#'   fit <- fit_glasso(obj = NULL, nfeatures = 100, rho = 0.2)
#'
#'   genes   <- fit$features
#'   target  <- genes[5]
#'   knocked <- genes[c(3, 7, 9)]
#'
#'   # 2) Pull the conditional distribution for target | knocked=0
#'   cond <- predict_knockout_from_fit(
#'     fit, target = target, knocked = knocked, use_original_units = TRUE
#'   )
#'   mu_c  <- as.numeric(cond$mean_cond[1])
#'   var_c <- as.numeric(cond$cov_cond[1, 1])
#'
#'   # 3) Simulate replicate measurements under that model
#'   # in practice, use knockout experiment 1D RNA counting vector
#'   set.seed(42)
#'   y_obs <- rnorm(40, mean = mu_c, sd = sqrt(max(var_c, 0)))
#'
#'   # 4) Validate (t-mean + chi-square)
#'   res <- validate_knockout_group_from_fit(
#'     fit, target = target, knocked = knocked, y = y_obs, use_original_units = TRUE
#'   )
#'   str(res)
#'   cat(sprintf("p_mean = %.3g, p_chisq = %.3g\n", res$p_mean, res$p_chisq))
#' }
#'
#' @seealso
#'   \code{\link{fit_glasso}}, \code{\link{auto_fit_glasso}},
#'   \code{\link{predict_knockout_from_fit}}, \code{\link{recover_covariance}}
#'
#' @references
#'   Friedman, J., Hastie, T., & Tibshirani, R. (2008).
#'   Sparse inverse covariance estimation with the graphical lasso.
#'   \emph{Biostatistics}, 9(3), 432–441. \doi{10.1093/biostatistics/kxm045}
#'
#'   Hao, Y., Hao, S., Andersen-Nissen, E., et al. (2021).
#'   Integrated analysis of multimodal single-cell data.
#'   \emph{Cell}, 184(13), 3573–3587. \doi{10.1016/j.cell.2021.04.048}
#'
#'   Bishop, C. M. (2006). \emph{Pattern Recognition and Machine Learning}.
#'   Springer. (Multivariate Normal conditioning identities)
#'
#'   Murphy, K. P. (2012). \emph{Machine Learning: A Probabilistic Perspective}.
#'   MIT Press. (Conditioning and marginalization in Gaussians)
#'
#'   Student (1908). "The probable error of a mean." \emph{Biometrika}, 6(1), 1–25.
#'
#'   Fisher, R.A. (1924). "On a distribution yielding the error functions of
#'   several well known statistics." \emph{Proceedings of the ICM}, 2, 805–813.
#'
#' @importFrom stats rnorm
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

# Gen Ai used for documentation and input verify