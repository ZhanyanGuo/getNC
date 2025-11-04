test_that("end-to-end workflow runs on packaged example", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("glasso")
  # plotly is optional; we guard those checks below

  # ---- Step 0: fit glasso on default packaged data ----
  set.seed(123)
  fit <- fit_glasso(
    obj       = NULL,   # uses the package's bundled example dataset
    nfeatures = 100,    # keep small for CI speed
    rho       = 0.2
  )

  # basic structure checks
  expect_true(is.list(fit))
  expect_true(all(c("omega","sigma","features") %in% names(fit)))
  expect_true(is.matrix(fit$omega) && is.matrix(fit$sigma))
  expect_true(length(fit$features) == nrow(fit$omega))
  expect_equal(dim(fit$omega), dim(fit$sigma))

  # ---- Step 2: conditional prediction for a target | knocked ----
  genes   <- fit$features
  expect_true(length(genes) >= 10)
  target  <- genes[5]
  knocked <- genes[c(3, 7, 9)]

  res <- predict_knockout_from_fit(
    fit,
    target  = target,
    knocked = knocked,
    use_original_units = TRUE
  )

  # Should be scalar target: 1 mean, 1x1 variance
  expect_true(is.numeric(res$mean_cond) && length(res$mean_cond) == 1L)
  expect_true(is.matrix(res$cov_cond) && all(dim(res$cov_cond) == c(1L, 1L)))
  expect_true(is.finite(res$mean_cond[1]))
  expect_true(is.finite(res$cov_cond[1,1]) && res$cov_cond[1,1] >= 0)

  # ---- Step 3: optional plotting of partner knockouts ----
  if (requireNamespace("plotly", quietly = TRUE)) {
    plots <- plot_partner_knockout_densities_dual(
      fit,
      target  = target,
      knocked = knocked,
      k = 10,
      use_original_units = TRUE
    )
    expect_true(is.list(plots))
    expect_true(all(c("mean_plot","var_plot") %in% names(plots)))
    # plotly objects usually carry class c("plotly","htmlwidget")
    expect_true(inherits(plots$mean_plot, "plotly"))
    expect_true(inherits(plots$var_plot, "plotly"))
  } else {
    skip("plotly not installed; skipping visualization part of integration test.")
  }

  # ---- Step 4: validation with synthetic observations ----
  mu_c  <- as.numeric(res$mean_cond[1])
  var_c <- as.numeric(res$cov_cond[1, 1])
  sd_c  <- sqrt(max(var_c, 0))

  set.seed(42)
  y_obs <- rnorm(40, mean = mu_c, sd = sd_c)

  val <- validate_knockout_group_from_fit(
    fit,
    target  = target,
    knocked = knocked,
    y       = y_obs,
    use_original_units = TRUE
  )

  # structure + bounds
  expect_true(is.list(val))
  expect_true(all(c("p_mean","p_chisq","n","mean_cond","var_cond") %in% names(val)))
  expect_true(is.finite(val$p_mean) && val$p_mean >= 0 && val$p_mean <= 1)
  expect_true(is.finite(val$p_chisq) && val$p_chisq >= 0 && val$p_chisq <= 1)
  expect_equal(val$n, length(y_obs))

  # Under true model, p-values should not be extremely small with high probability
  expect_gt(val$p_mean, 0.001)
  expect_gt(val$p_chisq, 0.001)
})