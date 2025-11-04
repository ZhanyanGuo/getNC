test_that("gaussian_1d_group_pvals: null behaves well (t-mean + chi-square)", {
  set.seed(1)
  y <- rnorm(400, mean = 0, sd = 1)
  pv <- gaussian_1d_group_pvals(y, mu = 0, var = 1)

  expect_true(is.numeric(pv) && length(pv) == 2L)
  expect_true(all(names(pv) == c("p_mean", "p_chisq")))
  # Under the null, p-values shouldn't be extreme
  expect_gt(pv["p_mean"], 0.01); expect_lt(pv["p_mean"], 0.99)
  expect_gt(pv["p_chisq"], 0.01); expect_lt(pv["p_chisq"], 0.99)
})

test_that("gaussian_1d_group_pvals: mean shift triggers t-mean (decoupled from var)", {
  set.seed(2)
  y <- rnorm(200, mean = 0.5, sd = 1)  # shift in mean only
  pv <- gaussian_1d_group_pvals(y, mu = 0, var = 1)
  expect_lt(pv["p_mean"], 0.01)       # t-mean should flag it
  # chi-square might also be small but we don't require it here
})

test_that("gaussian_1d_group_pvals: variance inflation triggers chi-square", {
  set.seed(3)
  y <- rnorm(200, mean = 0, sd = 1.7)  # variance change
  pv <- gaussian_1d_group_pvals(y, mu = 0, var = 1)
  expect_lt(pv["p_chisq"], 0.01)
})

test_that("gaussian_1d_group_pvals: degenerate variance var=0 handled", {
  y_equal    <- rep(2, 5)
  pv_equal   <- gaussian_1d_group_pvals(y_equal, mu = 2, var = 0)
  expect_equal(unname(pv_equal["p_chisq"]), 1)
  expect_equal(unname(pv_equal["p_mean"]), 1)

  y_notequal <- c(2, 2, 2, 3, 2)
  pv_noteq   <- gaussian_1d_group_pvals(y_notequal, mu = 2, var = 0)
  expect_equal(unname(pv_noteq["p_chisq"]), 0)
  expect_true(is.finite(pv_noteq["p_mean"]))
  expect_gte(unname(pv_noteq["p_mean"]), 0)
  expect_lte(unname(pv_noteq["p_mean"]), 1)
})
test_that("gaussian_1d_group_pvals: input validation", {
  expect_error(gaussian_1d_group_pvals(numeric(0), mu = 0, var = 1), "at least one")
  expect_error(gaussian_1d_group_pvals(1:3, mu = 0, var = -1), ">= 0")
})

test_that("validate_knockout_group_from_fit: works in original units (synthetic fit)", {
  # Build small positive-definite correlation, sd, covariance, mean
  R <- matrix(c(1,   0.25, -0.10,
                0.25, 1,    0.15,
               -0.10, 0.15, 1), 3, 3, byrow = TRUE)
  stopifnot(all(eigen(R, symmetric = TRUE)$values > 0))
  genes <- c("G1", "G2", "G3")
  sd    <- c(2, 3, 1.5)
  D     <- diag(sd, 3)
  Sigma <- D %*% R %*% D
  mu    <- c(0.3, -0.2, 0.1)

  # "fit" in the same style as run_glasso_seurat() (sigma on corr scale + sd)
  fit <- list(
    sigma    = R,
    sd       = sd,
    features = genes,
    mu       = mu
  )

  # Theoretical conditional for target=G1 | G2=0
  i1 <- 1; i2 <- 2
  mu_c  <- mu[i1] + Sigma[i1, i2] * (1 / Sigma[i2, i2]) * (0 - mu[i2])
  var_c <- Sigma[i1, i1] - Sigma[i1, i2] * (1 / Sigma[i2, i2]) * Sigma[i2, i1]

  set.seed(10)
  y <- rnorm(300, mean = mu_c, sd = sqrt(var_c))

  res <- validate_knockout_group_from_fit(
    fit,
    target  = "G1",
    knocked = "G2",
    y       = y,
    use_original_units = TRUE
  )

  expect_true(is.list(res))
  expect_equal(res$target, "G1")
  expect_equal(res$knocked, "G2")
  # Not extreme under the true model
  expect_gt(res$p_mean, 0.01); expect_lt(res$p_mean, 0.99)
  expect_gt(res$p_chisq, 0.01); expect_lt(res$p_chisq, 0.99)
  expect_equal(res$n, length(y))
})

test_that("validate_knockout_group_from_fit: correlation scale also OK", {
  R <- matrix(c(1, 0.2, 0.1,
                0.2, 1, 0.3,
                0.1, 0.3, 1), 3, 3, byrow = TRUE)
  genes <- c("A", "B", "C")

  fit <- list(
    sigma    = R,            # use correlation directly
    features = genes,
    mu       = rep(0, 3)     # zero mean on corr scale
  )

  # target A | B=0 : var = 1 - rho^2
  iA <- 1; iB <- 2
  mu_c  <- 0
  var_c <- 1 - R[iA, iB]^2

  set.seed(11)
  y <- rnorm(250, mean = mu_c, sd = sqrt(var_c))

  res <- validate_knockout_group_from_fit(
    fit,
    target  = "A",
    knocked = "B",
    y       = y,
    use_original_units = FALSE
  )
  expect_gt(res$p_mean, 0.01); expect_lt(res$p_mean, 0.99)
  expect_gt(res$p_chisq, 0.01); expect_lt(res$p_chisq, 0.99)
})

test_that("validate_knockout_group_from_fit: errors on multi-target", {
  R <- diag(1, 3)
  fit <- list(sigma = R, features = c("G1","G2","G3"), mu = c(0,0,0))
  expect_error(
    validate_knockout_group_from_fit(
      fit, target = c("G1","G2"), knocked = "G3", y = rnorm(10),
      use_original_units = FALSE
    ),
    "exactly one gene"
  )
})