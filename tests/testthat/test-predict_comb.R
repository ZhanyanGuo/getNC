test_that("predict_conditional_knockout works with names and indices", {
  # simple PD correlation matrix (3x3)
  R <- matrix(c(1, 0.2, -0.1,
                0.2, 1, 0.05,
               -0.1, 0.05, 1), nrow = 3, byrow = TRUE)
  stopifnot(all(eigen(R, symmetric = TRUE)$values > 0))

  genes <- c("G1", "G2", "G3")
  sd    <- c(2, 3, 4)
  D     <- diag(sd, 3)
  Sigma <- D %*% R %*% D  # covariance in original units

  # no mu: defaults to zero
  res1 <- predict_conditional_knockout(
    Sigma  = Sigma,
    genes  = genes,
    target = "G1",
    knocked = c("G2", "G3")
  )

  expect_type(res1, "list")
  expect_true(all(c("mean_cond", "cov_cond") %in% names(res1)))
  expect_equal(names(res1$mean_cond), "G1")
  expect_equal(dimnames(res1$cov_cond), list("G1", "G1"))
  expect_equal(length(res1$mean_cond), 1L)
  expect_equal(nrow(res1$cov_cond), 1L)
  expect_equal(ncol(res1$cov_cond), 1L)

  # with mu: should shift conditional mean appropriately
  mu <- c(1, 2, 3)
  res2 <- predict_conditional_knockout(
    Sigma  = Sigma,
    genes  = genes,
    target = "G1",
    knocked = c("G2", "G3"),
    mu     = mu
  )
  expect_false(isTRUE(all.equal(res1$mean_cond, res2$mean_cond)))

  # indices instead of names
  res3 <- predict_conditional_knockout(
    Sigma  = Sigma,
    genes  = genes,
    target = 1L,          # "G1"
    knocked = 2:3,        # "G2","G3"
    mu     = mu
  )
  expect_equal(res2$mean_cond, res3$mean_cond)
  expect_equal(res2$cov_cond,  res3$cov_cond)
})

test_that("predict_conditional_knockout errors on bad inputs", {
  R <- diag(1, 2)
  genes <- c("A", "B")
  Sigma <- R

  # non-square Sigma
  expect_error(
    predict_conditional_knockout(Sigma = matrix(1, 2, 3), genes = genes, target = "A", knocked = "B"),
    "square"
  )

  # genes length mismatch
  expect_error(
    predict_conditional_knockout(Sigma = Sigma, genes = "A", target = "A", knocked = "B"),
    "match nrow"
  )

  # missing gene names
  expect_error(
    predict_conditional_knockout(Sigma = Sigma, genes = genes, target = "Z", knocked = "B"),
    "not found"
  )
  expect_error(
    predict_conditional_knockout(Sigma = Sigma, genes = genes, target = "A", knocked = "Z"),
    "not found"
  )

  # mu length mismatch
  expect_error(
    predict_conditional_knockout(Sigma = Sigma, genes = genes, target = "A", knocked = "B", mu = 1:3),
    "mu.*length"
  )
})

test_that("recover_covariance returns D * R * D with proper dimnames", {
  # make a correlation and sd, then pack into a 'fit' list
  R <- matrix(c(1, 0.3, 0.1,
                0.3, 1, 0.2,
                0.1, 0.2, 1), nrow = 3, byrow = TRUE)
  stopifnot(all(eigen(R, symmetric = TRUE)$values > 0))

  genes <- c("G1","G2","G3")
  sd    <- c(2, 5, 3)
  fit <- list(
    sigma    = R,       # on correlation scale
    sd       = sd,
    features = genes
  )

  S_cov <- recover_covariance(fit, renormalize = TRUE)
  expect_equal(dim(S_cov), c(3L, 3L))
  expect_equal(rownames(S_cov), genes)
  expect_equal(colnames(S_cov), genes)

  # check elementwise: Sigma_ij = sd_i * R_ij * sd_j
  for (i in 1:3) for (j in 1:3) {
    expect_equal(S_cov[i, j], sd[i] * R[i, j] * sd[j], tolerance = 1e-12)
  }
})

test_that("predict_knockout_from_fit matches direct conditional on original units", {
  # build a small PD correlation, sd, and covariance
  R <- matrix(c(1, 0.25, -0.15,
                0.25, 1, 0.10,
               -0.15, 0.10, 1), nrow = 3, byrow = TRUE)
  stopifnot(all(eigen(R, symmetric = TRUE)$values > 0))

  genes <- c("G1","G2","G3")
  sd    <- c(2, 3, 4)
  D     <- diag(sd, 3)
  Sigma <- D %*% R %*% D

  mu <- c(1, 2, 0.5)

  # pretend this came from run_glasso_seurat()
  fit <- list(
    sigma    = R,        # correlation scale estimate
    sd       = sd,       # per-gene sd
    features = genes,
    mu       = mu        # original-units mean
  )

  # via wrapper (original units)
  out1 <- predict_knockout_from_fit(
    fit,
    target  = "G1",
    knocked = c("G2","G3"),
    use_original_units = TRUE,
    renormalize = TRUE
  )

  # direct call with Sigma + mu should match
  out2 <- predict_conditional_knockout(
    Sigma  = Sigma,
    genes  = genes,
    target = "G1",
    knocked = c("G2","G3"),
    mu     = mu
  )

  expect_equal(out1$mean_cond, out2$mean_cond, tolerance = 1e-10)
  expect_equal(out1$cov_cond,  out2$cov_cond,  tolerance = 1e-10)
})

test_that("predict_knockout_from_fit works on correlation scale when mu=0", {
  # correlation-scale test with zero mean
  R <- matrix(c(1, 0.2, 0.0,
                0.2, 1, 0.3,
                0.0, 0.3, 1), nrow = 3, byrow = TRUE)
  stopifnot(all(eigen(R, symmetric = TRUE)$values > 0))
  genes <- c("G1","G2","G3")

  fit <- list(
    sigma    = R,         # correlation
    features = genes,
    mu       = rep(0, 3)  # z-scored data has mean 0
    # no sd provided -> cannot recover original units; we won't request it
  )

  out <- predict_knockout_from_fit(
    fit,
    target  = "G1",
    knocked = c("G2","G3"),
    use_original_units = FALSE
  )

  # sanity: names and dims
  expect_equal(names(out$mean_cond), "G1")
  expect_equal(dim(out$cov_cond), c(1L, 1L))
})

test_that("predict_knockout_from_fit errors if original units requested without sd", {
  R <- diag(1, 2)
  fit <- list(
    sigma    = R,
    features = c("A","B")
    # sd missing
  )
  expect_error(
    predict_knockout_from_fit(fit, target = "A", knocked = "B",
                              use_original_units = TRUE),
    "fit\\$sd"
  )
})