##############################################
# Test: preprocess_seurat_data
##############################################

test_that("preprocess_seurat_data filters and normalizes", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  m <- Matrix::Matrix(rpois(2000, 5), nrow = 200, ncol = 10, sparse = TRUE)
  rownames(m) <- paste0("G", 1:200)
  colnames(m) <- paste0("cell", 1:10)

  obj  <- Seurat::CreateSeuratObject(m)
  obj2 <- preprocess_seurat_data(
    obj,
    min_genes = 0,
    max_genes = 5000,
    max_mt    = 100,
    normalize = TRUE
  )

  expect_s4_class(obj2, "Seurat")
  expect_true("percent.mt" %in% colnames(obj2[[]]))

  has_layer <- "layer" %in% names(formals(Seurat::GetAssayData))
  if (has_layer) {
    expect_true(nrow(Seurat::GetAssayData(obj2, layer = "data")) > 0)
  } else {
    expect_true(nrow(Seurat::GetAssayData(obj2, slot  = "data")) > 0)
  }
})



##############################################
# Test: run_glasso (Seurat path via dispatch)
##############################################

test_that("run_glasso correctly handles Seurat input", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("glasso")
  skip_if_not_installed("Matrix")

  m <- Matrix::Matrix(rpois(3000, 5), nrow = 300, ncol = 30, sparse = TRUE)
  rownames(m) <- paste0("G", 1:300)
  colnames(m) <- paste0("cell", 1:30)

  obj <- Seurat::CreateSeuratObject(m)
  obj <- Seurat::NormalizeData(obj)

  res <- run_glasso(obj, nfeatures = 40, rho = 0.1)

  expect_true(is.list(res))
  expect_equal(res$type, "seurat")
  expect_true(all(c("omega","sigma","mu","sd","features","glasso") %in% names(res)))

  p <- nrow(res$omega)
  expect_equal(ncol(res$omega), p)
  expect_equal(nrow(res$sigma),  p)
  expect_equal(ncol(res$sigma),  p)
  expect_equal(length(res$features), p)
  expect_equal(length(res$mu), p)
  expect_equal(length(res$sd), p)
})



##############################################
# Test: fit_glasso() (NULL → bundled example)
##############################################

test_that("fit_glasso runs on example data when obj = NULL", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("glasso")

  res <- fit_glasso(
    obj       = NULL,
    nfeatures = 20,
    rho       = 0.2
  )

  expect_true(is.list(res))
  expect_true(all(c("omega", "sigma", "mu", "sd", "features", "glasso") %in% names(res)))
})



##############################################
# Test: fit_glasso_raw()
##############################################

test_that("fit_glasso_raw works on raw matrix input", {
  skip_if_not_installed("glasso")

  set.seed(1)
  mat <- matrix(rpois(3000, 4), nrow = 100, ncol = 30)
  rownames(mat) <- paste0("cell", 1:100)
  colnames(mat) <- paste0("G", 1:30)

  res <- fit_glasso_raw(
    mat        = mat,
    nfeatures  = 10,
    rho        = 0.05,
    min_genes  = 0,
    max_genes  = 5000,
    max_mt     = 100,
    normalize  = TRUE
  )

  expect_true(is.list(res))
  expect_equal(res$type, "matrix")
  expect_true(all(c("omega","sigma","mu","sd","features","glasso") %in% names(res)))
})



##############################################
# Test: auto_fit_glasso() — Seurat dispatch
##############################################

test_that("auto_fit_glasso dispatches Seurat correctly", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("glasso")
  skip_if_not_installed("Matrix")

  m <- Matrix::Matrix(rpois(2000, 4), nrow = 100, ncol = 20, sparse = TRUE)
  rownames(m) <- paste0("G", 1:100)
  colnames(m) <- paste0("cell", 1:20)

  obj <- Seurat::CreateSeuratObject(m)
  obj <- Seurat::NormalizeData(obj)

  # relax QC for synthetic data
  res <- auto_fit_glasso(
    obj,
    min_genes = 0,
    max_genes = 1e6,
    max_mt    = 100,
    nfeatures = 10,
    rho       = 0.1
  )

  expect_equal(res$type, "seurat")
  expect_true(all(c("omega","sigma","mu","sd","features","glasso") %in% names(res)))
})



##############################################
# Test: auto_fit_glasso() — raw matrix dispatch
##############################################

test_that("auto_fit_glasso dispatches raw matrix correctly", {
  skip_if_not_installed("glasso")

  mat <- matrix(rpois(2000, 4), nrow = 100, ncol = 20)
  rownames(mat) <- paste0("cell", 1:100)
  colnames(mat) <- paste0("G", 1:20)

  res <- auto_fit_glasso(
    mat,
    min_genes = 0,
    max_genes = 1e6,
    max_mt    = 100,
    nfeatures = 5,   # pick < ncol(mat) = 20
    rho       = 0.1
  )

  expect_equal(res$type, "matrix")
  expect_true(all(c("omega","sigma","mu","sd","features","glasso") %in% names(res)))
})



##############################################
# Test: auto_fit_glasso() — NULL → fit_glasso()
##############################################

test_that("auto_fit_glasso handles NULL by calling fit_glasso", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("glasso")

  res <- auto_fit_glasso(NULL, nfeatures = 10, rho = 0.2)

  expect_equal(res$type, "seurat")   # fit_glasso() uses Seurat pathway
  expect_true(all(c("omega","sigma","mu","sd","features","glasso") %in% names(res)))
})
