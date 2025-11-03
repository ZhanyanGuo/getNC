test_that("preprocess_seurat_data filters and normalizes", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  # use sparse to avoid "Coercing to dgCMatrix" warnings
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

  # Seurat v5 uses layer=, v4 uses slot=
  has_layer <- "layer" %in% names(formals(Seurat::GetAssayData))
  if (has_layer) {
    expect_true(nrow(Seurat::GetAssayData(obj2, layer = "data")) > 0)
  } else {
    expect_true(nrow(Seurat::GetAssayData(obj2, slot  = "data")) > 0)
  }
})

test_that("zscore_seurat returns matrix with selected features", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  m <- Matrix::Matrix(rpois(3000, 5), nrow = 300, ncol = 20, sparse = TRUE)
  rownames(m) <- paste0("G", 1:300)
  colnames(m) <- paste0("cell", 1:20)

  obj <- Seurat::CreateSeuratObject(m)
  obj <- Seurat::NormalizeData(obj)

  zmat <- zscore_seurat(obj, nfeatures = 50)

  expect_true(is.matrix(zmat))
  expect_equal(nrow(zmat), 50)
  expect_equal(ncol(zmat), ncol(obj))
})

test_that("run_glasso_seurat returns glasso output incl. mu/sd", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("glasso")
  skip_if_not_installed("Matrix")

  m <- Matrix::Matrix(rpois(3000, 5), nrow = 300, ncol = 30, sparse = TRUE)
  rownames(m) <- paste0("G", 1:300)
  colnames(m) <- paste0("cell", 1:30)

  obj <- Seurat::CreateSeuratObject(m)
  obj <- Seurat::NormalizeData(obj)

  # This should pass now; if not enough features, the function will stop with a clear msg.
  res <- run_glasso_seurat(obj, nfeatures = 40, rho = 0.1)

  expect_true(is.list(res))
  expect_true(all(c("omega","sigma","mu","sd","features","glasso") %in% names(res)))
  p <- nrow(res$omega)
  expect_equal(ncol(res$omega), p)
  expect_equal(nrow(res$sigma),  p)
  expect_equal(ncol(res$sigma),  p)
  expect_equal(length(res$features), p)
  expect_equal(length(res$mu), p)
  expect_equal(length(res$sd), p)
})

test_that("fit_glasso runs on bundled example when obj is NULL", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("glasso")

  res <- fit_glasso(
    obj       = NULL,
    nfeatures = 20,   # small for speed in CI
    rho       = 0.2
  )

  expect_true(is.list(res))
  expect_true(all(c("omega", "sigma", "mu", "sd", "features", "glasso") %in% names(res)))
})