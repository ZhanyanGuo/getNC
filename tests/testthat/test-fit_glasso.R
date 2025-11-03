test_that("preprocess_seurat_data filters and normalizes", {
  skip_if_not_installed("Seurat")

  m <- matrix(rpois(2000, 5), nrow = 200, ncol = 10)
  rownames(m) <- paste0("G", 1:200)
  colnames(m) <- paste0("cell", 1:10)

  obj <- Seurat::CreateSeuratObject(m)

  obj2 <- preprocess_seurat_data(obj,
                                 min_genes = 0,
                                 max_genes = 5000,
                                 max_mt    = 100,
                                 normalize = TRUE)

  expect_true(inherits(obj2, "Seurat"))
  expect_true("percent.mt" %in% colnames(obj2[[]]))
  expect_true(nrow(Seurat::GetAssayData(obj2, slot = "data")) > 0)
})

test_that("zscore_seurat returns matrix with selected features", {
  skip_if_not_installed("Seurat")

  m <- matrix(rpois(3000, 5), nrow = 300, ncol = 20)
  rownames(m) <- paste0("G", 1:300)
  colnames(m) <- paste0("cell", 1:20)

  obj <- Seurat::CreateSeuratObject(m)
  obj <- Seurat::NormalizeData(obj)

  zmat <- zscore_seurat(obj, nfeatures = 50)

  expect_true(is.matrix(zmat))
  expect_equal(nrow(zmat), 50)
  expect_equal(ncol(zmat), ncol(obj))
})

test_that("run_glasso_seurat returns glasso output", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("glasso")

  m <- matrix(rpois(3000, 5), nrow = 300, ncol = 30)
  rownames(m) <- paste0("G", 1:300)
  colnames(m) <- paste0("cell", 1:30)

  obj <- Seurat::CreateSeuratObject(m)
  obj <- Seurat::NormalizeData(obj)

  res <- run_glasso_seurat(obj, nfeatures = 40, rho = 0.1)

  expect_true(is.list(res))
  expect_true(all(c("omega", "sigma", "features", "glasso") %in% names(res)))
  expect_equal(nrow(res$omega), ncol(res$omega))
  expect_equal(length(res$features), nrow(res$omega))
})

