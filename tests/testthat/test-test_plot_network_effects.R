test_that("new_GLnode creates correct structure", {
  node <- new_GLnode(index = 5, neighbours = c(2, 7), knocked = TRUE)

  expect_s3_class(node, "GLnode")
  expect_equal(node$index, 5L)
  expect_equal(node$neighbours, c(2L, 7L))
  expect_true(node$knocked)

  # defaults
  node2 <- new_GLnode(index = 10)
  expect_equal(node2$neighbours, integer(0))
  expect_false(node2$knocked)
})


test_that("new_GLgraph constructs graph properly without mocking", {

  # Create covariance matrix where gene 2 is strongly correlated with gene 1
  Sigma <- matrix(c(
    1,   0.9, 0.1,
    0.9, 1,   0.2,
    0.1, 0.2, 1
  ), nrow = 3, byrow = TRUE)

  # Precision matrix for testing (simple structure)
  Omega <- matrix(c(
    1,   0.5, 0,
    0.5, 1,   0.4,
    0,   0.4, 1
  ), nrow = 3, byrow = TRUE)

  fit <- list(
    sigma    = Sigma,
    omega    = Omega,
    features = c("G1", "G2", "G3")
  )

  # Build graph using real get_top_k_gene_indices()
  G <- new_GLgraph(fit, k = 2, target = "G1")

  # Should produce a GLgraph
  expect_s3_class(G, "GLgraph")

  # Top-k genes should be genes 1 and 2
  expect_equal(G$selected_idx, c(1, 2))

  # Should have 2 nodes
  expect_equal(G$k, 2)

  n1 <- G$nodes[[1]]
  n2 <- G$nodes[[2]]

  # Node structure should be GLnode
  expect_s3_class(n1, "GLnode")
  expect_s3_class(n2, "GLnode")

  # Check neighbors based on Omega_sub for indices (1,2):
  # 1 <-> 2 (since omega[1,2] = 0.5)
  expect_equal(n1$neighbours, 2L)
  expect_equal(n2$neighbours, 1L)
})
