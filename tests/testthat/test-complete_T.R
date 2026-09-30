test_that("complete_T recovers direct effects for any trait order", {
  set.seed(1)
  p <- 6
  for (rep in 1:20) {
    # Random DAG in a random (not lower triangular) trait order
    o <- sample(p)
    B <- matrix(0, p, p)
    B[lower.tri(B)] <- rbinom(p * (p - 1) / 2, 1, 0.5) * rnorm(p * (p - 1) / 2)
    B <- B[o, o]
    Tot <- solve(diag(p) - B) - diag(p)
    # Constrained: paths with no direct edge
    s <- which(Tot != 0 & B == 0, arr.ind = TRUE)
    colnames(s) <- c("row", "col")
    res <- complete_T(Tot * (B != 0), s)
    expect_equal(res$B, B)
    expect_equal(res$total_effects, Tot)
  }
})
