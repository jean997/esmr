# Setup graph
G <- matrix(
  c(0, 0, 0, 0,
    0.4, 0, 0, 0,
    0, 0, 0, 0,
    0, -0.5, 0.2, 0),
  nrow = 4,
  byrow = TRUE
)

B_lower <- lower.tri(G) + 0
B_correct <- (G != 0) + 0

test_that("complete_mvmr_cahced and non-cached versions produce similar results", {
  skip_if_not_installed("GWASBrewer")

  d <- ncol(G)

  # Simulate data
  dat <- GWASBrewer::sim_mv(
    G = G,
    N = 40000,
    J = 5000,
    h2 = 0.2,
    pi = 0.1,
    sporadic_pleiotropy = TRUE,
    est_s = TRUE
  )

  # Select variants
  Z <- dat$beta_hat / dat$s_estimate
  pval_select <- 2 * pnorm(-abs(Z))
  minp <- apply(pval_select, 1, min)
  ix <- which(minp < 5e-8)

  # Run NESMR with full graph
  nesmr_complete <- nesmr_complete(
    beta_hat = dat$beta_hat,
    se_beta_hat = dat$s_estimate
  )

  nesmr_complete_cached <- nesmr_complete_mvmr_cached_gfa(
    beta_hat = dat$beta_hat,
    se_beta_hat = dat$s_estimate
  )

  nesmr_complete$beta_hat - nesmr_complete_cached$beta_hat

  # Check output structure
  #expect_s3_class(nesmr_mod, "nesmr")
  #expect_true(!is.null(nesmr_mod$direct_effects))
  #expect_true(!is.null(nesmr_mod$standard_error_method))
})