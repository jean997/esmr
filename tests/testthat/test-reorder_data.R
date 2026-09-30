test_that("NESMR direct effects are invariant to trait relabeling", {
  skip_if_not_installed("GWASBrewer")

  # Same graph used in test-four-node-example.R.
  G <- matrix(
    c(0, 0, 0, 0,
      0.4, 0, 0, 0,
      0, 0, 0, 0,
      0, -0.5, 0.2, 0),
    nrow = 4,
    byrow = TRUE
  )
  d <- ncol(G)
  B_correct <- (G != 0) + 0
  B_lower <- lower.tri(B_correct) + 0

  dat <- GWASBrewer::sim_mv(
    G = G,
    N = 20000,
    J = 5000,
    h2 = 0.2,
    pi = 0.1,
    sporadic_pleiotropy = TRUE,
    est_s = TRUE
  )

  Z <- dat$beta_hat / dat$s_estimate
  pval_select <- 2 * pnorm(-abs(Z))
  minp <- apply(pval_select, 1, min)
  ix <- which(minp < 5e-8)

  # Reference fit: B_correct is lower triangular in this trait order.
  fit_ref <- nesmr(
    beta_hat_X = dat$beta_hat,
    se_X = dat$s_estimate,
    variant_ix = ix,
    direct_effect_template = B_lower,
    max_iter = 300,
    params = list(beta_prior_cov = 1)
  )

  # Relabel traits (a self-inverse permutation) so the identical graph now
  # has edges pointing from a higher index to a lower one.
  perm <- rev(seq_len(d))
  B_permuted <- B_lower[perm, perm]
  stopifnot(any(B_permuted[upper.tri(B_permuted)] != 0)) # confirm template is not lower triangular

  fit_perm <- nesmr(
    beta_hat_X = dat$beta_hat[, perm],
    se_X = dat$s_estimate[, perm],
    variant_ix = ix,
    direct_effect_template = B_permuted,
    max_iter = 300,
    params = list(beta_prior_cov = 1)
  )

  # Relabeling the input and then relabeling the output back should recover
  # the reference fit exactly, since this is the same likelihood under a
  # relabeling of the same data.
  expect_equal(
    fit_perm$direct_effects[perm, perm],
    fit_ref$direct_effects,
    tolerance = 1e-4
  )

  expect_equal(
    fit_perm$elbo,
    fit_ref$elbo,
    tolerance = 1e-4
  )
})

test_that("NESMR direct effects are invariant to trait relabeling - with complete_T function", {
  skip_if_not_installed("GWASBrewer")

  # Same graph used in test-four-node-example.R.
  G <- matrix(
    c(0, 0, 0, 0,
      0.4, 0, 0, 0,
      0, 0, 0, 0,
      0, -0.5, 0.2, 0),
    nrow = 4,
    byrow = TRUE
  )
  d <- ncol(G)
  B_correct <- (G != 0) + 0

  dat <- GWASBrewer::sim_mv(
    G = G,
    N = 20000,
    J = 5000,
    h2 = 0.2,
    pi = 0.1,
    sporadic_pleiotropy = TRUE,
    est_s = TRUE
  )

  Z <- dat$beta_hat / dat$s_estimate
  pval_select <- 2 * pnorm(-abs(Z))
  minp <- apply(pval_select, 1, min)
  ix <- which(minp < 5e-8)

  # Reference fit: B_correct is lower triangular in this trait order.
  fit_ref <- nesmr(
    beta_hat_X = dat$beta_hat,
    se_X = dat$s_estimate,
    variant_ix = ix,
    direct_effect_template = B_correct,
    max_iter = 300,
    params = list(beta_prior_cov = 1)
  )

  # Relabel traits (a self-inverse permutation) so the identical graph now
  # has edges pointing from a higher index to a lower one.
  perm <- rev(seq_len(d))
  B_permuted <- B_correct[perm, perm]
  stopifnot(any(B_permuted[upper.tri(B_permuted)] != 0)) # confirm template is not lower triangular

  fit_perm <- nesmr(
    beta_hat_X = dat$beta_hat[, perm],
    se_X = dat$s_estimate[, perm],
    variant_ix = ix,
    direct_effect_template = B_permuted,
    max_iter = 300,
    params = list(beta_prior_cov = 1)
  )

  # Relabeling the input and then relabeling the output back should recover
  # the reference fit exactly, since this is the same likelihood under a
  # relabeling of the same data.
  expect_equal(
    fit_perm$direct_effects[perm, perm],
    fit_ref$direct_effects,
    tolerance = 1e-3
  )

  expect_equal(
    fit_perm$elbo,
    fit_ref$elbo,
    tolerance = 1e-4
  )
})

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
