test_that("reorder_data() keeps G as identity after a permutation", {
  # dat$G <- diag(1, p) is how esmr_workhorse() initializes G for NESMR (it
  # has no separate factor basis, so G is just a placeholder identity over
  # the trait space). reorder_data() is called whenever order_upper_tri()
  # needs to topologically re-sort traits to make a direct_effect_template
  # lower triangular, and it must keep G an identity matrix no matter the
  # permutation - update_beta_joint()/update_beta_full_joint() use
  # dat$G %*% A %*% t(dat$G) to map between the fitting order and the
  # original trait basis, and a non-identity G there silently looks up the
  # wrong entries of A (or, in the worst case, an entry that is exactly
  # zero, which crashes Matrix::nearPD() with "Matrix seems negative
  # semi-definite").
  p <- 5
  dat <- list(
    Y = matrix(rnorm(p * 3), nrow = 3, ncol = p),
    S = matrix(1, nrow = 3, ncol = p),
    G = diag(1, p)
  )

  cols <- c(3, 1, 4, 2, 5) # arbitrary non-trivial permutation
  reordered <- esmr:::reorder_data(dat, cols)

  expect_equal(reordered$G, diag(1, p))
})

test_that("NESMR direct effects are invariant to trait relabeling that requires DAG reordering", {
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

  # Reference fit: B_correct is already lower triangular in this trait
  # order, so order_upper_tri() does not need to reorder anything.
  fit_ref <- nesmr(
    beta_hat_X = dat$beta_hat,
    se_X = dat$s_estimate,
    variant_ix = ix,
    direct_effect_template = B_correct,
    max_iter = 300,
    params = list(beta_prior_cov = 1)
  )

  # Relabel traits (a self-inverse permutation) so the identical graph now
  # has edges pointing from a higher index to a lower one, forcing
  # order_upper_tri()/reorder_data() to topologically re-sort internally.
  perm <- rev(seq_len(d))
  B_permuted <- B_correct[perm, perm]
  stopifnot(any(B_permuted[upper.tri(B_permuted)] != 0)) # confirm reordering really is required

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
  # relabeling of the same data. Before the reorder_data() fix, G's rows and
  # columns fell out of sync during the internal topo-sort, which silently
  # corrupted (or crashed) fits that actually required reordering.
  expect_equal(
    fit_perm$direct_effects[perm, perm],
    fit_ref$direct_effects,
    tolerance = 1e-4
  )
})
