test_that("NESMR runs with four-node DAG graph", {
  skip_if_not_installed("GWASBrewer")

  # Setup graph
  G <- matrix(
    c(0, 0, 0, 0,
      0.4, 0, 0, 0,
      0, 0, 0, 0,
      0, -0.5, 0.2, 0),
    nrow = 4,
    byrow = TRUE
  )

  d <- ncol(G)

  # Simulate data
  dat <- GWASBrewer::sim_mv(
    G = G,
    N = 20000,
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
  B_full <- matrix(1, d, d) - diag(d)
  nesmr_mod <- nesmr(
    beta_hat_X = dat$beta_hat,
    se_X = dat$s_estimate,
    direct_effect_template = B_full,
    variant_ix = ix,
    params = list(
      beta_prior_cov = 1,
      restrict_dag = FALSE
    )
  )

  # Check output structure
  expect_s3_class(nesmr_mod, "nesmr")
  expect_true(!is.null(nesmr_mod$direct_effects))
  expect_true(!is.null(nesmr_mod$standard_error_method))
})

test_that("NESMR runs with DAG restriction", {
  skip_if_not_installed("GWASBrewer")

  # Setup graph
  G <- matrix(
    c(0, 0, 0, 0,
      0.4, 0, 0, 0,
      0, 0, 0, 0,
      0, -0.5, 0.2, 0),
    nrow = 4,
    byrow = TRUE
  )

  d <- ncol(G)

  # Simulate data
  dat <- GWASBrewer::sim_mv(
    G = G,
    N = 20000,
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

  # Run NESMR with DAG restriction
  B_full <- matrix(1, d, d) - diag(d)
  # Error expected here due to non-DAG structure
  testthat::expect_error({
      nesmr_mod <- nesmr(
        beta_hat_X = dat$beta_hat,
        se_X = dat$s_estimate,
        direct_effect_template = B_full,
        variant_ix = ix,
        max_iter = 300,
        params = list(
          beta_prior_cov = 1
        )
      )
  })

  nesmr_mod <- nesmr(
    beta_hat_X = dat$beta_hat,
    se_X = dat$s_estimate,
    direct_effect_template = B_lower,
    variant_ix = ix,
    max_iter = 300,
    params = list(
      beta_prior_cov = 1
    )
  )


  # Check output structure
  expect_s3_class(nesmr_mod, "nesmr")
  expect_true(!is.null(nesmr_mod$direct_effects))
})

test_that("NESMR runs with correct DAG structure", {
  skip_if_not_installed("GWASBrewer")

  # Setup graph
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

  # Simulate data
  dat <- GWASBrewer::sim_mv(
    G = G,
    N = 20000,
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

  # Run NESMR with correct template
  nesmr_mod <- nesmr(
    beta_hat_X = dat$beta_hat,
    se_X = dat$s_estimate,
    variant_ix = ix,
    direct_effect_template = B_correct,
    max_iter = 300,
    params = list(
      beta_prior_cov = 1
    )
  )

  # Check output structure
  expect_s3_class(nesmr_mod, "nesmr")
  expect_true(!is.null(nesmr_mod$direct_effects))
})

test_that("MVMR complete works", {
  skip_if_not_installed("GWASBrewer")

  # Setup graph
  G <- matrix(
    c(0, 0, 0, 0,
      0.4, 0, 0, 0,
      0, 0, 0, 0,
      0, -0.5, 0.2, 0),
    nrow = 4,
    byrow = TRUE
  )

  d <- ncol(G)

  # Simulate data
  dat <- GWASBrewer::sim_mv(
    G = G,
    N = 20000,
    J = 5000,
    h2 = 0.2,
    pi = 0.1,
    sporadic_pleiotropy = TRUE,
    est_s = TRUE
  )

  # Select variants
  Z <- dat$beta_hat / dat$s_estimate
  pval_select <- 2 * pnorm(-abs(Z))

  # Run complete MVMR
  mvmr_res <- nesmr_complete_mvmr(
    beta_hat = dat$beta_hat,
    se_beta_hat = dat$s_estimate,
    pval_select = pval_select,
    alpha = 5e-8
  )

  # Check output structure
  expect_true(is.list(mvmr_res))
  expect_true(!is.null(mvmr_res))
  expect_equal(dim(mvmr_res$beta_hat), c(d, d))
  expect_equal(diag(mvmr_res$beta_hat), rep(0, d))
  expect_equal(dim(mvmr_res$se_beta_hat), c(d, d))
  expect_equal(diag(mvmr_res$se_beta_hat), rep(0, d))
})

test_that("Parametric bootstrap method works on nesmr object", {
  skip_if_not_installed("GWASBrewer")

  # Setup graph
  G <- matrix(
    c(0, 0, 0, 0,
      0.4, 0, 0, 0,
      0, 0, 0, 0,
      0, -0.5, 0.2, 0),
    nrow = 4,
    byrow = TRUE
  )

  d <- ncol(G)
  B_full <- matrix(1, d, d) - diag(d)

  # Simulate data
  dat <- GWASBrewer::sim_mv(
    G = G,
    N = 20000,
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

  # Run NESMR with non-DAG
  nesmr_mod <- nesmr(
    beta_hat_X = dat$beta_hat,
    se_X = dat$s_estimate,
    direct_effect_template = B_full,
    variant_ix = ix,
    params = list(
      beta_prior_cov = 1,
      restrict_dag = FALSE
    )
  )

  # Run parametric bootstrap on nesmr object
  pbstrap <- total_to_direct_parameteric_bootstrap(
    nesmr_mod,
    bootstrap_samples = 100
  )

  # Check output structure
  expect_true(is.list(pbstrap))
  expect_true(!is.null(pbstrap$direct_effects_se))
  expect_true(!is.null(pbstrap$direct_effects_bootstrap))
  expect_equal(nrow(pbstrap$direct_effects_se), d)
  expect_equal(ncol(pbstrap$direct_effects_se), d)
})

test_that("Parametric bootstrap method works with numeric inputs", {
  # Create simple test data
  beta <- c(0.1, 0.2, -0.15)
  beta_cov <- matrix(c(0.01, 0.002, -0.001,
                       0.002, 0.015, 0.003,
                       -0.001, 0.003, 0.012), nrow = 3)
  beta_ix <- cbind(c(1, 2, 3), c(2, 3, 1))
  d <- 3

  # Run parametric bootstrap with numeric inputs
  pbstrap <- total_to_direct_parameteric_bootstrap(
    x = beta,
    beta_cov = beta_cov,
    beta_ix = beta_ix,
    d = d,
    bootstrap_samples = 100
  )

  # Check output structure
  expect_true(is.list(pbstrap))
  expect_true(!is.null(pbstrap$direct_effects_se))
  expect_true(!is.null(pbstrap$direct_effects_bootstrap))
  expect_equal(nrow(pbstrap$direct_effects_se), d)
  expect_equal(ncol(pbstrap$direct_effects_se), d)
})
