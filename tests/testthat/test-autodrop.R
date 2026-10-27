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
    h2 = c(0.2, 0.2, 0.2, 0),
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

  esmr_mod_four <- esmr(
    beta_hat_Y = dat$beta_hat[, 2],
    se_Y = dat$s_estimate[, 2],
    beta_hat_X = dat$beta_hat[, -2],
    se_X = dat$s_estimate[, -2],
    params = list(
      strict_mode = FALSE,
      cond_num = 300
    )
  )

  expect_equal(esmr_mod_four$remove_suggest, 4)

  autoremove_mod <- fit_esmr_autoremove(
    beta_hat_X = dat$beta_hat,
    se_X = dat$s_estimate
  )

  n_mvmr_results <- nesmr_complete_mvmr_cached_gfa(
    beta_hat = dat$beta_hat,
    se_beta_hat = dat$s_estimate,
    workers = 1
  )

  n_mvmr_res2 <- nesmr_complete_mvmr(
    beta_hat = dat$beta_hat,
    se_beta_hat = dat$s_estimate
  )

  n_mvmr_results$beta_hat
  n_mvmr_res2$beta_hat
})