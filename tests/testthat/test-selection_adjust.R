test_that("selection_log_p1 exact matches Monte Carlo for two traits", {
  skip_if_not_installed("mvtnorm")
  dat <- list(p = 2, k = 2, S = matrix(c(0.01, 0.003), 200, 2, byrow = TRUE), s_equal = TRUE,
              omega = diag(1 / c(0.01, 0.003)^2))
  fgbar <- matrix(c(1, 0.3, 0, 1), 2, 2)
  sel <- selection_setup(dat, 5e-8)
  expect_true(sel$exact)
  pis <- c(0.1, 0.2); sg <- c(0.03, 0.02)
  lp_exact <- selection_log_p1(sel, fgbar, pis, sg)
  sel_mc <- sel
  sel_mc$exact <- FALSE
  set.seed(2)
  M <- 4e5
  sel_mc$u <- matrix(runif(M * 2), M, 2); sel_mc$xi <- matrix(rnorm(M * 2), M, 2)
  sel_mc$noise <- matrix(rnorm(M * 2), M, 2); sel_mc$pi_ref <- c(0.15, 0.15); sel_mc$mc_smooth <- 0.02
  lp_mc <- selection_log_p1(sel_mc, fgbar, pis, sg)
  expect_lt(lp_exact, 0)
  expect_equal(exp(lp_mc), exp(lp_exact), tolerance = 0.05)
  # more spread-out effects -> larger selection probability
  expect_gt(selection_log_p1(sel, fgbar, pis, sg * 2), lp_exact)
})

test_that("nesmr with selection_adjust runs and returns a larger objective than the unadjusted one", {
  skip_if_not_installed("GWASBrewer")
  set.seed(3)
  G <- matrix(0, 2, 2); G[1, 2] <- 0.4
  dat <- GWASBrewer::sim_mv(G = G, N = c(20000, 60000), J = 4000, h2 = 0.2, pi = 0.05, sporadic_pleiotropy = TRUE, est_s = TRUE)
  ix <- which(apply(2 * pnorm(-abs(dat$beta_hat / dat$s_estimate)), 1, min) < 5e-8)
  B <- matrix(0, 2, 2); B[1, 2] <- 1
  f0 <- suppressWarnings(nesmr(dat$beta_hat[ix, ], dat$s_estimate[ix, ], B, max_iter = 100,
                              params = list(beta_prior_cov = 1, strict_mode = FALSE)))
  f1 <- suppressWarnings(nesmr(dat$beta_hat[ix, ], dat$s_estimate[ix, ], B, max_iter = 100,
                              params = list(beta_prior_cov = 1, strict_mode = FALSE, selection_adjust = TRUE, selection_pval = 5e-8)))
  expect_true(is.finite(f1$elbo))
  expect_gt(f1$elbo, f0$elbo) # -n log P(S = 1) > 0
  expect_error(nesmr(dat$beta_hat[ix, ], dat$s_estimate[ix, ], B, max_iter = 5,
                     params = list(selection_adjust = TRUE, strict_mode = FALSE)), "selection")
})
