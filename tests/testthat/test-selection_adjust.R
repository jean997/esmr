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

# mh_graph_explore with the selection-aware likelihood on a two trait network where the cause (trait 1) has
# a 5 - 10 fold smaller GWAS than the downstream trait (trait 2). With the plain ELBO the edge is oriented
# 2 -> 1 by power for these seeds; with selection_adjust = TRUE the true direction 1 -> 2 is the best graph.
sim_imbalanced_pair <- function(ratio, seed) {
  set.seed(seed)
  G <- matrix(0, 2, 2)
  G[1, 2] <- 0.4 # true edge: trait 1 (small GWAS) -> trait 2 (large GWAS)
  suppressMessages(GWASBrewer::sim_mv(
    G = G, N = c(20000, 20000 * ratio), J = 20000, h2 = 0.2, pi = 0.05,
    sporadic_pleiotropy = TRUE, est_s = TRUE
  ))
}

run_mh_pair <- function(dat, adjust) {
  suppressWarnings(mh_graph_explore(
    dat, alpha = 5e-8, max_iter = 40, max_nesmr_fits = 10,
    graph_draw_mode = "fixed_prop", init_prob_threshold = 0.05,
    selection_adjust = adjust
  ))
}

test_that("mh_graph_explore(selection_adjust = TRUE) orients the edge correctly with a 10 fold smaller cause GWAS", {
  skip_on_cran()
  skip_if_not_installed("GWASBrewer")
  dat <- sim_imbalanced_pair(ratio = 10, seed = 1)

  res_plain <- run_mh_pair(dat, adjust = FALSE)
  res_adj <- run_mh_pair(dat, adjust = TRUE)

  B_adj <- res_adj$best_graph$B_template
  expect_equal(B_adj[1, 2], 1) # true direction 1 -> 2
  expect_equal(B_adj[2, 1], 0)
  expect_true(all(vapply(names(res_adj$visited_graphs), function(g) is_dag(flat_string_to_adj_mat(g)), logical(1))))

  # The adjustment adds -n log P(S = 1 | Theta) > 0 to the ELBO of every graph, so the same graph scores higher
  common <- intersect(names(res_plain$visited_graphs), names(res_adj$visited_graphs))
  expect_gt(length(common), 0)
  for (g in common) {
    expect_gt(res_adj$visited_graphs[[g]]$elbo, res_plain$visited_graphs[[g]]$elbo)
  }
})

test_that("mh_graph_explore(selection_adjust = TRUE) orients the edge correctly with a 5 fold smaller cause GWAS", {
  skip_on_cran()
  skip_if_not_installed("GWASBrewer")
  dat <- sim_imbalanced_pair(ratio = 5, seed = 6)
  res_adj <- run_mh_pair(dat, adjust = TRUE)
  B_adj <- res_adj$best_graph$B_template
  expect_equal(B_adj[1, 2], 1)
  expect_equal(B_adj[2, 1], 0)
})
