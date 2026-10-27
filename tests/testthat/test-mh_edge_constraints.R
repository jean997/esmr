test_that("edge_constraint_matrix accepts indices, names, lists and matrices", {
  nm <- c("a", "b", "c")
  m1 <- edge_constraint_matrix(rbind(c(1, 2), c(3, 1)), 3)
  expect_true(m1[1, 2] && m1[3, 1])
  expect_equal(sum(m1), 2)

  m2 <- edge_constraint_matrix(cbind("b", "a"), 3, nm)
  expect_true(m2[2, 1])
  expect_equal(sum(m2), 1)

  m3 <- edge_constraint_matrix(list(c("a", "c")), 3, nm)
  expect_true(m3[1, 3])

  adj <- matrix(0, 3, 3)
  adj[2, 3] <- 1
  expect_equal(edge_constraint_matrix(adj, 3), adj != 0)

  expect_equal(sum(edge_constraint_matrix(NULL, 3)), 0)
  expect_error(edge_constraint_matrix(cbind("x", "a"), 3, nm), "unknown trait")
  expect_error(edge_constraint_matrix(cbind(1, 4), 3), "out of range")
  expect_error(edge_constraint_matrix(cbind(2, 2), 3), "self edges")
})

test_that("enforce_edge_constraints respects forbidden and required edges and stays acyclic", {
  d <- 4
  forbidden <- edge_constraint_matrix(cbind(1, 2), d)
  required <- edge_constraint_matrix(cbind(2, 1), d)
  adj <- matrix(1, d, d)
  diag(adj) <- 0
  B <- enforce_edge_constraints(adj, forbidden, required, weight = matrix(runif(d^2), d, d))
  expect_true(is_dag(B))
  expect_equal(B[2, 1], 1)
  expect_equal(B[1, 2], 0)
  expect_false(violates_edge_constraints(B, forbidden, required))

  # Required edges forming a cycle are rejected
  cyc <- edge_constraint_matrix(rbind(c(1, 2), c(2, 1)), d)
  expect_error(check_edge_constraints(matrix(FALSE, d, d), cyc), "cycle")
  expect_error(check_edge_constraints(cyc, cyc), "both forbidden and required")
})

test_that("get_adjacent_graphs never proposes forbidden or required-edge changes", {
  d <- 4
  B <- matrix(0, d, d)
  B[1, 2] <- 1 # required
  B[2, 3] <- 1
  g <- igraph::graph_from_adjacency_matrix(B, mode = "directed")
  W <- matrix(3, d, d)

  forbidden <- edge_constraint_matrix(rbind(c(3, 4), c(3, 2)), d) # 3 -> 4 forbidden, flip of 2 -> 3 forbidden
  required <- edge_constraint_matrix(cbind(1, 2), d)
  adj <- get_adjacent_graphs(g, W, forbidden = forbidden, required = required)

  as_key <- function(m) if (nrow(m) == 0) character(0) else paste(m[, 1], m[, 2])
  expect_false("3 4" %in% as_key(adj$add_candidates))
  expect_true("4 3" %in% as_key(adj$add_candidates))
  expect_false("1 2" %in% as_key(adj$remove_candidates))
  expect_true("2 3" %in% as_key(adj$remove_candidates))
  expect_false("1 2" %in% as_key(adj$flip_candidates))
  expect_false("2 3" %in% as_key(adj$flip_candidates)) # would add 3 -> 2
  # lengths stay consistent with the candidate matrices
  expect_equal(length(adj$add_edge_prob), nrow(adj$add_candidates))
  expect_equal(length(adj$remove_edge_prob), nrow(adj$remove_candidates))
  expect_equal(length(adj$flip_edge_prob), nrow(adj$flip_candidates))
})

test_that("mh_graph_explore honors forbidden and required edges", {
  skip_if_not_installed("GWASBrewer")
  set.seed(11)
  # 4 traits: 1 -> 2 (strong), 3 -> 2, 2 -> 4
  G <- matrix(0, 4, 4)
  G[1, 2] <- 0.5
  G[3, 2] <- 0.4
  G[2, 4] <- -0.4
  dat <- GWASBrewer::sim_mv(G = G, N = 20000, J = 3000, h2 = 0.2, pi = 0.1,
                            sporadic_pleiotropy = TRUE, est_s = TRUE)
  d <- 4
  dat_in <- list(beta_hat = dat$beta_hat, s_estimate = dat$s_estimate)

  # Forbid the true direction 1 -> 2 and force the reverse 2 -> 1
  res <- mh_graph_explore(
    dat_in, alpha = 5e-8, max_iter = 40, max_nesmr_fits = 25,
    graph_draw_mode = "fixed_prop", init_prob_threshold = 0.05,
    forbidden_edges = cbind(1, 2), required_edges = cbind(2, 1)
  )
  graphs <- names(res$visited_graphs)
  mats <- lapply(graphs, flat_string_to_adj_mat)
  expect_true(all(vapply(mats, function(m) m[1, 2] == 0 && m[2, 1] == 1, logical(1))))
  expect_true(all(vapply(mats, is_dag, logical(1))))
  expect_true(res$forbidden_edges[1, 2])
  expect_true(res$required_edges[2, 1])
  expect_equal(res$best_graph$B_template[1, 2], 0)
  expect_equal(res$best_graph$B_template[2, 1], 1)

  # Direction constraint only: the edge is optional but can never point 1 -> 2
  res2 <- mh_graph_explore(
    dat_in, alpha = 5e-8, max_iter = 40, max_nesmr_fits = 25,
    graph_draw_mode = "fixed_prob", init_prob_threshold = 0.05,
    forbidden_edges = cbind(1, 2)
  )
  mats2 <- lapply(names(res2$visited_graphs), flat_string_to_adj_mat)
  expect_true(all(vapply(mats2, function(m) m[1, 2] == 0, logical(1))))

  # Continuing an unconstrained run with conflicting constraints is rejected
  res_free <- mh_graph_explore(
    dat_in, alpha = 5e-8, max_iter = 20, max_nesmr_fits = 10,
    graph_draw_mode = "fixed_prop", init_prob_threshold = 0.05
  )
  expect_error(
    mh_graph_explore(dat_in, alpha = 5e-8, max_iter = 5, max_nesmr_fits = 3, prev_run = res_free,
                     forbidden_edges = cbind(1, 2), required_edges = cbind(2, 1)),
    "violate"
  )
})
