library(testthat)
library(igraph)
library(GWASBrewer)

test_that("draw_graph adds edge correctly", {
  # Create a simple graph
  g <- graph_from_edgelist(matrix(c(1, 2), ncol = 2), directed = TRUE)
  # Ensure that all vertices exist
  g <- add_vertices(g, nv = 1, name = "3")

  # Create input list with add candidates only
  x <- list(
    add_candidates = matrix(c(2, 3), ncol = 2),
    add_edge_prob = c(0.8),
    remove_candidates = matrix(c(1, 2), ncol = 2),
    remove_edge_prob = c(0.2)
  )

  set.seed(123)
  result <- draw_graph(g, x)

  expect_is(result, "list")
  expect_named(result, c("g", "prob", "mod_edge", "insert_edge", "from", "to"))
  expect_is(result$g, "igraph")
  expect_is(result$prob, "numeric")
  expect_is(result$mod_edge, "character")
  expect_is(result$insert_edge, "logical")
})

test_that("draw_graph removes edge correctly", {
  # Create a graph with edges
  g <- graph_from_edgelist(matrix(c(1, 2, 2, 3), ncol = 2, byrow = TRUE), directed = TRUE)

  # Create input list forcing remove
  x <- list(
    add_candidates = matrix(numeric(0), ncol = 2),
    add_edge_prob = numeric(0),
    remove_candidates = matrix(c(1, 2), ncol = 2),
    remove_edge_prob = c(1.0)
  )

  set.seed(42)
  result <- draw_graph(g, x)

  expect_is(result, "list")
  expect_false(result$insert_edge)
  expect_equal(result$mod_edge, "1|2")
})

test_that("draw_graph returns correct from and to values", {
  g <- graph_from_edgelist(matrix(c(1, 2), ncol = 2), directed = TRUE)
  # Ensure that all vertices exist
  g <- add_vertices(g, nv = 2, name = c("3", "4"))

  x <- list(
    add_candidates = matrix(c(3, 4), ncol = 2),
    add_edge_prob = c(1.0),
    remove_candidates = matrix(numeric(0), ncol = 2),
    remove_edge_prob = numeric(0)
  )

  set.seed(99)
  result <- draw_graph(g, x)

  expect_equal(result$from, "3")
  expect_equal(result$to, "4")
})

test_that("draw_graph prob sums to 1", {
  g <- graph_from_edgelist(matrix(c(1, 2, 2, 3), ncol = 2, byrow = TRUE), directed = TRUE)
  # Ensure that all vertices exist
  g <- add_vertices(g, nv = 1, name = "3")

  x <- list(
    add_candidates = matrix(c(3, 4), ncol = 2),
    add_edge_prob = c(0.6),
    remove_candidates = matrix(c(1, 2, 2, 3), ncol = 2, byrow = TRUE),
    remove_edge_prob = c(0.2, 0.2)
  )

  set.seed(55)
  result <- draw_graph(g, x)

  expect_is(result$prob, "numeric")
  expect_gt(result$prob, 0)
  expect_lte(result$prob, 1)
})

test_that("get_adjacent_graphs filters out cycles", {
  g <- graph_from_edgelist(matrix(c(1, 2, 2, 3), ncol = 2, byrow = TRUE), directed = TRUE)
  weight_mat <- matrix(seq(0.1, 0.9, length.out = 9), nrow = 3)

  result <- get_adjacent_graphs(g, weight_mat)

  # Should not allow adding edge 3->1 (would create cycle)
  if (nrow(result$add_candidates) > 0) {
    edge_3_to_1 <- apply(result$add_candidates, 1, function(x) x[1] == 3 && x[2] == 1)
    expect_false(any(edge_3_to_1))
  }
})

test_that("get_adjacent_graphs validates weight matrix dimensions", {
  g <- graph_from_edgelist(matrix(c(1, 2), ncol = 2), directed = TRUE)
  # Create a weight matrix with wrong dimensions (3x3 instead of 2x2)
  weight_mat <- matrix(seq(0.1, 0.9, length.out = 9), nrow = 3)

  expect_error(get_adjacent_graphs(g, weight_mat),
               "weight_mat must be a square matrix with dimensions matching the number of vertices in the graph")
})

test_that("get_adjacent_graphs probabilities sum to 1", {
  g <- graph_from_edgelist(matrix(c(1, 2), ncol = 2), directed = TRUE)

  weight_mat <- matrix(c(0, 0.5, 0.3, 0), nrow = 2)

  result <- get_adjacent_graphs(g, weight_mat)

  total_prob <- sum(result$add_edge_prob) + sum(result$remove_edge_prob)
  expect_equal(total_prob, 1, tolerance = 1e-10)
})

test_that("get_adjacent_graphs respects logistic parameters", {
  g <- graph_from_edgelist(matrix(c(1, 2), ncol = 2), directed = TRUE)
  weight_mat <- matrix(c(0, 0.5, 0.3, 0), nrow = 2)

  result <- get_adjacent_graphs(g, weight_mat, logistic_scale = 2, logistic_location = 3)

  expect_equal(result$logistic_scale, 2)
  expect_equal(result$logistic_location, 3)
})

test_that("get_adjacent_graphs handles empty add candidates", {
  # Create a complete DAG - no valid edges to add
  g_complete <- graph_from_edgelist(matrix(c(1, 2, 1, 3, 2, 3), ncol = 2, byrow = TRUE), directed = TRUE)
  weight_mat <- matrix(seq(0.1, 0.9, length.out = 9), nrow = 3)

  result <- get_adjacent_graphs(g_complete, weight_mat)

  expect_equal(nrow(result$add_candidates), 0)
  expect_equal(length(result$add_edge_prob), 0)
})

test_that("get_adjacent_graphs remove candidates equal graph edges", {
  g <- graph_from_edgelist(matrix(c(1, 2, 2, 3, 1, 3), ncol = 2, byrow = TRUE), directed = TRUE)
  weight_mat <- matrix(seq(0.1, 0.9, length.out = 9), nrow = 3)

  result <- get_adjacent_graphs(g, weight_mat)

  expect_equal(nrow(result$remove_candidates), 3)
  expect_equal(length(result$remove_edge_prob), 3)
})

test_that("get_adjacent_graphs gives probability 1 to a single edge", {
  g <- graph_from_edgelist(matrix(c(1, 2), ncol = 2), directed = TRUE)
  # Ensure that all vertices exist
  #g <- add_vertices(g, nv = 1, name = "3")

  weight_mat <- matrix(c(0, 0, 1, 0), nrow = 2)

  result <- get_adjacent_graphs(g, weight_mat, logistic_scale = 1, logistic_location = 5)

  expect_equal(result$remove_edge_prob[1], 1, tolerance = 1e-6)
})

test_that("get_adjacent_graphs all probabilities are non-negative", {
  g <- graph_from_edgelist(matrix(c(1, 2), ncol = 2), directed = TRUE)
  # Ensure that all vertices exist
  g <- add_vertices(g, nv = 1, name = "3")
  weight_mat <- matrix(c(0, 0.5, 0.3, 0.2, 0.4, 0.6, 0.1, 0.7, 0), nrow = 3)

  result <- get_adjacent_graphs(g, weight_mat)

  expect_true(all(result$add_edge_prob >= 0))
  expect_true(all(result$remove_edge_prob >= 0))
})

### MH Discovery Test
test_that("mh_graph_explore function works", {
  set.seed(13)
  G <- matrix(
    c(0, 0, 0, 0,
      0.4, 0, 0, 0,
      0, 0, 0,  0,
      0, -0.5, 0.2, 0),
    nrow = 4,
    byrow = 4
  )

  d <- ncol(G)

  B_full <- matrix(1, d, d) - diag(d)
  B_lower <- lower.tri(G) + 0
  B_correct <- (G != 0) + 0

  h2 <- 0.3
  J <- 5000
  N <- 40000
  pi_J <- 0.1
  alpha <- 5e-8

  dat <- sim_mv(
      G = G,
      N = N,
      J = J,
      h2 = h2,
      pi = pi_J,
      sporadic_pleiotropy = TRUE,
      est_s = TRUE
  )

  Z <- dat$beta_hat / dat$s_estimate
  pval_select <- 2*pnorm(-abs(Z))
  minp <- apply(pval_select, 1, min)
  ix <- which(minp < alpha)

  complete_mvmr <- nesmr_complete_mvmr(
      beta_hat = dat$beta_hat,
      se_beta_hat = dat$s_estimate,
      alpha = alpha
  )

  nesmr_discovery <- mh_graph_explore(
        dat, pval_select = pval_select, verbose = FALSE, debug = FALSE,
        n_mvmr_res = complete_mvmr, max_nesmr_fits = 20)

  best_graph <- top_i_graph(nesmr_discovery) %>%
    as.matrix(value = "direct_effect") %>% unname()

  expect_equal(best_graph, G, tolerance = 0.1)
  expect_equal((best_graph != 0) + 0, B_correct)
})
