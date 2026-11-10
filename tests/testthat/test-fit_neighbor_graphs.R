random_dag <- function(d, p = 0.35, seed = 1) {
  set.seed(seed)
  B <- matrix(0, d, d)
  perm <- sample(d)
  for (i in seq_len(d)) for (j in seq_len(d)) if (i < j && runif(1) < p) B[perm[i], perm[j]] <- 1
  B
}

test_that("graph_neighbors(distance = 1) matches the MH proposal neighbourhood", {
  B <- random_dag(7)
  a <- get_adjacent_graphs(igraph::graph_from_adjacency_matrix(B), matrix(1, 7, 7))
  key <- function(M) paste0(M, collapse = "")
  expected <- c(
    apply(a$add_candidates, 1, function(e) { M <- B; M[e[1], e[2]] <- 1; key(M) }),
    apply(a$remove_candidates, 1, function(e) { M <- B; M[e[1], e[2]] <- 0; key(M) }),
    apply(a$flip_candidates, 1, function(e) { M <- B; M[e[1], e[2]] <- 0; M[e[2], e[1]] <- 1; key(M) })
  )
  nb <- graph_neighbors(B, 1)
  expect_setequal(nb$graph, expected)
  expect_equal(anyDuplicated(nb$graph), 0)
  expect_true(all(nb$distance == 1))
})

test_that("graph_neighbors(distance = 2) returns distinct DAGs at minimal distance", {
  B <- random_dag(6, seed = 2)
  nb <- graph_neighbors(B, 2)
  expect_equal(anyDuplicated(nb$graph), 0)
  expect_false(paste0(B, collapse = "") %in% nb$graph)
  expect_true(all(vapply(nb$graph, function(s) is_dag(flat_string_to_adj_mat(s)), logical(1))))
  # Every distance-1 graph keeps distance 1 even though it is also reachable in 2 moves
  expect_setequal(nb$graph[nb$distance == 1], graph_neighbors(B, 1)$graph)
  expect_true(all(nb$distance %in% 1:2))
})

test_that("moves argument restricts the neighbourhood", {
  B <- random_dag(6, seed = 3)
  expect_true(all(graph_neighbors(B, 1, moves = "remove")$n_edges == sum(B) - 1))
  expect_true(all(graph_neighbors(B, 1, moves = "add")$n_edges == sum(B) + 1))
  expect_true(all(graph_neighbors(B, 1, moves = "flip")$n_edges == sum(B)))
})

test_that("fit_neighbor_graphs dry_run and max_fits do not fit anything", {
  d <- 5
  B <- random_dag(d, seed = 4)
  key <- paste0(B, collapse = "")
  nb1 <- graph_neighbors(B, 1)
  explore <- structure(
    list(visited_graphs = setNames(
      list(list(elbo = 10), list(elbo = 5)), c(key, nb1$graph[1])
    )),
    class = "nesmr_mh_graph_explore"
  )
  dat <- list(beta_hat = matrix(0, 10, d), s_estimate = matrix(1, 10, d))
  dr <- fit_neighbor_graphs(explore, dat, edit_distance = 1, dry_run = TRUE)
  expect_s3_class(dr, "nesmr_neighbor_fits")
  expect_false(nb1$graph[1] %in% dr$graphs$graph) # already visited is skipped
  expect_equal(nrow(dr$graphs), nrow(nb1[nb1$n_edges > 0, ]) - 1)
  expect_false(any(dr$graphs$fitted))
  expect_error(
    fit_neighbor_graphs(explore, dat, edit_distance = 1, max_fits = 1, check_consistency = FALSE),
    "max_fits"
  )
})
