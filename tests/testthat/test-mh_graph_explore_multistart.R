key <- function(M) paste0((M != 0) + 0, collapse = "")

test_that(".sample_start_graphs returns unique acyclic graphs and respects max_tries", {
  set.seed(1)
  d <- 6
  z <- matrix(rnorm(d * d, sd = 3), d, d)
  diag(z) <- 0
  thr <- c(0.5, 5e-3)
  s <- .sample_start_graphs(z, thr, n_graphs = 15, max_tries = 200)
  expect_equal(anyDuplicated(s$starts$graph), 0)
  expect_true(all(vapply(s$starts$graph, function(g) is_dag(flat_string_to_adj_mat(g)), NA)))
  expect_true(all(s$starts$n_edges > 0))
  expect_true(all(s$tries <= 200))
  expect_true(all(c("mvmr", "draw") %in% s$starts$source))
  # at most n_graphs random draws per threshold (+ the deterministic graph)
  expect_true(all(table(s$starts$threshold) <= 16))

  # a tiny max_tries gives a warning instead of looping
  expect_warning(.sample_start_graphs(z, 5e-7, n_graphs = 50, max_tries = 3))
})

test_that("edit distance counts a flip as one edit", {
  A <- matrix(0, 3, 3); A[1, 2] <- 1
  B <- matrix(0, 3, 3); B[2, 1] <- 1 # flip of A
  C <- A; C[2, 3] <- 1 # one added edge
  D <- .edit_distance_matrix(c(key(A), key(B), key(C)))
  expect_equal(D[1, 2], 1)
  expect_equal(D[1, 3], 1)
  expect_equal(D[2, 3], 2)
  expect_equal(diag(D), c(0, 0, 0))
})

test_that(".select_chain_starts keeps the best graph of each group within merge_distance", {
  # A-B, A-C and B-E, C-E are one edit apart (one connected group); F is far below the window
  A <- matrix(0, 4, 4); A[1, 2] <- 1; A[2, 3] <- 1
  B <- A; B[3, 4] <- 1
  C <- A; C[1, 4] <- 1
  E <- A; E[3, 4] <- 1; E[1, 4] <- 1
  F <- matrix(0, 4, 4); F[4, 1] <- 1
  starts <- data.frame(
    graph = vapply(list(A, B, C, E, F), key, ""),
    elbo = c(-10, -5, -8, -12, -100), stringsAsFactors = FALSE
  )
  # B-C are 2 edits apart but both within 1 of A -> same connected group
  out <- .select_chain_starts(starts, window = 20, merge_distance = 1)
  expect_equal(out$in_window, c(TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_equal(out$kept, c(FALSE, TRUE, FALSE, FALSE, FALSE)) # E joins the group through B / C
  expect_equal(length(unique(out$group[out$in_window])), 1)

  # merge_distance = 0 keeps every start inside the window
  out0 <- .select_chain_starts(starts, window = 20, merge_distance = 0)
  expect_equal(out0$kept, c(TRUE, TRUE, TRUE, TRUE, FALSE))
})

test_that(".select_chain_starts: three graphs pairwise within one edit leave only the best", {
  A <- matrix(0, 4, 4); A[1, 2] <- 1
  B <- A; B[2, 3] <- 1
  C <- A; C[3, 4] <- 1
  starts <- data.frame(graph = vapply(list(A, B, C), key, ""), elbo = c(-3, -1, -2), stringsAsFactors = FALSE)
  out <- .select_chain_starts(starts, window = 20, merge_distance = 1)
  # B-C differ by two edits, but both are one edit from A so the group is connected
  expect_equal(sum(out$kept), 1)
  expect_true(out$kept[2])
})

test_that(".seen_graphs collects proposed and current graphs", {
  info <- data.frame(
    curr_graph = c(NA, "a", "b"), prop_graph = c("a", "b", "c"), stringsAsFactors = FALSE
  )
  expect_setequal(.seen_graphs(info), c("a", "b", "c"))
  expect_equal(.seen_graphs(NULL), character(0))
})

# ---- end to end on a small simulation ----------------------------------------

sim_dat <- function() {
  G <- matrix(c(0, 0, 0, 0, 0.4, 0, 0, 0, 0, 0, 0, 0, 0, -0.5, 0.2, 0), nrow = 4, byrow = TRUE)
  set.seed(11)
  GWASBrewer::sim_mv(G = G, N = 20000, J = 2000, h2 = 0.2, pi = 0.1, sporadic_pleiotropy = TRUE, est_s = TRUE)
}

test_that("kill_no_improve_models stops a chain and restarts the counter for each chain", {
  skip_if_not_installed("GWASBrewer")
  dat <- sim_dat()
  set.seed(2)
  res <- mh_graph_explore(
    dat, max_iter = 300, max_nesmr_fits = 300, kill_no_improve_models = 5,
    init_prob_threshold = 0.05, mh_chain_init = list(), random_starts = 2
  )
  cs <- res$chain_status
  expect_equal(nrow(cs), 3) # best_approx + 2 random starts
  # Every chain that stopped early needed exactly >= 5 new models since its last improvement
  early <- cs[cs$stop_reason == "early_stop", ]
  expect_true(all(early$models_since_improve >= 5))
  # A later chain still gets to fit models rather than stopping right after the first iteration
  expect_true(all(cs$n_nesmr_fits >= 1))
  info <- res$mh_chain_info
  iters <- tapply(info$iter, info$chain, max)
  expect_true(all(iters >= 1))
})

test_that("mh_graph_explore_multistart runs chains in ELBO order into one object", {
  skip_if_not_installed("GWASBrewer")
  dat <- sim_dat()
  set.seed(3)
  # the 4-node problem has few distinct graphs at strict thresholds, which warns
  res <- suppressWarnings(mh_graph_explore_multistart(
    dat,
    thresholds = c(0.5, 5e-3), n_graphs = 6, max_tries = 100,
    start_elbo_window = 1e6, merge_distance = 1, early_stop_iter = 4,
    max_iter = 60, max_nesmr_fits = 60, max_chains = 3
  ))
  expect_s3_class(res, "nesmr_mh_graph_explore")
  ms <- res$multistart
  expect_equal(anyDuplicated(ms$starts$graph), 0)
  # chains are ordered by decreasing start ELBO
  expect_false(is.unsorted(rev(ms$chains$start_elbo)))
  ran <- ms$chains[ms$chains$status == "ran", ]
  expect_lte(nrow(ran), 3)
  expect_gte(nrow(ran), 1)
  # separate, named chains in a single object
  expect_setequal(unique(res$mh_chain_info$chain_name), ran$chain_name)
  expect_true(all(grepl("^start_", ran$chain_name)))
  expect_equal(nrow(res$chain_status), nrow(ran))
  # the start fits are in visited_graphs
  expect_true(all(ms$starts$graph %in% names(res$visited_graphs)))
  # starts within one edit of a better start are not used as chain starts
  kept <- ms$starts[ms$starts$kept, ]
  D <- .edit_distance_matrix(kept$graph)
  expect_true(all(D[upper.tri(D)] > 1))
  # a graph already seen by a chain is never started again
  seen_before <- character(0)
  for (r in seq_len(nrow(ran))) {
    expect_false(ran$graph[r] %in% seen_before)
    info <- res$mh_chain_info[res$mh_chain_info$chain_name == ran$chain_name[r], ]
    seen_before <- union(seen_before, .seen_graphs(info))
  }
})
