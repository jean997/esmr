make_mh_fixture <- function() {
  B1 <- matrix(0, 4, 4); B1[1, 2] <- .5; B1[2, 3] <- -.4; B1[1, 4] <- .3
  B2 <- B1; B2[3, 4] <- .2
  gs <- list(B1, B2)
  elbos <- c(-100, -90)
  vg <- lapply(1:2, function(i) list(
    elbo = elbos[i], elbo_without_prior = elbos[i], beta_hat = gs[[i]],
    se_beta_hat = gs[[i]] * 0 + .1, visited_count = 1))
  names(vg) <- vapply(gs, function(B) paste((B != 0) * 1, collapse = ""), "")
  ci <- data.frame(
    chain = 1, iter = 0:3,
    curr_graph = c(NA, names(vg)[1], names(vg)[2], names(vg)[2]),
    prop_graph = names(vg)[c(1, 2, 2, 1)],
    accepted = c(NA, 1, 0, 1),
    move_type = c("init", "add", NA, "remove"),
    curr_elbo = c(NA, -100, -90, -90),
    prop_elbo = c(-100, -90, -90, -100),
    mod_edge = c(NA, "3|4", "3|4", "3|4"),
    insert_edge = c(NA, TRUE, TRUE, FALSE)
  )
  structure(list(visited_graphs = vg, mh_chain_info = ci), class = "nesmr_mh_graph_explore")
}

test_that("tile plots work for every type and with labels", {
  x <- make_mh_fixture()
  for (type in c("edge_inc_prob", "best", "total", "indirect", "start", "proposal_freq")) {
    p <- plot_nesmr_tile(x, type = type, labels = LETTERS[1:4])
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  p <- plot_nesmr_tile(top_i_graph(x), labels = LETTERS[1:4], node_order = 4:1)
  expect_equal(levels(p$data$from), c("D", "C", "B", "A"))
})

test_that("MH diagnostic plots build", {
  x <- make_mh_fixture()
  plots <- c(
    lapply(c("elbo", "exp_elbo", "cum_best_elbo", "cum_total_elbo", "new_best"),
           function(t) plot_mh_trace(x, t)),
    lapply(c("edit_distance", "elbo_diff", "frob_norm"),
           function(y) plot_mh_best_path(x, y, label_gain = NULL)),
    list(plot_mh_accept_rate(x, window = 2))
  )
  for (p in plots) expect_no_error(ggplot2::ggplot_build(p))
})

test_that("graph differences accept labels", {
  x <- make_mh_fixture()
  p <- plot_graph_differences(top_i_graph(x), discovery_start_graph(x), labels = LETTERS[1:4])
  expect_no_error(ggplot2::ggplot_build(p))
})
