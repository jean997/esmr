#' @importFrom tidygraph as_tbl_graph tbl_graph activate as_tibble
#' @export
as.matrix.nesmr_tbl_graph <- function(x, value, ...) {
  graph_names <- tidygraph::activate(x, nodes) |> pull(name)

  edgelist <- x |>
    tidygraph::activate(edges) |>
    tidygraph::as_tibble()

  edgelist <- if (ncol(edgelist) > 3) {
    edgelist %>% dplyr::select(from, to, !!value) %>% as.matrix()
  } else {
    edgelist %>% as.matrix()
  }

  rtn_mat <- matrix(0, nrow = length(graph_names), ncol = length(graph_names))
  rtn_mat[edgelist[, 1:2]] <- edgelist[,3]
  colnames(rtn_mat) <- rownames(rtn_mat) <- graph_names
  return(rtn_mat)
}

#' @export
as_matrix.nesmr_tbl_graph <- function(x, ...) {
  as.matrix.nesmr_tbl_graph(x, ...)
}

## Build a nesmr_tbl_graph from a matrix, with one edge per nonzero entry
matrix_to_nesmr_tbl_graph <- function(mat, weight = "direct_effect") {
  edges <- which(mat != 0, arr.ind = TRUE)
  edges <- data.frame(from = edges[, 1], to = edges[, 2], value = mat[edges])
  names(edges)[3] <- weight
  tg <- tidygraph::tbl_graph(
    nodes = data.frame(name = seq_len(ncol(mat))),
    edges = edges,
    directed = TRUE
  )
  class(tg) <- c("nesmr_tbl_graph", "esmr_tbl_graph", class(tg))
  tg
}

#' Starting graph of an MH graph exploration
#'
#' @param x A `nesmr_mh_graph_explore` object.
#' @return A `nesmr_tbl_graph` with the estimated direct effects of the graph the first
#'   chain started from.
#' @export
discovery_start_graph <- function(x) {
  info <- x$mh_chain_info
  start <- if (is.na(info$curr_graph[[1]])) info$prop_graph[[1]] else info$curr_graph[[1]]
  matrix_to_nesmr_tbl_graph(x$visited_graphs[[start]]$beta_hat)
}

#' Frequency of each edge among the proposed graphs of an MH graph exploration
#'
#' @param x A `nesmr_mh_graph_explore` object.
#' @return A `nesmr_tbl_graph` with edge column `proposal_freq`.
#' @export
edge_proposal_freq <- function(x) {
  n_proposed <- table(x$mh_chain_info$prop_graph)
  freq <- Reduce("+", Map(
    function(g, n) flat_string_to_adj_mat(g) * n,
    names(n_proposed), as.integer(n_proposed)
  )) / sum(n_proposed)
  matrix_to_nesmr_tbl_graph(freq, weight = "proposal_freq")
}
