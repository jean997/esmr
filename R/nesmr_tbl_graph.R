#' @importFrom tidygraph as_tbl_graph tbl_graph
#' @export
as.matrix.nesmr_tbl_graph <- function(x, value = c("direct_effect", "total_effect"), ...) {
  graph_names <- activate(x, nodes) |> pull(name)

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