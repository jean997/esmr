topo_sort_mat <- function(x) {
  ig <- igraph::graph_from_adjacency_matrix(
      t(x) != 0,
      mode = 'directed'
    )
  as.numeric(igraph::topo_sort(ig))
}
