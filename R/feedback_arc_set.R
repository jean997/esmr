maximal_acyclic_subgraph <- function(W) {
  g <- igraph::graph_from_adjacency_matrix(W, diag = FALSE, weighted = TRUE)
  adj_mat <- igraph::as_adjacency_matrix(g - igraph::feedback_arc_set(g), sparse = FALSE)
  return(W * adj_mat)
}
