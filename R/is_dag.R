# Checks if a matrix is a DAG
# Uses Proposition 1 from arXiv:1803.01422v2
is_dag_matrix <- function(x) {
  d <- ncol(x)
  stopifnot(nrow(x) == d)
  tryCatch({
    all(max(abs(eigen(x, only.values = TRUE)$values)) < 1) &&
      sum(diag(solve(diag(d) - x))) == d
  }, error = function(e) {
    FALSE
  })
}

# Faster to topological sort than to check for cycles
is_dag <- function(x) {
  tryCatch({
    !is.null(topo_sort_mat(x))
  }, error = function(e) {
    FALSE
  })
}

