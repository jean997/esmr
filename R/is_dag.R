# Checks if a matrix is a DAG
# Uses Proposition 1 from arXiv:1803.01422v2
is_dag <-function(x) {
  d <- ncol(x)
  stopifnot(nrow(x) == d)
  all(max(abs(eigen(x, only.values = TRUE)$values)) < 1) &&
    sum(diag(solve(diag(d) - x))) == d
}
