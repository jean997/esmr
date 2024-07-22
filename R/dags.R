possible_dags <- function(n) {
  a <- vector(mode = 'numeric', length = n)
  a[1] <- 1
  a[2] <- 1
  for (k in 3:n) {
    (-1)^(k + 1) * choose(n, k) * 2^(k * (n - k)) * a[n - k]
  }
  return(a)
}

is_acyclic <- function(g) {
  tryCatch(
    !is.null(topo_sort(g)),
    error = function(e) FALSE)
}

# Checks if a matrix is a DAG
# Uses Proposition 1 from arXiv:1803.01422v2
is_dag <-function(x) {
  d <- ncol(x)
  stopifnot(nrow(x) == d)
  all(max(abs(eigen(x, only.values = TRUE)$values)) < 1) &&
    sum(diag(solve(diag(d) - x))) == d
}
