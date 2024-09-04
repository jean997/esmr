threshold_beta <- function(W, threshold = 0.01) {
  W[abs(W) < threshold_max] <- 0
  W
}


# Thresholds the DAG at the smallest value that gives a DAG
threshold_DAG <- function(W, return_threshold = TRUE) {
  # Max dag zeros is d
  d <- ncol(W)
  stopifnot(nrow(W) == d)
  non_zero_vals <- sort(abs(matrix_to_edgelist(W, remove_diag = TRUE)$value))
  max_i <- d * (d - 1)
  start_i <- max_i / 2
  is_W_dag <- is_dag(W)
  while (! is_dag(W) && start_i <= max_i) {
    threshold <- non_zero_vals[start_i]
    W <- threshold_beta(W, threshold = threshold)
    cat("start_i", start_i, "threshold", threshold, "\n")
    start_i <- start_i + 1
  }
  if (start_i > max_i) {
    warning("Could not find a threshold that gives a DAG")
    return(NULL)
  }
  if (return_threshold) {
    return(list(W = W, threshold = threshold))
  } else {
    return(W)
  }
}
