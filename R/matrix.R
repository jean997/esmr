#' Checks if matrix is positive semi-definite
#' @param X matrix to test
#' @param tol tolerance to allow for zero eigenvalues
is_psd <- function(X, tol = 1e-8) {
  eigen_vals <- eigen(X, only.values = TRUE)$values
  # First check is to ensure non-complex eigenvalues
  return(is.numeric(eigen_vals) &&
      all(eigen(X, only.values = TRUE)$values > tol))
}

#' Computes a faster version of a diagonal matrix times a square matrix
#' times the same diagonal matrix
#'
#' @param X square matrix with dimensions d x d
#' @param s vector of length d
solve_diag_psd_diag <- function(X, s, check_psd = FALSE) {
  if (!is.null(dim(s))) s <- diag(s)
  stopifnot(length(s) == ncol(X))
  # Same as diag(s) %*% X %*% diag(s)
  # But ~4x faster
  # s * X does column-wise multiplication: s * X[, 1], s * X[, 2], ...
  res <- t(t(s * X) * s)

  # Assume PSD if no check requested
  is_psd <- ! check_psd || is_psd(res, tol = 0)
  if (! is_psd) {
    solve(res)
  } else {
    # Assumes PSD
    chol2inv(chol(res))
  }
}
