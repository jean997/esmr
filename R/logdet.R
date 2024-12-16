#' LogDeterminant acyclicity functions
#'
#' h_det(X) and h_det_grad(X) are zero exactly when X is an adjacency matrix that is a DAG
#'
#' From: Bello, K., Aragam, B., & Ravikumar, P. (2023)
#' DAGMA: Learning DAGs via M-matrices and a Log-Determinant Acyclicity Characterization (arXiv:2209.08037).
#' arXiv. https://doi.org/10.48550/arXiv.2209.08037
#' @export
#' @rdname logdet
#' @param X adjacency matrix
#' @param x vector representation of the adjacency matrix
#' @param s h_det parameter. Should be greater than the largest eigenvalues of X^2
#' @param d number of rows/columns in the matrix for vector versions of functions
#' @param ix index of the non-zero elements of the matrix in the vector versions of functions
h_det <- function(X, s = 1) {
  d <- ncol(X)
  -log(det(s * diag(d) - X * X)) + d * log(s)
}

#' @export
#' @rdname logdet
h_det_grad <- function(X, s = 1) {
  d <- ncol(X)
  2 * solve(t(s * diag(d) - X * X)) * X
}

#' @export
#' @rdname logdet
h_det_vec <- function(x, s = 1, d, ix = NULL) {
  if (is.null(ix)) {
    X <- matrix(x, ncol = d)
    h_det(X, s)
  } else {
    X <- matrix(0, ncol = d, nrow = d)
    X[ix] <- x
    h_det(X, s)
  }
}

#' @export
#' @rdname logdet
h_det_grad_vec <- function(x, s = 1, d, ix = NULL) {
  if (is.null(ix)) {
    X <- matrix(x, ncol = d)
    grad <- h_det_grad(X, s)
    as.vector(grad)
  } else {
    X <- matrix(0, ncol = d, nrow = d)
    X[ix] <- x
    grad <- h_det_grad(X, s)
    as.vector(grad[ix])
  }
}

# Sparse matrix version
# From https://en.wikipedia.org/wiki/Commutation_matrix#R
comm_mat <- function(m, n){
  i = seq_len(m * n)
  j = NULL
  for (k in 1:m) {
    j = c(j, m * 0:(n-1) + k)
  }
  Matrix::sparseMatrix(
    i = i, j = j, x = 1
  )
}

#' @export
#' @rdname logdet
h_det_hessian <- function(X, s = 1) {
  d <- nrow(X)
  N <- solve(diag(s, nrow = d) - X * X)
  W_vec <- as.numeric(X)
  W_t_vec <- as.numeric(t(X))
  N_t_vec <- as.numeric(t(X))
  K_comm <- comm_mat(d, d)
  rtn <- 4 * diag(W_vec) %*% kronecker(X, t(X))
  rtn <- rtn %*% diag(W_t_vec) %*% K_comm
  rtn <- rtn + 2 * diag(N_t_vec)
  rtn
}

#' @export
#' @param W_vec vector representation of a matrix (as.vector or as.numeric(W))
#' @param d number of rows/columns in the matrix
#' @param s h_det parameter. Should be greater than the highest eigenvalues of W^2
#' @rdname logdet
h_det_hessian_vec <- function(x, d = sqrt(length(x)), s = 1) {
  X <- matrix(x, nrow = d)
  h_det_hessian(X, s)
}

#' @export
spectral_radius <- function(A) {
  max(abs(eigen(A)$values))
}
