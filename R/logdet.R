# From:
# Bello, K., Aragam, B., & Ravikumar, P. (2023).
# DAGMA: Learning DAGs via M-matrices and a Log-Determinant Acyclicity Characterization (arXiv:2209.08037).
# arXiv. https://doi.org/10.48550/arXiv.2209.08037
h_det <- function(X, s = 1) {
  d <- ncol(X)
  -log(det(s * diag(d) - X * X)) + d * log(s)
}

h_det_grad <- function(X, s = 1) {
  d <- ncol(X)
  2 * solve(t(s * diag(d) - X * X)) * X
}

# Wrapper function to convert matrix to vector and vice versa for optim
h_det_vector <- function(x, s = 1, d) {
  X <- matrix(x, ncol = d)
  h_det(X, s)
}

h_det_ell <- function(x, ell_func, d, s = 1, phi = 1) {
  X <- matrix(x, ncol = d)
  fg <- X + diag(d)

  ell_func(fg) + phi * h_det(X, s)
}

h_det_ell_grad <- function(x, ell_func, d, s = 1, phi = 1) {
  X <- matrix(x, ncol = d)
  fg <- X + diag(d)

  numDeriv::grad(ell_func, x = fg) + phi * h_det_grad(X, s)
}

h_det_grad_vector <- function(x, s = 1, d) {
  X <- matrix(x, ncol = d)
  grad <- h_det_grad(X, s)
  as.vector(grad)
}

spectral_radius <- function(A) {
  max(abs(eigen(A)$values))
}
