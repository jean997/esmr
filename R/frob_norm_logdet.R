frob_logdet_grad_vec <- function(
    x, W, d, s = 1, lambda = 1, ix = NULL,
    penalty = c("L2", "L1")) {
  penalty <- match.arg(penalty)
  if (penalty == "L2") {
    penalty_grad <- function(.X) { 2 * (.X - W) }
  } else if (penalty == "L1") {
    penalty_grad <- function(.X) { sign(.X - W) }
  }

  if (is.null(ix)) {
    X <- matrix(x, ncol = d, nrow = d)
    grad <- h_det_grad(X, s) + penalty_grad(X) * lambda
    as.vector(grad)
  } else {
    X <- matrix(0, ncol = d, nrow = d)
    X[ix] <- x
    grad <- h_det_grad(X, s) + penalty_grad(X) * lambda
    as.vector(grad[ix])
  }
}

frob_logdet_vec <- function(
    x, W, d, s = 1, lambda = 1, ix = NULL, penalty = c("L1", "L2")) {
  penalty <- match.arg(penalty)
  if (is.null(ix)) {
    X <- matrix(x, ncol = d)
  } else {
    X <- matrix(0, ncol = d, nrow = d)
    X[ix] <- x
  }
  score <- h_det(X, s)
  if (penalty == "L1") {
    score <- score + sum(abs(X - W)) * lambda
  } else if (penalty == "L2") {
    score <- score + sum((X - W)^2) * lambda
  }
  return(score)
}
