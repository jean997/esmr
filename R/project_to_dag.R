#' @export
project_to_DAG <- function(
    X, s = 1,
    method = c("BFGS","NR")) {
  # TODO: Maybe something more efficient than this?
  method <- match.arg(method)
  d <- nrow(X)
  logdet_ix <- which(X != 0, arr.ind = TRUE)
  X_dag <- matrix(0, ncol = d, nrow = d)
  if (method == "BFGS") {
    result <- optim(
      par = X[logdet_ix],
      fn = h_det_vec,
      gr = h_det_grad_vec,
      s = s,
      ix = logdet_ix,
      d = d,
      method = "BFGS",
      control = list(maxit = 2000, trace = 0)
    )

    if (result$convergence != 0) {
      warning("Optimization did not converge")
      }

    X_dag[logdet_ix] <- result$par
    } else if (method == "NR") {
      result <- lava::NR(as.numeric(X),
               objective = h_det_vec,
               gradient = h_det_grad_vec,
               hessian = h_det_hessian_vec,
               args = list(
                 s = s,
                 d = d
               ),
               control = list(
                 backtrack = TRUE
               )
      )
      X_dag <- matrix(result$par, ncol = d, nrow = d)
    } else {
      stop("Invalid method")
    }
  return(X_dag)
}
