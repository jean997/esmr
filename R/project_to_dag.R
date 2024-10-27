#' Project a matrix to a DAG using the Frobenius norm and logdet penalty (DAGMA - Bello 2022)
#'
#' @param X a square matrix, not necessarily a DAG
#' @param s log-det penalty parameter
#' @param lambda vector for path finding parameter
#' @param penalty L2 or L1 penalty (default and more stable: L2)
#' @param threshold_to_DAG should the matrix be thresholded to a true DAG?
#' Otherwise, the matrix will have near zeros but not exact zeros.
#' Thresholding is done via \code{\link{threshold_DAG}} which iterativly finds the minimal threshold that ensures acyclicity.
#' @param ... additional arguments to the `control` argument of \code{\link{optim}}
#'
#' @return
#' @export
#'
#' @examples
project_to_DAG <- function(
    X, s = max(max(X^2) + 0.1, 1), lambda = c(1, 10^-seq(1,5), 0),
    penalty = c("L2", "L1"),
    threshold_to_DAG = FALSE,
    maxit = 100
    ) {
  penalty <- match.arg(penalty)
  d <- nrow(X)
  logdet_ix <- which(X != 0, arr.ind = TRUE)
  results <- list()
  # Initialize to 0 matrix rather than non-DAG starting matrix
  curr_pars <- rep(0, nrow(logdet_ix))
  for (i in seq_along(lambda)) {
    results[[i]] <- optim(
      par = curr_pars,
      fn = frob_logdet_vec,
      gr = frob_logdet_grad_vec,
      s = s,
      ix = logdet_ix,
      W = X,
      d = d,
      penalty = penalty,
      lambda = lambda[i],
      method = "BFGS",
      control = list(maxit = maxit)
    )

    if (results[[i]]$convergence != 0) {
      warning("Optimization did not converge at lambda: ", lambda[i])
    }

    curr_pars <- results[[i]]$par
  }

  X_dag <- matrix(0, ncol = d, nrow = d)
  X_dag[logdet_ix] <- results[[i]]$par

  if (threshold_to_DAG) {
    X_dag <- threshold_DAG(X_dag)
  }
  X_dag
}

#' Bootstrap LogDet
#'
#' Start with an effect matrix and standard errors and sample from an (independent) MV normal distribution
#' and project to DAG using LogDet characterization
#'
#' @param total_est total or direct effects matrix
#' @param total_est_se standard errors of total or direct effects matrix
#' @param reps total bootstrap samples
#' @param s logdet parameter
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
project_to_DAG_bootstrap <- function(
    total_est, total_est_se, reps = 100,
    s = 1.1,
    maxit = 500) {
  d <- ncol(total_est)
  non_diag_i <- -seq(1, d^2, by = d + 1)
  replicate(reps, {
    bootstrap_tot_est <- matrix(0, nrow = d, ncol = d)
    bootstrap_tot_est[non_diag_i] <- total_est[non_diag_i] + rnorm(d * (d - 1), mean = 0, sd = total_est_se[non_diag_i])
    project_to_DAG(
      X = bootstrap_tot_est,
      penalty = "L2",
      threshold_to_DAG = TRUE,
      lambda = c(1, 10^-seq(1,5), 0),
      s = s, # TODO: Why do we need s > 1? Always fails when s = 1
      maxit = maxit
    )
  }, simplify = FALSE)
  # TODO: Standard error for each configuration ?
}
