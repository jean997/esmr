#' @export
project_to_DAG <- function(
    X, s = max(max(X^2) + 0.1, 1), lambda = c(1, 10^-seq(1,5), 0), ix = NULL,
    penalty = c("L2", "L1"),
    threshold_to_DAG = FALSE,
    ...
    ) {
  penalty <- match.arg(penalty)
  d <- nrow(X)
  logdet_ix <- which(X != 0, arr.ind = TRUE)
  results <- list()
  curr_pars <- X[logdet_ix]
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
      control = list(...)
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

project_to_DAG_bootstrap <- function(
    total_est, total_est_se, reps = 100,
    s = 1.1,
    ...) {
  d <- ncol(total_est)
  non_diag_i <- -seq(1, d^2, by = d + 1)
  replicate(reps, {
    bootstrap_tot_est <- matrix(0, nrow = d, ncol = d)
    bootstrap_tot_est[non_diag_i] <- total_est[non_diag_i] + rnorm(d * (d - 1), mean = 0, sd = total_est_se[non_diag_i])
    project_to_DAG(
      bootstrap_tot_est,
      threshold_to_DAG = TRUE,
      maxit = 2000,
      trace = 5,
      s = s, # TODO: Why do we need s > 1? Always fails when s = 1
      ...
    )
  }, simplify = FALSE)
  # TODO: Standard error for each configuration ?
}
