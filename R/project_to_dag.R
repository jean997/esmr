project_to_DAG <- function(
    X, s = 1, lambda = c(1, 10^-seq(1,5), 0), ix = NULL,
    penalty = c("L1", "L2")) {
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
      control = list(maxit = 2000, trace = 0)
    )

    if (results[[i]]$convergence != 0) {
      warning("Optimization did not converge at lambda: ", lambda[i])
    }

    curr_pars <- results[[i]]$par
    print(curr_pars)
  }

  X_dag <- matrix(0, ncol = d, nrow = d)
  X_dag[logdet_ix] <- results[[i]]$par
  X_dag
}
