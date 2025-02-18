# This is the true KL (positive), but throughout NESMR KL is actually -KL
kl_mvn <- function(
    mu_q, Sigma_q,
    mu_g = 0,
    Sigma_g) {
  d <- length(mu_q)

  Sigma_g_inv <- solve(Sigma_g)
  # TODO: Error handling here ?
  log_det_ratio <- determinant(Sigma_g, logarithm = TRUE)$modulus - determinant(Sigma_q, logarithm = TRUE)$modulus

  # Equivalent but faster version of
  # trace_term <- sum(diag(Sigma_g_inv %*% Sigma_q))
  trace_term <- sum(t(Sigma_g_inv) * Sigma_q)

  mean_diff <- mu_g - mu_q
  quadratic_term <- t(mean_diff) %*% Sigma_g_inv %*% mean_diff

  kl_div <- 0.5 * (trace_term + quadratic_term - d + log_det_ratio)

  return(as.numeric(kl_div))
}
