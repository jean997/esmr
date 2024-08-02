nesmr_all_permn <- function(
    beta_hat,
    se_beta_hat,
    B_templates = NULL,
    posterior_probs = TRUE,
    ...
    ) {
  if (is.null(B_templates)) {
    d <- ncol(beta_hat)
    B_lower <- matrix(0, nrow = d, ncol = d)
    B_lower[lower.tri(B_lower)] <- 1
    if (d >= 9) {
      warning(
      sprintf(
        "With %d traits there are %d! = %d permutations. This will be slow...consider supplying only the permutations with parameter `B_templates`",
        d, d, factorial(d))
      )
    }
    all_perms <- combinat::permn(seq_len(d))
    B_templates <- lapply(all_perms, function(perm) {
      B_lower[perm, perm]
    })
  }

  nesmr_models <- lapply(B_templates, function(B) {
    res <- esmr(
      beta_hat_X = beta_hat,
      se_X = se_beta_hat,
      G = diag(d),
      direct_effect_template = B,
      ...
      )
    res$log_lik <- log_py(res)
    res
  })

  if (posterior_probs) {
    # This is basically just softmax
    mod_log_lik <- sapply(nesmr_models, function(x) x$log_lik)
    log_denom <- log_sum_exp(mod_log_lik)
    posterior_probs <- exp(mod_log_lik - log_denom)
    return(
      list(
        nesmr_models = nesmr_models,
        posterior_probs = posterior_probs
      )
    )
  }

  return(nesmr_models)
}
