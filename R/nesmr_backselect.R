nesmr_backselect <- function(
    mod, beta_hat, se_beta_hat, Z_true = NULL,
    aic_cutoff = 2, pvalue_cutoff = 0.05,
    alpha = 5e-8) {
  # While there is still a significant edge
  # Take the lowest edge out of the model
  # Refit the model, keep model if stop rule is not hit
  # This is kind of a shortest path? Breadth first search?
  if (! "log_lik" %in% names(mod)) {
    mod$log_lik <- log_py(mod)
  }
  if (! "aic" %in% names(mod)) {
    mod$num_params <- sum(mod$B_template)
    mod$aic <- -2 * mod$log_lik + 2 * mod$num_params
  }

  if (is.null(Z_true)) {
    Z <- beta_hat/se_beta_hat
  } else {
    Z <- Z_true
  }
  select_pval <- 2*pnorm(-abs(Z))

  # Extract current p-values
  B_template <- mod$B_template
  d <- nrow(B_template)
  # Do this to get only the maximum non-zero p-value
  edge_ix <- which(mod$B_template != 0, arr.ind = TRUE)
  log_pvalues <- mod$pvals_dm[edge_ix]
  max_pvalue_ix <- which.max(log_pvalues)
  return_mods <- list()
  # TODO: Do we want to put the edge back on when we're done?
  while(log_pvalues[max_pvalue_ix] > log(pvalue_cutoff)) {
    print(B_template)
    curr_edge <- edge_ix[max_pvalue_ix,, drop = FALSE]
    # Remove the edge with the highest p-value
    B_template[curr_edge] <- 0
    # Refit the model
    trait_ix <- which(rowSums(B_template) > 0)
    if (length(trait_ix) <= 0) {
      stop("No traits left in model. Stopping.")
    }

    minp <- apply(select_pval[, trait_ix, drop = FALSE], 1, min)
    variant_ix <- which(minp < alpha)

    new_mod <- esmr(
      beta_hat_X = beta_hat,
      se_X = se_beta_hat,
      variant_ix = variant_ix,
      G = diag(d),
      direct_effect_template = B_template)

    new_mod$log_lik <- log_py(new_mod)
    new_mod$num_params <- sum(B_template)
    new_mod$aic <- -2 * new_mod$log_lik + 2 * new_mod$num_params

    if (abs(new_mod$aic - mod$aic) <= aic_cutoff) {
      # Keep this model and recursively call
      # ...
      return_mods <- append(return_mods, list(new_mod))
      return_mods <- append(
        return_mods,
        # Recursively call... should turn this into non-recursive algo at some point
        nesmr_backselect(new_mod, beta_hat, se_beta_hat, Z_true)
        )
    }

    # TODO: Should we select the next highest p-value? Put the edge back in the graph?
    edge_ix <- which(B_template != 0, arr.ind = TRUE)
    log_pvalues <- mod$pvals_dm[edge_ix]
    max_pvalue_ix <- which.max(log_pvalues)
    # Put the edge back into the model
    # B_template[curr_edge] <- 1
  }

  return(return_mods)
}
