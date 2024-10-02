nesmr_backselect <- function(
    mod, beta_hat, se_beta_hat, Z_true = NULL,
    aic_cutoff = 2, pvalue_cutoff = 0.05,
    alpha = 5e-8, visited = list()) {
  print(match.call())
  print(mod$B_template)
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
  minp <- apply(select_pval, 1, min)
  variant_ix <- which(minp < alpha)

  # Extract current p-values
  B_template <- mod$B_template

  B_template_chr <- paste0(as.character(B_template), collapse = "")
  if (B_template_chr %in% visited) {
    return(list())
  }

  visited <- append(visited, B_template_chr)
  d <- nrow(B_template)
  # Do this to get only the maximum non-zero p-value
  edge_ix <- which(mod$B_template != 0, arr.ind = TRUE)
  if (nrow(edge_ix) <= 1) {
    return(list(mod))
  }

  log_pvalues <- mod$pvals_dm[edge_ix]
  # Filter out the ones < log(pvalue_cutoff)
  pvalue_order <- order(log_pvalues, decreasing = TRUE)
  log_pvalues <- log_pvalues[pvalue_order]
  edge_ix <- edge_ix[pvalue_order,, drop = FALSE]
  non_sig <- log_pvalues > log(pvalue_cutoff)
  n_non_sig <- sum(non_sig)
  return_mods <- list()

  if (n_non_sig <= 1 || all(! non_sig)) {
    return(list(mod))
  }

  # TODO: Do we want to put the edge back on when we're done?
  for (i in seq_len(n_non_sig)) {
#   while(log_pvalues[curr_pvalue_ix] > log(pvalue_cutoff)) {
    curr_edge <- edge_ix[i,, drop = FALSE]
    # Remove the edge with the highest p-value
    B_template[curr_edge] <- 0
    if (sum(B_template) <= 1) {
      next
    }
    # print(B_template)

    # Skip if we have visited node
    B_template_chr <- paste0(as.character(B_template), collapse = "")
    if (B_template_chr %in% visited) {
      next
    }
#    visited <- append(visited, B_template_chr)

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
      return_mods <- append(return_mods, list(new_mod))
      return_mods <- append(
        return_mods,
        # Recursively call... should turn this into non-recursive algo at some point
        nesmr_backselect(
          new_mod,
          beta_hat = beta_hat,
          se_beta_hat = se_beta_hat,
          Z_true = Z_true,
          aic_cutoff = aic_cutoff, pvalue_cutoff = pvalue_cutoff,
          alpha = alpha,
          visited = visited)
        )
    } else {
      # TODO: Do we put the edge back into the model?
      # Put the edge back into the model
      # B_template[curr_edge] <- 1
      # The way to explore everything would be to return the models with the edge,
      # and the models without the edge.
      # If we don't put the edge back then we are only exploring the greedy path
      # I think we instead want to break this branch entirely
      return_mods <- append(
        return_mods,
        # Recursively call... should turn this into non-recursive algo at some point
        nesmr_backselect(
          mod,
          beta_hat = beta_hat,
          se_beta_hat = se_beta_hat,
          Z_true = Z_true,
          aic_cutoff = aic_cutoff, pvalue_cutoff = pvalue_cutoff,
          alpha = alpha,
          visited = visited))

      B_template[curr_edge] <- 1
#      break
    }
    B_template_chr <- paste0(as.character(B_template), collapse = "")
    # if (B_template_chr %in% visited) {
    #   next
    # }
    visited <- append(visited, B_template_chr)
  }

  return(return_mods)
}
