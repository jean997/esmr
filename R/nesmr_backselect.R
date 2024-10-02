nesmr_backselect <- function(
    mod_list,
    beta_hat, se_beta_hat,
    mod_log_liks = NULL,
    Z_true = NULL,
    aic_cutoff = 2, pvalue_cutoff = 0.05,
    alpha = 5e-8) {
  # While there is still a significant edge
  # Take the lowest edge out of the model
  # Refit the model, keep model if stop rule is not hit
  # This is kind of a shortest path? Breadth first search?
  n_params <- sapply(mod_list, function(x) { sum(x$B_template) })

  if (is.null(mod_log_liks)) {
    mod_log_liks <- sapply(mod_list, function(x) {
      if ("log_lik" %in% names(x)) x$log_lik
      else log_py(x)
      })
  }

  mod_aic <- -2 * mod_log_liks + 2 * n_params

  best_mod_aic <- min(mod_aic)

  if (is.null(Z_true)) {
    Z <- beta_hat/se_beta_hat
  } else {
    Z <- Z_true
  }

  select_pval <- 2*pnorm(-abs(Z))
  minp <- apply(select_pval, 1, min)
  variant_ix <- which(minp < alpha)

  # Extract current p-values
  # B_template <- mod$B_template


  # Treat as a stack,

  queue <- rpqueue()

  for (mod in mod_list) {
    queue <- queue %>% insert_back(mod)
  }

  visited <- list()
  return_mods <- mod_list

  while(! rstackdeque::empty(queue)) {
    curr_mod <- peek_front(queue)
    queue <- without_front(queue)

    B_template <- curr_mod$B_template
    # B_template_chr <- paste(B_template, collapse = "")
    # if (B_template_chr %in% visited) {
    #   next
    # }

    # Do this to get only the maximum non-zero p-value
    edge_ix <- which(B_template != 0, arr.ind = TRUE)
    # if (nrow(edge_ix) <= 1) {
    #   print("TODO...>>")
    #   # return(list(mod))
    # }

    log_pvalues <- curr_mod$pvals_dm[edge_ix]
    # Filter out the ones < log(pvalue_cutoff)
    pvalue_order <- order(log_pvalues, decreasing = TRUE)
    log_pvalues <- log_pvalues[pvalue_order]
    edge_ix <- edge_ix[pvalue_order,, drop = FALSE]
    non_sig <- log_pvalues > log(pvalue_cutoff)
    n_non_sig <- sum(non_sig)

    # Might need to update this?
    if (n_non_sig <= 0 || all(! non_sig)) {
      #return_mods <- append(return_mods, list(new_mod))
      next
    }

    for (i in seq_len(n_non_sig)) {
      # TODO: Check here if we have already visited this configuration
      curr_edge <- edge_ix[i,, drop = FALSE]
      # Remove the edge with the highest p-value
      B_template[curr_edge] <- 0

      B_template_chr <- paste(B_template, collapse = "")
      if (B_template_chr %in% visited) {
        next
      } else {
        visited <- append(visited, B_template_chr)
      }

      new_mod <- esmr(
        beta_hat_X = beta_hat,
        se_X = se_beta_hat,
        variant_ix = variant_ix,
        G = diag(d),
        direct_effect_template = B_template,
        direct_effect_init = B_template * mod$direct_effects
      )

      new_mod$log_lik <- log_py(new_mod)
      new_mod$num_params <- sum(B_template)
      new_mod$aic <- -2 * new_mod$log_lik + 2 * new_mod$num_params

      # TODO: Should this really be from the original starting model or the best one?
      if (abs(best_mod_aic - new_mod$aic) <= aic_cutoff) {
        best_mod_aic <- min(best_mod_aic, new_mod$aic)
        return_mods <- append(return_mods, list(new_mod))
        queue <- queue %>% insert_back(new_mod)
      } else {
        # Stop looking at edges with lower likelihood
        #break
      }

      # Put the edge back
      B_template[curr_edge] <- 1
    }
  }

  return(return_mods)
}
