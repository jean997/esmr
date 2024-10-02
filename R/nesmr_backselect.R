#' NESMR backselect algorithm
#'
#' @param mod_list List of initial models to search
#' @param beta_hat Estimated beta coefficients from GWAS
#' @param se_beta_hat Standard errors of beta coefficients from GWAS
#' @param mod_log_liks Optional log-likelihoods of the model. Otherwise, look at model$log_lik and then attempt to call log_py(model)
#' @param Z_true Optional true Z-scores to use for selection
#' @param method Currently only does AIC based select within `aic_cutoff` of the best current model
#' @param aic_cutoff AIC difference cutoff for selecting models
#' @param pvalue_cutoff P-value cutoff for selecting edges
#' @param alpha Significance level for selecting variants
#'
#' @return List of models selected by NESMR backselect
#' @export
#'
#' @examples
#' set.seed(13)
#' d <- 4
#' G <- matrix(0, nrow = d, ncol = d)
#' G_low <- lower.tri(G)
#' B_lower <- G_low + 0
#' G[d, 1] <- 0.1
#'
#' B_true <- (G != 0) + 0
#'
#' h2 <- 0.3
#' J <- 5000
#' N <- 30000
#' pi_J <- 0.1
#' alpha <- 5e-8
#' lambda <- qnorm(1 - alpha / 2)
#' dat <- GWASBrewer::sim_mv(
#'   G = G,
#'   N = N,
#'   J = J,
#'   h2 = h2,
#'   pi = pi_J,
#'   sporadic_pleiotropy = TRUE,
#'   est_s = TRUE
#' )
#'
#' Ztrue <- with(dat, beta_marg/se_beta_hat)
#' pval_true <- 2*pnorm(-abs(Ztrue))
#' minp <- apply(pval_true, 1, min)
#' ix <- which(minp < alpha)
#'
#' backselect_results <- nesmr_backselect(
#'   list(full_mod),
#'   beta_hat = dat$beta_hat,
#'   se_beta_hat = dat$s_estimate, Z_true = Ztrue,
#'   aic_cutoff = 10, pvalue_cutoff = 0.05 / sum(B_lower),
#'   alpha = 5e-8)
#'
#' # Print the AIC of the models and indicate the true generating model
#' back_aic <- sapply(backselect_results, function(x) x$aic)
#' min_aic <- which.min(back_aic)
#' true_ix <- which(sapply(backselect_results, function(x) all(x$B_template == B_true)))
#' plot(min(back_aic) - back_aic,
#'     pch = ifelse(1:length(back_aic) == true_ix, 8, 1)
#' )
#' legend("topleft", legend = c("True model"),
#'       pch = c(8, 1))
nesmr_backselect <- function(
    mod_list,
    beta_hat, se_beta_hat,
    mod_log_liks = NULL,
    Z_true = NULL,
    method = c("aic", "posterior_prob"),
    aic_cutoff = 2,
    pvalue_cutoff = 0.05,
    alpha = 5e-8) {
  method <- match.arg(method)
  n_params <- sapply(mod_list, function(x) { sum(x$B_template) })

  if (is.null(mod_log_liks)) {
    mod_list <- lapply(mod_list, function(x) {
      if (!"log_lik" %in% names(x)) x$log_lik <- log_py(x)
      if (! "aic" %in% names(x)) x$aic <- -2 * x$log_lik + 2 * sum(x$B_template)

      return(x)
    })

    mod_log_liks <- sapply(mod_list, function(x) x$log_lik)
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

  # Use a queue to simulate breadth first search
  # First we remove all possible edges from the starting model and check stopping criterion
  # Then, we look at all possible edges from each of these subsets
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

    # Do this to get only the maximum non-zero p-value
    edge_ix <- which(B_template != 0, arr.ind = TRUE)

    log_pvalues <- curr_mod$pvals_dm[edge_ix]
    # Filter out the ones < log(pvalue_cutoff)
    pvalue_order <- order(log_pvalues, decreasing = TRUE)
    log_pvalues <- log_pvalues[pvalue_order]
    edge_ix <- edge_ix[pvalue_order,, drop = FALSE]
    non_sig <- log_pvalues > log(pvalue_cutoff)
    n_non_sig <- sum(non_sig)

    # If all edges are significant, skip this model
    # The model will already be in the return_mods list from the previous iteration
    if (n_non_sig <= 0 || all(! non_sig)) {
      next
    }

    for (i in seq_len(n_non_sig)) {
      curr_edge <- edge_ix[i,, drop = FALSE]
      # Remove the edge with the highest p-value
      B_template[curr_edge] <- 0

      # Check if we have already visited this configuration
      # If so, skip all checks as it is either in the model list or the queue
      B_template_chr <- paste(B_template, collapse = "")
      if (B_template_chr %in% visited) {
        next
      } else {
        visited <- append(visited, B_template_chr)
      }

      # Fit the new model
      # Currently use the same variants across all models but this could change
      new_mod <- esmr(
        beta_hat_X = beta_hat,
        se_X = se_beta_hat,
        variant_ix = variant_ix,
        G = diag(d),
        direct_effect_template = B_template,
        direct_effect_init = B_template * mod$direct_effects
      )

      # TODO: Eventually want to do something faster than computing likelihood each time
      new_mod$log_lik <- log_py(new_mod)
      new_mod$num_params <- sum(B_template)
      new_mod$aic <- -2 * new_mod$log_lik + 2 * new_mod$num_params

      # TODO: Generalize this to "stopping criterion" function
      # Either aic difference, or posterior probability with different priors
      if (abs(best_mod_aic - new_mod$aic) <= aic_cutoff) {
        best_mod_aic <- min(best_mod_aic, new_mod$aic)
        return_mods <- append(return_mods, list(new_mod))
        queue <- queue %>% insert_back(new_mod)
      } else {
        # TODO: Can we stop looking at edges with lower likelihood?
        # I think we can make the claim that if p-value for an edge is higher then
        # the likelihood will be lower but this is not guaranteed to be true
        # for all subsets of the graph
        # break
      }

      # Put the edge back
      B_template[curr_edge] <- 1
    }
  }

  return(return_mods)
}
