#' @export
nesmr_all_permn <- function(
    beta_hat,
    se_beta_hat,
    B_templates = NULL,
    all_DAGs = FALSE,
    posterior_probs = TRUE,
    return_model = FALSE,
    catch_errors = TRUE,
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

    all_perms <- combinat::permn(d)
    if (! all_DAGs) {
      B_templates <- lapply(all_perms, function(perm) {
        B_lower[perm, perm]
      })
    } else {

      B_templates <- unique(unlist(lapply(seq_len(d), function(i) {
        keep_indices <- combinat::combn(seq_len(d), i, simplify = F)
        curr_dags <- unlist(lapply(keep_indices, function(one_index) {
          B_lower_combn <- B_lower
          B_lower_combn[lower.tri(B_lower_combn)][-one_index] <- 0
          unique(lapply(all_perms, function(perm) {
            B_lower_combn[perm, perm]
          }))
        }), recursive = FALSE)
      }), recursive = FALSE))

      # Add the null templates
      # B_templates <- c(B_templates, list(matrix(0, nrow = d, ncol = d)))
    }

    B_templates
  }

  nesmr_models <- lapply(B_templates, function(B) {
    res <- tryCatch({
      z <- esmr(
        beta_hat_X = beta_hat,
        se_X = se_beta_hat,
        G = diag(d),
        direct_effect_template = B,
        ...
      )
      z$log_lik <- log_py(z)
      z
    }, error = function(e) {
      warning(e)
      list(
        error = e,
        log_lik = NA
      )
    })
    res
  })

  rtn <- list()

  if (return_model) {
    rtn$nesmr_models <- nesmr_models
  } else {
    rtn$direct_effects <- lapply(nesmr_models, function(x) x$direct_effects)
    rtn$se_beta_hat <- lapply(nesmr_models, function(x) x$se_dm)
    rtn$B_template <- B_templates
    rtn$beta_hat <- lapply(nesmr_models, function(x) x$direct_effects)
    rtn$pvals_dm <- lapply(nesmr_models, function(x) x$pvals_dm)
    rtn$log_lik <- sapply(nesmr_models, function(x) x$log_lik)
  }

  if (posterior_probs) {
    mod_log_lik <- sapply(nesmr_models, function(x) x$log_lik)
    log_denom <- log_sum_exp(mod_log_lik)
    posterior_probs <- exp(mod_log_lik - log_denom)
    rtn$posterior_probs <- posterior_probs
  }

  return(rtn)
}
