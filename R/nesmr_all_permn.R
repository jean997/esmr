# From chatGPT: https://chatgpt.com/c/35ed431d-dea5-4dcb-bdb3-7425395d83cf
# Function to generate all possible DAGs for a given matrix size n
#'@export
generate_dags_with_permutations <- function(n, remove_empty = TRUE) {
  # Generate all permutations of node indices (i.e., topological sorts)
  node_permutations <- gtools::permutations(n, n)  # All permutations of length n

  # Get the number of possible edges in the upper triangular part (excluding diagonal)
  num_edges <- choose(n, 2)

  # Generate all possible combinations of edges (0 or 1) for the upper triangular part
  all_combinations <- expand.grid(replicate(num_edges, c(0, 1), simplify = FALSE))

  # List to store all generated DAGs
  dag_list <- list()

  # Loop over each permutation of nodes
  for (perm in 1:nrow(node_permutations)) {
    node_perm <- node_permutations[perm, ]

    # Loop over each combination of edges in the upper triangular matrix
    for (i in 1:nrow(all_combinations)) {
      # Create an empty adjacency matrix
      adj_matrix <- matrix(0, n, n)

      # Fill the upper triangular part with the current combination of edges
      upper_tri_indices <- which(upper.tri(adj_matrix))
      adj_matrix[upper_tri_indices] <- as.numeric(all_combinations[i, ])

      # Apply the permutation of nodes
      permuted_matrix <- adj_matrix[node_perm, node_perm]

      # Add the generated DAG (permuted adjacency matrix) to the list
      dag_list[[length(dag_list) + 1]] <- permuted_matrix
    }
  }

  dag_list <- unique(dag_list)

  if (remove_empty) {
    dag_list <- dag_list[sapply(dag_list, function(x) sum(x) > 0)]
  }

  # TODO: Optimize this since we generate duplicates
  return(dag_list)
}

#' Run NESMR on all permutation of the traits or all DAGs
#'
#'
#' @param beta_hat Matrix of effect estimates
#' @param se_beta_hat
#' @param B_templates
#' @param all_DAGs Should all DAGs be generated? FALSE means all permutations of the traits
#' @param posterior_probs Should posterior probabilities be computed?
#' @param return_model Should the full model object be returned?
#' @param direct_effect_init Initial direct effect matrix for the current order of the traits
#' @param variant_ix Index of variants to subset
#' @param alpha Significance level for p-value cutoff
#' @param ... Additional arguments to pass to `esmr`
#'
#' @return List of results
#' @export
nesmr_all_permn <- function(
    beta_hat,
    se_beta_hat,
    B_templates = NULL,
    all_DAGs = FALSE,
    posterior_probs = TRUE,
    return_model = FALSE,
    direct_effect_init = NULL,
    variant_ix = NULL,
    alpha = 5e-8,
    ...
) {
  d <- ncol(beta_hat)
  if (!is.null(direct_effect_init)) {
    stopifnot(all(dim(direct_effect_init) == c(d, d)))
  } else {
    direct_effect_init <- matrix(0, nrow = d, ncol = d)
  }

  if (is.null(B_templates)) {
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
      B_templates <- generate_dags_with_permutations(d)
    }

    B_templates
  }

  nesmr_models <- lapply(B_templates, function(B) {
    if (sum(B) == 0) {
      return(NULL)
    }
    res <- tryCatch({
      B_direct_effect_init <- direct_effect_init * B

      z <- esmr(
        beta_hat_X = beta_hat,
        se_X = se_beta_hat,
        G = diag(d),
        direct_effect_template = B,
        direct_effect_init = B_direct_effect_init,
        variant_ix = variant_ix,
        ...
      )
      z$num_params <- sum(B)
      z$log_lik <- log_py(z)
      z$aic <- -2 * z$log_lik + 2 * z$num_params
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

  nesmr_models <- nesmr_models[sapply(nesmr_models, Negate(is.null))]

  rtn <- list()

  if (return_model) {
    rtn$nesmr_models <- nesmr_models
  }
  rtn$direct_effects <- lapply(nesmr_models, function(x) x$direct_effects)
  rtn$se_beta_hat <- lapply(nesmr_models, function(x) x$se_dm)
  rtn$B_template <- B_templates
  rtn$beta_hat <- lapply(nesmr_models, function(x) x$direct_effects)
  rtn$pvals_dm <- lapply(nesmr_models, function(x) x$pvals_dm)
  rtn$log_lik <- sapply(nesmr_models, function(x) x$log_lik)

  if (posterior_probs) {
    mod_log_lik <- sapply(nesmr_models, function(x) x$log_lik)
    mod_aic <- sapply(nesmr_models, function(x) x$aic)
    mod_params <- sapply(nesmr_models, function(x) x$num_params)
    log_denom <- log_sum_exp(mod_log_lik)
    posterior_probs <- exp(mod_log_lik - log_denom)
    rtn$posterior_probs <- posterior_probs

    # Posterior probs weighted by number of edges
    log_denom_n_edge <- log_sum_exp(mod_log_lik - log(mod_params))
    rtn$posterior_probs_n_edge <- exp((mod_log_lik - log(mod_params)) - log_denom_n_edge)

    #rtn$posterior_probs_sparse <-
    # TODO: Maybe change posterior_probs argument to different name...
    rtn$aic <- mod_aic
  }

  return(rtn)
}
