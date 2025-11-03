nesmr_explore <- function(
  beta_hat, se_beta_hat, pval_select = NULL,
  zscore_filter = qnorm(0.975), chains = 5, elbo_threshold = 0.05,
  graph_prior_pi = 0.5,
  alpha = 5e-8,
  R = NULL) {
    d <- ncol(beta_hat)
    max_edges <- d * (d - 1) / 2
    if (is.null(pval_select)) {
      Z_cursed <- beta_hat/se_beta_hat
      pval_cursed <- 2 * pnorm(-abs(Z_cursed))
      pval_select <- pval_cursed
    }
    minp <- apply(pval_select, 1, min)
    ix <- which(minp < alpha)

    all_mvmr_mod <- esmr::nesmr_complete_mvmr(
        beta_hat = beta_hat,
        se_beta_hat = se_beta_hat,
        pval_select = pval_select,
        R = R
    )

    # TODO: Replace this with a true hash
    visited_graphs <- list()
    elbo_denom <- NA

    .draw_new_edges <- function(g, weight_mat = NULL) {
        g_comp <- igraph::complementer(g, loops = FALSE)

        # Get all candidate edges as pairs (as a vector: from, to)
        candidates <- igraph::as_edgelist(g_comp, names = FALSE)
        if (!is.null(weight_mat)) {
            # Get the weights of the candidates
            candidate_weights <- weight_mat[candidates]
            candidate_ordering <- order(abs(candidate_weights), decreasing = TRUE)
            candidate_ordering <- candidate_ordering[abs(candidate_weights[candidate_ordering]) > zscore_filter]
        } else {
            candidate_ordering <- seq_len(nrow(candidates))
        }

        candidate_graphs <- Filter(Negate(is.null), lapply(candidate_ordering, function(i) {
            from <- candidates[i, 1]
            to <- candidates[i, 2]

            # Try adding the edge
            g_test <- igraph::add_edges(ig, c(from, to))

            g_B <- igraph::as_adjacency_matrix(g_test, sparse = FALSE)
            g_B_str <- paste0(g_B, collapse = "")
            if (g_B_str %in% names(visited_graphs)) {
                return(NULL)
            }

            is_dag_test <- igraph::is_dag(g_test)
            if (is_dag_test) {
                return(g_test)
            } else {
              return(NULL)
            }
        }))

        return(candidate_graphs)
    }

    queue <- rstackdeque::rpqueue()

    full_graph_zscores <- all_mvmr_mod$beta_hat / all_mvmr_mod$se_beta_hat
    diag(full_graph_zscores) <- 0
    non_diag_i <- -seq(1, d^2, by = d + 1)
    for (chain_i in seq_len(chains)) {
      print(sprintf("Starting chain: %s", chain_i))
      # Note: This could be outside of the while or inside..
      mat_init <- matrix(0, nrow = d, ncol = d)
      # TODO: Should we just start from the "best" graph? no noise?
      mat_init[non_diag_i] <- full_graph_zscores[non_diag_i]
      # If we have more than one chain, then we add noise to each chain
      if (chain_i > 1) {
        mat_init[non_diag_i] <- mat_init[non_diag_i] + rnorm(max_edges)
      }

      init_filter_zscore <- mat_init * (abs(mat_init) > zscore_filter)
      if (sum(init_filter_zscore != 0) == 0) {
        warning("No edges found in initial graph. Skipping this chain.")
        next
      }
      # Note: Better to do L2 or L1 norm?
      curr_adj_mat <- sqrt(maximal_acyclic_subgraph((init_filter_zscore)^2)) * sign(init_filter_zscore)
      # TODO: Should we fit this initial model ? Probably...
      curr_B <- (curr_adj_mat != 0) + 0

      # Now we expore the graph starting from curr_adj_mat
      # Note: Not sure if it really matter if we have weighted or not...
      ig <- igraph::graph_from_adjacency_matrix(
          curr_adj_mat != 0,
          mode = 'directed'
          )

      queue <- queue %>% rstackdeque::insert_back(ig)
      while(! is.null(ig)) {
        curr_B <- igraph::as_adjacency_matrix(ig, sparse = FALSE)
        curr_B_str <- paste0(curr_B, collapse = "")
        if (curr_B_str %in% names(visited_graphs)) {
            visited_graphs[[curr_B_str]]$visited_count <- visited_graphs[[curr_B_str]]$visited_count + 1
        } else {
          capture.output({
            new_mod <- esmr::esmr(
                beta_hat_X = beta_hat,
                se_X = se_beta_hat,
                variant_ix = ix,
                G = diag(d),
                direct_effect_template = curr_B,
                max_iter = 300,
                restrict_dag = T,
                R = R,
                beta_prior_cov = 1 # TODO: Make these parameters?
                )
          }, file = nullfile())
            k <- sum(curr_B != 0)
            new_elbo <- new_mod$elbo + log_graph_prior(k, d, pi_0 = graph_prior_pi)
            visited_graphs[[curr_B_str]] <- list(
                elbo = new_elbo,
                beta_hat = new_mod$beta_mat$beta_hat,
                se_beta_hat = new_mod$beta_mat$beta_se,
                visited_count = 1,
                chain = chain_i
                )

            old_elbo_denom <- elbo_denom
            elbo_denom <- matrixStats::logSumExp(c(elbo_denom, new_elbo), na.rm = TRUE)

            new_elbo_prop <- exp(new_elbo - elbo_denom)
            print(sprintf("New graph elbo: %s old elbo denom: %s", round(new_elbo, 4), round(old_elbo_denom, 4)))
            print(sprintf("New graph elbo is %s proportion of total", round(new_elbo_prop, 4)))

            if (new_elbo_prop > elbo_threshold) {
              for (e in igraph::E(ig)) {
                new_dg <- igraph::delete_edges(ig, e)
                new_B <- igraph::as_adjacency_matrix(new_dg, sparse = FALSE)
                new_B_str <- paste0(new_B, collapse = "")
                if (! new_B_str %in% names(visited_graphs) && sum(new_B) > 0) {
                    queue <- queue %>% rstackdeque::insert_back(new_dg)
                }
              }
            }

            canidate_add_edges <- .draw_new_edges(ig, weight_mat = mat_init)
            # TODO: Does it make sense to add all these edges?
            # I think we should change this to breadth/best first search:
            # Here only add the best edge (highest abs Z-score)
            # If we hit a graph that does not add lower elbo, then do not search any more of the subgraphs
            # E.g. might need to add these to the "visited graph" and remove them from possible candidates
            for (g in canidate_add_edges) {
                queue <- queue %>% rstackdeque::insert_back(g)
            }
        }

        if (rstackdeque::empty(queue)) {
            ig <- NULL
        } else {
          print(sprintf("Queue size: %s", length(queue)))
          ig <- rstackdeque::peek_front(queue)
          queue <- rstackdeque::without_front(queue)
        }
      }
    }
    return(visited_graphs)
}