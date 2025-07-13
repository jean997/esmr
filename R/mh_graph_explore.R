mh_graph_explore <- function(
  dat,
  n_mvmr_res = NULL,
  pval_select = NULL,
  R = NULL,
  alpha = 5e-8,
  max_Z = 5,
  max_prob = 0.7,
  min_prob = 0.01,
  init_prob_threshold = 0.1,
  sparse_chain = FALSE,
  dense_chain = FALSE,
  max_iter = 1000,
  max_nesmr_fits = 100,
  visited_graphs = list(),
  checkpoint_file = NULL,
  checkpoint_every = 0
  ) {
  d <- ncol(dat$beta_hat)
  max_edges <- d * (d - 1)
  if (is.null(pval_select)) {
      dat_Z <- dat$beta_hat / dat$s_estimate
      pval_select <- 2 * pnorm(-abs(dat_Z))
  }

  minp <- apply(pval_select, 1, min)
  ix <- which(minp < alpha)

  if (is.null(n_mvmr_res)) {
    ## Start function here
    n_mvmr_res <- esmr::nesmr_complete_mvmr(
        beta_hat = dat$beta_hat,
        se_beta_hat = dat$s_estimate,
        pval_select = pval_select,
        R = R
    )
  }

  # TODO: Replace this with a true hash
  #visited_graphs <- list()

  ## For MH algorithm we need to compute g(M|M') and g(M'|M)

  full_graph_zscores <- n_mvmr_res$beta_hat / n_mvmr_res$se_beta_hat
  stopifnot(ncol(n_mvmr_res$beta_hat) == d)
  non_diag_i <- -seq(1, d^2, by = d + 1)
  diag(full_graph_zscores) <- 0

  # TODO: Move this into own function
  # Compute the Z -> probability based on the MVMR full graph
  z_max <- max(abs(full_graph_zscores[non_diag_i]), na.rm = TRUE)
  #    q_min <- qlogis(min_prob)
  #    q_max <- qlogis(max_prob)
  #    beta_z <- (q_min - q_max) / (z_min - z_max)
  #    beta_0 <- q_max - beta_z * z_max
  z <- c(0, min(z_max, max_Z))
  .y <- c(min_prob, max_prob)
  # Note: this closure keeps mod_coefs in scope but does not same GLM object itself
  Z_to_prob <- esmr:::get_Z_to_prob(z, .y)

  # z_range <- seq(z_min, z_max, length.out = 1000)
  # plot(z_range, Z_to_prob(z_range), type = "l")
  # Z-score filter:

  ## Rather than randomly sampling from Z-scores, instead start from:
  # Maximal acyclic subgraph (without filtering); Should have as many edges as possible
  # Maximal acyclic subgraph (with filtering); Should be close to the true graph
  # A graph with single best Z-score

  mh_chain_init <- list(
    best_approx = {
      mat_init <- matrix(0, nrow = d, ncol = d)

      mat_init[non_diag_i] <- full_graph_zscores[non_diag_i]
      edge_prob_matrix <- Z_to_prob(abs(mat_init))
      init_filter_zscore <- mat_init * (edge_prob_matrix > init_prob_threshold)
      sqrt(esmr:::maximal_acyclic_subgraph((init_filter_zscore)^2)) * sign(init_filter_zscore)
    }
  )

  if (sparse_chain) {
    mh_chain_init$min_graph <- {
        mat_init <- matrix(0, nrow = d, ncol = d)
        mat_init[non_diag_i] <- 0
        max_index <- which(abs(full_graph_zscores[non_diag_i]) == max(abs(full_graph_zscores[non_diag_i]), na.rm = TRUE))
        mat_init[non_diag_i][max_index] <- full_graph_zscores[non_diag_i][max_index]
        mat_init
    }
  }

  if (dense_chain) {
    mh_chain_init$max_graph <- {
        mat_init <- matrix(0, nrow = d, ncol = d)
        mat_init[non_diag_i] <- full_graph_zscores[non_diag_i]
        sqrt(esmr:::maximal_acyclic_subgraph((mat_init)^2)) * sign(mat_init)
    }
  }

  mh_chain <- list()
  mh_accept <- list()
  mh_elbo_chain <- list()
  for (i in seq_along(mh_chain_init)) {
      # Note: This could be outside of the while or inside..
    curr_adj_mat <- mh_chain_init[[i]]
#       print(sprintf("Initial graph: correct = %s", B_correct_str == paste0((curr_adj_mat != 0) + 0, collapse = "")))
      print(curr_adj_mat)
      # Collect the graphs as flattened strings as we go
      # TODO: Should we fit this initial model ? Probably...
      curr_B <- (curr_adj_mat != 0) + 0
      curr_B_str <- paste0(curr_B, collapse = "")

      if (is.null(visited_graphs[[curr_B_str]])) {

        # Initial NESMR fit
#        capture.output({
        start_time <- Sys.time()
        init_mod <- esmr::esmr(
            beta_hat_X = dat$beta_hat,
            se_X = dat$s_estimate,
            variant_ix = ix,
            G = diag(d),
            direct_effect_template = curr_B,
            max_iter = 300,
            restrict_dag = T,
            beta_prior_cov = 1,
            R = R
        )
        end_time <- Sys.time()
        init_fit_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
        print(sprintf("Initial fit time: %.2f seconds", init_fit_time))
        print(sprintf("Expect the total time to be around %.2f minutes", init_fit_time * max_nesmr_fits / 60))

#        }, file = nullfile())

        visited_graphs[[curr_B_str]]$elbo <- init_mod$elbo
        visited_graphs[[curr_B_str]]$beta_hat <- init_mod$direct_effects
        visited_graphs[[curr_B_str]]$se_beta_hat <- init_mod$se_dm
        visited_graphs[[curr_B_str]]$proposed <- (visited_graphs[[curr_B_str]]$proposed %||% 0) + 1
      }

      elbo_denom <- visited_graphs[[curr_B_str]]$elbo

      mh_chain[[i]] <- list(curr_B_str)
      mh_accept[[i]] <- 1
      mh_elbo_chain[[i]] <- visited_graphs[[curr_B_str]]$elbo
      # Now we expore the graph starting from curr_adj_mat
      # Note: Not sure if it really matter if we have weighted or not...
      ig <- igraph::graph_from_adjacency_matrix(
          curr_adj_mat != 0,
          mode = "directed"
      )
      iter <- 1
      nesmr_fits <- 1
      #while (hit_old_graph <= no_new_graph_limit && iter < max_iter) {
      while(iter < max_iter && nesmr_fits < max_nesmr_fits) {
          print(curr_B)
          if (curr_B_str %in% names(visited_graphs) && !is.null(visited_graphs[[curr_B_str]]$adj_graph_info)) {
              adj_graph_info <- visited_graphs[[curr_B_str]]$adj_graph_info
          } else {
              visited_graphs[[curr_B_str]]$adj_graph_info <- esmr:::get_adjacent_graphs(ig, edge_prob_matrix)
              adj_graph_info <- visited_graphs[[curr_B_str]]$adj_graph_info
          }

          print(adj_graph_info)
          candidate_draw <- esmr:::draw_graph(ig, adj_graph_info)
          tmp_ig <- candidate_draw$g
          # Denominator: h(G'|G)
          prop_denom <- candidate_draw$prob

          # Same procedure for if we are removing or taking away:
          prop_B <- igraph::as_adjacency_matrix(tmp_ig, sparse = FALSE)
          prop_B_str <- paste0(prop_B, collapse = "")
          proposal_graph_info <- visited_graphs[[prop_B_str]]

          if (is.null(proposal_graph_info)) {
              proposal_graph_info <- list()
          }

          # Check if we have neighboring graph information
          if (is.null(proposal_graph_info$adj_graph_info)) {
              proposal_graph_info$adj_graph_info <- esmr:::get_adjacent_graphs(tmp_ig, edge_prob_matrix)
              # Get the remove_candidate probability that we are removing
          }

          # Get h(G|G') - Reverse direction
          # If we added the edge: Check the prob for removing the edge
          # If we removed the edge: Check the prob for adding the edge
          if (candidate_draw$insert_edge) {
              # We added the edge
              remove_candidates <- proposal_graph_info$adj_graph_info$remove_candidates
              remove_edge_prob <- proposal_graph_info$adj_graph_info$remove_edge_prob
              remove_edge_ix <- which(apply(remove_candidates, 1, function(x) {
                  paste0(x, collapse = "|")
              }) == candidate_draw$mod_edge)
              # Numerator: h(G|G')
              prop_num <- remove_edge_prob[remove_edge_ix]
          } else {
              # We removed the edge
              add_candidates <- proposal_graph_info$adj_graph_info$add_candidates
              add_candidate_prob <- proposal_graph_info$adj_graph_info$add_candidate_prob
              add_candidate_ix <- which(apply(add_candidates, 1, function(x) {
                  paste0(x, collapse = "|")
              }) == candidate_draw$mod_edge)
              # Numerator: h(G|G')
              prop_num <- add_candidate_prob[add_candidate_ix]
          }

          if (is.null(proposal_graph_info$elbo)) {
              # If we have zero edges; continue
              # Eventually esmr should support having zero edges
              if (sum(prop_B) == 0) {
                  print("Zero edges; continue")
                  mh_chain[[i]] <- append(mh_chain[[i]], curr_B_str)
                  mh_accept[[i]] <- append(mh_accept[[i]], 0)
                  mh_elbo_chain[[i]] <- append(mh_elbo_chain[[i]], visited_graphs[[curr_B_str]]$elbo)
                  iter <- iter + 1
                  next
              }

              # TODO: Probably want to just re-fit from initial/previous chain
#              capture.output(
#                  {
                      new_mod <- esmr::esmr(
                          beta_hat_X = dat$beta_hat,
                          se_X = dat$s_estimate,
                          variant_ix = ix,
                          G = diag(d),
                          direct_effect_template = prop_B,
                          max_iter = 300,
                          restrict_dag = T,
                          beta_prior_cov = 1,
                          R = R
                      )
                      nesmr_fits <- nesmr_fits + 1
                      elbo_denom <- matrixStats::logSumExp(c(elbo_denom, new_mod$elbo), na.rm = TRUE)
#                  },
#                  file = nullfile()
#              )

              visited_graphs[[prop_B_str]]$elbo <- new_mod$elbo
              visited_graphs[[prop_B_str]]$beta_hat <- new_mod$direct_effects
              visited_graphs[[prop_B_str]]$se_beta_hat <- new_mod$se_dm
              proposal_graph_info <- visited_graphs[[prop_B_str]]
          }

          visited_graphs[[prop_B_str]] <- proposal_graph_info
          visited_graphs[[prop_B_str]]$proposed <- (visited_graphs[[prop_B_str]]$proposed %||% 0) + 1

          elbo_diff <- proposal_graph_info$elbo - visited_graphs[[curr_B_str]]$elbo
          print(sprintf("ELBO diff: %s", round(elbo_diff, 4)))
          print(sprintf("ELBO curr: %s", round(visited_graphs[[curr_B_str]]$elbo, 4)))
          print(sprintf("ELBO denom: %s", round(elbo_denom, 4)))
          print(sprintf("exp(ELBO curr - ELBO denom): %s", round(exp(visited_graphs[[curr_B_str]]$elbo - elbo_denom), 4)))


          # TODO: Switch to log scale for everything
          prop_ratio <- prop_num / prop_denom
          print(sprintf("Proposal ratio: %s", round(prop_ratio, 4)))

          # Check accept/reject
          # max(1, exp(elbo(tmp_ig) - elbo(ig)))
          accept_prob <- min(
              1, exp(elbo_diff) * prop_ratio
          )
          print(sprintf("Total proposal ratio: %s", round(exp(elbo_diff) * prop_ratio, 4)))


        # Note: This is not really the "chain elbo" but rather the elbo of the proposal at each step
          mh_elbo_chain[[i]] <- append(mh_elbo_chain[[i]], proposal_graph_info$elbo)

          mh_chain[[i]] <- append(mh_chain[[i]], prop_B_str)

          if (accept_prob == 1 || runif(1) < accept_prob) {
              curr_B <- prop_B
              curr_B_str <- prop_B_str
              ig <- tmp_ig
              visited_graphs[[prop_B_str]]$visited_count <- (visited_graphs[[prop_B_str]]$visited_count %||% 0) + 1

              mh_accept[[i]] <- append(mh_accept[[i]], 1)
              #mh_elbo_chain[[i]] <- append(mh_elbo_chain[[i]], proposal_graph_info$elbo)
          } else {
              #mh_chain[[i]] <- append(mh_chain[[i]], curr_B_str)
              mh_accept[[i]] <- append(mh_accept[[i]], 0)
              #mh_elbo_chain[[i]] <- append(mh_elbo_chain[[i]], visited_graphs[[curr_B_str]]$elbo)
              visited_graphs[[curr_B_str]]$visited_count <- (visited_graphs[[curr_B_str]]$visited_count %||% 0 ) + 1
          }
          print(sprintf("Accept/Reject ratio: %.2f", mean(unlist(mh_accept[[i]]))))

          # print(sprintf("Diff from true elbo: %.2f", true_mod$elbo - proposal_graph_info$elbo))

          iter <- iter + 1
          print(sprintf("Starting next iteration: %d", iter))
          print(sprintf("Number of nesmr fits: %d", nesmr_fits))
            if (!is.null(checkpoint_file) && nesmr_fits %% checkpoint_every == 0) {
                saveRDS(
                    list(
                        visited_graphs = visited_graphs,
                        mh_chain = mh_chain,
                        mh_accept = mh_accept,
                        mh_elbo_chain = mh_elbo_chain,
                        iter = iter,
                        nesmr_fits = nesmr_fits,
                        elbo_denom = elbo_denom
                    ),
                    file = checkpoint_file
                )
                print(sprintf("Checkpoint saved at iteration %d", iter))
            }
      }
  }

  rtn <- list(
      visited_graphs = visited_graphs,
      mh_chain = mh_chain,
      mh_accept = mh_accept,
      mh_accept_ratio = mean(unlist(mh_accept[[i]])),
      elbo_chain = mh_elbo_chain,
      iter = iter,
      nesmr_fits = nesmr_fits,
      elbo_denom = elbo_denom,
      mvmr_all = n_mvmr_res
  )
  class(rtn) <- "nesmr_mh_graph_explore"
  return(rtn)
}

draw_graph <- function(g, x) {
  add_candidates <- x$add_candidates
  add_candidate_prob <- x$add_candidate_prob
  total_add_prob <- x$total_add_prob
  total_remove_prob <- x$total_remove_prob
  remove_candidates <- x$remove_candidates
  remove_edge_prob <- x$remove_edge_prob

  total_prob <- total_add_prob + total_remove_prob

  insert_edge <- runif(1) < total_add_prob / total_prob
  if (insert_edge) {
      # Draw from the add candidates
      cond_prob <- add_candidate_prob / total_add_prob
      add_candidate_ix <- sample(seq_along(cond_prob), 1, prob = cond_prob)
      new_graph <- igraph::add_edges(g, add_candidates[add_candidate_ix, ])
      prob <- cond_prob[add_candidate_ix] / total_prob
      mod_edge <- paste0(add_candidates[add_candidate_ix, ], collapse = "|")
  } else {
      # Draw from the remove edges
      cond_prob <- remove_edge_prob / total_remove_prob
      remove_edge_ix <- sample(seq_along(remove_edge_prob), 1, prob = cond_prob)
      mod_edge <- paste0(remove_candidates[remove_edge_ix, ], collapse = "|")
      new_graph <- igraph::delete_edges(g, mod_edge)
      prob <- cond_prob[remove_edge_ix] / total_prob
  }
  return(
      list(
          g = new_graph, prob = prob, mod_edge = mod_edge, insert_edge = insert_edge))
}

# TODO: Import igraph
get_adjacent_graphs <- function(g, weight_mat) {
    g_comp <- igraph::complementer(g, loops = FALSE)

#    # TODO: Is it from here that we ned
    add_candidates <- igraph::as_edgelist(g_comp, names = FALSE)
    add_candidate_g <- apply(add_candidates, 1, function(x) {
        from <- x[1]
        to <- x[2]
        # Try adding the edge
        g_test <- igraph::add_edges(g, c(from, to))
        is_dag_test <- igraph::is_dag(g_test)
        if (is_dag_test) {
            return(g_test)
        } else {
            return(NULL)
        }
    })
    keep_graphs <- sapply(add_candidate_g, Negate(is.null))

    if (length(keep_graphs) == 0) {
        add_candidate_prob <- 0
        add_candidates <- 0
    } else {
        add_candidate_g <- add_candidate_g[keep_graphs]
        add_candidates <- add_candidates[keep_graphs,, drop = FALSE]
        add_candidate_prob <- weight_mat[add_candidates]
    }

    total_add_prob <- sum(add_candidate_prob)

    # Remove candidates : all edges
    remove_candidates <- igraph::as_edgelist(g, names = FALSE)
    remove_edge_prob <- 1 - weight_mat[remove_candidates]

    total_remove_prob <- sum(remove_edge_prob)

    total_prob <- total_add_prob + total_remove_prob
    # Normalize both the probabilities
    add_candidate_prob <- add_candidate_prob / total_prob
    remove_edge_prob <- remove_edge_prob / total_prob

    # First draw: Add w prob total_add_prob / (total_add_prob + total_remove_prob)
    # Second if add: Draw from one of the add candidates w weights in add_candidate_prob
    # Third if remove:
    #   - Draw from one of the remove candidates w weights in remove_edge_prob
    #   - Fit the graph that is removed
    list(
        add_candidates = add_candidates,
        add_candidate_prob = add_candidate_prob,
        total_add_prob = total_add_prob,
        total_remove_prob = total_remove_prob,
        remove_candidates = remove_candidates,
        remove_edge_prob = remove_edge_prob
    )
}

get_Z_to_prob <- function(x, y) {
    suppressWarnings(mod_coefs <- unname(glm(y ~ x, family = binomial)$coef))
    return(function(x) plogis(mod_coefs[1] + mod_coefs[2] * x))
}

summary.nesmr_mh_graph_explore <- discovery_summary <- function(x) {
    # Want:
    # Table of graphs with number of edges, elbo, norm_elbo, and visited count
    n <- length(x$visited_graphs)
    all_elbos <- sapply(x$visited_graphs, function(g) g$elbo)
    norm_elbo <- exp(all_elbos - x$elbo_denom)
    flat_graphs <- names(x$visited_graphs)
    num_edges <- sapply(flat_graphs, function(s) {
        sum(as.numeric(stringr::str_split(s, "")[[1]]))
    })
    visit_count <- sapply(x$visited_graphs, function(g) g$visited_count %||% 0)

    graph_summary <- data.frame(
        graph = flat_graphs,
        num_edges = num_edges,
        elbo = all_elbos,
        norm_elbo = norm_elbo,
        visit_count = visit_count
    )
    rownames(graph_summary) <- NULL
    graph_summary <- graph_summary[order(graph_summary$norm_elbo, decreasing = TRUE), ]

    return(graph_summary)
}