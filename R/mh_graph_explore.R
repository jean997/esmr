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

get_adjacent_graphs <- function(g, weight_mat) {
    g_comp <- complementer(g, loops = FALSE)

    add_candidates <- as_edgelist(g_comp, names = FALSE)
    add_candidate_g <- apply(add_candidates, 1, function(x) {
        from <- x[1]
        to <- x[2]
        # Try adding the edge
        g_test <- add_edges(g, c(from, to))
        is_dag_test <- igraph::is_dag(g_test)
        if (is_dag_test) {
            return(g_test)
        } else {
            return(NULL)
        }
    })
    keep_graphs <- sapply(add_candidate_g, Negate(is.null))
    add_candidate_g <- add_candidate_g[keep_graphs]
    add_candidates <- add_candidates[keep_graphs, ]

    add_candidate_prob <- weight_mat[add_candidates]

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