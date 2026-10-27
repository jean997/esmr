# ============================================================
# TOPO-style ordering search for lower-triangular SEMs
#
# Assumption:
# You already have a function:
#
#   fit <- nesmr(Y, order = perm, ...)
#
# returning at least:
#
#   fit$elbo        : scalar ELBO
#   fit$grad_upper  : K x K matrix
#
# where grad_upper[i,j] gives the gradient signal for a
# forbidden edge under the current ordering
# (i < j in ordering coordinates).
#
# If grad_upper is unavailable, the code falls back to
# adjacent swap search using ELBO only.
#
# ------------------------------------------------------------
# Main function:
#
# topo_search(...)
#
# returns:
#
#   best_order
#   best_fit
#   history
#
# ============================================================


# ------------------------------------------------------------
# Utility: adjacent swap
# ------------------------------------------------------------

swap_adjacent <- function(perm, i) {
  # swap positions i and i+1
  stopifnot(i < length(perm))

  new_perm <- perm
  tmp <- new_perm[i]
  new_perm[i] <- new_perm[i + 1]
  new_perm[i + 1] <- tmp

  new_perm
}


# ------------------------------------------------------------
# Utility: arbitrary swap
# ------------------------------------------------------------

swap_positions <- function(perm, i, j) {
  stopifnot(i != j)

  new_perm <- perm
  tmp <- new_perm[i]
  new_perm[i] <- new_perm[j]
  new_perm[j] <- tmp

  new_perm
}


# ------------------------------------------------------------
# Extract candidate swaps from gradient matrix
# ------------------------------------------------------------

get_candidate_swaps <- function(grad_upper,
                                top_n = 10,
                                adjacent_only = FALSE) {

  K <- nrow(grad_upper)

  candidates <- list()

  if (adjacent_only) {
    for (i in 1:(K - 1)) {
      score <- abs(grad_upper[i, i + 1])

      candidates[[length(candidates) + 1]] <- data.frame(
        i = i,
        j = i + 1,
        score = score
      )
    }
  } else {
    for (i in 1:(K - 1)) {
      for (j in (i + 1):K) {
        score <- abs(grad_upper[i, j])

        candidates[[length(candidates) + 1]] <- data.frame(
          i = i,
          j = j,
          score = score
        )
      }
    }
  }

  cand_df <- do.call(rbind, candidates)

  cand_df <- cand_df[order(-cand_df$score), , drop = FALSE]

  head(cand_df, top_n)
}


# ------------------------------------------------------------
# Evaluate one swap by re-fitting model
# ------------------------------------------------------------

# evaluate_swap <- function(current_perm,
#                           i,
#                           j,
#                           ...) {

#   fit <- nesmr(
#     order = new_perm,
#     ...
#   )

#   list(
#     order = new_perm,
#     fit = fit,
#     elbo = fit$elbo,
#     i = i,
#     j = j
#   )
# }


# ------------------------------------------------------------
# Main TOPO search
# ------------------------------------------------------------

topo_search <- function(
    beta_hat,
    beta_hat_se,
    initial_order = seq_len(ncol(beta_hat)),
    max_iter = 50,
    top_n_swaps = 10,
    adjacent_only = FALSE,
    min_improvement = 1e-6,
    verbose = TRUE,
    beta_prior_cov = 1,
    ...) {

  B_lower <- matrix(0, nrow = ncol(beta_hat), ncol = ncol(beta_hat))
  B_lower[lower.tri(B_lower)] <- 1

  current_order <- initial_order

  if (verbose) {
    cat("Initial fit...\n")
  }

  current_fit <- nesmr(
    beta_hat_X = beta_hat,
    se_X = beta_hat_se,
    direct_effect_template = B_lower[current_order, current_order],
    max_iter = 300,
    params = list(
      beta_prior_cov = beta_prior_cov
      )
    )

  current_fit$grad_upper <- nesmr_with_grad(current_fit)$grad_upper

  current_elbo <- current_fit$elbo

  history <- list()

  for (iter in 1:max_iter) {

    if (verbose) {
      cat(
        "\nIteration:", iter,
        "| Current ELBO:", round(current_elbo, 4),
        "\n"
      )
    }

    # --------------------------------------------------------
    # Candidate generation
    # --------------------------------------------------------

    if (!is.null(current_fit$grad_upper)) {

      candidates <- get_candidate_swaps(
        grad_upper = current_fit$grad_upper,
        top_n = top_n_swaps,
        adjacent_only = adjacent_only
      )

      if (verbose) {
        cat("Using gradient-based candidate swaps\n")
        print(candidates)
      }

    } else {

      if (verbose) {
        cat("No gradient supplied; using adjacent swaps\n")
      }

      candidates <- data.frame(
        i = 1:(length(current_order) - 1),
        j = 2:length(current_order),
        score = NA_real_
      )
    }

    # --------------------------------------------------------
    # Evaluate candidates (greedy: stop at first improvement)
    # --------------------------------------------------------

    best_candidate <- NULL
    best_candidate_elbo <- current_elbo
    best_swap_i <- NULL
    best_swap_j <- NULL
    improvement <- 0

    for (r in 1:nrow(candidates)) {

      i <- candidates$i[r]
      j <- candidates$j[r]

      if (verbose) {
        cat("Trying swap:", i, "<->", j, " (score:", round(candidates$score[r], 6), ")\n")
      }

      new_perm <- swap_positions(current_order, i, j)

      candidate <- nesmr(
        beta_hat_X = beta_hat,
        se_X = beta_hat_se,
        direct_effect_template = B_lower[new_perm, new_perm],
        max_iter = 300,
        params = list(
          beta_prior_cov = beta_prior_cov
        )
      )

      improvement <- candidate$elbo - current_elbo

      if (improvement >= min_improvement) {
        # Greedy: accept first improving swap
        best_candidate <- candidate
        best_candidate_elbo <- candidate$elbo
        best_swap_i <- i
        best_swap_j <- j

        if (verbose) {
          cat("Found improving swap! ELBO improvement:", round(improvement, 6), "\n")
        }

        break
      } else {
        if (verbose) {
          cat("No improvement (delta:", round(improvement, 6), ")\n")
        }
      }
    }

    # --------------------------------------------------------
    # Stopping criterion
    # --------------------------------------------------------

    if (is.null(best_candidate)) {

      if (verbose) {
        cat("Converged: no improving swap found after evaluating", nrow(candidates), "candidates.\n")
      }

      break
    }

    # --------------------------------------------------------
    # Accept move
    # --------------------------------------------------------

    new_order <- swap_positions(current_order, best_swap_i, best_swap_j)
    current_order <- new_order
    current_fit <- best_candidate
    current_elbo <- best_candidate_elbo
    current_fit$grad_upper <- nesmr_with_grad(current_fit)$grad_upper

    history[[iter]] <- list(
      iter = iter,
      order = current_order,
      elbo = current_elbo,
      swap = c(best_swap_i, best_swap_j)
    )

    if (verbose) {
      cat(
        "Accepted swap:",
        best_swap_i, "<->", best_swap_j,
        "| New ELBO:",
        round(current_elbo, 4),
        "\n"
      )
    }
  }

  list(
    best_order = current_order,
    best_fit = current_fit,
    history = history
  )
}
