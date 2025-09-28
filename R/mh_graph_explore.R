#' Metropolis-Hastings Graph Exploration for NESMR
#'
#' This function performs Metropolis-Hastings (MH) sampling to explore the space of possible directed acyclic graphs (DAGs) which are then fit using NESMR.
#' It starts from one or more initial graphs ("best", "sparse", "dense") and iteratively proposes new graphs by adding or removing edges, accepting or rejecting proposals based on the change in model evidence (ELBO) and a proposal probability.
#'
#' Edges are proposed based on a logistic transformation of the Z-scores. The shape of the logistic function is determined from the points: (0, min_prob) and (max_Z, max_prob).
#' A lower max_prob will flatten the curve and make the difference between higher and lower Z-scores less. A higher max_Z will shift the curve to the right and make it more flat.
#'
#' @param dat A list containing GWAS summary statistics which must include `beta_hat` and `s_estimate` (matrices of effect sizes and standard errors).
#' @param n_mvmr_res Optional. Precomputed NESMR MVMR results. If NULL, will be computed internally. This is computed from the nesmr_complete_mvmr and is beneficial pass if graph explorations are run multiple times.
#' @param pval_select Optional. Matrix of p-values for variant selection. If NULL, will be computed from `dat` as 2 * pnorm(-abs( dat$beta_hat / dat$s_estimate))
#' @param R Optional. LD matrix or NULL. See `?esmr::esmr` for more details.
#' @param alpha P-value threshold for variant selection. Default is 5e-8. Ignored if `pval_select` is provided.
#' @param max_Z Maximum Z-score for edge proposal probability mapping. Default is 5.
#' @param max_prob Maximum probability for edge proposal mapping. Default is 0.7.
#' @param min_prob Minimum probability for edge proposal mapping. Default is 0.01.
#' @param init_prob_threshold Initial probability threshold for edge inclusion based on Z-score p-values. Default is 0.05.
#' @param init_prob_method Method for determining initial edge inclusion threshold. Options are "pvalue", "fdr". If "pvalue", then directly use the p-values from Z-score if "fdr", then use FDR adjusted values.
#' @param mh_chain_init A named list of initial adjacency matrices to start MH chains from. If empty, will start one chain from the "best_approx" graph derived from the full MVMR results.
#' @param sparse_chain Logical. Whether to initialize a sparse chain (single best edge). Default is FALSE.
#' @param dense_chain Logical. Whether to initialize a dense chain (all possible edges). Default is FALSE.
#' @param max_iter Maximum number of MH iterations per chain. Default is 1000.
#' @param max_nesmr_fits Maximum number of NESMR model fits per chain. Default is 100.
#' @param visited_graphs List of previously visited graphs (for warm start or checkpointing). Default is empty list.
#' @param checkpoint_file Optional. File path to save checkpoints. Default is NULL.
#' @param checkpoint_every Integer. Save checkpoint every N NESMR fits. Default is 0 (no checkpointing).
#' @param verbose Logical. Print progress and debug information. Default is FALSE.
#' @param debug Logical. If TRUE, collects additional debug information during the MH sampling. Default is FALSE.
#'
#' @return An object of class `nesmr_mh_graph_explore`, a list containing:
#'   - visited_graphs: List of all visited graphs and their ELBOs
#'   - norm_elbo: Normalized ELBOs for all graphs
#'   - mh_chain: List of graph chains (one per initialization)
#'   - mh_accept: List of accept/reject indicators per chain
#'   - mh_accept_ratio: Mean acceptance ratio per chain
#'   - elbo_chain: List of ELBOs per chain
#'   - iter: Number of iterations per chain
#'   - nesmr_fits: Number of NESMR fits performed
#'   - elbo_denom: Log-sum-exp denominator for normalization
#'   - mvmr_all: The full MVMR NESMR result
#'
#' @examples
#'
#' library(GWASBrewer)
#' library(esmr)
#' # Generate a simple 4-node DAG and simulate using GWASBrewer
#' G <- matrix(
#'   c(0, 0, 0, 0,
#'     0.25, 0, 0, 0,
#'     0, 0, 0,  0,
#'     0, -0.15, 0.2, 0),
#'   nrow = 4,
#'   byrow = 4
#' )
#' d <- ncol(G)
#' h2 <- 0.2
#' J <- 5000
#' N <- 20000
#' pi_J <- 0.1
#' alpha <- 5e-8
#'
#' dat <- sim_mv(
#'     G = G,
#'     N = N,
#'     J = J,
#'     h2 = h2,
#'     pi = pi_J,
#'     sporadic_pleiotropy = TRUE,
#'     est_s = TRUE
#' )
#'
#' Z <- dat$beta_hat / dat$s_estimate
#' pval_select <- 2*pnorm(-abs(Z))
#' minp <- apply(pval_select, 1, min)
#' ix <- which(minp < alpha)
#' discovery_results <- mh_graph_explore(
#'   dat, pval_select = pval_select, verbose = TRUE)
#' print(discovery_results)
#'
#' @export
mh_graph_explore <- function(
  dat,
  n_mvmr_res = NULL,
  pval_select = NULL,
  R = NULL,
  alpha = 5e-8,
  max_Z = 5,
  max_prob = 0.7,
  min_prob = 0.01,
  init_prob_threshold = 0.05, # TODO: Change this to Z-threshold/p-value threshold
  init_prob_method = c("pvalue", "fdr"),
  mh_chain_init = list(), # Default is empty list which is one chain at "best_approx"
  sparse_chain = FALSE, # TODO: Remove this ?
  dense_chain = FALSE, # TODO: Remove this?
  max_iter = 1000,
  max_nesmr_fits = 100,
  visited_graphs = list(),
  checkpoint_file = NULL,
  checkpoint_every = 0,
  temperature = FALSE,
  burnin = round(max_iter / 10),
  max_heat = 5,
  verbose = FALSE,
  debug = FALSE
  ) {
  # Create a logging function based on verbose parameter
  log_msg <- function(...) {
    if (verbose) {
      message("mh_graph_explore: ", ...)
    }
  }

  init_prob_method <- match.arg(init_prob_method)

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
    capture.output({n_mvmr_res <- esmr::nesmr_complete_mvmr(
        beta_hat = dat$beta_hat,
        se_beta_hat = dat$s_estimate,
        pval_select = pval_select,
        R = R
    )}, file = nullfile()) # if (verbose) "" else nullfile())
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

  edge_prob_matrix <- Z_to_prob(abs(full_graph_zscores))

  # Edge filter matrix
  init_filter_matrix <- if (init_prob_method == "pvalue") {
      (2 * pnorm(-abs(full_graph_zscores)) < init_prob_threshold) * 1
  } else if (init_prob_method == "fdr") {
      pval_matrix <- 2 * pnorm(-abs(full_graph_zscores))
      # pval_matrix[diag(d)] <- NA
      pval_vector <- pval_matrix[non_diag_i]
      fdr_vector <- p.adjust(pval_vector, method = "fdr")
      fdr_matrix <- matrix(0, nrow = d, ncol = d)
      fdr_matrix[non_diag_i] <- fdr_vector
      (fdr_matrix < init_prob_threshold) * 1 - diag(d)
  } else {
      stop("Unknown init_prob_method")
  }

  mh_chain_init$best_approx = {
      mat_init <- matrix(0, nrow = d, ncol = d)
      mat_init[non_diag_i] <- full_graph_zscores[non_diag_i]
      init_filter_zscore <- mat_init * init_filter_matrix
      sqrt(esmr:::maximal_acyclic_subgraph((init_filter_zscore)^2)) * sign(init_filter_zscore)
    }

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
  elbo_denom <- -Inf
  if (debug) {
    mh_prop_denom <- list()
    mh_prop_num <- list()
    mh_prop_ratio <- list()
    mh_exp_elbo_diff <- list()
    mh_accept_prob <- list()
    mh_heat <- list()
    mh_insert_edge <- list()
    mh_mod_edge <- list()
    mh_prop_edge_Z_score <- list()
  }

  for (i in seq_along(mh_chain_init)) {
    if (debug) {
        mh_prop_denom[[i]] <- list()
        mh_prop_num[[i]] <- list()
        mh_prop_ratio[[i]] <- list()
        mh_exp_elbo_diff[[i]] <- list()
        mh_accept_prob[[i]] <- list()
        mh_heat[[i]] <- list()
        mh_insert_edge[[i]] <- list()
        mh_mod_edge[[i]] <- list()
        mh_prop_edge_Z_score[[i]] <- list()
    }
      # Note: This could be outside of the while or inside..
    curr_adj_mat <- mh_chain_init[[i]]
#       log_msg(sprintf("Initial graph: correct = %s", B_correct_str == paste0((curr_adj_mat != 0) + 0, collapse = "")))
      #log_msg("Initial adjacency matrix:")
      # if (verbose) print(curr_adj_mat)
      # Collect the graphs as flattened strings as we go
      curr_B <- (curr_adj_mat != 0) + 0
      curr_B_str <- paste0(curr_B, collapse = "")

      if (is.null(visited_graphs[[curr_B_str]])) {

        # Initial NESMR fit
#        capture.output({
        start_time <- Sys.time()
        capture.output({
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
        )}, file = nullfile()) # if (verbose) stdout() else nullfile())
        end_time <- Sys.time()
        init_fit_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
        log_msg(sprintf("Initial fit time: %.2f seconds", init_fit_time))
        log_msg(sprintf("Expect the total time to be around %.2f minutes", init_fit_time * max_nesmr_fits / 60))

#        }, file = nullfile())

        visited_graphs[[curr_B_str]]$elbo <- init_mod$elbo
        visited_graphs[[curr_B_str]]$beta_hat <- init_mod$direct_effects
        visited_graphs[[curr_B_str]]$se_beta_hat <- init_mod$se_dm
        visited_graphs[[curr_B_str]]$proposed <- (visited_graphs[[curr_B_str]]$proposed %||% 0) + 1

        elbo_denom <- matrixStats::logSumExp(c(elbo_denom, init_mod$elbo), na.rm = TRUE)
      }


      if (i == 1) elbo_denom <- visited_graphs[[curr_B_str]]$elbo

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
      heat_param_func <- approxfun(
        x = c(1, burnin),
        y = c(max_heat, 1),
      )

      #while (hit_old_graph <= no_new_graph_limit && iter < max_iter) {
      while(iter < max_iter && nesmr_fits < max_nesmr_fits) {
        log_msg("========================")
        log_msg(sprintf("Chain %s, Iteration %d of max %d", names(mh_chain_init)[i], iter, max_iter))
        log_msg(sprintf("Current number of unique graphs: %d", length(visited_graphs)))
        log_msg(sprintf("Current number of NESMR fits: %d of max %d", nesmr_fits, max_nesmr_fits))
          #log_msg("Current graph B matrix:")
          #if (verbose) print(curr_B)
          if (curr_B_str %in% names(visited_graphs) && !is.null(visited_graphs[[curr_B_str]]$adj_graph_info)) {
              adj_graph_info <- visited_graphs[[curr_B_str]]$adj_graph_info
          } else {
              visited_graphs[[curr_B_str]]$adj_graph_info <- esmr:::get_adjacent_graphs(ig, edge_prob_matrix)
              adj_graph_info <- visited_graphs[[curr_B_str]]$adj_graph_info
          }

          # log_msg("Adjacent graph information:")
          # if (verbose) print(adj_graph_info)
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
              add_edge_prob <- proposal_graph_info$adj_graph_info$add_edge_prob
              add_candidate_ix <- which(apply(add_candidates, 1, function(x) {
                  paste0(x, collapse = "|")
              }) == candidate_draw$mod_edge)
              # Numerator: h(G|G')
              prop_num <- add_edge_prob[add_candidate_ix]
          }

          if (is.null(proposal_graph_info$elbo)) {
              # If we have zero edges; continue
              # Eventually esmr should support having zero edges
              if (sum(prop_B) == 0) {
                  log_msg("Zero edges; continue")
                  mh_chain[[i]] <- append(mh_chain[[i]], curr_B_str)
                  mh_accept[[i]] <- append(mh_accept[[i]], 0)
                  mh_elbo_chain[[i]] <- append(mh_elbo_chain[[i]], visited_graphs[[curr_B_str]]$elbo)
                  iter <- iter + 1
                  next
              }

              # TODO: Probably want to just re-fit from initial/previous chain
              capture.output(
                  {
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
                  },
                  file = nullfile() #if (verbose) stdout() else nullfile()
              )

              visited_graphs[[prop_B_str]]$elbo <- new_mod$elbo
              visited_graphs[[prop_B_str]]$beta_hat <- new_mod$direct_effects
              visited_graphs[[prop_B_str]]$se_beta_hat <- new_mod$se_dm
              proposal_graph_info <- visited_graphs[[prop_B_str]]
          }

          visited_graphs[[prop_B_str]] <- proposal_graph_info
          visited_graphs[[prop_B_str]]$proposed <- (visited_graphs[[prop_B_str]]$proposed %||% 0) + 1

          elbo_diff <- proposal_graph_info$elbo - visited_graphs[[curr_B_str]]$elbo
          #log_msg(sprintf("\tELBO diff: %s", round(elbo_diff, 4)))
          #log_msg(sprintf("\tELBO curr: %s", round(visited_graphs[[curr_B_str]]$elbo, 4)))
          #log_msg(sprintf("\tELBO denom: %s", round(elbo_denom, 4)))
          # log_msg(sprintf("\texp(ELBO curr - ELBO denom): %s", round(exp(visited_graphs[[curr_B_str]]$elbo - elbo_denom), 4)))


          # TODO: Switch to log scale for everything
          prop_ratio <- prop_num / prop_denom

          log_msg(sprintf("\tProposal ratio: %s", round(prop_ratio, 4)))

          heat_param <- if (temperature && iter <= burnin) heat_param_func(iter) else 1
          # This should flatten the distribution so that more weight is on the proposal
          mh_ratio <- prop_ratio * exp(elbo_diff / heat_param)
          # Check accept/reject
          # max(1, exp(elbo(tmp_ig) - elbo(ig)))
          accept_prob <- min(
              1, mh_ratio
          )

          if (debug) {
            mh_prop_denom[[i]] <- append(mh_prop_denom[[i]], prop_denom)
            mh_prop_num[[i]] <- append(mh_prop_num[[i]], prop_num)
            mh_prop_ratio[[i]] <- append(mh_prop_ratio[[i]], prop_ratio)
            mh_exp_elbo_diff[[i]] <- append(mh_exp_elbo_diff[[i]], exp(elbo_diff))
            mh_accept_prob[[i]] <- append(mh_accept_prob[[i]], mh_ratio)
            mh_heat[[i]] <- append(mh_heat[[i]], heat_param)
            mh_insert_edge[[i]] <- append(mh_insert_edge[[i]], candidate_draw$insert_edge)
            mh_mod_edge[[i]] <- append(mh_mod_edge[[i]], candidate_draw$mod_edge)

            mh_prop_edge_Z_score[[i]] <- append(
                mh_prop_edge_Z_score[[i]],
                full_graph_zscores[
                    as.numeric(candidate_draw$from), as.numeric(candidate_draw$to)
                ])
          }

          log_msg(sprintf("\tTotal proposal ratio: %s (Heat parameter: %s)", round(mh_ratio, 4), round(heat_param, 4)))

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
          log_msg(sprintf("\tAccept/Reject ratio: %.2f", mean(unlist(mh_accept[[i]]))))

          iter <- iter + 1
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
                log_msg(sprintf("Checkpoint saved at iteration %d", iter))
            }
      }
  }

  all_elbos <- sapply(visited_graphs, function(g) g$elbo)
  norm_elbo <- exp(all_elbos - elbo_denom)

  rtn <- list(
      visited_graphs = visited_graphs,
      norm_elbo = norm_elbo,
      mh_chain = mh_chain,
      mh_accept = mh_accept,
      mh_accept_ratio = mean(unlist(mh_accept[[i]])),
      elbo_chain = mh_elbo_chain,
      iter = iter,
      nesmr_fits = nesmr_fits,
      elbo_denom = elbo_denom,
      mvmr_all = n_mvmr_res,
      mh_chain_init = mh_chain_init
  )

  if (debug) {
    rtn$mh_prop_denom <- mh_prop_denom
    rtn$mh_prop_num <- mh_prop_num
    rtn$mh_prop_ratio <- mh_prop_ratio
    rtn$mh_exp_elbo_diff <- mh_exp_elbo_diff
    rtn$mh_accept_prob <- mh_accept_prob
    rtn$mh_mod_edge <- mh_mod_edge
    rtn$mh_insert_edge <- mh_insert_edge
    rtn$mh_heat <- mh_heat
    rtn$mh_prop_edge_Z_score <- mh_prop_edge_Z_score
  }
  class(rtn) <- "nesmr_mh_graph_explore"
  return(rtn)
}

draw_graph <- function(g, x) {
  add_candidates <- x$add_candidates
  add_edge_prob <- x$add_edge_prob
  total_add_prob <- x$total_add_prob
  total_remove_prob <- x$total_remove_prob
  remove_candidates <- x$remove_candidates
  remove_edge_prob <- x$remove_edge_prob

  total_prob <- total_add_prob + total_remove_prob

  insert_edge <- runif(1) < total_add_prob / total_prob
  if (insert_edge) {
      # Draw from the add candidates
      cond_prob <- add_edge_prob / total_add_prob
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

  mod_edge_ix <- strsplit(mod_edge, "\\|")

  return(
      list(
          g = new_graph, prob = prob, mod_edge = mod_edge, insert_edge = insert_edge,
          from = mod_edge_ix[[1]][1], to = mod_edge_ix[[1]][2]
          )
)
}

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
        add_edge_prob <- numeric(0)
        add_candidates <- matrix(numeric(0), ncol = 2)
    } else {
        add_candidate_g <- add_candidate_g[keep_graphs]
        add_candidates <- add_candidates[keep_graphs,, drop = FALSE]
        add_edge_prob <- weight_mat[add_candidates]
    }

    total_add_prob <- sum(add_edge_prob)

    # Remove candidates : all edges
    remove_candidates <- igraph::as_edgelist(g, names = FALSE)
    remove_edge_prob <- 1 - weight_mat[remove_candidates]

    total_remove_prob <- sum(remove_edge_prob)

    total_prob <- total_add_prob + total_remove_prob
    # Normalize both the probabilities
    add_edge_prob <- add_edge_prob / total_prob
    remove_edge_prob <- remove_edge_prob / total_prob

    # First draw: Add w prob total_add_prob / (total_add_prob + total_remove_prob)
    # Second if add: Draw from one of the add candidates w weights in add_edge_prob
    # Third if remove:
    #   - Draw from one of the remove candidates w weights in remove_edge_prob
    #   - Fit the graph that is removed
    list(
        add_candidates = add_candidates,
        add_edge_prob = add_edge_prob,
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