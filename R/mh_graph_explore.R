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
#' @param init_prob_threshold Initial probability threshold for edge inclusion based on Z-score p-values. Default is 0.05.
#' @param init_prob_method Method for determining initial edge inclusion threshold. Options are "pvalue", "fdr". If "pvalue", then directly use the p-values from Z-score if "fdr", then use FDR adjusted values.
#' @param logistic_location Location parameter for logistic scaling. Default is 5 which centers the logistic function at the midpoint of the logistic function such that Z-score of 5 maps to 0.5.
#' @param logistic_scale_range A numeric vector of length 2 indicating the range of logistic scaling. Default c(3, 0.0001) which linearly decrease the scale from 3 to 0.0001 over [1, max_iter]. This progressively puts more weight on the ELBO difference and more on the initial Z-score MVMR estimates as the chain progresses.
#' @param mh_chain_init A named list of initial adjacency matrices to start MH chains from. If empty, will start one chain from the "best_approx" graph derived from the full MVMR results.
#' @param mh_chain_params A list of additional parameters for the MH chains. If the chain has an entry in this list, it will use the parameters (logistic_location, logistic_scale) from this list instead of the global parameters. Some of the names of chains provided are "best_approx", "min_graph", "max_graph", "random_start_1", "random_start_2", etc.
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
#'     c(
#'         0, 0, 0, 0,
#'         0.25, 0, 0, 0,
#'         0, 0, 0, 0,
#'         0, -0.15, 0.2, 0
#'     ),
#'     nrow = 4,
#'     byrow = 4
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
#' pval_select <- 2 * pnorm(-abs(Z))
#' minp <- apply(pval_select, 1, min)
#' ix <- which(minp < alpha)
#' discovery_results <- mh_graph_explore(
#'     dat,
#'     pval_select = pval_select, verbose = TRUE
#' )
#' print(discovery_results)
#'
#' @export
mh_graph_explore <- function(
    dat,
    n_mvmr_res = NULL,
    pval_select = NULL,
    R = NULL,
    alpha = 5e-8,
    init_prob_threshold = 0.05, # TODO: Change this to Z-threshold/p-value threshold
    init_prob_method = c("pvalue", "fdr"),
    logistic_location = 4,
    logistic_scale = 1, # c(0.25, 3),
    mh_chain_init = list(), # Default is empty list which is one chain at "best_approx"
    mh_chain_params = list(),
    sparse_chain = FALSE, # TODO: Remove this ?
    dense_chain = FALSE, # TODO: Remove this?
    random_starts = 0,
    max_iter = 1000,
    max_nesmr_fits = 100,
    visited_graphs = list(),
    checkpoint_file = NULL,
    checkpoint_every = 0,
    temperature = FALSE,
    burnin = round(max_iter / 10),
    max_heat = 5,
    verbose = FALSE,
    debug = FALSE) {
    # Create a logging function based on verbose parameter
    log_msg <- function(...) {
        if (verbose) {
            message("mh_graph_explore: ", ...)
        }
    }

    if (!length(logistic_scale) %in% c(1, 2)) {
        stop("logistic_scale must be a numeric vector of length 1 or 2")
    } else if (length(logistic_scale) == 1) {
        logistic_scale <- rep(logistic_scale, 2)
    }
    logistic_scale_range <- logistic_scale

    if (!length(logistic_location) %in% c(1, 2)) {
        stop("logistic_location must be a numeric vector of length 1 or 2")
    } else if (length(logistic_location) == 1) {
        logistic_location <- rep(logistic_location, 2)
    }
    logistic_location_range <- logistic_location

    if (diff(logistic_scale_range) < 0) {
        warning("logistic_scale_range should be increasing. Consider reversing the order.")
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
        capture.output(
            {
                n_mvmr_res <- esmr::nesmr_complete_mvmr(
                    beta_hat = dat$beta_hat,
                    se_beta_hat = dat$s_estimate,
                    pval_select = pval_select,
                    R = R
                )
            },
            file = nullfile()
        )
    }

    full_graph_zscores <- n_mvmr_res$beta_hat / n_mvmr_res$se_beta_hat
    stopifnot(ncol(n_mvmr_res$beta_hat) == d)
    non_diag_i <- -seq(1, d^2, by = d + 1)
    diag(full_graph_zscores) <- 0

    MVMR_abs_Z_scores <- abs(full_graph_zscores)
    # edge_prob_matrix <- Z_to_prob(abs(full_graph_zscores))

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

    mh_chain_init$best_approx <- {
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

    if (random_starts > 0) {
        # Note: Question here about unique or not.
        # If we are fitting with a single set of parameters, only makes sense to fit unique graphs
        # If we are fitting with different parameters, then makes sense to fit non-unique graphs
        random_start_graphs <- unique(map(seq_len(random_starts), function(i) {
            noisy_zscores <- matrix(0, nrow = K, ncol = K)
            noisy_zscores[non_diag_i] <- map(full_graph_zscores[non_diag_i], ~ rnorm(1, mean = .x, sd = 1)) %>% unlist()
            diag(noisy_zscores) <- 0
            initial_filter <- (abs(noisy_zscores) > qnorm(init_prob_threshold / 2, lower.tail = FALSE)) + 0
            noisy_zscores <- noisy_zscores * initial_filter
            (maximal_acyclic_subgraph(noisy_zscores^2) != 0) + 0
            }))
        mh_chain_init <- append(mh_chain_init, random_start_graphs %>% setNames(paste0("random_start_", seq_along(.))))
    }

    mh_chain_info <- lapply(seq_along(mh_chain_init), function(x) vector("list", length = max_iter))
    elbo_denom <- -Inf

    for (i in seq_along(mh_chain_init)) {
        accept_count <- 0
        curr_adj_mat <- mh_chain_init[[i]]
        curr_B <- (curr_adj_mat != 0) + 0
        curr_B_str <- paste0(curr_B, collapse = "")
        chain_name <- names(mh_chain_init)[i]

        if (is.null(visited_graphs[[curr_B_str]])) {
            # Initial NESMR fit
            #        capture.output({
            start_time <- Sys.time()
            capture.output(
                {
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
                },
                file = nullfile()
            )
            end_time <- Sys.time()
            init_fit_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
            log_msg(sprintf("Initial fit time: %.2f seconds", init_fit_time))
            log_msg(sprintf("Expect the total time to be around %.2f minutes", init_fit_time * max_nesmr_fits / 60))

            visited_graphs[[curr_B_str]]$elbo <- init_mod$elbo
            visited_graphs[[curr_B_str]]$beta_hat <- init_mod$direct_effects
            visited_graphs[[curr_B_str]]$se_beta_hat <- init_mod$se_dm
            visited_graphs[[curr_B_str]]$proposed <- (visited_graphs[[curr_B_str]]$proposed %||% 0) + 1

            elbo_denom <- matrixStats::logSumExp(c(elbo_denom, init_mod$elbo), na.rm = TRUE)
        }

        if (i == 1) elbo_denom <- visited_graphs[[curr_B_str]]$elbo

        ig <- igraph::graph_from_adjacency_matrix(
            curr_adj_mat != 0,
            mode = "directed"
        )

        iter <- 1
        nesmr_fits <- 1
        heat_param_func <- approxfun(
            x = c(1, burnin),
            y = c(max_heat, 1),
            rule = 2
        )


        if (chain_name %in% names(mh_chain_params) && !is.null(mh_chain_params[[chain_name]]$logistic_scale)) {
            logistic_scale_func <- function(i) mh_chain_params[[chain_name]]$logistic_scale
        } else {
            logistic_scale_func <- approxfun(
                x = c(1, max_iter),
                y = logistic_scale_range,
                method = "linear",
                rule = 2
            )
        }

        if (chain_name %in% names(mh_chain_params) && !is.null(mh_chain_params[[chain_name]]$logistic_location)) {
            logistic_location_func <- function(i) mh_chain_params[[chain_name]]$logistic_location
        } else {
            logistic_location_func <- approxfun(
                x = c(1, max_iter),
                y = logistic_location_range,
                method = "linear",
                rule = 2
            )
        }

        while (iter <= max_iter && nesmr_fits <= max_nesmr_fits) {
            logistic_scale <- logistic_scale_func(iter)
            logistic_location <- logistic_location_func(iter)
            log_msg("========================")
            log_msg(sprintf("Chain %s, Iteration %d of max %d", names(mh_chain_init)[i], iter, max_iter))
            log_msg(sprintf("Current number of unique graphs: %d", length(visited_graphs)))
            log_msg(sprintf("Current number of NESMR fits: %d of max %d", nesmr_fits, max_nesmr_fits))
            # log_msg("Current graph B matrix:")
            # if (verbose) print(curr_B)
            if (curr_B_str %in% names(visited_graphs) && !is.null(visited_graphs[[curr_B_str]]$adj_graph_info)) {
                adj_graph_info <- visited_graphs[[curr_B_str]]$adj_graph_info
            } else {
                visited_graphs[[curr_B_str]]$adj_graph_info <- esmr:::get_adjacent_graphs(
                    ig, MVMR_abs_Z_scores,
                    logistic_scale = logistic_scale
                )
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
                proposal_graph_info$adj_graph_info <- esmr:::get_adjacent_graphs(
                    tmp_ig, MVMR_abs_Z_scores,
                    logistic_location = logistic_location, # Use logistic_location from previous graph as need same proposal dist
                    logistic_scale = logistic_scale
                )
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
                                # Record the chain info
                    mh_chain_info[[i]][[iter]] <- list(
                        chain = i,
                        iter = iter,
                        curr_graph = curr_B_str,
                        prop_graph = prop_B_str,
                        accepted = 0,
                        curr_elbo = visited_graphs[[curr_B_str]]$elbo,
                        prop_elbo = -Inf,
                        elbo_diff = -Inf,
                        prop_num = 0,
                        prop_denom = 0,
                        accept_prob = 0,
                        mod_edge = candidate_draw$mod_edge,
                        insert_edge = candidate_draw$insert_edge,
                        heat_param = heat_param,
                        logistic_scale = logistic_scale,
                        logistic_location = logistic_location
                    )
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
                    file = nullfile() # if (verbose) stdout() else nullfile()
                )

                visited_graphs[[prop_B_str]]$elbo <- new_mod$elbo
                visited_graphs[[prop_B_str]]$beta_hat <- new_mod$direct_effects
                visited_graphs[[prop_B_str]]$se_beta_hat <- new_mod$se_dm
                proposal_graph_info <- visited_graphs[[prop_B_str]]
            }

            visited_graphs[[prop_B_str]] <- proposal_graph_info

            visited_graphs[[prop_B_str]]$adj_graph_info <- NULL
            visited_graphs[[prop_B_str]]$proposed <- (visited_graphs[[prop_B_str]]$proposed %||% 0) + 1

            elbo_diff <- proposal_graph_info$elbo - visited_graphs[[curr_B_str]]$elbo

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

            log_msg(sprintf(
                "\tTotal proposal ratio: %s (%sLogistic scale: %s Logistic location: %s)",
                round(mh_ratio, 4),
                if (temperature) paste0("Heat parameter: ", round(heat_param, 4), " ") else "",
                round(logistic_scale, 4),
                round(logistic_location, 4))
            )

            # Record all info for this step
            accepted <- as.integer(accept_prob == 1 || runif(1) < accept_prob)

            # Record the chain info
            mh_chain_info[[i]][[iter]] <- list(
                chain = i,
                iter = iter,
                curr_graph = curr_B_str,
                prop_graph = prop_B_str,
                accepted = accepted,
                curr_elbo = visited_graphs[[curr_B_str]]$elbo,
                prop_elbo = proposal_graph_info$elbo,
                elbo_diff = elbo_diff,
                prop_num = prop_num,
                prop_denom = prop_denom,
                accept_prob = accept_prob,
                mod_edge = candidate_draw$mod_edge,
                insert_edge = candidate_draw$insert_edge,
                heat_param = heat_param,
                logistic_scale = logistic_scale,
                logistic_location = logistic_location
            )


            if (accepted) {
                curr_B <- prop_B
                curr_B_str <- prop_B_str
                ig <- tmp_ig
                visited_graphs[[prop_B_str]]$visited_count <- (visited_graphs[[prop_B_str]]$visited_count %||% 0) + 1
            } else {
                visited_graphs[[curr_B_str]]$visited_count <- (visited_graphs[[curr_B_str]]$visited_count %||% 0) + 1
            }

            accept_count <- accept_count + accepted
            accept_ratio <- accept_count / iter

            log_msg(sprintf("\tAccept/Reject ratio: %.2f", accept_ratio))
            iter <- iter + 1
            if (!is.null(checkpoint_file) && nesmr_fits %% checkpoint_every == 0) {
                saveRDS(
                    # TODO: Need to update this for the resume mh_graph_explore
                    list(
                        visited_graphs = visited_graphs,
                        mh_chain_info = dplyr::bind_rows(mh_chain_info),
                        n_nesmr_fits = nesmr_fits,
                        mvmr_all = n_mvmr_res,
                        mh_chain_init = mh_chain_init,
                        temperature = temperature,
                        burnin = burnin,
                        max_heat = max_heat,
                        logistic_location = logistic_location_range,
                        logistic_scale = logistic_scale_range
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
        mh_chain_info = dplyr::bind_rows(mh_chain_info),
        n_nesmr_fits = nesmr_fits,
        mvmr_all = n_mvmr_res,
        mh_chain_init = mh_chain_init,
        temperature = temperature,
        burnin = burnin,
        max_heat = max_heat,
        logistic_location = logistic_location_range,
        logistic_scale = logistic_scale_range
    )

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

# If missing location, then use the maximum value in the current adjacency matrix
# TODO: May be an issue for the proposal distribution for the inverse?
get_adjacent_graphs <- function(
    g, weight_mat, logistic_scale = 1,
    logistic_location = 5) {
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
        add_edge_weight <- add_edge_prob <- numeric(0)
        add_candidates <- matrix(numeric(0), ncol = 2)
    } else {
        add_candidate_g <- add_candidate_g[keep_graphs]
        add_candidates <- add_candidates[keep_graphs, , drop = FALSE]
        add_edge_weight <- weight_mat[add_candidates]
        # add_edge_prob <- weight_mat[add_candidates]
    }

    # Remove candidates : all edges
    remove_candidates <- igraph::as_edgelist(g, names = FALSE)
    remove_edge_weight <- weight_mat[remove_candidates]

    # if (is.null(location)) {
    #     if (length(add_edge_weight) > 0) {
    #         location <- max(add_edge_weight, na.rm = TRUE)
    #         if (is.na(location) || is.infinite(location)) {
    #             location <- max(weight_mat, na.rm = TRUE)
    #         }
    #     } else if (length(add_edge_weight) == 0) {
    #         # Use the global max if no add candidates
    #         location <- max(weight_mat, na.rm = TRUE)
    #     }
    # }
    add_edge_prob <- plogis(weight_mat[add_candidates], location = logistic_location, scale = logistic_scale)
    remove_edge_prob <- 1 - plogis(weight_mat[remove_candidates], location = logistic_location, scale = logistic_scale)

    total_add_prob <- sum(add_edge_prob)
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
    # TODO: Put this into a single dataframe/matrix instead?
    list(
        add_candidates = add_candidates,
        add_edge_prob = add_edge_prob,
        total_add_prob = total_add_prob,
        total_remove_prob = total_remove_prob,
        remove_candidates = remove_candidates,
        remove_edge_prob = remove_edge_prob,
        logistic_scale = logistic_scale,
        logistic_location = logistic_location
    )
}

get_Z_to_prob <- function(x, y) {
    suppressWarnings(mod_coefs <- unname(glm(y ~ x, family = binomial)$coef))
    return(function(x) plogis(mod_coefs[1] + mod_coefs[2] * x))
}
