#' Multi-start MH graph exploration from randomized MVMR graphs
#'
#' Random MVMR draws -> maximal acyclic subgraph give many different starting
#' topologies. This function fits those starts, keeps the good and distinct ones,
#' and runs one [mh_graph_explore()] chain from each, all collected in a single
#' `nesmr_mh_graph_explore` object.
#'
#' Procedure:
#' 1. **Sample starts.** For every p-value threshold, draw `z + N(0, 1)` from the MVMR
#'    z-scores, keep edges with `|z|` above the threshold, and take the maximal acyclic
#'    subgraph. Draw until `n_graphs` graphs that were not seen before are found (never
#'    a duplicate template, also across thresholds) or `max_tries` draws are used. The
#'    deterministic (noise free) MVMR graph of each threshold is added as well.
#' 2. **Fit every start once** with [nesmr()] (ELBO includes the graph prior).
#' 3. **Choose chain starts.** Keep the starts within `start_elbo_window` ELBO of the
#'    best start. Group the kept graphs: two graphs within `merge_distance` edits
#'    (add / remove / flip each count as 1) are in the same group, and groups are the
#'    connected components of that relation. Only the highest-ELBO graph of each group
#'    is kept.
#' 4. **Run chains in order of ELBO.** One [mh_graph_explore()] call per start, each
#'    with `prev_run` set to the previous call. A start is skipped if its graph was
#'    already proposed or visited by an earlier chain. Each chain stops after
#'    `early_stop_iter` newly fitted models without beating the overall best ELBO
#'    (the best over the starts and all earlier chains).
#'
#' Chains are named `start_<rank>_p<threshold>` in the `chain_name` column of
#' `mh_chain_info`.
#'
#' @inheritParams mh_graph_explore
#' @param thresholds P-value thresholds for an MVMR edge to be included in a start graph.
#' @param n_graphs Number of unique random start graphs per threshold (the noise free
#'   MVMR graph comes on top of this).
#' @param max_tries Maximum number of draws per threshold.
#' @param start_elbo_window Only starts with ELBO within this distance of the best start
#'   are considered for chains.
#' @param merge_distance Starts within this many edits of each other are merged (the best
#'   ELBO represents the group). Use 0 to disable.
#' @param early_stop_iter Per-chain early stop: number of newly fitted models without a new
#'   overall best ELBO.
#' @param max_iter,max_nesmr_fits Per-chain limits passed to [mh_graph_explore()].
#' @param max_chains Run at most this many chains (default all).
#' @param prev_run Optional `nesmr_mh_graph_explore` (or RDS path). Its fitted graphs are
#'   reused, and its chains count as "seen".
#' @param graph_edge_prior Edge prior used by [log_graph_prior()].
#' @param nesmr_max_iter `max_iter` for the [nesmr()] fits of the start graphs
#'   ([mh_graph_explore()] uses 300).
#' @param cores Cores for fitting the start graphs (unix only). The chains themselves run
#'   sequentially because each depends on the graphs seen by the earlier ones.
#' @param checkpoint_file Optional path. The result so far is saved here after the
#'   start fits and after each chain.
#' @param ... Further arguments passed to [mh_graph_explore()] (e.g. `logistic_location`,
#'   `move_weight`, `graph_draw_mode`).
#'
#' @return A `nesmr_mh_graph_explore` object with an extra element `multistart`:
#'   - `starts`: one row per fitted start (`graph`, `threshold`, `source` = "mvmr" or
#'     "draw", `n_edges`, `elbo`, `in_window`, `group`, `kept`, `rank`)
#'   - `chains`: one row per kept start with `status` ("ran", "skipped_seen"),
#'     `chain_name`, `stop_reason`, `n_nesmr_fits`
#'   - `tries`: draws used per threshold, and the arguments used.
#' @export
mh_graph_explore_multistart <- function(
    dat,
    n_mvmr_res = NULL,
    pval_select = NULL,
    R = NULL,
    alpha = 5e-8,
    thresholds = c(5e-1, 5e-3, 5e-5, 5e-7),
    n_graphs = 50,
    max_tries = 1000,
    start_elbo_window = 20,
    merge_distance = 1,
    early_stop_iter = 100,
    max_iter = 5000,
    max_nesmr_fits = 1000,
    max_chains = Inf,
    prev_run = NULL,
    nesmr_params = list(beta_prior_cov = 1, strict_mode = FALSE),
    graph_edge_prior = 0.5,
    nesmr_max_iter = 300,
    cores = 1L,
    checkpoint_file = NULL,
    verbose = FALSE,
    ...) {
    log_msg <- function(...) if (verbose) message("mh_graph_explore_multistart: ", ...)
    d <- ncol(dat$beta_hat)
    if (is.character(prev_run) && length(prev_run) == 1) prev_run <- readRDS(prev_run)

    if (is.null(pval_select)) pval_select <- 2 * pnorm(-abs(dat$beta_hat / dat$s_estimate))
    ix <- which(apply(pval_select, 1, min) < alpha)
    if (is.null(n_mvmr_res)) {
        capture.output(
            {
                n_mvmr_res <- nesmr_complete_mvmr(
                    beta_hat = dat$beta_hat, se_beta_hat = dat$s_estimate,
                    pval_select = pval_select, R = R
                )
            },
            file = nullfile()
        )
    }
    z <- n_mvmr_res$beta_hat / n_mvmr_res$se_beta_hat
    z[!is.finite(z)] <- 0
    diag(z) <- 0

    # ---- 1. sample unique starts ----
    samp <- .sample_start_graphs(z, thresholds, n_graphs, max_tries)
    starts <- samp$starts
    log_msg(sprintf("%d unique start graphs (draws per threshold: %s)", nrow(starts), paste(samp$tries, collapse = ", ")))
    if (nrow(starts) == 0) stop("No start graphs were generated; try larger thresholds.")

    # ---- 2. fit every start once (reuse graphs already fitted in prev_run) ----
    old_visited <- if (!is.null(prev_run)) prev_run$visited_graphs else list()
    todo <- which(!(starts$graph %in% names(old_visited)))
    log_msg(sprintf("Fitting %d start graphs (%d already fitted in prev_run)", length(todo), nrow(starts) - length(todo)))
    fit_one <- function(g) {
        M <- flat_string_to_adj_mat(g)
        tryCatch(.fit_graph(M, dat, ix, R, nesmr_params, nesmr_max_iter, graph_edge_prior), error = function(e) e)
    }
    fits <- .lapply_cores(starts$graph[todo], fit_one, cores = cores)
    names(fits) <- starts$graph[todo]
    failed <- vapply(fits, inherits, NA, "error")
    if (any(failed)) warning(sprintf("%d start graphs failed to fit and are dropped", sum(failed)))
    fits <- fits[!failed]

    visited <- old_visited
    for (g in names(fits)) {
        visited[[g]] <- fits[[g]][c("elbo_without_prior", "elbo", "beta_hat", "se_beta_hat")]
    }
    starts <- starts[starts$graph %in% names(visited), , drop = FALSE]
    starts$elbo <- vapply(starts$graph, function(g) visited[[g]]$elbo, numeric(1))

    # best model overall (mh_graph_explore needs the model object in `best_graph`)
    prev_best_elbo <- if (!is.null(prev_run) && length(old_visited)) max(vapply(old_visited, `[[`, 1, "elbo")) else -Inf
    best_graph <- if (!is.null(prev_run)) prev_run$best_graph else NULL
    fit_elbo <- vapply(fits, `[[`, 1, "elbo")
    if (length(fits) && max(fit_elbo) > prev_best_elbo) {
        best_graph <- fits[[which.max(fit_elbo)]]$model
    }
    state <- list(
        visited_graphs = visited,
        mh_chain_info = if (!is.null(prev_run)) prev_run$mh_chain_info else NULL,
        mvmr_all = n_mvmr_res, best_graph = best_graph,
        chain_status = if (!is.null(prev_run)) prev_run$chain_status else NULL
    )
    class(state) <- "nesmr_mh_graph_explore"

    # ---- 3. choose chain starts ----
    starts <- .select_chain_starts(starts, start_elbo_window, merge_distance)
    chain_starts <- starts[starts$kept, , drop = FALSE]
    chain_starts <- chain_starts[order(-chain_starts$elbo), , drop = FALSE]
    chain_starts$rank <- seq_len(nrow(chain_starts))
    starts$rank <- chain_starts$rank[match(starts$graph, chain_starts$graph)]
    log_msg(sprintf(
        "%d starts fitted; %d within %g ELBO of the best; %d chain starts after merging graphs within %d edit(s)",
        nrow(starts), sum(starts$in_window), start_elbo_window, nrow(chain_starts), merge_distance
    ))

    result_so_far <- function(chains) {
        out <- state
        out$multistart <- list(
            starts = starts, chains = chains,
            tries = list(
                draws = samp$tries,
                args = list(
                    thresholds = thresholds, n_graphs = n_graphs, max_tries = max_tries,
                    start_elbo_window = start_elbo_window, merge_distance = merge_distance,
                    early_stop_iter = early_stop_iter
                )
            )
        )
        out
    }
    chains <- data.frame(
        rank = chain_starts$rank, graph = chain_starts$graph, threshold = chain_starts$threshold,
        start_elbo = chain_starts$elbo, status = NA_character_, chain_name = NA_character_,
        stop_reason = NA_character_, n_nesmr_fits = NA_integer_, stringsAsFactors = FALSE
    )
    if (!is.null(checkpoint_file)) saveRDS(result_so_far(chains), checkpoint_file)

    # ---- 4. one chain per start ----
    seen <- .seen_graphs(state$mh_chain_info)
    n_run <- 0
    for (r in seq_len(nrow(chains))) {
        g <- chains$graph[r]
        if (g %in% seen) {
            chains$status[r] <- "skipped_seen"
            log_msg(sprintf("Start %d (ELBO %.1f): already seen by an earlier chain; skipped", r, chains$start_elbo[r]))
            next
        }
        if (n_run >= max_chains) {
            chains$status[r] <- "not_run_max_chains"
            next
        }
        nm <- sprintf("start_%02d_p%s", chains$rank[r], format(chains$threshold[r], scientific = TRUE, digits = 2))
        log_msg(sprintf("Start %d of %d (ELBO %.1f, %s)", r, nrow(chains), chains$start_elbo[r], nm))
        res <- mh_graph_explore(
            dat = dat, n_mvmr_res = n_mvmr_res, pval_select = pval_select, R = R, alpha = alpha,
            max_iter = max_iter, max_nesmr_fits = max_nesmr_fits,
            kill_no_improve_models = early_stop_iter,
            prev_run = state, starting_graph = flat_string_to_adj_mat(g), chain_name = nm,
            nesmr_params = nesmr_params, graph_edge_prior = graph_edge_prior,
            verbose = FALSE, ...
        )
        n_run <- n_run + 1
        cs <- res$chain_status[nrow(res$chain_status), ]
        chains$status[r] <- "ran"
        chains$chain_name[r] <- cs$chain_name
        chains$stop_reason[r] <- cs$stop_reason
        chains$n_nesmr_fits[r] <- cs$n_nesmr_fits
        info <- res$mh_chain_info
        seen <- union(seen, .seen_graphs(info[info$chain_name == cs$chain_name, , drop = FALSE]))
        res$chain_status <- rbind(state$chain_status, res$chain_status)
        state <- res
        log_msg(sprintf(
            "  chain stopped (%s) after %d fits; best ELBO so far %.1f",
            cs$stop_reason, cs$n_nesmr_fits, max(vapply(res$visited_graphs, `[[`, 1, "elbo"))
        ))
        if (!is.null(checkpoint_file)) saveRDS(result_so_far(chains), checkpoint_file)
    }
    result_so_far(chains)
}

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

.graph_key <- function(M) paste0((M != 0) + 0, collapse = "")

# Unique start graphs per threshold. Returns list(starts = data.frame, tries = named int).
.sample_start_graphs <- function(z, thresholds, n_graphs, max_tries) {
    d <- ncol(z)
    mas <- function(zz, thr) {
        zz <- zz * (abs(zz) > qnorm(thr / 2, lower.tail = FALSE))
        if (all(zz == 0)) {
            return(NULL)
        }
        (maximal_acyclic_subgraph(zz^2) != 0) + 0
    }
    seen <- character(0)
    rows <- list()
    tries <- integer(length(thresholds))
    names(tries) <- format(thresholds, scientific = TRUE, digits = 2)
    add <- function(M, thr, source) {
        if (is.null(M) || sum(M) == 0) {
            return(FALSE)
        }
        k <- .graph_key(M)
        if (k %in% seen) {
            return(FALSE)
        }
        seen <<- c(seen, k)
        rows[[length(rows) + 1]] <<- data.frame(
            graph = k, threshold = thr, source = source, n_edges = sum(M), stringsAsFactors = FALSE
        )
        TRUE
    }
    for (t in seq_along(thresholds)) {
        thr <- thresholds[t]
        add(mas(z, thr), thr, "mvmr")
        n_new <- 0
        while (n_new < n_graphs && tries[t] < max_tries) {
            zn <- z + matrix(rnorm(d * d), d, d)
            diag(zn) <- 0
            tries[t] <- tries[t] + 1L
            if (add(mas(zn, thr), thr, "draw")) n_new <- n_new + 1
        }
        if (n_new < n_graphs) {
            warning(sprintf("threshold %g: only %d of %d unique graphs after %d draws", thr, n_new, n_graphs, tries[t]))
        }
    }
    list(starts = do.call(rbind, rows), tries = tries)
}

# nesmr fit of a template; mirrors the fits done inside mh_graph_explore()
.fit_graph <- function(M, dat, ix, R, nesmr_params, nesmr_max_iter, graph_edge_prior) {
    capture.output(
        {
            mod <- nesmr(
                beta_hat_X = dat$beta_hat, se_X = dat$s_estimate, variant_ix = ix,
                direct_effect_template = M, max_iter = nesmr_max_iter, R = R, params = nesmr_params
            )
        },
        file = nullfile()
    )
    list(
        elbo_without_prior = mod$elbo,
        elbo = mod$elbo + log_graph_prior(sum(M), ncol(M), pi_0 = graph_edge_prior),
        beta_hat = mod$direct_effects, se_beta_hat = mod$direct_effects_se, model = mod
    )
}

# Number of edits between graphs: one per unordered pair whose state
# (none / i->j / j->i / both) differs, so a flip counts 1.
.edit_distance_matrix <- function(keys) {
    X <- t(vapply(keys, function(k) {
        M <- flat_string_to_adj_mat(k)
        d <- nrow(M)
        up <- which(upper.tri(M), arr.ind = TRUE)
        M[up] + 2 * M[up[, 2:1, drop = FALSE]]
    }, numeric(choose(sqrt(nchar(keys[1])), 2))))
    X <- matrix(X, nrow = length(keys))
    n <- length(keys)
    D <- matrix(0, n, n)
    for (a in seq_len(n)) D[a, ] <- rowSums(X != matrix(X[a, ], n, ncol(X), byrow = TRUE))
    D
}

# Adds columns in_window, group and kept to `starts` (which has graph and elbo).
.select_chain_starts <- function(starts, window, merge_distance) {
    starts$in_window <- starts$elbo > max(starts$elbo) - window
    starts$group <- NA_integer_
    starts$kept <- FALSE
    w <- which(starts$in_window)
    if (length(w) == 1 || merge_distance <= 0) {
        starts$group[w] <- seq_along(w)
        starts$kept[w] <- TRUE
        return(starts)
    }
    D <- .edit_distance_matrix(starts$graph[w])
    A <- (D <= merge_distance) & (D > 0)
    memb <- igraph::components(igraph::graph_from_adjacency_matrix(A * 1, mode = "undirected"))$membership
    starts$group[w] <- memb
    for (gr in unique(memb)) {
        members <- w[memb == gr]
        starts$kept[members[which.max(starts$elbo[members])]] <- TRUE
    }
    starts
}

# Graph strings that appeared (proposed or current) in chain info
.seen_graphs <- function(info) {
    if (is.null(info) || nrow(info) == 0) {
        return(character(0))
    }
    s <- unique(c(info$curr_graph, info$prop_graph))
    s[!is.na(s)]
}
