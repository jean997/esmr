#' Enumerate all graphs within a given edit distance of a DAG
#'
#' Breadth-first enumeration of every DAG reachable from `B` using at most
#' `edit_distance` single-edge moves. The moves are the same ones used by
#' [mh_graph_explore()]: add an edge, remove an edge, or reverse an edge
#' (flip). A move is only allowed if the result is still a DAG. The distance
#' reported for a graph is the *minimum* number of moves needed to reach it.
#'
#' No models are fit here, so this is cheap and can be used to see how many fits
#' [fit_neighbor_graphs()] would need before committing to them.
#'
#' @param B A square 0/1 (or non-zero) adjacency matrix. `B[i, j] != 0` means an
#'   edge `i -> j`. Must be a DAG.
#' @param edit_distance Maximum number of moves from `B`. Default 1.
#' @param moves Which moves are allowed. Any of `"add"`, `"remove"`, `"flip"`.
#'
#' @return A data frame with one row per graph at distance `1..edit_distance`
#'   (`B` itself is not included) and columns `graph` (flat string key, the same
#'   format as the names of `visited_graphs`), `distance` and `n_edges`.
#' @export
graph_neighbors <- function(B, edit_distance = 1, moves = c("add", "remove", "flip")) {
    moves <- match.arg(moves, several.ok = TRUE)
    stopifnot(is.matrix(B), nrow(B) == ncol(B))
    stopifnot(is.numeric(edit_distance), length(edit_distance) == 1, edit_distance >= 1)
    B <- (B != 0) + 0
    diag(B) <- 0
    if (!is_dag(B)) stop("B must be a DAG")

    key <- function(M) paste0(M, collapse = "")
    seen <- key(B)
    frontier <- list(B)
    out <- list()

    for (dist in seq_len(edit_distance)) {
        next_frontier <- list()
        for (G in frontier) {
            for (M in .one_edit_neighbors(G, moves)) {
                k <- key(M)
                if (!(k %in% seen)) {
                    seen <- c(seen, k)
                    next_frontier[[length(next_frontier) + 1]] <- M
                    out[[length(out) + 1]] <- data.frame(
                        graph = k, distance = dist, n_edges = sum(M),
                        stringsAsFactors = FALSE
                    )
                }
            }
        }
        frontier <- next_frontier
        if (length(frontier) == 0) break
    }
    if (length(out) == 0) {
        return(data.frame(graph = character(0), distance = numeric(0), n_edges = numeric(0)))
    }
    do.call(rbind, out)
}

# Transitive closure of a 0/1 adjacency matrix: reach[a, b] is TRUE iff there is a
# directed path a -> ... -> b of length >= 1.
.reachability <- function(A) {
    reach <- A != 0
    repeat {
        new <- reach | ((reach + 0) %*% (A != 0) > 0)
        if (identical(new, reach)) break
        reach <- new
    }
    reach
}

# All DAGs one move away from G (list of 0/1 matrices).
.one_edit_neighbors <- function(G, moves) {
    d <- nrow(G)
    out <- list()
    if ("add" %in% moves) {
        reach <- .reachability(G)
        for (i in seq_len(d)) for (j in seq_len(d)) {
            # Adding i -> j closes a cycle iff j already reaches i. This also rules out
            # adding i -> j when j -> i exists (a 2-cycle).
            if (i != j && G[i, j] == 0 && !reach[j, i]) {
                M <- G
                M[i, j] <- 1
                out[[length(out) + 1]] <- M
            }
        }
    }
    edges <- which(G != 0, arr.ind = TRUE)
    for (e in seq_len(nrow(edges))) {
        i <- edges[e, 1]
        j <- edges[e, 2]
        M <- G
        M[i, j] <- 0
        if ("remove" %in% moves) out[[length(out) + 1]] <- M
        if ("flip" %in% moves && G[j, i] == 0) {
            # Reverse i -> j to j -> i. Valid iff i does not still reach j without i -> j.
            if (!.reachability(M)[i, j]) {
                M2 <- M
                M2[j, i] <- 1
                out[[length(out) + 1]] <- M2
            }
        }
    }
    out
}

#' Fit all unvisited graphs near a center graph from an MH graph exploration
#'
#' Given a `nesmr_mh_graph_explore` object (the result of [mh_graph_explore()]),
#' enumerate every DAG within `edit_distance` single-edge moves (add / remove /
#' flip) of a center graph (by default the best visited graph), drop the ones that
#' were already visited, and fit the rest with [nesmr()]. This is a deterministic
#' complement to the MH chain: it answers "is the best graph a local optimum?" and
#' fills in the neighborhood of the top graphs, which improves edge inclusion
#' weights.
#'
#' The `explore` object is **not modified**. Fits are returned in a separate
#' `nesmr_neighbor_fits` object whose `fits` element has the same structure as
#' `explore$visited_graphs`, so it can be merged later if desired.
#'
#' ## Matching the original run
#' The `explore` object does not store the data or fitting options, so they must be
#' supplied again and must match the call that produced `explore`: `dat`, the
#' variant selection (`alpha` / `pval_select`), `R`, `nesmr_params` (including any
#' `fix_g` / `g_init`) and `graph_edge_prior`. Otherwise ELBOs are not comparable.
#' With `check_consistency = TRUE` (default) one graph that was already fit is
#' re-fit and its ELBO is compared with the stored value, which catches most
#' mismatches.
#'
#' @param explore A `nesmr_mh_graph_explore` object.
#' @param dat List with `beta_hat` and `s_estimate`, as given to [mh_graph_explore()].
#' @param edit_distance Maximum number of moves from the center. Default 1. The number
#'   of graphs grows quickly with this; use `dry_run = TRUE` to see the count first.
#' @param center Center graph. `NULL` (default) uses the best visited graph. May also
#'   be a flat string key from `explore$visited_graphs`, an adjacency matrix, or an
#'   integer `i` for the `i`-th best visited graph.
#' @param moves Allowed moves; any of `"add"`, `"remove"`, `"flip"`.
#' @param alpha,pval_select Variant selection, as in [mh_graph_explore()].
#' @param R LD / residual correlation matrix, as in [mh_graph_explore()].
#' @param nesmr_params List passed to [nesmr()] as `params`.
#' @param graph_edge_prior Edge prior used by [log_graph_prior()]. Default 0.5.
#' @param nesmr_max_iter `max_iter` passed to [nesmr()]. [mh_graph_explore()] uses 300.
#' @param max_fits Refuse to fit more than this many graphs (guards against an
#'   accidental large `edit_distance`). Default 500.
#' @param dry_run If `TRUE`, only enumerate and return the planned graphs.
#' @param check_consistency Re-fit one already-visited graph and compare with the
#'   stored ELBO before fitting anything else. Default `TRUE`.
#' @param consistency_tol Absolute ELBO tolerance for the consistency check.
#' @param cores Number of cores for [parallel::mclapply()] (Unix only). Default 1.
#' @param checkpoint_file Optional path. Partial results are saved here with
#'   `saveRDS()` after each chunk of `checkpoint_every` fits.
#' @param checkpoint_every Chunk size for fitting/checkpointing. Default 50.
#' @param verbose Print progress.
#'
#' @return An object of class `nesmr_neighbor_fits`, a list with
#'   - `center`: flat string of the center graph and its stored `center_elbo`
#'   - `graphs`: data frame (`graph`, `distance`, `n_edges`, `elbo`, `elbo_gain`
#'     relative to the center, `fitted`, `error`)
#'   - `fits`: named list (by graph string) of `elbo`, `elbo_without_prior`,
#'     `beta_hat`, `se_beta_hat`, same layout as `explore$visited_graphs`
#'   - `config`, `consistency`, and counts of graphs enumerated / already visited.
#' @export
fit_neighbor_graphs <- function(
    explore,
    dat,
    edit_distance = 1,
    center = NULL,
    moves = c("add", "remove", "flip"),
    alpha = 5e-8,
    pval_select = NULL,
    R = NULL,
    nesmr_params = list(beta_prior_cov = 1, strict_mode = FALSE),
    graph_edge_prior = 0.5,
    nesmr_max_iter = 300,
    max_fits = 500,
    dry_run = FALSE,
    check_consistency = TRUE,
    consistency_tol = 1e-3,
    cores = 1L,
    checkpoint_file = NULL,
    checkpoint_every = 50,
    verbose = FALSE) {
    stopifnot(inherits(explore, "nesmr_mh_graph_explore"))
    log_msg <- function(...) if (verbose) message("fit_neighbor_graphs: ", ...)
    moves <- match.arg(moves, several.ok = TRUE)
    visited <- explore$visited_graphs
    d <- ncol(dat$beta_hat)
    visited_elbo <- vapply(visited, function(g) g$elbo %||% NA_real_, numeric(1))

    # ---- center graph ----
    if (is.null(center)) center <- 1L
    if (is.matrix(center)) {
        center_str <- paste0((center != 0) + 0, collapse = "")
    } else if (is.character(center)) {
        center_str <- center
    } else if (is.numeric(center) && length(center) == 1) {
        center_str <- names(sort(visited_elbo, decreasing = TRUE))[center]
    } else {
        stop("center must be NULL, an adjacency matrix, a graph string or a rank")
    }
    if (is.na(center_str) || nchar(center_str) != d^2) stop("center is not a valid graph for dat")
    B <- flat_string_to_adj_mat(center_str)

    # ---- candidates ----
    cand <- graph_neighbors(B, edit_distance = edit_distance, moves = moves)
    cand <- cand[cand$n_edges > 0, , drop = FALSE] # nesmr cannot fit an empty graph
    n_enum <- nrow(cand)
    already <- cand$graph %in% names(visited)
    todo <- cand[!already, , drop = FALSE]
    log_msg(sprintf(
        "%d graphs within distance %d of the center; %d already visited; %d to fit",
        n_enum, edit_distance, sum(already), nrow(todo)
    ))
    center_elbo <- visited_elbo[center_str]
    if (is.na(center_elbo)) center_elbo <- NA_real_

    config <- list(
        alpha = alpha, graph_edge_prior = graph_edge_prior, nesmr_params = nesmr_params,
        edit_distance = edit_distance, moves = moves, nesmr_max_iter = nesmr_max_iter
    )
    new_result <- function(graphs, fits, consistency) {
        structure(list(
            center = center_str, center_elbo = unname(center_elbo),
            graphs = graphs, fits = fits, config = config, consistency = consistency,
            n_enumerated = n_enum, n_already_visited = sum(already)
        ), class = "nesmr_neighbor_fits")
    }
    empty_graphs <- function(df) {
        df$elbo <- NA_real_
        df$elbo_gain <- NA_real_
        df$fitted <- FALSE
        df$error <- NA_character_
        rownames(df) <- NULL
        df
    }

    if (dry_run) {
        return(new_result(empty_graphs(todo), list(), list(checked = FALSE)))
    }
    if (nrow(todo) > max_fits) {
        stop(sprintf(
            "%d graphs need fitting, which exceeds max_fits = %d. Lower edit_distance, restrict `moves`, or raise max_fits (use dry_run = TRUE to inspect).",
            nrow(todo), max_fits
        ))
    }

    # ---- variant selection and fitter: identical to mh_graph_explore ----
    if (is.null(pval_select)) {
        pval_select <- 2 * pnorm(-abs(dat$beta_hat / dat$s_estimate))
    }
    ix <- which(apply(pval_select, 1, min) < alpha)
    config$n_variants <- length(ix)

    fit_one <- function(M) {
        capture.output(
            {
                mod <- nesmr(
                    beta_hat_X = dat$beta_hat,
                    se_X = dat$s_estimate,
                    variant_ix = ix,
                    direct_effect_template = M,
                    max_iter = nesmr_max_iter,
                    R = R,
                    params = nesmr_params
                )
            },
            file = nullfile()
        )
        list(
            elbo_without_prior = mod$elbo,
            elbo = mod$elbo + log_graph_prior(sum(M), d, pi_0 = graph_edge_prior),
            beta_hat = mod$direct_effects,
            se_beta_hat = mod$direct_effects_se
        )
    }

    # ---- consistency check: refit an already-visited graph ----
    consistency <- list(checked = FALSE)
    if (check_consistency && !is.na(center_elbo)) {
        log_msg("Refitting the center graph to check that dat / options match the original run")
        refit <- fit_one(B)
        diff <- refit$elbo - center_elbo
        consistency <- list(checked = TRUE, stored_elbo = unname(center_elbo), refit_elbo = refit$elbo, diff = unname(diff))
        if (!is.finite(diff) || abs(diff) > consistency_tol) {
            stop(sprintf(
                "Refit ELBO of the center graph (%.4f) differs from the stored ELBO (%.4f) by %.4g. `dat`, alpha/pval_select, R, nesmr_params or graph_edge_prior probably do not match the original mh_graph_explore() call.",
                refit$elbo, center_elbo, diff
            ))
        }
    }

    # ---- fit in chunks (parallel within a chunk, checkpoint between chunks) ----
    graphs <- empty_graphs(todo)
    fits <- list()
    chunk_size <- max(1L, as.integer(if (checkpoint_every > 0) checkpoint_every else max(nrow(todo), 1L)))
    chunks <- split(seq_len(nrow(todo)), ceiling(seq_len(nrow(todo)) / chunk_size))
    for (ch in chunks) {
        res <- .lapply_cores(ch, function(r) {
            M <- flat_string_to_adj_mat(graphs$graph[r])
            tryCatch(fit_one(M), error = function(e) e)
        }, cores = cores)
        for (k in seq_along(ch)) {
            r <- ch[k]
            if (inherits(res[[k]], "error") || inherits(res[[k]], "try-error") || is.null(res[[k]])) {
                graphs$error[r] <- if (inherits(res[[k]], "error")) conditionMessage(res[[k]]) else as.character(res[[k]])
            } else {
                fits[[graphs$graph[r]]] <- res[[k]]
                graphs$elbo[r] <- res[[k]]$elbo
                graphs$elbo_gain[r] <- res[[k]]$elbo - center_elbo
                graphs$fitted[r] <- TRUE
            }
        }
        log_msg(sprintf("Fitted %d / %d graphs", sum(graphs$fitted), nrow(graphs)))
        if (!is.null(checkpoint_file)) {
            saveRDS(new_result(graphs, fits, consistency), checkpoint_file)
        }
    }
    new_result(graphs, fits, consistency)
}

.lapply_cores <- function(x, f, cores = 1L) {
    if (cores > 1L && .Platform$OS.type != "windows") {
        parallel::mclapply(x, f, mc.cores = cores)
    } else {
        lapply(x, f)
    }
}

#' Merge several neighbor-fit results into one cumulative object
#'
#' Combines `nesmr_neighbor_fits` objects (for example one per round of
#' [fit_neighbor_graphs()], each centered on the best graph found so far) into a
#' single object. The merged object has the same structure, so it can be passed to
#' [merge_neighbor_fits()] again or used as `prior_fits` in the pipeline.
#'
#' A graph fitted in more than one object keeps the earliest fit (fits are
#' deterministic for identical data and options). In `graphs`, `distance` and
#' `elbo_gain` stay relative to the center of the run that fitted the graph; the
#' columns `source` and `center` identify that run, and `elbo_vs_best` is the ELBO
#' minus the best ELBO over the merged fits. `runs` has one row per input.
#'
#' @param ... `nesmr_neighbor_fits` objects, oldest first. Already merged objects are
#'   allowed.
#' @param sources Optional character vector with one label per object in `...`
#'   (default `run1`, `run2`, ...). Ignored for rows of an already merged object, which
#'   keep their `source`.
#' @return An object of class `nesmr_neighbor_fits` with the latest object's `center`,
#'   `center_elbo`, `config` and `consistency`, merged `graphs` and `fits`, a `runs` data
#'   frame, and `best_graph` / `best_elbo` over all fitted graphs.
#' @export
merge_neighbor_fits <- function(..., sources = NULL) {
    objs <- list(...)
    stopifnot(length(objs) >= 1, all(vapply(objs, inherits, NA, "nesmr_neighbor_fits")))
    if (is.null(sources)) sources <- paste0("run", seq_along(objs))
    stopifnot(length(sources) == length(objs))

    graphs <- list()
    runs <- list()
    fits <- list()
    for (k in seq_along(objs)) {
        o <- objs[[k]]
        g <- o$graphs
        g$elbo_vs_best <- NULL # recomputed after merging
        if (is.null(g$source)) g$source <- sources[k]
        if (is.null(g$center)) g$center <- o$center
        graphs[[k]] <- g
        runs[[k]] <- if (!is.null(o$runs)) o$runs else data.frame(
            source = sources[k], center = o$center, center_elbo = o$center_elbo,
            edit_distance = o$config$edit_distance, moves = paste(o$config$moves, collapse = ","),
            n_fitted = sum(o$graphs$fitted), stringsAsFactors = FALSE
        )
        fits <- c(fits, o$fits[setdiff(names(o$fits), names(fits))])
    }
    graphs <- do.call(rbind, graphs)
    graphs <- graphs[order(!graphs$fitted), , drop = FALSE] # fitted rows win the de-duplication
    graphs <- graphs[!duplicated(graphs$graph), , drop = FALSE]
    graphs <- graphs[order(match(graphs$source, unique(unlist(lapply(runs, `[[`, "source"))))), , drop = FALSE]
    graphs$elbo_vs_best <- graphs$elbo - suppressWarnings(max(graphs$elbo, na.rm = TRUE))
    rownames(graphs) <- NULL
    latest <- objs[[length(objs)]]
    best_ix <- if (any(graphs$fitted)) which.max(graphs$elbo) else NA_integer_
    structure(list(
        center = latest$center, center_elbo = latest$center_elbo,
        graphs = graphs, fits = fits, config = latest$config, consistency = latest$consistency,
        n_enumerated = latest$n_enumerated, n_already_visited = latest$n_already_visited,
        runs = do.call(rbind, runs),
        best_graph = if (is.na(best_ix)) NA_character_ else graphs$graph[best_ix],
        best_elbo = if (is.na(best_ix)) NA_real_ else graphs$elbo[best_ix]
    ), class = "nesmr_neighbor_fits")
}

#' @export
print.nesmr_neighbor_fits <- function(x, ...) {
    g <- x$graphs
    merged <- !is.null(x$runs)
    latest <- if (merged) g[g$source == x$runs$source[nrow(x$runs)], , drop = FALSE] else g
    cat("Neighbor-graph fits around the center graph\n")
    cat(sprintf(
        "  %sedit distance <= %s | %d graphs enumerated, %d already visited, %d to fit, %d fitted, %d failed\n",
        if (merged) "latest run: " else "", x$config$edit_distance, x$n_enumerated, x$n_already_visited,
        nrow(latest), sum(latest$fitted), sum(!latest$fitted & !is.na(latest$error))
    ))
    if (isTRUE(x$consistency$checked)) {
        cat(sprintf("  consistency check: refit ELBO differs by %.2g\n", x$consistency$diff))
    }
    if (any(latest$fitted)) {
        best <- latest[which.max(latest$elbo), ]
        cat(sprintf(
            "  center ELBO %.3f | best new graph ELBO %.3f (gain %+.3f, distance %d)\n",
            x$center_elbo, best$elbo, best$elbo_gain, best$distance
        ))
        cat(sprintf("  new graphs better than center: %d\n", sum(latest$elbo_gain > 0, na.rm = TRUE)))
    }
    if (merged) {
        cat(sprintf(
            "  merged: %d runs, %d graphs fitted in total | best ELBO %.3f\n",
            nrow(x$runs), sum(g$fitted), x$best_elbo
        ))
    }
    invisible(x)
}
