# Backfill missing flip candidate information for saved visited_graph objects.
backfill_visited_graph_flips <- function(visited_graphs, weight_mat, verbose = FALSE) {
    if (!is.list(visited_graphs)) {
        stop("visited_graphs must be a list")
    }
    if (!is.matrix(weight_mat) || nrow(weight_mat) != ncol(weight_mat)) {
        stop("weight_mat must be a square matrix")
    }

    d <- nrow(weight_mat)
    graph_names <- names(visited_graphs)
    if (is.null(graph_names)) {
        stop("visited_graphs must be a named list where names are flattened adjacency matrices")
    }

    n_updated <- 0L
    for (graph_name in graph_names) {
        graph_info <- visited_graphs[[graph_name]]

        # Nothing to backfill unless adjacency proposal metadata exists.
        if (is.null(graph_info$adj_graph_info)) {
            next
        }

        adj_graph_info <- graph_info$adj_graph_info
        # has_flip <- !is.null(adj_graph_info$flip_candidates) ||
        #     !is.null(adj_graph_info$flip_canidates)
        # has_flip_prob <- !is.null(adj_graph_info$flip_edge_prob)
        # if (has_flip && has_flip_prob) {
        #     next
        # }

        graph_vec <- as.numeric(strsplit(graph_name, "")[[1]])
        if (length(graph_vec) != d * d) {
            warning(sprintf("Skipping graph '%s': name length does not match d^2", graph_name))
            next
        }

        adj_mat <- matrix(graph_vec, nrow = d, ncol = d)
        g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "directed")

        logistic_scale <- adj_graph_info$logistic_scale
        if (is.null(logistic_scale)) logistic_scale <- 1

        logistic_location <- adj_graph_info$logistic_location
        if (is.null(logistic_location)) logistic_location <- 5

        flip_logistic_scale <- adj_graph_info$flip_logistic_scale
        if (is.null(flip_logistic_scale)) flip_logistic_scale <- 3

        recomputed_adj <- get_adjacent_graphs(
            g = g,
            weight_mat = weight_mat,
            logistic_scale = logistic_scale,
            logistic_location = logistic_location,
            flip_logistic_scale = flip_logistic_scale
        )

        visited_graphs[[graph_name]]$adj_graph_info <- recomputed_adj
    }

    if (isTRUE(verbose)) {
        message(sprintf("backfill_visited_graph_flips: updated %d graph(s)", n_updated))
    }

    visited_graphs
}