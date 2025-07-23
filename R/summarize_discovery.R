#' @export
discovery_summary <- function(x) {
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

    # Get the order that the models were visited
    # Could either/both do the actual index time as well as the model visit index
    # lapply(seq_along(x$mh_chain), function(i) {
    #    match(flat_graphs, x$mh_chain[[i]])
    # })

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

#' @export
summary.nesmr_mh_graph_explore <- discovery_summary

# TODO: Add variable names here or read from the results object
#' @export
edge_inclusion_probs <- function(x, min_prob_threshold = 0) {
    discovery_table <- discovery_summary(x)

    discovery_table <- discovery_table[discovery_table$norm_elbo > min_prob_threshold, ]

    inclusion_prob_long <- purrr::map_dfr(seq_along(discovery_table$graph), function(i) {
        g <- discovery_table$graph[i]
        cbind(matrix_to_edgelist(flat_string_to_adj_mat(g)),
              norm_elbo = discovery_table$norm_elbo[i], graph_i = i) |>
        dplyr::filter(value > 0)
    }) |>
    dplyr::group_by(from, to) |>
    dplyr::summarise(
        inclusion_prob = sum(norm_elbo)
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(desc(inclusion_prob)) |>
    as.data.frame()

    return(inclusion_prob_long)
}

#' @export
top_i_graph <- function(x, i = 1) {
    all_elbos <- sapply(x$visited_graphs, function(g) g$elbo)
    norm_elbo <- exp(all_elbos - x$elbo_denom)
    # Get the min index of the graph
    graph_idx <- order(norm_elbo, decreasing = TRUE)[i]

    top_graph <- x$visited_graphs[[graph_idx]]
    return(top_graph)
}