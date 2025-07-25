#' @export
discovery_summary <- function(x) {
    # Want:
    # Table of graphs with number of edges, elbo, norm_elbo, and visited count
    n <- length(x$visited_graphs)
    all_elbos <- sapply(x$visited_graphs, function(g) g$elbo)
    norm_elbo <- x$norm_elbo
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

#' @export
print.nesmr_mh_graph_explore <- function(x, ...) {
    cat("Summary of discovery from MH graph exploration:\n")
    summary_table <- discovery_summary(x)
    print(summary_table)
    invisible(x)
}

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

    # TODO: Remove this check if we make the graph a tidygraph
    if (!inherits(top_graph, "tidygraph")) {
        edgelist <- which(top_graph$beta_hat != 0, arr.ind = TRUE)

        # edgelist <- which(x$direct_effect != 0, arr.ind = T)
        edgelist <- data.frame(from = edgelist[, 1],
                            to = edgelist[, 2],
                            direct_effect = top_graph$beta_hat[edgelist],
                            direct_effect_se = top_graph$se_beta_hat[edgelist])

        # TODO: Add node names from the object
        nodes <- data.frame(name = 1:ncol(top_graph$beta_hat))

        tg <- tidygraph::tbl_graph(nodes = nodes, edges = edgelist)
        tg$elbo <- top_graph$elbo
        tg$norm_elbo <- norm_elbo[graph_idx]
        class(tg) <- c("discovery_tbl_graph", "nesmr_tbl_graph", "esmr_tbl_graph", class(tg))
        top_graph <- tg
    }
    return(top_graph)
}

#' @export
plot.nesmr_mh_graph_explore <- function(
    x,
    max_graphs = 20,
    plot_type = c("norm_elbo", "cum_norm_elbo", "both"),
    ...) {
    plot_type <- match.arg(plot_type)
    if (max_graphs > length(x$norm_elbo)) {
        max_graphs <- length(x$norm_elbo)
    }
    if (max_graphs < 1) {
        stop("max_graphs must be at least 1.")
    }

    top_norm_elbo <- sort(x$norm_elbo, decreasing = TRUE)[1:max_graphs]
    cum_norm_elbo <- cumsum(top_norm_elbo)
    plot_df <- data.frame(
        Index = 1:max_graphs,
        norm_elbo = top_norm_elbo,
        cum_norm_elbo = cum_norm_elbo
    )

    if (plot_type == "norm_elbo") {
        ggplot(plot_df, aes(x = Index, y = norm_elbo)) +
            geom_point(size = 3) +
            scale_y_continuous(limits = c(0, 1.05)) +
            labs(
                title = "Normalized ELBO for Top Graphs",
                x = "Graph Index",
                y = "Normalized ELBO"
            ) +
            theme_classic(base_size = 20)
    } else if (plot_type == "cum_norm_elbo") {
        ggplot(plot_df, aes(x = Index, y = cum_norm_elbo)) +
            geom_point(size = 3) +
            geom_line(size = 1) +
            scale_y_continuous(limits = c(0, 1.05)) +
            labs(
                title = "Cumulative Normalized ELBO for Top Graphs",
                x = "Graph Index",
                y = "Cumulative Normalized ELBO"
            ) +
            theme_classic(base_size = 20)
    } else if (plot_type == "both") {
        ggplot(plot_df, aes(x = Index)) +
            geom_point(aes(y = cum_norm_elbo, color = "Cumulative Norm ELBO"), size = 3) +
            geom_line(aes(y = cum_norm_elbo, color = "Cumulative Norm ELBO"), size = 1) +
            geom_point(aes(y = norm_elbo, color = "Norm ELBO"), size = 3) +
            scale_y_continuous(limits = c(0, 1.05)) +
            labs(
                title = "Normalized ELBO and Cumulative ELBO for Top Graphs",
                x = "Graph Index",
                y = "Normalized ELBO",
                color = "Legend"
            ) +
            scale_color_manual(values = c("Norm ELBO" = "blue", "Cumulative Norm ELBO" = "orange")) +
            theme_classic(base_size = 20)
    }
}