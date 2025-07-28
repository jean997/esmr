#' @export
plot_nesmr_tile.nesmr_tbl_graph <- function(
  x, weight = c("direct_effect", "total_effect"),
  neg_color = "#1f77b4",
  pos_color = "#d62728",
  ...) {
  # If we are plotting a graph, use the weights from direct/total effects
  weight <- match.arg(weight)

  if (!inherits(x, "nesmr_tbl_graph")) {
    x <- tidygraph::as_tbl_graph(x)
  }

  node_order <- tryCatch(
    unlist(layered_topological_sort(x)),
    error = function(e) {
      # TODO: Might be a better way to handle this?
      # This does maximal feedback arc set and topo sort from there
      ig <- igraph::as.igraph(x)
      ts <- as.integer(igraph::topo_sort(ig - igraph::feedback_arc_set(ig)))
      c(ts, setdiff(seq_len(nrow(x)), ts))  # Ensure all nodes are included
    }
  )

  # Get the weight values for scale limits
  weight_values <- x %>%
    tidygraph::activate(edges) %>%
    pull(!!sym(weight))

  scale_limit <- max(abs(weight_values), na.rm = TRUE)

  # Create nice title from weight name
  weight_title <- tools::toTitleCase(gsub("_", " ", weight))

  # Get all node names for diagonal elements
  all_nodes <- x %>% tidygraph::activate(nodes) %>% pull(name)
  diag_df <- data.frame(
    from = all_nodes,
    to = all_nodes
  )

  x %>%
    tidygraph::activate(edges) %>%
    tidygraph::as_tibble() %>%
    full_join(
      diag_df,
      by = c("from", "to")
    ) %>%
    mutate(
      from = factor(from, levels = node_order),
      to = factor(to, levels = rev(node_order)),
      !!sym(weight) := ifelse(from == to, NA, !!sym(weight))
    ) %>%
    ggplot(., aes(x = to, y = from, fill = !!sym(weight))) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = ifelse(abs(!!sym(weight)) > 0.01, sprintf("%.2f", !!sym(weight)), "")),
      size = 5, color = "black") +
    scale_fill_gradient2(
      low = neg_color,
      mid = "white",
      high = pos_color,
      midpoint = 0,
      limits = c(-scale_limit, scale_limit),
      name = weight_title,
      na.value = "lightgrey"
    ) +
    labs(title = paste(weight_title, "Matrix"),
         x = "To",
         y = "From") +
    scale_y_discrete(drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_equal() +
    geom_vline(xintercept = seq(0.5, length(node_order) + 0.5, by = 1), color = "grey80") +
    geom_hline(yintercept = seq(0.5, length(node_order) + 0.5, by = 1), color = "grey80")
}

#' @export
plot_nesmr_tile <- function(x) {
  UseMethod("plot_nesmr_tile")
}

#' @export
edge_inclusion_tile_plot <- function(x, ...) {
  UseMethod("edge_inclusion_tile_plot")
}

#' @export
edge_inclusion_tile_plot.nesmr_mh_graph_explore <- function(x, ...) {
  eip <- edge_inclusion_probs(x)
  best_graph <- top_i_graph(x, i = 1)
  node_order <- unlist(layered_topological_sort(best_graph))
  edge_inclusion_tile_plot(eip, node_order = node_order, ...)
}

#' @export
edge_inclusion_tile_plot.nesmr_tbl_graph <- function(
  x,
  weight = c("direct_effect", "total_effect"),
  plot_type = c("adj", "beta"),
  pos_color = "#d62728",
  neg_color = "#1f77b4",
  scale_limit = NULL,
  scale_factor = 1,
  node_order = NULL) {

  weight <- match.arg(weight)
  plot_type <- match.arg(plot_type)

  if (!inherits(x, "nesmr_tbl_graph")) {
    x <- tidygraph::as_tbl_graph(x)
  }

  if (is.null(node_order)) {
    # Get the node order from the layered topological sort
    node_order <- tryCatch(
      unlist(layered_topological_sort(x)),
      error = function(e) {
        # TODO: Might be a better way to handle this?
        # This does maximal feedback arc set and topo sort from there
        ts <- as.integer(igraph::topo_sort(x - igraph::feedback_arc_set(x)))
        c(ts, setdiff(ts, tidygraph::activate(x, nodes) %>% pull(name)))  # Ensure all nodes are included
      })
  }

  # Get all node names for diagonal elements
  all_nodes <- x %>% tidygraph::activate(nodes) %>% pull(name)
  diag_df <- data.frame(
    from = all_nodes,
    to = all_nodes
  )

  # TODO: Add a extra row of tiles (like BPG plots) for node layer
  x %>%
    tidygraph::activate(edges) %>%
    tidygraph::as_tibble() %>%
    full_join(
      diag_df,
      by = c("from", "to")
    ) %>%
    mutate(
      from = factor(from, levels = node_order),
      to = factor(to, levels = rev(node_order)),
      inclusion_prob = ifelse(from == to, NA, inclusion_prob)
    ) %>%
  ggplot(., aes(x = to, y = from, fill = round(inclusion_prob, 4))) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = ifelse(inclusion_prob > 0.01, round(inclusion_prob, 2), "")), size = 5, color = "black") +
    scale_fill_gradient(
      low = "white", high = "orange",
      name = "Edge Inclusion",
      limits = c(0, 1),
      na.value = "lightgrey"
    ) +
    labs(title = "Edge Inclusion Probability",
        x = "To",
        y = "From") +
    scale_y_discrete(drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_equal() +
    geom_vline(xintercept = seq(0.5, length(node_order) + 0.5, by = 1), color = "grey80") +
    geom_hline(yintercept = seq(0.5, length(node_order) + 0.5, by = 1), color = "grey80")
}

#' @export
plot_nesmr_tile.nesmr_mh_graph_explore <- function(x, type = c("edge_inc_prob", "best")) {
  # Plot the edge inclusion probability or the best graph
  type <- match.arg(type)
  best_graph <- top_i_graph(x, i = 1)
  if (type == "edge_inc_prob") {
    # First get the ordering from th best graph
    node_order <- unlist(layered_topological_sort(best_graph))
    edge_inclusion_tile_plot.nesmr_tbl_graph(edge_inclusion_probs(x), node_order = node_order)
  } else if (type == "best") {
    # Get the best graph
    plot_nesmr_tile.nesmr_tbl_graph(best_graph, weight = "direct_effect", plot_type = "beta")
  } else {
    stop("Unknown type for nesmr_mh_graph_explore.")
  }
}