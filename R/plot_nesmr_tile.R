#' @export
plot_nesmr_tile.nesmr_tbl_graph <- function(x, weight = c("direct_effect", "total_effect"), ...) {
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
      as.integer(igraph::topo_sort(ig - igraph::feedback_arc_set(ig)))
    }
  )

  # Get the weight values for scale limits
  weight_values <- x %>%
    tidygraph::activate(edges) %>%
    pull(!!sym(weight))

  scale_limit <- max(abs(weight_values), na.rm = TRUE)

  # Create nice title from weight name
  weight_title <- tools::toTitleCase(gsub("_", " ", weight))

  x %>%
    tidygraph::activate(edges) %>%
    mutate(
      from = factor(from, levels = node_order),
      to = factor(to, levels = node_order)
    ) %>%
    tidygraph::as_tibble() %>%
    ggplot(., aes(x = to, y = from, fill = !!sym(weight))) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = ifelse(abs(!!sym(weight)) > 0.01, sprintf("%.2f", !!sym(weight)), "")),
      size = 5, color = "black") +
    scale_fill_gradient2(
      low = "#1f77b4",
      mid = "white",
      high = "#d62728",
      midpoint = 0,
      limits = c(-scale_limit, scale_limit),
      name = weight_title
    ) +
    scale_y_discrete(limits = factor(rev(node_order))) +  # Reverse y-axis for matrix ordering
    scale_x_discrete(limits = factor(node_order)) +
    labs(title = paste(weight_title, "Matrix"),
         x = "To",
         y = "From") +
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
  edge_inclusion_tile_plot(eip, ...)
}

#' @export
edge_inclusion_tile_plot.nesmr_tbl_graph <- function(
  x,
  weight = c("direct_effect", "total_effect"),
  plot_type = c("adj", "beta"),
  pos_color = "#d62728",
  neg_color = "#1f77b4",
  scale_limit = NULL,
  scale_factor = 1) {

  weight <- match.arg(weight)
  plot_type <- match.arg(plot_type)

  if (!inherits(x, "nesmr_tbl_graph")) {
    x <- tidygraph::as_tbl_graph(x)
  }

  node_order <- tryCatch(
    unlist(layered_topological_sort(x)),
    error = function(e) {
      # TODO: Might be a better way to handle this?
      # This does maximal feedback arc set and topo sort from there
      as.integer(igraph::topo_sort(x - igraph::feedback_arc_set(x)))
    }
  )
  # TODO: Add a extra row of tiles (like BPG plots) for node layer
  x %>%
    tidygraph::activate(edges) %>%
    mutate(
      from = factor(from, levels = node_order),
      to = factor(to, levels = node_order)
    ) %>%
    tidygraph::as_tibble() %>%
  ggplot(., aes(x = to, y = from, fill = round(inclusion_prob, 4))) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = ifelse(inclusion_prob > 0.01, round(inclusion_prob, 2), "")), size = 5, color = "black") +
    scale_fill_gradient(
      low = "white", high = "orange",
      name = "Edge Inclusion",
      limits = c(0, 1)
    ) +
    scale_y_discrete(limits = factor(rev(node_order))) +  # Reverse y-axis for matrix ordering
    scale_x_discrete(limits = factor(node_order)) +
    labs(title = "Edge Inclusion Probability",
        x = "To",
        y = "From") +
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
  if (type == "edge_inc_prob") {
    edge_inclusion_tile_plot.nesmr_tbl_graph(edge_inclusion_probs(x))
  } else if (type == "best") {
    # Get the best graph
    best_graph <- top_i_graph(x, i = 1)
    plot_nesmr_tile.nesmr_tbl_graph(best_graph, weight = "direct_effect", plot_type = "beta")
  } else {
    stop("Unknown type for nesmr_mh_graph_explore.")
  }
}