#' Tile plots of NESMR graphs
#'
#' Matrix-style tile plots of the edges of a graph. `plot_nesmr_tile()` shows effect sizes
#' of a `nesmr_tbl_graph` (or one of several summaries of a `nesmr_mh_graph_explore`
#' object); `edge_inclusion_tile_plot()` shows edge inclusion probabilities.
#'
#' @param x A `nesmr_tbl_graph` or `nesmr_mh_graph_explore` object.
#' @param weight Edge column to plot. `"indirect_effect"` is `total_effect - direct_effect`
#'   and requires a graph with total effects.
#' @param type For a `nesmr_mh_graph_explore`, what to plot: the edge inclusion
#'   probabilities (`"edge_inc_prob"`), direct (`"best"`), total (`"total"`) or indirect
#'   (`"indirect"`) effects of the best graph, the starting graph (`"start"`), or the
#'   frequency of each edge among proposed graphs (`"proposal_freq"`).
#' @param node_order Order of the nodes, either positions in the graph or node names (after
#'   `labels` are applied). Defaults to the layered topological order, of the best graph for
#'   a `nesmr_mh_graph_explore`.
#' @param labels Optional labels for the nodes, one per node in graph order or named by the
#'   current node names.
#' @param pos_color,neg_color Colors for positive and negative effects.
#' @param x_axis_position `"top"` or `"bottom"`: orientation of the matrix.
#' @param ... Further arguments passed to `create_tile_plot()`, e.g. `lower_tri_only`,
#'   `x_axis_at` or `x_text_angle`.
#' @return A ggplot object.
#' @name nesmr_tile_plots
NULL

#' @rdname nesmr_tile_plots
#' @export
plot_nesmr_tile <- function(x, ...) {
  UseMethod("plot_nesmr_tile")
}

#' @rdname nesmr_tile_plots
#' @export
plot_nesmr_tile.nesmr_tbl_graph <- function(
  x, weight = c("direct_effect", "total_effect", "indirect_effect"),
  neg_color = "#1f77b4",
  pos_color = "#d62728",
  x_axis_position = c("top", "bottom"),
  node_order = NULL,
  labels = NULL,
  ...) {
  weight <- match.arg(weight)
  x_axis_position <- match.arg(x_axis_position)

  x <- relabel_nodes(x, labels)
  node_order <- resolve_node_order(x, node_order)

  if (weight == "indirect_effect") {
    x <- x |>
      tidygraph::activate(edges) |>
      mutate(indirect_effect = total_effect - direct_effect)
  }

  weight_values <- x |>
    tidygraph::activate(edges) |>
    pull(!!sym(weight))
  scale_limit <- max(abs(weight_values), na.rm = TRUE)
  weight_title <- tools::toTitleCase(gsub("_", " ", weight))

  create_tile_plot(
    data = prepare_tile_plot_data(x, node_order, weight, x_axis_position),
    fill_var = weight,
    scale_type = "gradient2",
    scale_params = list(
      low = neg_color,
      mid = "white",
      high = pos_color,
      midpoint = 0,
      limits = c(-scale_limit, scale_limit),
      name = weight_title
    ),
    title = paste(weight_title, "Matrix"),
    x_axis_position = x_axis_position,
    node_order = node_order,
    ...
  )
}

#' @rdname nesmr_tile_plots
#' @export
plot_nesmr_tile.nesmr_mh_graph_explore <- function(
  x,
  type = c("edge_inc_prob", "best", "total", "indirect", "start", "proposal_freq"),
  node_order = NULL, ...) {
  type <- match.arg(type)
  if (is.null(node_order)) node_order <- best_node_order(x)

  switch(type,
    edge_inc_prob = edge_inclusion_tile_plot(x, node_order = node_order, ...),
    proposal_freq = plot_edge_prob_tile(
      edge_proposal_freq(x), "proposal_freq", "Proposal Frequency",
      "Frequency of Edges in Proposed Graphs", node_order = node_order, ...),
    start = plot_nesmr_tile(discovery_start_graph(x), node_order = node_order, ...),
    best = plot_nesmr_tile(top_i_graph(x, i = 1), node_order = node_order, ...),
    total = plot_nesmr_tile(
      top_i_graph(x, i = 1, include_total_effect = TRUE),
      weight = "total_effect", node_order = node_order, ...),
    indirect = plot_nesmr_tile(
      top_i_graph(x, i = 1, include_total_effect = TRUE),
      weight = "indirect_effect", node_order = node_order, ...)
  )
}

#' @rdname nesmr_tile_plots
#' @export
edge_inclusion_tile_plot <- function(x, ...) {
  UseMethod("edge_inclusion_tile_plot")
}

#' @rdname nesmr_tile_plots
#' @export
edge_inclusion_tile_plot.nesmr_mh_graph_explore <- function(x, node_order = NULL, ...) {
  if (is.null(node_order)) node_order <- best_node_order(x)
  edge_inclusion_tile_plot(edge_inclusion_probs(x), node_order = node_order, ...)
}

#' @rdname nesmr_tile_plots
#' @export
edge_inclusion_tile_plot.nesmr_tbl_graph <- function(x, ...) {
  plot_edge_prob_tile(
    x, "inclusion_prob", "Edge Inclusion", "Edge Inclusion Probability", ...)
}

## Tile plot of an edge column on a [0, 1] scale, shared by the edge inclusion and
## proposal frequency plots.
plot_edge_prob_tile <- function(
  x, weight, fill_name, title,
  high_color = "orange",
  node_order = NULL,
  labels = NULL,
  x_axis_position = c("top", "bottom"),
  ...) {
  x_axis_position <- match.arg(x_axis_position)

  x <- relabel_nodes(x, labels)
  node_order <- resolve_node_order(x, node_order)

  create_tile_plot(
    data = prepare_tile_plot_data(x, node_order, weight, x_axis_position),
    fill_var = weight,
    scale_type = "gradient",
    scale_params = list(
      low = "white",
      high = high_color,
      name = fill_name,
      limits = c(0, 1.0001)
    ),
    title = title,
    x_axis_position = x_axis_position,
    node_order = node_order,
    text_format = "%.2f",
    ...
  )
}

#' Prepare data for tile plots
#'
#' Internal function to prepare edge data with diagonal elements for tile plotting
#'
#' @param x A tidygraph object
#' @param node_order Character vector specifying node order
#' @param weight_var Character name of the weight variable to handle
#' @param x_axis_position Position of x-axis ("top" or "bottom") for factor level ordering
#'
#' @return Data frame prepared for tile plotting
prepare_tile_plot_data <- function(x, node_order, weight_var, x_axis_position = "top") {
  all_nodes <- as.character(tbl_node_names(x))
  diag_df <- data.frame(
    from = seq_along(all_nodes),
    to = seq_along(all_nodes)
  )

  x |>
    tidygraph::activate(edges) |>
    tidygraph::as_tibble() |>
    full_join(diag_df, by = c("from", "to")) |>
    mutate(
      from = factor(all_nodes[from], levels = node_order),
      to = factor(all_nodes[to], levels = if (x_axis_position == "top") rev(node_order) else node_order),
      # Set diagonal elements to NA for the weight variable
      !!weight_var := ifelse(from == to, NA, !!sym(weight_var))
    )
}

#' Common tile plot function
#'
#' Internal function to create consistent tile plots across all plot types
#'
#' @param data Data frame with from, to, and fill columns
#' @param fill_var Character name of the fill variable
#' @param scale_type One of "gradient2", "gradient", or "manual"
#' @param scale_params List of parameters for the scale function
#' @param title Plot title
#' @param x_label X-axis label (default "To")
#' @param y_label Y-axis label (default "From")
#' @param x_axis_position Orientation of the matrix: the order of the columns (`"top"`
#'   reverses them relative to the rows).
#' @param x_axis_at Side of the plot where the x-axis is drawn, defaults to
#'   `x_axis_position`.
#' @param x_text_angle Angle of the x-axis labels.
#' @param text_size Size of text labels on tiles
#' @param text_threshold Threshold for showing text labels
#' @param text_format sprintf format for text labels
#' @param base_size Base font size for theme
#' @param node_order Character vector of node order for grid lines
#' @param lower_tri_only Only draw grid lines for the lower triangle
#'
#' @return ggplot object
create_tile_plot <- function(
  data,
  fill_var,
  scale_type = c("gradient2", "gradient", "manual"),
  scale_params = list(),
  title = "",
  x_label = "To",
  y_label = "From",
  x_axis_position = c("top", "bottom"),
  x_axis_at = x_axis_position,
  x_text_angle = 0,
  text_size = 5,
  text_threshold = 0.01,
  text_format = "%.2f",
  base_size = 16,
  node_order = NULL,
  lower_tri_only = FALSE
) {
  scale_type <- match.arg(scale_type)
  x_axis_position <- match.arg(x_axis_position)
  x_axis_at <- match.arg(x_axis_at, c("top", "bottom"))

  # Create base plot
  p <- ggplot(data, aes(x = to, y = from, fill = !!sym(fill_var))) +
    geom_tile(color = "white", show.legend = TRUE)

  # Add text only if threshold is not Inf
  if (fill_var != "diff_sign_factor" && text_threshold < Inf) {
    p <- p + geom_text(
      data = data |>
        filter(!is.na(!!sym(fill_var))),
      aes(label = ifelse(
        abs(!!sym(fill_var)) > text_threshold,
        sprintf(text_format, !!sym(fill_var)), "")),
      size = text_size, color = "black"
    )
  }

  # Add appropriate scale
  scale_fn <- switch(scale_type,
    gradient2 = scale_fill_gradient2,
    gradient = scale_fill_gradient,
    manual = scale_fill_manual
  )
  p <- p + do.call(scale_fn, c(scale_params, list(na.value = "lightgrey")))

  # Add labels and scales
  p <- p +
    labs(title = title, x = x_label, y = y_label) +
    scale_y_discrete(drop = FALSE) +
    scale_x_discrete(drop = FALSE, position = x_axis_at) +
    theme_classic(base_size = base_size) +
    theme(
      axis.text.x = element_text(
        angle = x_text_angle,
        hjust = if (x_axis_at == "top") 0 else 1,
        vjust = if (x_text_angle == 0) NULL else if (x_axis_at == "top") 0 else 1
      )
    ) +
    coord_equal()

  if (is.null(node_order)) node_order <- levels(data$from)
  n_nodes <- length(node_order)

  if (lower_tri_only) {
    # Add grid lines, but only for the lower triangle
    grid_lines <- bind_rows(
      data.frame(
        x = seq_len(n_nodes) + 0.5,
        xend = seq_len(n_nodes) + 0.5,
        y = n_nodes - seq_len(n_nodes) + 0.5,
        yend = 0.5
      ),
      data.frame(
        x = 0.5,
        xend = seq_len(n_nodes) + 0.5,
        y = n_nodes - seq_len(n_nodes) + 0.5,
        yend = n_nodes - seq_len(n_nodes) + 0.5
      )
    )

    p <- p +
      geom_segment(
        data = grid_lines,
        inherit.aes = FALSE,
        aes(x = x, xend = xend, y = y, yend = yend),
        color = "grey80"
      ) +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.4)
      )
  } else {
    p <- p +
      geom_vline(xintercept = seq(0.5, n_nodes + 0.5, by = 1), color = "grey80") +
      geom_hline(yintercept = seq(0.5, n_nodes + 0.5, by = 1), color = "grey80")
  }

  p
}
