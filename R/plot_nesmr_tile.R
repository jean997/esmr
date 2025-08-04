#' @export
plot_nesmr_tile.nesmr_tbl_graph <- function(
  x, weight = c("direct_effect", "total_effect"),
  neg_color = "#1f77b4",
  pos_color = "#d62728",
  x_axis_position = c("top", "bottom"),
  ...) {
  # If we are plotting a graph, use the weights from direct/total effects
  weight <- match.arg(weight)
  x_axis_position <- match.arg(x_axis_position)

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

  plot_data <- prepare_tile_plot_data(x, node_order, weight, x_axis_position)

  # Create scale parameters
  scale_params <- list(
    low = neg_color,
    mid = "white",
    high = pos_color,
    midpoint = 0,
    limits = c(-scale_limit, scale_limit),
    name = weight_title
  )

  # Use common tile plot function
  create_tile_plot(
    data = plot_data,
    fill_var = weight,
    scale_type = "gradient2",
    scale_params = scale_params,
    title = paste(weight_title, "Matrix"),
    x_axis_position = x_axis_position,
    node_order = node_order
  )
}

#' @export
plot_nesmr_tile <- function(x, ...) {
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
  node_order = NULL,
  x_axis_position = c("top", "bottom")) {

  weight <- match.arg(weight)
  plot_type <- match.arg(plot_type)
  x_axis_position <- match.arg(x_axis_position)

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


  # TODO: Add a extra row of tiles (like BPG plots) for node layer
  plot_data <- prepare_tile_plot_data(x, node_order, "inclusion_prob", x_axis_position)

  # Create scale parameters
  scale_params <- list(
    low = "white",
    high = "orange",
    name = "Edge Inclusion",
    limits = c(0, 1)
  )

  # Use common tile plot function
  create_tile_plot(
    data = plot_data,
    fill_var = "inclusion_prob",
    scale_type = "gradient",
    scale_params = scale_params,
    title = "Edge Inclusion Probability",
    x_axis_position = x_axis_position,
    node_order = node_order,
    text_format = "%.2f"
  )
}

#' @export
plot_nesmr_tile.nesmr_mh_graph_explore <- function(x, type = c("edge_inc_prob", "best"), ...) {
  # Plot the edge inclusion probability or the best graph
  type <- match.arg(type)
  best_graph <- top_i_graph(x, i = 1)
  if (type == "edge_inc_prob") {
    # First get the ordering from th best graph
    node_order <- unlist(layered_topological_sort(best_graph))
    edge_inclusion_tile_plot.nesmr_tbl_graph(edge_inclusion_probs(x), node_order = node_order)
  } else if (type == "best") {
    # Get the best graph
    plot_nesmr_tile.nesmr_tbl_graph(best_graph, weight = "direct_effect", plot_type = "beta", ...)
  } else {
    stop("Unknown type for nesmr_mh_graph_explore.")
  }
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
  # Get all node names for diagonal elements
  all_nodes <- x %>% tidygraph::activate(nodes) %>% pull(name)
  diag_df <- data.frame(
    from = all_nodes,
    to = all_nodes
  )

  # Prepare plot data with diagonal elements
  plot_data <- x %>%
    tidygraph::activate(edges) %>%
    tidygraph::as_tibble() %>%
    full_join(
      diag_df,
      by = c("from", "to")
    ) %>%
    mutate(
      from = factor(from, levels = node_order),
      to = factor(to, levels = if (x_axis_position == "top") rev(node_order) else node_order),
      # Set diagonal elements to NA for the weight variable
      !!weight_var := ifelse(from == to, NA, !!sym(weight_var))
    )

  return(plot_data)
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
#' @param x_axis_position Position of x-axis ("top" or "bottom")
#' @param text_size Size of text labels on tiles
#' @param text_threshold Threshold for showing text labels
#' @param text_format sprintf format for text labels
#' @param base_size Base font size for theme
#' @param node_order Character vector of node order for grid lines
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
  text_size = 5,
  text_threshold = 0.01,
  text_format = "%.2f",
  base_size = 16,
  node_order = NULL
) {
  scale_type <- match.arg(scale_type)
  x_axis_position <- match.arg(x_axis_position)

  # Create base plot
  p <- ggplot(data, aes(x = to, y = from, fill = !!sym(fill_var))) +
    geom_tile(color = "white", show.legend = TRUE)

  # Add text only if threshold is not Inf
  if (fill_var != "diff_sign_factor" && text_threshold < Inf) {
    p <- p + geom_text(
      aes(label = ifelse(
        abs(!!sym(fill_var)) > text_threshold,
        sprintf(text_format, !!sym(fill_var)), "")),
      size = text_size, color = "black"
    )
  }

  # Add appropriate scale
  if (scale_type == "gradient2") {
    p <- p + do.call(scale_fill_gradient2, c(scale_params, list(na.value = "lightgrey")))
  } else if (scale_type == "gradient") {
    p <- p + do.call(scale_fill_gradient, c(scale_params, list(na.value = "lightgrey")))
  } else if (scale_type == "manual") {
    p <- p + do.call(scale_fill_manual, c(scale_params, list(na.value = "lightgrey")))
  }

  # Add labels and scales
  p <- p +
    labs(title = title, x = x_label, y = y_label) +
    scale_y_discrete(drop = FALSE) +
    scale_x_discrete(drop = FALSE, position = x_axis_position) +
    theme_classic(base_size = base_size) +
    theme(
      axis.text.x = element_text(hjust = if(x_axis_position == "top") 0 else 1)
    ) +
    coord_equal()

  # Add grid lines if node_order is provided
  p <- p +
    geom_vline(xintercept = seq(0.5, length(node_order) + 0.5, by = 1), color = "grey80") +
    geom_hline(yintercept = seq(0.5, length(node_order) + 0.5, by = 1), color = "grey80")

  return(p)
}