#' Plot Differences Between Two Graphs
#'
#' This function visualizes the differences between two graphs, either by adjacency matrices or beta coefficients.
#'
#' @param tg1 first graph object
#' @param tg2 second graph object
#' @param diff_type Character vector specifying the type of difference to plot. Options are \code{"adj"} for adjacency matrix differences or \code{"beta"} for beta coefficient differences.
#'
#' @return A plot visualizing the differences between the two graphs.
#' @export
plot_graph_differences <- function(
  tg1, tg2, diff_type = c("adj", "beta"),
  effect = c("direct_effect", "total_effect"),
  pos_color = "#d62728",
  neg_color = "#1f77b4",
  scale_colors = c(
    `-2` = "#D7191C", `-1` = "#FDAE61",
    `0` = "white", `1` = "#ABDDA4", `2` = "#2B83BA"),
  node_order = NULL,
  x_axis_position = c("top", "bottom")) {
  effect <- match.arg(effect)
  diff_type <- match.arg(diff_type)
  x_axis_position <- match.arg(x_axis_position)
  names1 <- tg1 %>% tidygraph::activate(nodes) %>% pull(name)
  names2 <- tg2 %>% tidygraph::activate(nodes) %>% pull(name)
  stopifnot(identical(names1, names2))

  if (is.null(node_order)) {
    node_order <- as.integer(igraph::topo_sort(tg1))
  }

  diag_df <- data.frame(
    from = names1,
    to = names1
  )

  graph_diffs <- tg1 %>%
    tidygraph::activate(edges) %>%
    tidygraph::as_tibble() %>%
    select(from, to, effect = !!effect) %>%
    full_join(
      tg2 %>% tidygraph::activate(edges) %>% tidygraph::as_tibble() %>%
        select(from, to, effect = !!effect),
      by = c("from", "to"),
      suffix = c("_tg1", "_tg2")
    ) %>%
    mutate(
      across(ends_with("_tg1"), ~tidyr::replace_na(., 0)),
      across(ends_with("_tg2"), ~tidyr::replace_na(., 0))
    ) %>%
    full_join(
      diag_df,
      by = c("from", "to")
    ) %>%
    mutate(
      effect_diff = effect_tg1 - effect_tg2,
      sign_diff = sign(effect_diff),
      diff_sign = sign(effect_tg1) - sign(effect_tg2),
      inc_sign = (effect_tg1 != 0) - (effect_tg2 != 0),
      from = factor(from, levels = node_order),
      to = factor(to, levels = rev(node_order)),
      effect_diff = ifelse(from == to, NA, effect_diff),
      diff_sign = ifelse(from == to, NA, diff_sign)
    )

  # TODO: Different plots for different diff types
  # tg_diff <- tbl_graph(
  #   edges = graph_diffs,
  #   nodes = data.frame(name = names1, ix = seq_along(names1)))

  if (diff_type == "beta") {
    # For beta differences, show the actual effect differences with gradient
    max_diff <- max(abs(graph_diffs$effect_diff), na.rm = TRUE)

    # Create scale parameters for gradient2
    scale_params <- list(
      low = neg_color,
      mid = "white",
      high = pos_color,
      midpoint = 0,
      limits = c(-max_diff, max_diff),
      name = paste(tools::toTitleCase(gsub("_", " ", effect)), "Difference")
    )

    # Use common tile plot function
    p <- create_tile_plot(
      data = graph_diffs,
      fill_var = "effect_diff",
      scale_type = "gradient2",
      scale_params = scale_params,
      title = paste("Effect Differences Between Graphs (", tools::toTitleCase(gsub("_", " ", effect)), ")"),
      x_axis_position = x_axis_position,
      node_order = names1
    )
  } else {
    # For adjacency differences, show discrete edge changes
    all_levels <- c("-2", "-1", "0", "1", "2")
    scale_params <- list(
      values = scale_colors,
      name = "Edge Change",
      labels = c(
        "-2" = "Sign Flip: - → +",
        "-1" = "Edge Removed",
        "0" = "No Change",
        "1" = "Edge Added",
        "2" = "Sign Flip: + → -"
      ),
      breaks = all_levels,
      limits = all_levels,
      drop = FALSE,
      guide = guide_legend(override.aes = list(color = "black", size = 1))
    )

    # Prepare data with factor levels - ensure all levels are represented
    plot_data <- graph_diffs %>%
      mutate(diff_sign_factor = factor(as.character(inc_sign), levels = all_levels))    # Use common tile plot function

    p <- create_tile_plot(
      data = plot_data,
      fill_var = "diff_sign_factor",
      scale_type = "manual",
      scale_params = scale_params,
      title = "Edge Differences Between Graphs",
      x_axis_position = x_axis_position,
      node_order = names1,
      text_threshold = Inf  # Don't show text for manual scale
    )
  }

  return(p)
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