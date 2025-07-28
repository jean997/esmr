#' Layered topological sort
#' @param adj A square adjacency matrix, where [i,j] = 1 means an edge from i -> j
#' @return A list of layers, where each layer is a vector of node names (or indices) that can be processed in parallel
layered_topological_sort <- function(g, names = NULL) {
  if (is.matrix(g)) {
    g <- igraph::graph_from_adjacency_matrix(g, mode = "directed")
  } else if (inherits(g, "tbl_graph")) {
    g <- igraph::as.igraph(g)
  } else if (!inherits(g, "igraph")) {
    stop("Input must be a matrix, tbl_graph, or igraph object.")
  }

  layers <- list()
  remaining <- g

  while (igraph::vcount(remaining) > 0) {
    in_deg <- igraph::degree(remaining, mode = "in")
    # TODO: This will fail if no names
    current_layer <- igraph::V(remaining)[in_deg == 0]$name

    if (length(current_layer) == 0) {
      stop("Cycle detected: topological sort not possible.")
    }

    layers[[length(layers) + 1]] <- current_layer
    remaining <- igraph::delete_vertices(remaining, current_layer)
  }

  return(layers)
}

coordinates_from_layers <- function(layers, x_spacing = 1, y_spacing = 1) {
  coords <- data.frame(name = character(), x = numeric(), y = numeric(), stringsAsFactors = FALSE)

  for (i in seq_along(layers)) {
    layer <- layers[[i]]
    n <- length(layer)

    even_shift_x <- (i %% 2) * 0.25
    even_shift_y <- (n %% 2) * 0.25

    # Spread nodes in layer evenly on y-axis, center around y = 0
    y_positions <- seq(from = -(n - 1) / 2 + even_shift_x - even_shift_y, to = (n - 1) / 2 - even_shift_x + even_shift_y, length.out = n) * y_spacing
    layer_df <- data.frame(
      name = layer,
      x = rep(i * x_spacing, n),
      y = y_positions,
      stringsAsFactors = FALSE
    )

    coords <- rbind(coords, layer_df)
  }

  rownames(coords) <- NULL
  return(coords)
}

layered_topo_with_edges <- function(adj, x_spacing = 1, y_spacing = 1, nice_names = NULL) {
  if (is.null(nice_names)) {
    nice_names <- colnames(adj) %||% as.character(seq_len(ncol(adj)))
  }
  # Step 1: Topo sort to layers
  layers <- layered_topological_sort(adj != 0)
  coords <- coordinates_from_layers(layers, x_spacing = x_spacing, y_spacing = y_spacing)

  coords$name <- nice_names[as.integer(coords$name)]

  # Step 2: Build igraph and extract edges
  #g <- igraph::graph_from_adjacency_matrix(adj, mode = "directed")
  edgelist <- which(adj != 0, arr.ind = T)
  edgelist <- data.frame(name = nice_names[edgelist[, 1]],
                          to = nice_names[edgelist[, 2]],
                          value = adj[edgelist])
  #edges <- igraph::as_data_frame(g, what = "edges") %>%
  #  mutate(from = as.character(from), to = as.character(to))

  # Step 3: Join coordinates for both ends
  edges_coords <- edgelist %>%
    full_join(coords, by = "name") %>%
    left_join(
      rename(coords, xend = x, yend = y),
      by = c("to" = "name")
    )

  list(
    layers = layers,
    node_coords = coords,
    edge_coords = edges_coords
  )
}


#' @export
#' @importFrom ggraph guide_edge_colorbar
plot_layered_topo <- function(
  tg,
  weight = c("direct_effect", "total_effect"),
  plot_type = c("adj", "beta"),
  pos_color = "#d62728",
  neg_color = "#1f77b4",
  scale_limit = NULL,
  scale_factor = 1) {
  weight <- match.arg(weight)
  plot_type <- match.arg(plot_type)
  if (!inherits(tg, "nesmr_tbl_graph")) {
    tg <- tidygraph::as_tbl_graph(tg)
  }
  ts_graph <- layered_topological_sort(tg)
  coords <- coordinates_from_layers(ts_graph)

  tg <- tg %>%
    tidygraph::activate(nodes) %>%
    left_join(coords, by = "name")

  tg <- tg %>%
    tidygraph::activate(edges) %>%
    left_join(coords, by = c("from" = "name"), suffix = c("_from", "_to")) %>%
    left_join(coords, by = c("to" = "name"), suffix = c("_from", "_to"))

  # Calculate scale limits for beta plots
  if (plot_type == "beta" && is.null(scale_limit)) {
    weight_values <- tg %>%
      tidygraph::activate(edges) %>%
      pull(!!sym(weight))
    scale_limit <- max(abs(weight_values), na.rm = TRUE)
  }

  # Create appropriate color variable based on plot type
  if (plot_type == "adj") {
    tg <- tg %>%
      tidygraph::activate(edges) %>%
      mutate(
        colour = factor(
          ifelse(!!sym(weight) > 0, "Positive", "Negative"),
          levels = c("Negative", "Positive")
        )
      )
  } else if (plot_type == "beta") {
    tg <- tg %>%
      tidygraph::activate(edges) %>%
      mutate(colour = !!sym(weight))
  }

  g <- ggraph::ggraph(
    tg, layout = "manual",
    x = x,
    y = y
    ) +
    ggraph::geom_edge_arc(
    data = ~filter(ggraph::get_edges()(.x), abs(x_from - x_to) > 2 | abs(y_from - y_to) > 0.5),
      aes(
        edge_color = colour,
        label = if (plot_type == "beta") sprintf("%.2f", round(!!sym(weight), 2)) else NULL
        ),
      strength = 0.05,
      arrow = grid::arrow(length = grid::unit(5 * scale_factor, "pt"), type = "closed"),
      edge_width = 1.5 * scale_factor,
      label_size = 8 * scale_factor,
      start_cap = ggraph::circle(1 * scale_factor, 'cm'),
      end_cap = ggraph::circle(1 * scale_factor, 'cm'),
      angle_calc = "along",
      force_flip = F,
      check_overlap = T
    ) +
  ggraph::geom_edge_link(
    data = ~filter(ggraph::get_edges()(.x), abs(x_from - x_to) <= 2 & abs(y_from - y_to) <= 0.5),
      aes(
        edge_color = colour,
        label = if (plot_type == "beta") sprintf("%.2f", round(!!sym(weight), 2)) else NULL),
      arrow = grid::arrow(length = grid::unit(5 * scale_factor, "pt"), type = "closed"),
      edge_width = 1.5 * scale_factor,
      label_size = 10 * scale_factor,
      start_cap = ggraph::circle(1 * scale_factor, 'cm'),
      end_cap = ggraph::circle(1 * scale_factor, 'cm'),
      check_overlap = T) +
    ggraph::geom_node_point(size = 3 * scale_factor) +
    ggraph::geom_node_label(aes(label = name), size = 8 * scale_factor) +
    theme(legend.position = c(0.98, 0.02), legend.justification = c("right", "bottom")) +
    theme_void(base_size = 30 * scale_factor)

  # Add appropriate color scale based on plot type
  if (plot_type == "adj") {
    g <- g + ggraph::scale_edge_color_manual(
      values = c("Negative" = neg_color, "Positive" = pos_color),
      name = "Effect Sign"
    )
  } else if (plot_type == "beta") {
    g <- g + ggraph::scale_edge_color_gradient2(
      low = neg_color,
      mid = "white",
      high = pos_color,
      midpoint = 0,
      limits = c(-scale_limit, scale_limit),
      name = "Beta"
    ) +
    guides(edge_color = guide_edge_colorbar(barheight = 10 * scale_factor, barwidth = 0.5 * scale_factor))
  }
  return(g)
}