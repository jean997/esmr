#' Generate coordinates for layered topological sort layout
#'
#' @export
#' @param g A graph object (matrix, tbl_graph, or igraph)
#' @param names Optional character vector of node names (if g is a matrix)
layered_topo_sort_coords <- function(g, names = NULL, x_spacing = 1, y_spacing = 1) {
  layers <- layered_topological_sort(g, names)
  coords <- coordinates_from_layers(layers, x_spacing = x_spacing, y_spacing = y_spacing)
  return(coords)
}

#' Layered topological sort
#' @param adj A square adjacency matrix, where [i,j] = 1 means an edge from i -> j
#' @return A list of layers, where each layer is a vector of node names (or indices) that can be processed in parallel
layered_topological_sort <- function(g, names = NULL) {
  # Accept adjacency matrix, tbl_graph, igraph, or logical adjacency (matrix != 0)
  if (is.matrix(g)) {
    adj <- (g != 0) * 1
    node_names <- colnames(g)
    if (is.null(node_names)) node_names <- as.character(seq_len(ncol(g)))
  } else if (inherits(g, "tbl_graph")) {
    g_ig <- igraph::as.igraph(g)
    adj <- as.matrix(igraph::as_adj(g_ig, sparse = FALSE))
    node_names <- igraph::V(g_ig)$name
    if (is.null(node_names)) node_names <- as.character(seq_len(ncol(adj)))
  } else if (inherits(g, "igraph")) {
    adj <- as.matrix(igraph::as_adj(g, sparse = FALSE))
    node_names <- igraph::V(g)$name
    if (is.null(node_names)) node_names <- as.character(seq_len(ncol(adj)))
  } else {
    stop("Input must be a matrix, tbl_graph, or igraph object.")
  }

  # Kahn's algorithm for layered/topological sort
  n <- ncol(adj)
  colnames(adj) <- node_names
  rownames(adj) <- node_names

  in_deg <- colSums(adj != 0)
  remaining <- rep(TRUE, n)
  names(remaining) <- node_names

  layers <- list()
  while (any(remaining)) {
    zero_in <- names(which(in_deg == 0 & remaining))
    if (length(zero_in) == 0) {
      # cycle detected; return partial layers and then remaining nodes as a final layer
      warning("Cycle detected: returning partial order with remaining nodes grouped in final layer.")
      rem_nodes <- names(remaining[remaining])
      layers[[length(layers) + 1]] <- rem_nodes
      break
    }
    # add current layer
    layers[[length(layers) + 1]] <- zero_in
    # remove these nodes from remaining and decrement in-degrees of their targets
    for (v in zero_in) {
      remaining[v] <- FALSE
      # for each outgoing edge v -> u, decrement in_deg[u]
      outs <- which(adj[v, ] != 0)
      if (length(outs) > 0) {
        out_names <- colnames(adj)[outs]
        in_deg[out_names] <- in_deg[out_names] - 1
      }
    }
  }

  # If names parameter was NULL and node_names are numeric-like, convert layers to numeric indices
  maybe_numeric <- suppressWarnings(as.numeric(node_names))
  if (!is.null(names) || (!any(is.na(maybe_numeric)) && all(node_names == as.character(maybe_numeric)))) {
    layers <- lapply(layers, function(x) as.integer(x))
  }

  return(layers)
}

coordinates_from_layers <- function(layers, x_spacing = 1, y_spacing = 1) {
  coords <- data.frame(name = character(), x = numeric(), y = numeric(), stringsAsFactors = FALSE)

  # Determine maximum layer size to allow consistent vertical spacing
  max_n <- if (length(layers) > 0) max(vapply(layers, length, integer(1))) else 0

  for (i in seq_along(layers)) {
    layer <- layers[[i]]
    n <- length(layer)
    if (n == 0) next

    # Base horizontal position for this layer (0-indexed)
    x_base <- (i - 1) * x_spacing

    # Spread nodes within layer horizontally to reduce parallel edge overlap
    intra_frac <- 0.6
    if (n == 1) {
      x_offsets <- 0
    } else {
      # offsets centered around 0, scaled by intra_frac * x_spacing
      x_offsets <- seq(from = -(n - 1) / 2, to = (n - 1) / 2, length.out = n) / max(1, n - 1) * (x_spacing * intra_frac)
    }

    # Vertical positions: center nodes around y = 0 with spacing y_spacing
    if (n == 1) {
      y_positions <- 0
    } else {
      y_positions <- seq(from = (n - 1) / 2, to = - (n - 1) / 2, length.out = n) * y_spacing
    }

    layer_df <- data.frame(
      name = layer,
      x = x_base + x_offsets,
      y = y_positions,
      stringsAsFactors = FALSE
    )

    coords <- rbind(coords, layer_df)
  }

  # Center horizontally so graph is symmetric around x = 0
  if (nrow(coords) > 0) {
    x_range <- range(coords$x)
    x_center <- mean(x_range)
    coords$x <- coords$x - x_center
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
plot_layered_topo <- function(
  tg,
  weight = c("direct_effect", "total_effect"),
  plot_type = c("adj", "beta"),
  pos_color = "#d62728",
  neg_color = "#1f77b4",
  scale_limit = NULL,
  scale_factor = 1,
  coords = NULL,
  ...) {
    #ts_graph <- layered_topological_sort(tg)
    #coords <- coordinates_from_layers(ts_graph)

    plot_nesmr_graph(
      tg,
      weight = weight,
      plot_type = plot_type,
      pos_color = pos_color,
      neg_color = neg_color,
      scale_limit = scale_limit,
      scale_factor = scale_factor,
      coords = coords
      #coords = coords
    )
}
#' @export
#' @importFrom ggraph guide_edge_colorbar
plot_nesmr_graph <- function(
  tg,
  weight = c("direct_effect", "total_effect"),
  plot_type = c("adj", "beta"),
  pos_color = "#d62728",
  neg_color = "#1f77b4",
  scale_limit = NULL,
  scale_factor = 1,
  coords = NULL) {
  weight <- match.arg(weight)
  plot_type <- match.arg(plot_type)
  if (!inherits(tg, "nesmr_tbl_graph")) {
    tg <- tidygraph::as_tbl_graph(tg)
  }

  if (is.null(coords)) {
    ts_graph <- layered_topological_sort(tg)
    coords <- coordinates_from_layers(ts_graph)
  }

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
  ggraph::geom_edge_fan(
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
#    ggraph::geom_node_point(size = 3 * scale_factor) +
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