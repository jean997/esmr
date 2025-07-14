plot_nesmr_network <- function(x) {
# Probably want to time order using:
# https://github.com/r-causal/ggdag/blob/535386358b86db3e713eacd7dbbd99e8498acd67/R/layouts.R#L69
}

#' Layered topological sort
#' @param adj A square adjacency matrix, where [i,j] = 1 means an edge from i -> j
#' @return A list of layers, where each layer is a vector of node names (or indices) that can be processed in parallel
layered_topological_sort <- function(adj) {
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "directed")
  igraph::V(g)$name <- if (is.null(colnames(adj))) as.character(1:nrow(adj)) else colnames(adj)

  layers <- list()
  remaining <- g

  while (igraph::vcount(remaining) > 0) {
    in_deg <- igraph::degree(remaining, mode = "in")
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

    even_shift <- (i %% 2) * 0.5

    # Spread nodes in layer evenly on y-axis, center around y = 0
    y_positions <- seq(from = -floor((n - 1) / 2) + even_shift, to = ceil((n - 1) / 2) - even_shift, length.out = n) * y_spacing
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