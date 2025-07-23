#' @export
plot_nesmr_network <- function(
  adj_mat, nice_names = NULL, layout = "stress",
  pos_col = "#EE220C",
  neg_col = "#00A2FF",
  EIP = NULL) {
# Probably want to time order using:
# https://github.com/r-causal/ggdag/blob/535386358b86db3e713eacd7dbbd99e8498acd67/R/layouts.R#L69
  if (is.null(nice_names)) {
    nice_names <- colnames(adj_mat) %||% as.character(seq_len(ncol(adj_mat)))
  }

  edgelist <- which(adj_mat != 0, arr.ind = T)
  edgelist <- data.frame(name = edgelist[, 1],
                          to = edgelist[, 2],
                          value = adj_mat[edgelist])
  if (!is.null(EIP)) {
    edgelist <- edgelist %>%
      left_join(EIP, by = c("name" = "from", "to" = "to"))
  }

  edgelist$name <- nice_names[as.integer(edgelist$name)]
  edgelist$to <- nice_names[as.integer(edgelist$to)]

  g <- edgelist %>%
    as_tidy_dagitty(layout = layout) %>%
    mutate(color = ifelse(value > 0, "pos", "neg")) %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_dag_point(col = "white") +
      geom_dag_edges(
        aes(edge_colour = color, edge_alpha = inclusion_prob, label = round(inclusion_prob, 2)), show.legend = TRUE) +
      geom_dag_label() +
      ggraph::theme_graph() + ggraph::scale_edge_colour_manual(
        name = "Edge sign",
        values = c("pos" = pos_col, "neg" = neg_col),
        breaks = c("pos", "neg"),
        labels = c("Positive", "Negative"),
        guide = guide_legend(override.aes = list(edge_width = 2))
      )

  g
}



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

plot_layered_topo <- function(tg, plot_type = "adj") {
  ts_graph <- layered_topological_sort(tg)
  coords <- coordinates_from_layers(ts_graph)

  if (is.matrix(tg)) {
    edgelist <- which(g_adj != 0, arr.ind = T)
    edgelist <- data.frame(from = edgelist[, 1],
                          to = edgelist[, 2],
                          from_name = nice_names[edgelist[, 1]],
                          to_name = nice_names[edgelist[, 2]])
    edgelist[[plot_type]] <- g_adj[edgelist]

    tg <- tbl_graph(edges = edgelist, nodes = data.frame(name = nice_names, ix = seq_along(nice_names)))
  }

  tg <- tg %>%
    activate(nodes) %>%
    left_join(coords)

  tg <- tg %>%
    activate(edges) %>%
    mutate(
      `adj` = sign(zapsmall(beta)),
    ) %>%
    left_join(coords, by = c("from_name" = "name"), suffix = c("_from", "_to")) %>%
    left_join(coords, by = c("to_name" = "name"), suffix = c("_from", "_to"))

  ggraph(
    tg, layout = "manual",
    x = x,
    y = y
    ) +
    geom_edge_arc(
      aes(
        filter = abs(x_from - x_to) > 2 | abs(y_from - y_to) > 0.5,
        colour = get(plot_type) > 0),
      strength = 0.05,
      arrow = grid::arrow(length = grid::unit(5, "pt"), type = "closed"),
      start_cap = circle(1, 'cm'),
      end_cap = circle(1, 'cm'),
      angle_calc = "along",
      force_flip = F,
      check_overlap = T
    ) +
  geom_edge_link(
      aes(
        filter = abs(x_from - x_to) <= 2 & abs(y_from - y_to) <= 0.5,
        colour = get(plot_type) > 0),
      arrow = grid::arrow(length = grid::unit(5, "pt"), type = "closed"),
      start_cap = circle(1, 'cm'),
      end_cap = circle(1, 'cm'),
      check_overlap = T) +
      geom_node_point() +
    geom_node_label(aes(label = name)) +
    theme_dag()
}


plot_graph_differences <- function(tg1, tg2, diff_type = c("adj", "beta")) {
  names1 <- tg1 %>% tidygraph::activate(nodes) %>% pull(name)
  names2 <- tg2 %>% tidygraph::activate(nodes) %>% pull(name)
  stopifnot(identical(names1, names2))
  #print(tidygraph::activate(tg1, edges))
  #print(tidygraph::activate(tg2, edges))

  graph_diffs <- tg1 %>%
    activate(edges) %>%
    as_tibble() %>%
    select(from_name, to_name, beta) %>%
    full_join(
      tg2 %>% activate(edges) %>% as_tibble(),
      by = c("from_name", "to_name"),
      suffix = c("_tg1", "_tg2")
    ) %>%
    mutate(
      across(ends_with("_tg1"), ~replace_na(., 0)),
      across(ends_with("_tg2"), ~replace_na(., 0))
    ) %>%
    mutate(
      beta_diff = beta_tg1 - beta_tg2,
      sign_diff = sign(beta_diff),
      diff_sign = sign(beta_tg1) - sign(beta_tg2),
      inc_sign = (beta_tg1 != 0) - (beta_tg2 != 0)
    )

  # TODO: Different plots for different diff types
  # tg_diff <- tbl_graph(
  #   edges = graph_diffs,
  #   nodes = data.frame(name = names1, ix = seq_along(names1)))

  ggplot(graph_diffs, aes(x = to_name, y = from_name, fill = factor(inc_sign))) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c("-1" = "red", "0" = "white", "1" = "blue"),
    name = "Difference",
    labels = c("-1" = "Removed", "0" = "No Change", "1" = "Added")
  ) +
  scale_x_discrete(limits = names1) +  # Ensure x-axis is in the same order as names1
  scale_y_discrete(limits = rev(names1)) +  # Reverse y-axis for matrix ordering
  labs(title = "Edge Differences Between Top Two Graphs",
       x = "To",
       y = "From") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_equal() +
  geom_vline(xintercept = seq(0.5, length(nice_names) + 0.5, by = 1), color = "grey80") +
  geom_hline(yintercept = seq(0.5, length(nice_names) + 0.5, by = 1), color = "grey80")
}