make_lower_tri_gridlines <- function(n_nodes) {
  n_lines <- n_nodes + 1
    grid_lines <- bind_rows(
      data.frame(
        x = (seq_len(n_lines) + 0.5) - 1,
        xend = (seq_len(n_lines) + 0.5) - 1,
        y = (n_nodes - seq_len(n_lines) + 0.5) + 1,
        yend = 0.5
      ),
      data.frame(
        x = 0.5,
        xend = seq_len(n_nodes) + 0.5,
        y = n_nodes - seq_len(n_nodes) + 0.5,
        yend = n_nodes - seq_len(n_nodes) + 0.5
      )
    )

    geom_segment(
        data = grid_lines,
        inherit.aes = FALSE,
        aes(x = x, xend = xend, y = y, yend = yend),
        color = "grey80"
    )
}

make_gridlines_full <- function(n_nodes) {
  n_lines <- n_nodes + 1
  grid_lines <- bind_rows(
      data.frame(
        x = (seq_len(n_lines) + 0.5) - 1,
        xend = (seq_len(n_lines) + 0.5) - 1,
        y = n_lines - 0.5,
        yend = 0.5
      ),
      data.frame(
        x = 0.5,
        xend = n_nodes + 0.5,
        y = (seq_len(n_lines) + 0.5) - 1,
        yend = (seq_len(n_lines) + 0.5) - 1
      )
    )
  geom_segment(
        data = grid_lines,
        inherit.aes = FALSE,
        aes(x = x, xend = xend, y = y, yend = yend),
        color = "grey80"
    )
}

#' @importFrom ggplot2 ggplot aes geom_tile geom_text geom_line geom_point geom_step
#' @importFrom ggplot2 geom_segment geom_vline geom_hline facet_grid labs theme theme_classic
#' @importFrom ggplot2 element_text scale_fill_gradient scale_fill_gradient2 scale_fill_manual
#' @importFrom ggplot2 scale_color_manual scale_x_discrete scale_y_discrete scale_y_continuous
#' @importFrom ggplot2 coord_equal guide_legend
NULL

## Helpers shared by the tile, graph, and MH-chain plotting functions -----------------------

tbl_node_names <- function(x) {
  tidygraph::activate(x, nodes) |> pull(name)
}

## Replace the node names of a tbl_graph. `labels` is either positional (one per node) or
## named by the current node names.
relabel_nodes <- function(x, labels = NULL) {
  if (is.null(labels)) return(x)
  nm <- as.character(tbl_node_names(x))
  stopifnot(length(labels) == length(nm))
  if (!is.null(names(labels)) && all(nm %in% names(labels))) labels <- labels[nm]
  new_names <- unname(labels)
  tidygraph::activate(x, nodes) |> mutate(name = new_names)
}

## Node order as a character vector of node names. A numeric `node_order` is taken as
## positions in the graph; the default is the layered topological order.
resolve_node_order <- function(x, node_order = NULL) {
  nm <- as.character(tbl_node_names(x))
  if (is.null(node_order)) return(as.character(unlist(layered_topological_sort(x))))
  if (is.numeric(node_order)) nm[node_order] else as.character(node_order)
}

## Positions of the layered topological order of `x`, used to order every plot of an
## object by its best graph.
best_node_order <- function(x) {
  best <- top_i_graph(x, i = 1)
  match(resolve_node_order(best), as.character(tbl_node_names(best)))
}

move_type_colors <- c(add = "blue", remove = "orange", flip = "green", init = "grey40")

## Chain info collapsed to the first visit of each proposed graph, in the order they were
## fit, with running best / total ELBO.
##   iter      - row of the proposal in mh_chain_info (across all chains)
##   graph_fit - index of the proposal among unique graphs
##   new_best  - whether the graph had the highest ELBO so far
mh_new_graph_chain <- function(x, labels = NULL) {
  info <- x$mh_chain_info
  if (is.null(labels)) labels <- seq_len(sqrt(nchar(info$prop_graph[[1]])))
  edge_node <- function(k) {
    labels[as.integer(vapply(strsplit(info$mod_edge, "|", fixed = TRUE), `[`, character(1), k))]
  }
  info |>
    mutate(
      iter = row_number(),
      from = edge_node(1),
      to = edge_node(2),
      move_type = case_when(
        is.na(move_type) & insert_edge ~ "add",
        is.na(move_type) & !insert_edge ~ "remove",
        TRUE ~ move_type
      )
    ) |>
    filter(!duplicated(prop_graph), is.finite(prop_elbo)) |>
    mutate(
      graph = prop_graph,
      elbo = prop_elbo,
      graph_fit = row_number(),
      cum_best_elbo = cummax(elbo),
      new_best = elbo == cum_best_elbo,
      cum_total_elbo = max(elbo) + log(cumsum(exp(elbo - max(elbo))))
    )
}
