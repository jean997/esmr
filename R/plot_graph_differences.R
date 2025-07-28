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
  node_order = NULL) {
  effect <- match.arg(effect)
  diff_type <- match.arg(diff_type)
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

    p <- ggplot(graph_diffs, aes(x = to, y = from, fill = effect_diff)) +
      geom_tile(color = "white") +
      geom_text(
        aes(label = ifelse(abs(effect_diff) > 0.01, sprintf("%.2f", effect_diff), "")),
        size = 5, color = "black") +
      scale_fill_gradient2(
        low = neg_color,
        mid = "white",
        high = pos_color,
        midpoint = 0,
        limits = c(-max_diff, max_diff),
        name = paste(tools::toTitleCase(gsub("_", " ", effect)), "Difference"),
        na.value = "lightgrey"
      ) +
      labs(title = paste("Effect Differences Between Graphs (", tools::toTitleCase(gsub("_", " ", effect)), ")"),
           x = "To",
           y = "From")
  } else {
    # For adjacency differences, show discrete edge changes
    p <- ggplot(graph_diffs,
      aes(x = to, y = from, fill = factor(as.character(diff_sign), levels = c("-2", "-1", "0", "1", "2")))) +
      geom_tile(color = "white", show.legend = TRUE) +
      geom_tile(data = subset(graph_diffs, is.na(diff_sign)), fill = "grey") +
      scale_fill_manual(
        values = scale_colors,
        name = "Edge Change",
        labels = c(
          "-2" = "Sign Flip: - → +",
          "-1" = "Edge Added",
          "0" = "No Change",
          "1" = "Edge Removed",
          "2" = "Sign Flip: + → -"
        ),
        drop = FALSE,
        na.value = "lightgrey",
        na.translate = FALSE,
        guide = guide_legend(
          override.aes = list(color = "black", size = 1)
        )
      ) +
      labs(title = "Edge Differences Between Graphs",
           x = "To",
           y = "From")
  }

  # Add common styling elements
  p <- p +
    scale_y_discrete(drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_equal() +
    geom_vline(xintercept = seq(0.5, length(names1) + 0.5, by = 1), color = "grey80") +
    geom_hline(yintercept = seq(0.5, length(names1) + 0.5, by = 1), color = "grey80")

  return(p)
}