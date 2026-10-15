#' Diagnostic plots for an MH graph exploration
#'
#' `plot_mh_trace()` plots a summary of the unique graphs fit over the course of the chain.
#' `plot_mh_best_path()` plots how each new best graph compares to the overall best graph.
#' `plot_mh_accept_rate()` plots the acceptance and new graph rates over the chain.
#' See also `plot()` for the normalized ELBO of the top graphs and
#' `plot_nesmr_tile(type = "proposal_freq")` for edge proposal frequencies.
#'
#' @param x A `nesmr_mh_graph_explore` object.
#' @param type What to plot against the chain: the ELBO (`"elbo"`), exp(ELBO)
#'   (`"exp_elbo"`) relative to the best graph, the cumulative best ELBO
#'   (`"cum_best_elbo"`), the cumulative total exp(ELBO) explored (`"cum_total_elbo"`),
#'   or the iterations that found a new best graph (`"new_best"`).
#' @param x_axis Iteration of the chain (`"iter"`) or number of unique graphs fit
#'   (`"graph_fit"`).
#' @param y For `plot_mh_best_path()`, the distance of each new best graph to the overall best:
#'   the edit distance of the graphs, the ELBO difference, or the Frobenius norm of the
#'   difference in estimated direct effects.
#' @param labels Optional labels for the nodes in the move annotations.
#' @param label_gain Annotate the moves that improved the ELBO by more than this. Requires
#'   ggrepel; `NULL` for no annotations.
#' @param window Window of the rolling rates.
#' @return A ggplot object.
#' @name plot_mh_chain
NULL

#' @rdname plot_mh_chain
#' @export
plot_mh_trace <- function(
  x,
  type = c("elbo", "exp_elbo", "cum_best_elbo", "cum_total_elbo", "new_best"),
  x_axis = c("iter", "graph_fit")) {
  type <- match.arg(type)
  x_axis <- match.arg(x_axis)

  d <- mh_new_graph_chain(x)
  d$x_val <- d[[x_axis]]
  d$y_val <- switch(type,
    elbo = d$elbo - max(d$elbo),
    exp_elbo = exp(d$elbo - max(d$elbo)),
    cum_best_elbo = d$cum_best_elbo - max(d$cum_best_elbo),
    cum_total_elbo = exp(d$cum_total_elbo - max(d$cum_total_elbo)),
    new_best = as.numeric(d$new_best)
  )
  y_label <- c(
    elbo = "ELBO - max(ELBO)",
    exp_elbo = "exp(ELBO - max(ELBO))",
    cum_best_elbo = "Cumulative best ELBO - max(ELBO)",
    cum_total_elbo = "exp(cum-ELBO - max(cum-ELBO))",
    new_best = "Found new best graph?"
  )[[type]]
  title <- c(
    elbo = "ELBO chain",
    exp_elbo = "ELBO chain",
    cum_best_elbo = "Cumulative best ELBO discovered",
    cum_total_elbo = "Cumulative total exp(ELBO) explored",
    new_best = "New best graph discovered"
  )[[type]]
  x_label <- c(iter = "MH Iteration", graph_fit = "Graph Fit")[[x_axis]]

  p <- ggplot(d, aes(x = x_val, y = y_val))
  p <- switch(type,
    elbo = {
      p <- p + geom_line()
      if (x_axis == "iter") {
        # Accepted state at each iteration, behind the proposed graphs
        accepted <- x$mh_chain_info |>
          mutate(iter = row_number(), elbo = ifelse(accepted %in% 1, prop_elbo, curr_elbo)) |>
          filter(is.finite(elbo))
        p <- p + geom_point(data = accepted, aes(x = iter, y = elbo - max(d$elbo)), alpha = 0.3)
      }
      p + geom_point(aes(color = move_type))
    },
    exp_elbo = p + geom_point(aes(color = move_type)),
    cum_best_elbo = p + geom_step() +
      geom_point(data = d[d$new_best, ], aes(color = move_type)),
    cum_total_elbo = p + geom_step() +
      geom_point(data = d[d$new_best, ], aes(color = move_type)),
    new_best = p + geom_segment(aes(xend = x_val, y = 0, yend = y_val)) +
      geom_point() +
      scale_y_continuous(breaks = c(0, 1), labels = c("No", "Yes"))
  )

  if (type != "new_best") p <- p + scale_color_manual(values = move_type_colors, name = "Move type")
  p +
    labs(title = title, x = x_label, y = y_label) +
    theme_classic() +
    theme(legend.position = "bottom")
}

#' @rdname plot_mh_chain
#' @export
plot_mh_best_path <- function(
  x,
  y = c("edit_distance", "elbo_diff", "frob_norm"),
  labels = NULL,
  label_gain = 10) {
  y <- match.arg(y)

  chain <- mh_new_graph_chain(x, labels)
  best_graph <- chain$graph[which.max(chain$elbo)]
  best_beta <- x$visited_graphs[[best_graph]]$beta_hat

  d <- chain |>
    filter(new_best) |>
    mutate(
      edit_distance = as.numeric(utils::adist(graph, best_graph)),
      elbo_diff = elbo - max(elbo),
      frob_norm = vapply(
        graph, function(g) sqrt(sum((x$visited_graphs[[g]]$beta_hat - best_beta)^2)),
        numeric(1)),
      elbo_gain = elbo - lag(elbo)
    )
  d$y_val <- d[[y]]
  last_fit <- max(chain$graph_fit)
  y_label <- c(
    edit_distance = "Edit distance",
    elbo_diff = "ELBO - max(ELBO)",
    frob_norm = "Frobenius norm to overall best graph"
  )[[y]]

  p <- ggplot(d, aes(x = graph_fit, y = y_val)) +
    geom_step() +
    # The last best graph stays the best for the rest of the chain
    geom_segment(
      data = d[nrow(d), ],
      aes(xend = last_fit, yend = y_val),
      linetype = "dashed"
    ) +
    geom_point(aes(color = move_type)) +
    scale_color_manual(values = move_type_colors, name = "New best graph move type") +
    labs(
      title = "Current best to overall best graph",
      x = "Graph Fit",
      y = y_label
    ) +
    theme_classic() +
    theme(legend.position = "bottom")

  if (!is.null(label_gain) && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = d[!is.na(d$elbo_gain) & d$elbo_gain > label_gain, ],
      aes(label = paste0(move_type, ": ", from, "->", to, "\nELBO gain: ", round(elbo_gain))),
      nudge_y = 1
    )
  }
  p
}

#' @rdname plot_mh_chain
#' @export
plot_mh_accept_rate <- function(x, window = 50) {
  info <- x$mh_chain_info |>
    mutate(
      iter = row_number(),
      new_graph = !duplicated(prop_graph),
      chain_start = chain != lag(chain, default = chain[[1]])
    ) |>
    filter(!is.na(accepted))

  rolling_mean <- function(v) as.numeric(stats::filter(v, rep(1 / window, window), sides = 1))
  metrics <- c(
    "Cumulative\nAcceptance\nRate", "Cumulative\nNew Graph\nRate",
    sprintf("Rolling\nAcceptance\nRate (window=%d)", window),
    sprintf("Rolling\nNew Graph\nRate (window=%d)", window)
  )
  rates <- bind_rows(
    data.frame(iter = info$iter, metric = metrics[1], value = cumsum(info$accepted) / seq_len(nrow(info))),
    data.frame(iter = info$iter, metric = metrics[2], value = cumsum(info$new_graph) / seq_len(nrow(info))),
    data.frame(iter = info$iter, metric = metrics[3], value = rolling_mean(info$accepted)),
    data.frame(iter = info$iter, metric = metrics[4], value = rolling_mean(info$new_graph))
  )
  rates$metric <- factor(rates$metric, levels = metrics)

  ggplot(rates, aes(x = iter, y = value, color = metric)) +
    facet_grid(metric ~ ., scales = "free_y") +
    geom_vline(
      data = info[info$chain_start, ], aes(xintercept = iter),
      linetype = "dashed", color = "grey"
    ) +
    geom_line() +
    labs(
      title = "Acceptance and new graph discovery rates over MH iterations",
      x = "MH Iteration",
      y = "Rates"
    ) +
    theme_classic() +
    theme(legend.position = "none")
}
