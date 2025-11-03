#' @title Convert esmr object to tidygraph
#' @description
#' Converts an esmr object to a tidygraph representation, including total and direct effects as well as standard errors for edges.
#'
#' @param param1 An esmr object to be converted.
#' @param param2 Additional parameters for conversion (if applicable).
#' @return A tidygraph object with edge attributes for total/direct effects and their standard errors.
#' @export
as_tbl_graph.nesmr <- function(x, ...) {
  # Check if x is an esmr object or nesmr object
  edgelist <- as.data.frame(
    x$beta[c("beta_k", "beta_j", "fix_beta", "beta_m", "beta_s")]
  )
  colnames(edgelist) <- c("from", "to", "fix_beta", "total_effect", "total_effect_se")

  if (!is.null(x$direct_effect_template)) {
    # If direct effects are available, use them
    #dir_edgelist <- cbind(edgelist, direct_effect = dat$direct_effects[edgelist$beta_j, edgelist$beta_k])
    dir_edgelist <- which(x$direct_effect_template != 0, arr.ind = TRUE)

    # edgelist <- which(x$direct_effect != 0, arr.ind = T)
    dir_edgelist <- data.frame(from = dir_edgelist[, 1],
                        to = dir_edgelist[, 2],
                        direct_effect = x$direct_effects[dir_edgelist],
                        direct_effect_se = x$se_dm[dir_edgelist])
    edgelist <- dplyr::full_join(edgelist, dir_edgelist, by = c("from", "to"))
  }

  nodes = data.frame(name = 1:ncol(x$beta_mat$beta_hat))

  tg <- tidygraph::tbl_graph(nodes = nodes, edges = edgelist)
  tg$elbo <- x$elbo
  class(tg) <- c("nesmr_tbl_graph", "esmr_tbl_graph", class(tg))
  return(tg)
}

#' @export
print.esmr_tbl_graph <- function(x, ...) {
  cat("esmr tidygraph object\n")
  if (!is.null(x$elbo)) {
    cat("Evidence lower bound (ELBO) objective function:", x$elbo, "\n")
  }
  if (!is.null(x$norm_elbo)) {
    cat("Normalized ELBO (Approximate posterior probability):", round(x$norm_elbo, 2), "\n")
  }
  NextMethod("print", x)
}

#' @export
print.nesmr_diff_tbl_graph <- function(x, ...) {
  cat("esmr tidygraph object graph differences\n")
  cat("Difference in evidence lower bound (ELBO) objective function:", x$elbo_diff, "\n")
  if (x$elbo1 > x$elbo2) {
    cat("\t(Graph 1 is better than Graph 2)\n")
  } else if (x$elbo1 < x$elbo2) {
    cat("\t(Graph 2 is better than Graph 1)\n")
  } else {
    cat("\t(Both graphs have the same ELBO)\n")
  }
  NextMethod("print", x)
}

#' @export
`+.esmr_tbl_graph` <- function(g1, g2) {
  stopifnot(inherits(g2, "esmr_tbl_graph"))
  edges1 <- as.data.frame(tidygraph::as_tibble(tidygraph::activate(g1, "edges")))
  edges2 <- as.data.frame(tidygraph::as_tibble(tidygraph::activate(g2, "edges")))
  merged <- dplyr::full_join(edges1, edges2, by = c("from", "to"), suffix = c(".1", ".2"))

  if ("total_effect.1" %in% names(merged) && "total_effect.2" %in% names(merged)) {
    # Sum total effects and standard errors
    merged$total_effect <- rowSums(merged[, c("total_effect.1", "total_effect.2")], na.rm = TRUE)
    merged$total_effect_se <- sqrt(rowSums(merged[, c("total_effect_se.1", "total_effect_se.2")]^2, na.rm = TRUE))
  }

  if ("direct_effect.1" %in% names(merged) && "direct_effect.2" %in% names(merged)) {
    merged$direct_effect <- rowSums(merged[, c("direct_effect.1", "direct_effect.2")], na.rm = TRUE)
    merged$direct_effect_se <- sqrt(rowSums(merged[, c("direct_effect_se.1", "direct_effect_se.2")]^2, na.rm = TRUE))
  }

  # Keep only relevant columns
  edge_cols <- c("from", "to", "direct_effect", "direct_effect_se")
  if ("total_effect" %in% names(merged)) edge_cols <- c(edge_cols, "total_effect", "total_effect_se")
  new_edges <- merged[, edge_cols, drop = FALSE]

  nodes <- as.data.frame(tidygraph::as_tibble(tidygraph::activate(g1, "nodes")))
  tg <- tidygraph::tbl_graph(nodes = nodes, edges = new_edges)
   class(tg) <- c("nesmr_tbl_graph", "esmr_tbl_graph", "tbl_graph", "igraph")
  return(tg)
}

#' @export
`-.esmr_tbl_graph` <- function(g1, g2) {
  stopifnot(inherits(g2, "esmr_tbl_graph"))
  edges1 <- as.data.frame(tidygraph::as_tibble(tidygraph::activate(g1, "edges")))
  edges2 <- as.data.frame(tidygraph::as_tibble(tidygraph::activate(g2, "edges"))) %>%
    dplyr::mutate(
      direct_effect = -direct_effect
    )
  if ("total_effect" %in% names(edges2)) {
    edges2$total_effect <- -edges2$total_effect
  }
  merged <- dplyr::full_join(edges1, edges2, by = c("from", "to"), suffix = c(".1", ".2")) %>%
    dplyr::mutate(
      across(ends_with(".1"), ~tidyr::replace_na(., 0)),
      across(ends_with(".2"), ~tidyr::replace_na(., 0))
    )

  # adj difference + sign difference + weight difference
  if ("direct_effect.1" %in% names(merged) && "direct_effect.2" %in% names(merged)) {
    merged$adj_diff <- as.integer((merged$direct_effect.1 != 0) - (merged$direct_effect.2 != 0))
    merged$sign_diff <- sign(merged$direct_effect.1) - sign(merged$direct_effect.2)

    merged$direct_effect <- rowSums(merged[, c("direct_effect.1", "direct_effect.2")], na.rm = TRUE)
    merged$direct_effect_se <- sqrt(rowSums(merged[, c("direct_effect_se.1", "direct_effect_se.2")]^2, na.rm = TRUE))
  }

  if ("total_effect.1" %in% names(merged) && "total_effect.2" %in% names(merged)) {
    # Sum total effects and standard errors
    merged$total_effect <- rowSums(merged[, c("total_effect.1", "total_effect.2")], na.rm = TRUE)
    merged$total_effect_se <- sqrt(rowSums(merged[, c("total_effect_se.1", "total_effect_se.2")]^2, na.rm = TRUE))
  }

  # Keep only relevant columns
  edge_cols <- c("from", "to", "adj_diff", "sign_diff", "direct_effect", "direct_effect_se")
  if ("total_effect" %in% names(merged)) edge_cols <- c(edge_cols, "total_effect", "total_effect_se")
  new_edges <- merged[, edge_cols, drop = FALSE]

  nodes <- as.data.frame(tidygraph::as_tibble(tidygraph::activate(g1, "nodes")))
  tg <- tidygraph::tbl_graph(nodes = nodes, edges = new_edges)
  tg$elbo1 <- g1$elbo
  tg$elbo2 <- g2$elbo
  tg$elbo_diff <- tg$elbo1 - tg$elbo2
  class(tg) <- c("nesmr_diff_tbl_graph", "nesmr_tbl_graph", "esmr_tbl_graph", "tbl_graph", "igraph")
  return(tg)
}