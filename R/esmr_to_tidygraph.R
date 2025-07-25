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
  cat("Evidence lower bound (ELBO) objective function:", x$elbo, "\n")
  if (!is.null(x$norm_elbo)) {
    cat("Normalized ELBO (Approximate posterior probability):", round(x$norm_elbo, 2), "\n")
  }
  NextMethod("print", x)
}
