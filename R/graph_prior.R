#' Graph prior
#' @k Number of edges
#' @n Total number of traits
#' @pi_0 Parameter to penalize the number of edges. Uniform prior is 0.5, AIC is ~0.73.
#' @return prior distribution on the number of edges
#' @examples
#' # Prior weights for all edges on 5 traits
#' n <- 5
#' total_edges <- n*(n-1)/2
#' k <- seq(0, total_edges)
#' prior <- log_graph_prior(k, n, pi = 0.6)
log_graph_prior <- function(k, n, pi_0 = 0.5, normalize = TRUE) {
  M <- n * (n - 1) / 2
  .f <- function(.k) {
    .k * (log(1 - pi_0) - log(pi_0)) + M * log(pi_0)
  }
  res <- .f(k)
  if (normalize) {
    # TODO: Norm term is not correct but something like this
    # norm_term <- log(pi_0^(M + 1) - (1 - pi_0)^(M + 1)) - log(2*pi - 1)
    # Slower version sums all other values
    all_k <- seq(0, M)
    res <- res - log_sum_exp(.f(all_k))
  }
  return(res)
}
