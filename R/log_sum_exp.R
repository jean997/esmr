#' Compute log(exp(a) + exp(b) + ...) with more numerical stability
#' @param x numeric vector or first element in the sum
#' @param ... additional arguments to sum
#' @return log(exp(a) + exp(b) + ...)
#' @examples
#' log_sum_exp(1, 2, 3)
#' log_sum_exp(c(1, 2, 3))
#' log_sum_exp(c(1,2), 3)
#' # Too big to compute directly, exp(log_liks) = c(Inf, Inf, Inf)
#' log_liks <- c(3000.5, 3001, 3000)
#' log_denom <- log_sum_exp(log_liks)
#' # Posterior probability with uniform prior on the three discrete states
#' exp(log_liks - log_denom)
#'
log_sum_exp <- function(x, ...) {
  x <- unlist(c(x, list(...)))
  x_max <- max(x)
  x_max + log(sum(exp(x - x_max)))
}
