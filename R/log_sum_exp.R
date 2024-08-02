log_sum_exp <- function(x, ...) {
  x <- unlist(c(x, list(...)))
  x_max <- max(x)
  x_max + log(sum(exp(x - x_max)))
}
