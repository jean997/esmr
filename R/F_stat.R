#' F-statistic for instrument strength
#' @param x An esmr object
#' @export
F_stat <- function(x) {
  if (!is(x, "esmr")) {
    stop("x must be an esmr object")
  }

  colSums(x$l$lbar^2)
}