#' @export
nesmr_aic <- function(fit, ...) {
  if ("B_template" %in% names(fit)) {
    -2 * log_py(fit, ...) + 2 * sum(fit$B_template)
  }
}
