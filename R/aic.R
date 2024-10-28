#' @export
nesmr_aic <- function(fit, ...) {
  if (fit$is_nesmr) {
    -2 * log_py(fit, ...) + 2 * sum(fit$direct_effect_template)
  }
}
