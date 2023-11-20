#' @export
logLik.esmr <- function(x) {
  with(x, log_py(Y, l$g_hat, f$fgbar, omega))
}
