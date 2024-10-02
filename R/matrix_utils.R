#' Coverts a matrix into a string with newlines at the rows
#' @B A matrix
#' @collapse A string to collapse the columns
#' @examples
#' B <- matrix(0, nrow = 4, ncol = 4)
#' B[lower.tri(B)] <- 1
#' cat(matrix_to_str(B), "\n")
matrix_to_str <- function(B, collapse = "") {
  paste0(
    apply(B, 1, function(x) {
      paste0(x, collapse = collapse)
    }),
    collapse = "\n")
}
