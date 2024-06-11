#' Computes direct effects and fills in constrained total effects
#'
#' O(n^3)
#'
#' @param total_effects Total effect matrix (lower triangular).
#' @param s A matrix of indices with columns: row, col such that the direct effects B
#' of B[row, col] = 0
#'
#' @return list of direct effects B and total effects filled in
complete_T <- function(total_effects, s) {
  total_effects_complete <- total_effects
  n <- nrow(total_effects)
  X <- diag(n)
  X_inv <- (diag(n) + total_effects)
  # Order in order that we are iterating
  s <- s[order(s[,'row'], -s[,'col']),, drop = FALSE]
  s_idx <- 1
  s_max <- nrow(s)

  # Iterate rows
  for (i in 2:n) {
    # Iterate cols right to left
    for (j in (i - 1):1) {
      if ((j + 1) <= (i - 1)) {
        partial_sum <- - sum(
          X_inv[i, (j + 1):(i - 1)] *
            X[(j + 1):(i - 1), j]
        )
      } else {
        partial_sum <- 0
      }

      if (s_idx <= s_max && all(s[s_idx, ] == c(i,j))) {
        # Fill in t_ij
        s_idx <- s_idx + 1
        X[i,j] <- 0
        X_inv[i,j] <- partial_sum
      } else {
        X[i,j] = - X_inv[i,j] + partial_sum
      }
    }
  }

  return(
    list(
      B = diag(n) - X,
      total_effects = X_inv - diag(n)
    )
  )
}
