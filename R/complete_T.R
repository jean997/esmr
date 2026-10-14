#' Computes direct effects and fills in constrained total effects
#'
#' Does not depend on the ordering of the traits. Writing `T = (I - B)^-1 - I`,
#' the direct effects `B` are supported on every off-diagonal position not in
#' `s` and satisfy `B = T - (T - B)` on that support, where `T - B` is the
#' contribution of paths of length >= 2, which only depends on `B`. Iterating
#' `B <- (total_effects - (T(B) - B))` on the support adds one path length
#' per pass, so for a DAG it converges exactly in at most `n - 1` passes.
#' O(n^4) worst case, but it stops as soon as `B` stops changing.
#'
#' @param total_effects Total effect matrix. Entries in `s` are ignored, all
#' other off-diagonal entries are treated as known.
#' @param s A matrix of indices with columns: row, col such that the direct effects B
#' of B[row, col] = 0
#'
#' @return list of direct effects B and total effects filled in
complete_T <- function(total_effects, s) {
  total_effects_complete <- total_effects

  # Reorder total effects to be lower triangular
  topo_order <- esmr:::topo_sort_mat(total_effects_complete)
  topo_order_inv <- order(topo_order)
  total_effects_complete <- total_effects_complete[topo_order, topo_order]
  s <- s[topo_order_inv[s[, 'row']], topo_order_inv[s[, 'col']], drop = FALSE]

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
      B = (diag(n) - X)[topo_order_inv, topo_order_inv],
      total_effects = (X_inv - diag(n))[topo_order_inv, topo_order_inv]
    )
  )
}
