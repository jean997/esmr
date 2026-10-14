#' Computes direct effects and fills in constrained total effects for a DAG
#'
#' Infers a topological order and applies the triangular recurrence.
#' O(n^3) time and O(n^2) memory; independent of input trait ordering.
#'
#' @param total_effects Total effect matrix. Entries in s and the diagonal
#'   are ignored; all other off-diagonal entries are treated as known.
#' @param s Matrix with named columns row and col specifying zero direct
#'   effects. NULL means no constraints.
#' @return List containing B and completed total_effects in input order.
complete_T <- function(total_effects, s = NULL) {
  if (is.null(s)) {
    s <- matrix(integer(), 0, 2, dimnames = list(NULL, c("row", "col")))
  }
  s <- s[, c("row", "col"), drop = FALSE]
  total_effects_complete <- total_effects
  diag(total_effects_complete) <- 0
  total_effects_complete[s] <- 0

  # Reorder total effects to be lower triangular.
  # Inside esmr, call its internal helper directly.
  topo_order <- topo_sort_mat(total_effects_complete)
  topo_order_inv <- order(topo_order)
  total_effects_complete <- total_effects_complete[
    topo_order, topo_order, drop = FALSE
  ]
  s[, "row"] <- topo_order_inv[s[, "row"]]
  s[, "col"] <- topo_order_inv[s[, "col"]]
  # The recurrence visits only lower-triangular positions, once each.
  s <- unique(s[s[, "row"] > s[, "col"], , drop = FALSE])

  n <- nrow(total_effects)
  X <- diag(n)
  X_inv <- diag(n) + total_effects_complete
  s <- s[order(s[, "row"], -s[, "col"]), , drop = FALSE]
  s_idx <- 1
  s_max <- nrow(s)

  if (n > 1L) {
    for (i in 2:n) {
      for (j in (i - 1):1) {
        if ((j + 1) <= (i - 1)) {
          partial_sum <- -sum(
            X_inv[i, (j + 1):(i - 1)] *
              X[(j + 1):(i - 1), j]
          )
        } else {
          partial_sum <- 0
        }

        if (s_idx <= s_max && all(s[s_idx, ] == c(i, j))) {
          s_idx <- s_idx + 1
          X[i, j] <- 0
          X_inv[i, j] <- partial_sum
        } else {
          X[i, j] <- -X_inv[i, j] + partial_sum
        }
      }
    }
  }

  B <- (diag(n) - X)[topo_order_inv, topo_order_inv, drop = FALSE]
  T <- (X_inv - diag(n))[topo_order_inv, topo_order_inv, drop = FALSE]
  dimnames(B) <- dimnames(T) <- dimnames(total_effects)
  list(B = B, total_effects = T)
}
