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
  n <- nrow(total_effects)
  I <- diag(n)
  free <- matrix(TRUE, n, n)
  diag(free) <- FALSE
  if (nrow(s) > 0) free[cbind(s[, "row"], s[, "col"])] <- FALSE

  tgt <- total_effects * free
  B <- tgt
  for (it in seq_len(max(n - 1, 1))) {
    higher <- solve(I - B) - I - B
    B_new <- (tgt - higher) * free
    done <- all(B_new == B)
    B <- B_new
    if (done) break
  }

  list(
    B = B,
    total_effects = solve(I - B) - I
  )
}
