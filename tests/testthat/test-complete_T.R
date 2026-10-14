# Package use: place in tests/testthat/test-complete_T.R.
# Standalone use, with complete_T.R beside this file:
# topo_sort_mat <- esmr:::topo_sort_mat
# source("complete_T.R")
# testthat::test_file("test_complete_T.R")
# Package testing supplies complete_T through the package test environment.

local_edition(3)

# User's original function, unchanged. Its tests use n >= 2 and lower s.
complete_T_lower_reference <- function(total_effects, s) {
  total_effects_complete <- total_effects
  n <- nrow(total_effects)
  X <- diag(n)
  X_inv <- (diag(n) + total_effects)
  s <- s[order(s[, 'row'], -s[, 'col']), , drop = FALSE]
  s_idx <- 1
  s_max <- nrow(s)
  for (i in 2:n) {
    for (j in (i - 1):1) {
      if ((j + 1) <= (i - 1)) {
        partial_sum <- -sum(
          X_inv[i, (j + 1):(i - 1)] * X[(j + 1):(i - 1), j]
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
  list(B = diag(n) - X, total_effects = X_inv - diag(n))
}

indices <- function(mask) {
  s <- which(mask, arr.ind = TRUE)
  colnames(s) <- c("row", "col")
  s
}
totals <- function(B) solve(diag(nrow(B)) - B) - diag(nrow(B))
verify <- function(fit, input, s, true_B = NULL) {
  n <- nrow(input)
  known <- matrix(TRUE, n, n)
  diag(known) <- FALSE
  if (nrow(s)) known[s] <- FALSE
  expect_equal(unname(fit$total_effects[known]), unname(input[known]), tolerance = 1e-10)
  expect_equal(fit$total_effects, totals(fit$B), tolerance = 1e-10)
  if (nrow(s)) expect_equal(unname(fit$B[s]), rep(0, nrow(s)), tolerance = 1e-10)
  expect_equal(unname(diag(fit$B)), rep(0, n), tolerance = 1e-10)
  if (!is.null(true_B)) expect_equal(unname(fit$B), unname(true_B), tolerance = 1e-10)
}

test_that("complete_T matches the original triangular recurrence before and after permutation", {
  withr::local_seed(20261014)

  # Direct comparisons with the original on random ordered DAGs. Only some
  # zero edges are constrained: other zero edges remain known zero totals.
  for (n in c(2L, 3L, 5L, 10L, 25L)) {
    for (replicate_id in seq_len(40L)) {
      B <- matrix(0, n, n)
      B[lower.tri(B)] <- rnorm(n * (n - 1L) / 2L, sd = 0.2)
      B[lower.tri(B) & matrix(runif(n * n) < 0.65, n, n)] <- 0
      T <- totals(B)
      constrained <- lower.tri(B) & B == 0 &
        matrix(runif(n * n) < 0.8, n, n)
      s <- indices(constrained)
      if (nrow(s)) s <- s[sample.int(nrow(s)), , drop = FALSE]
      input <- T
      input[s] <- rnorm(nrow(s), sd = 100) # Ignored placeholders.
      ref <- complete_T_lower_reference(input, s)
      fit <- complete_T(input, s)
      expect_equal(fit$B, ref$B, tolerance = 1e-10)
      expect_equal(fit$total_effects, ref$total_effects, tolerance = 1e-10)
      verify(fit, input, s, B)

      # Permute input and constraints, then compare to the original's result.
      p <- sample.int(n)
      inverse_p <- order(p)
      sp <- s
      if (nrow(sp)) {
        sp[, "row"] <- inverse_p[s[, "row"]]
        sp[, "col"] <- inverse_p[s[, "col"]]
      }
      permuted <- input[p, p, drop = FALSE]
      permuted[sp] <- NA_real_ # Ignored missing values must not propagate.
      fp <- complete_T(permuted, sp)
      expect_equal(fp$B[inverse_p, inverse_p], ref$B, tolerance = 1e-10)
      expect_equal(fp$total_effects[inverse_p, inverse_p], ref$total_effects, tolerance = 1e-10)
      verify(fp, permuted, sp, B[p, p, drop = FALSE])
    }
  }
})

test_that("complete_T recovers nonzero direct effects when total effects cancel", {
  # Exact cancellation: a nonzero direct effect has a zero known total effect.
  # This tests why merely treating nonzero totals as the DIRECT support is wrong.
  B <- matrix(0, 4, 4)
  B[2, 1] <- 2
  B[3, 2] <- 3
  B[3, 1] <- -6
  B[4, 3] <- 0.5
  T <- totals(B)
  expect_equal(T[3, 1], 0)
  expect_true(B[3, 1] != 0)
  s <- indices(lower.tri(B) & B == 0)
  for (p in list(1:4, c(4, 2, 1, 3), c(3, 1, 4, 2))) {
    ip <- order(p)
    sp <- s
    sp[, "row"] <- ip[s[, "row"]]
    sp[, "col"] <- ip[s[, "col"]]
    verify(complete_T(T[p, p], sp), T[p, p], sp, B[p, p])
  }
})

test_that("complete_T completes constrained effects along long mediated paths", {
  # Long mediated path: constrained total effects must be completed recursively.
  n <- 12L
  B <- matrix(0, n, n)
  B[cbind(2:n, 1:(n - 1L))] <- seq(0.2, 0.7, length.out = n - 1L)
  s <- indices(B == 0 & row(B) != col(B))
  T <- totals(B)
  T[s] <- NA_real_
  verify(complete_T(T, s), T, s, B)
})

test_that("complete_T handles disconnected graphs and preserves trait names", {
  # No constraints, disconnected components, negative effects, and named axes.
  B <- matrix(0, 6, 6)
  B[2, 1] <- -0.7
  B[6, 4] <- 0.8
  B[3, 2] <- 0.3
  p <- c(3, 6, 1, 5, 2, 4)
  T <- totals(B)[p, p]
  dimnames(T) <- list(paste0("trait", p), paste0("trait", p))
  fit <- complete_T(T)
  verify(fit, T, indices(matrix(FALSE, 6, 6)), B[p, p])
  expect_identical(dimnames(fit$B), dimnames(T))
  expect_identical(dimnames(fit$total_effects), dimnames(T))
})

test_that("complete_T ignores missing placeholders when every edge is constrained", {
  # Every edge constrained; even ignored diagonal values may be missing.
  T <- matrix(NA_real_, 5, 5)
  s <- indices(row(T) != col(T))
  fit <- complete_T(T, s)
  expect_equal(fit$B, matrix(0, 5, 5))
  expect_equal(fit$total_effects, matrix(0, 5, 5))
})

test_that("complete_T handles singleton inputs and empty constraint matrices", {
  # Singleton and empty constraint matrices.
  fit <- complete_T(matrix(NA_real_, 1, 1))
  expect_equal(fit$B, matrix(0, 1, 1))
  expect_equal(fit$total_effects, matrix(0, 1, 1))
  expect_equal(complete_T(matrix(0, 2, 2), indices(matrix(FALSE, 2, 2)))$B,
               matrix(0, 2, 2))
})

test_that("complete_T accepts duplicate constraints and reversed index columns", {
  # Duplicate constraints and reversed column order are harmless.
  s <- rbind(c(row = 3, col = 1), c(row = 3, col = 1))
  T <- matrix(0, 3, 3)
  T[2, 1] <- 0.3
  T[3, 2] <- -0.4
  fit <- complete_T(T, s[, c("col", "row"), drop = FALSE])
  expect_equal(fit$total_effects[3, 1], -0.12, tolerance = 1e-14)
})

test_that("complete_T preserves tiny nonzero ordering information", {
  # Do not threshold away tiny but genuine ordering information.
  T <- matrix(0, 3, 3)
  T[3, 2] <- 1e-14
  T[1, 3] <- 0.4
  s <- rbind(c(row = 1, col = 2))
  fit <- complete_T(T, s)
  expect_gt(fit$total_effects[1, 2], 0)
  expect_lt(abs(fit$total_effects[1, 2] - 4e-15), 1e-29)
})

test_that("complete_T rejects cycles in known total effects", {
  # Cyclic systems are outside the DAG contract, even when I - B is invertible.
  expect_error(complete_T(matrix(c(0, 0.2, 0.3, 0), 2, 2)), "cycle")
  T <- matrix(0, 3, 3)
  T[cbind(c(2, 3, 1), c(1, 2, 3))] <- 0.1
  expect_error(complete_T(T), "cycle")
})

