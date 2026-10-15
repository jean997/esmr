# Requires subset_gfa_traits() from subset_gfa_traits.R.
# Fits GFA once on all traits, then runs ESMR with each trait as outcome.
# Set workers > 1 to use future multisession workers (also works on Windows).
# seed controls GFA and parallel RNG streams; the previous future plan is restored.
nesmr_complete_mvmr_cached_gfa <- function(
    beta_hat, se_beta_hat, pval_select = NULL, alpha = 5e-8,
    lower_tri = FALSE, R = NULL, gfa_args = list(),
    params = list(strict_mode = FALSE),
    return_models = FALSE, return_gfa_fit = FALSE,
    workers = 1L, seed = NULL, ...) {

  if (!exists("subset_gfa_traits", mode = "function")) {
    stop("Source subset_gfa_traits.R before calling this function.")
  }
  if (!is.matrix(beta_hat) || !is.numeric(beta_hat) ||
      !is.matrix(se_beta_hat) || !is.numeric(se_beta_hat) ||
      !identical(dim(beta_hat), dim(se_beta_hat))) {
    stop("beta_hat and se_beta_hat must be numeric SNP-by-trait matrices of equal dimension.")
  }
  d <- ncol(beta_hat)
  if (d < 2L || nrow(beta_hat) < 1L) stop("At least two traits and one SNP are required.")
  if (any(!is.finite(se_beta_hat) | se_beta_hat <= 0)) {
    stop("se_beta_hat must contain positive, finite values.")
  }
  if (is.null(pval_select)) {
    pval_select <- 2 * stats::pnorm(-abs(beta_hat / se_beta_hat))
  }
  if (!is.matrix(pval_select) || !identical(dim(pval_select), dim(beta_hat))) {
    stop("pval_select must have the same dimensions as beta_hat.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) stop("alpha must lie in (0, 1).")
  if (!is.null(R)) {
    if (is.null(dim(R)) || !identical(dim(R), c(d, d))) {
      stop("R must be a trait-by-trait matrix in beta_hat column order.")
    }
  }
  if (!is.list(gfa_args) || any(names(gfa_args) %in% c("B_hat", "S", "R"))) {
    stop("gfa_args must be a list without B_hat, S, or R.")
  }
  if (!is.list(params)) stop("params must be a list of esmr::esmr parameters.")
  if (is.null(params$strict_mode)) params$strict_mode <- FALSE
  augment <- if (is.null(params$augment_G)) TRUE else params$augment_G
  if (!is.logical(augment) || length(augment) != 1L || is.na(augment)) {
    stop("params$augment_G must be TRUE or FALSE.")
  }

  if (!is.numeric(workers) || length(workers) != 1L || !is.finite(workers) ||
      workers < 1 || workers != floor(workers)) stop("workers must be a positive integer.")
  workers <- as.integer(min(workers, d))
  if (workers > 1L && (!requireNamespace("future", quietly = TRUE) ||
                       !requireNamespace("future.apply", quietly = TRUE))) {
    stop("Parallel execution requires packages future and future.apply.")
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
      stop("seed must be NULL or a nonnegative integer.")
    }
    set.seed(as.integer(seed))
  }
  esmr_args <- list(...)
  reserved <- c("beta_hat_Y", "se_Y", "beta_hat_X", "se_X", "G", "R", "variant_ix", "params")
  if (length(esmr_args) && (is.null(names(esmr_args)) || any(!nzchar(names(esmr_args))) ||
      any(names(esmr_args) %in% reserved) || anyDuplicated(names(esmr_args)))) {
    stop("Additional ESMR arguments must be uniquely named and must not override wrapper inputs.")
  }

  # No column of the outcome's associations is used for IV selection.
  eligible <- which(colSums(pval_select < alpha, na.rm = TRUE) > 0L)

  # Exactly one GFA call. In contrast to esmr::estimate_G(), this uses ALL
  # traits, including those later treated as outcomes; see interpretation note.
  gfa <- do.call(GFA::gfa_fit, c(
    list(B_hat = beta_hat, S = se_beta_hat, R = R), gfa_args
  ))
  if (is.null(gfa$F_hat) || !is.matrix(gfa$F_hat) || nrow(gfa$F_hat) != d) {
    stop("GFA did not return a wrapped F_hat with one row per trait.")
  }

  labels <- colnames(beta_hat)
  if (is.null(labels)) labels <- paste0("trait_", seq_len(d))
  if (anyDuplicated(labels)) labels <- make.unique(labels)
  effect <- se <- matrix(0, d, d, dimnames = list(labels, labels))
  models <- if (return_models) vector("list", d) else NULL
  if (return_models) names(models) <- labels
  failures <- vector("list", d)
  names(failures) <- labels
  exposure_sets <- vector("list", d)
  names(exposure_sets) <- labels

  run_outcome <- function(i) {
    result <- list(failure = NULL, exposures = integer(0), beta = NULL, se = NULL, model = NULL)
    tryCatch({
      candidates <- if (lower_tri) which(seq_len(d) > i) else setdiff(seq_len(d), i)
      x_idx <- intersect(candidates, eligible)
      if (!length(x_idx)) {
        result$failure <- "No eligible exposures."
        return(result)
      }
      iv_idx <- which(rowSums(pval_select[, x_idx, drop = FALSE] < alpha,
                              na.rm = TRUE) > 0L)
      if (!length(iv_idx)) {
        result$failure <- "No instruments selected from exposures."
        return(result)
      }

      repeat {
        # The factor projection renormalizes F and its paired L. ESMR needs F.
        projected <- subset_fn(projection_gfa, keep = x_idx)
        Fx <- projected$F_hat
        if (length(x_idx) == 1L) {
          # ESMR's built-in two-trait special case uses the identity basis.
          Gi <- diag(2L)
        } else {
          # Matches estimate_G(..., add_trait1 = TRUE), including augmentation.
          basis <- if (augment) cbind(Fx, diag(length(x_idx))) else Fx
          Gi <- rbind(c(1, rep(0, ncol(basis))),
                      cbind(rep(0, length(x_idx)), basis))
        }
        Ri <- if (is.null(R)) NULL else R[c(i, x_idx), c(i, x_idx), drop = FALSE]
        fit <- tryCatch(
          do.call(esmr::esmr, c(list(
            beta_hat_Y = beta_hat[, i], se_Y = se_beta_hat[, i],
            beta_hat_X = beta_hat[, x_idx, drop = FALSE],
            se_X = se_beta_hat[, x_idx, drop = FALSE],
            G = Gi, R = Ri, variant_ix = iv_idx, params = params
          ), esmr_args)),
          error = function(e) e
        )
        if (inherits(fit, "error")) {
          result$failure <- conditionMessage(fit)
          break
        }
        # ESMR numbers its outcome 1 and exposures 2, ..., length(x_idx)+1.
        rs <- fit$remove_suggest
        if (is.null(rs)) {
          if (length(fit$beta$beta_m) != length(x_idx) ||
              length(fit$beta$beta_s) != length(x_idx)) {
            result$failure <- "Unexpected number of ESMR coefficients."
            break
          }
          result$beta <- fit$beta$beta_m
          result$se <- fit$beta$beta_s
          result$exposures <- x_idx
          if (return_models) result$model <- fit
          break
        }
        rs <- as.integer(rs)
        if (length(rs) != 1L || is.na(rs) || rs < 2L || rs > length(x_idx) + 1L) {
          result$failure <- "Invalid remove_suggest index returned by ESMR."
          break
        }
        x_idx <- x_idx[-(rs - 1L)]
        if (!length(x_idx)) {
          result$failure <- "All exposures were suggested for removal."
          break
        }
        # Keep iv_idx fixed, matching nesmr_complete_mvmr's retry behavior.
      }
      result
    }, error = function(e) {
      result$failure <- conditionMessage(e)
      result
    })
  }

  environment(run_outcome) <- list2env(list(
    d = d, lower_tri = lower_tri, eligible = eligible,
    beta_hat = beta_hat, se_beta_hat = se_beta_hat,
    pval_select = pval_select, alpha = alpha, R = R, augment = augment,
    params = params, esmr_args = esmr_args, return_models = return_models,
    subset_fn = get("subset_gfa_traits", mode = "function"),
    projection_gfa = gfa[c("F_hat", "L_hat", "scale", "R", "method", "params", "gfa_pve")]
  ), parent = baseenv())

  if (workers == 1L) {
    results <- lapply(seq_len(d), run_outcome)
  } else {
    previous_plan <- future::plan()
    on.exit(future::plan(previous_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)
    results <- future.apply::future_lapply(
      seq_len(d), run_outcome, future.seed = TRUE,
      future.packages = "esmr", future.scheduling = Inf
    )
  }
  # Workers return one result per outcome; assemble matrices in original order.
  for (i in seq_len(d)) {
    result <- results[[i]]
    failures[i] <- list(result$failure)
    exposure_sets[i] <- list(result$exposures)
    if (length(result$exposures)) {
      effect[result$exposures, i] <- result$beta
      se[result$exposures, i] <- result$se
    }
    if (return_models) models[i] <- list(result$model)
  }

  out <- list(beta_hat = effect, se_beta_hat = se,
              failures = failures, exposures = exposure_sets,
              gfa_F_hat = gfa$F_hat)
  if (return_models) out$models <- models
  if (return_gfa_fit) out$gfa_fit <- gfa
  out
}
