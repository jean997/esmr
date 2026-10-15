# Project a fitted jean997/GFA::gfa_fit result onto a subset of traits.
# This is a projection of an existing fit, not a refit on reduced data.
subset_gfa_traits <- function(gfa, keep = NULL, drop = NULL,
                              single_trait_thresh = 0.95,
                              zero_tol = 0) {
  if (is.null(gfa$F_hat) || is.null(gfa$L_hat) || is.null(gfa$scale)) {
    stop("Expected a wrapped GFA::gfa_fit object with F_hat, L_hat, and scale.")
  }
  F <- gfa$F_hat
  L <- gfa$L_hat
  if (!is.matrix(F) || !is.matrix(L) || ncol(F) != ncol(L)) {
    stop("F_hat and L_hat must be matrices with the same number of factors.")
  }
  p <- nrow(F)
  if (length(gfa$scale) != p) stop("scale does not match F_hat.")
  if (xor(is.null(keep), is.null(drop)) == FALSE) {
    stop("Supply exactly one of keep or drop.")
  }
  if (!is.numeric(single_trait_thresh) || length(single_trait_thresh) != 1L ||
      is.na(single_trait_thresh) || single_trait_thresh <= 0 ||
      single_trait_thresh >= 1) stop("single_trait_thresh must lie in (0, 1).")
  if (!is.numeric(zero_tol) || length(zero_tol) != 1L ||
      is.na(zero_tol) || zero_tol < 0) stop("zero_tol must be nonnegative.")

  trait_names <- rownames(F)
  resolve <- function(x) {
    if (is.character(x)) {
      if (is.null(trait_names) || anyDuplicated(trait_names)) {
        stop("Character trait selection requires unique rownames(F_hat).")
      }
      ix <- match(x, trait_names)
      if (anyNA(ix)) stop("Unknown trait: ", paste(x[is.na(ix)], collapse = ", "))
    } else if (is.numeric(x) && all(is.finite(x)) &&
               all(x == as.integer(x)) && all(x >= 1L & x <= p)) {
      ix <- as.integer(x)
    } else {
      stop("Trait selection must contain valid names or 1-based indices.")
    }
    if (anyDuplicated(ix)) stop("Repeated traits in selection.")
    ix
  }
  ix <- if (!is.null(keep)) resolve(keep) else setdiff(seq_len(p), resolve(drop))
  if (!length(ix)) stop("At least one trait must remain.")

  F_sub <- F[ix, , drop = FALSE]
  # The original F_hat columns have unit norm. Re-normalize on retained traits
  # and absorb the removed norm into L_hat, preserving L %*% t(F_sub).
  norm <- sqrt(colSums(F_sub^2))
  active <- which(norm > zero_tol)
  F_new <- F_sub[, active, drop = FALSE]
  L_new <- L[, active, drop = FALSE]
  if (length(active)) {
    F_new <- sweep(F_new, 2L, norm[active], "/")
    L_new <- sweep(L_new, 2L, norm[active], "*")
  }
  singleton <- if (length(active)) {
    which(apply(abs(F_new), 2L, max) > single_trait_thresh)
  } else integer(0)

  pve <- gfa$gfa_pve
  if (!is.null(pve)) {
    if (!is.null(pve$genet_var)) pve$genet_var <- pve$genet_var[ix]
    if (!is.null(pve$pve)) {
      if (!is.matrix(pve$pve) || nrow(pve$pve) != p ||
          ncol(pve$pve) != ncol(F)) {
        stop("gfa_pve$pve dimensions do not match F_hat.")
      }
      pve$pve <- pve$pve[ix, active, drop = FALSE]
    }
  }
  R_new <- gfa$R
  if (!is.null(R_new)) {
    if (is.null(dim(R_new)) || length(dim(R_new)) != 2L ||
        any(dim(R_new) != p)) stop("R dimensions do not match F_hat.")
    R_new <- R_new[ix, ix, drop = FALSE] # covariance marginal, not Schur complement
  }
  out <- list(
    L_hat = L_new, F_hat = F_new,
    F_hat_single = if (length(singleton)) F_new[, singleton, drop = FALSE] else NULL,
    F_hat_multi = F_new[, setdiff(seq_len(ncol(F_new)), singleton), drop = FALSE],
    gfa_pve = pve, scale = gfa$scale[ix], R = R_new,
    method = gfa$method, params = gfa$params,
    trait_index = ix, factor_index = active,
    # The original flashier fit, posterior sampler, and likelihood refer to ALL
    # original traits, including fixed factors from the eigendecomposition of R.
    fit = NULL
  )
  class(out) <- "gfa_trait_projection"
  out
}