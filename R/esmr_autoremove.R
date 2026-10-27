#' Fit ESMR, automatically removing low-information exposures
#'
#' Estimates G once, then refits after each removal suggestion.
#' Arguments follow esmr(). R and RE include the outcome first.
#' The returned exposure_ix and removed_exposure_ix refer to the original X.
#' The removal suggestion is read from fit$remove_suggest (row index in G/R,
#' with the outcome in row 1). Errors from esmr() are propagated.
#'
#' Note: Should probably be used with a low cond_num in parmas as
#' @export
fit_esmr_autoremove <- function(beta_hat_X, se_X,
                                beta_hat_Y = NULL, se_Y = NULL,
                                G = NULL, R = NULL,
                                pval_thresh = NULL, variant_ix = NULL,
                                ld_scores = NULL, RE = NULL,
                                max_iter = 100, params = list()){
  beta_hat_X <- as.matrix(beta_hat_X)
  se_X <- as.matrix(se_X)
  if(is.null(beta_hat_Y)){
    beta_hat_Y <- beta_hat_X[,1]
    se_Y <- se_X[,1]
    beta_hat_X <- beta_hat_X[,-1,drop = FALSE]
    se_X <- se_X[,-1,drop = FALSE]
  }
  exposure_ix <- seq_len(ncol(beta_hat_X))
  removed_exposure_ix <- integer(0)
  params$strict_mode <- FALSE

  augment <- if(is.null(params$augment_G)) TRUE else params$augment_G
  R_X <- if(is.null(R)) NULL else R[-1,-1,drop = FALSE]
  # G is rebuilt from the exposure loadings F after each removal. F is
  # projected onto the retained exposures (columns renormalized, empty columns
  # dropped) rather than just deleting rows/columns of G.
  gfa <- NULL
  if(is.null(G)){
    if(ncol(beta_hat_X) > 1){
      gfa <- GFA::gfa_fit(B_hat = beta_hat_X, S = se_X, R = R_X)
      Fx <- gfa$F_hat
    }
  }else{
    # User-supplied G: its exposure block is projected as is (any augmentation
    # columns it contains are projected along with the rest).
    Fx <- as.matrix(G)[-1,-1,drop = FALSE]
    augment <- FALSE
  }

  repeat{
    G <- if(ncol(beta_hat_X) == 1) diag(2) else build_autoremove_G(Fx, augment)
    fit <- esmr::esmr(beta_hat_X = beta_hat_X, se_X = se_X,
                      beta_hat_Y = beta_hat_Y, se_Y = se_Y,
                      G = G, R = R, pval_thresh = pval_thresh,
                      variant_ix = variant_ix, ld_scores = ld_scores,
                      RE = RE, max_iter = max_iter, params = params)

    # remove_suggest indexes rows of G/R/RE, where row 1 is the outcome,
    # so exposure k corresponds to remove_suggest = k + 1.
    remove_suggest <- fit$remove_suggest
    if(is.null(remove_suggest)){
      fit$exposure_ix <- exposure_ix
      fit$removed_exposure_ix <- removed_exposure_ix
      return(fit)
    }

    remove <- remove_suggest - 1L
    if(!remove %in% seq_along(exposure_ix)) stop("Invalid removal suggestion.")
    if(length(exposure_ix) == 1) stop("No exposures remain after removal.")
    removed_exposure_ix <- c(removed_exposure_ix, exposure_ix[remove])
    exposure_ix <- exposure_ix[-remove]
    beta_hat_X <- beta_hat_X[,-remove,drop = FALSE]
    se_X <- se_X[,-remove,drop = FALSE]
    if(ncol(beta_hat_X) > 1){
      Fx <- if(!is.null(gfa)){
        subset_gfa_traits(gfa, keep = exposure_ix)$F_hat
      }else{
        project_autoremove_F(Fx, remove)
      }
    }
    if(!is.null(R)) R <- R[-remove_suggest,-remove_suggest,drop = FALSE]
    if(!is.null(RE)) RE <- RE[-remove_suggest,-remove_suggest,drop = FALSE]
    bpc <- params$beta_prior_cov
    if(is.matrix(bpc)){
      params$beta_prior_cov <- bpc[-remove,-remove,drop = FALSE]
    }else if(length(bpc) > 1){
      params$beta_prior_cov <- bpc[-remove]
    }
  }
}

# Rows are traits, columns are factors; outcome row/column added here.
build_autoremove_G <- function(Fx, augment){
  p <- nrow(Fx)
  basis <- if(augment) cbind(Fx, diag(p)) else Fx
  rbind(c(1, rep(0, ncol(basis))),
        cbind(rep(0, p), basis))
}

# Drop exposure `remove` from F, renormalize columns, drop empty columns.
# Matches the F_hat produced by subset_gfa_traits().
project_autoremove_F <- function(Fx, remove){
  Fs <- Fx[-remove,,drop = FALSE]
  nrm <- sqrt(colSums(Fs^2))
  active <- which(nrm > 0)
  sweep(Fs[,active,drop = FALSE], 2L, nrm[active], "/")
}
