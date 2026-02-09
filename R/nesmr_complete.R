#' Perform leave-one-out MVMR, treating each trait as the outcome in turn and all other traits as exposures.
#'
#' This can be used to get an initial estimate of the graph structure for NESMR, or as a standalone method for performing network MR with a large number of traits.
#'
#' @param beta_hat Matrix of SNP-trait associations (n by p)
#' @param se_beta_hat Matrix of standard errors of beta_hat
#' @param pval_select Matrix of p-values for variant selection. If NULL, p-values will be calculated from beta_hat and se_beta_hat.
#' @param alpha P-value threshold for variant selection.
#' @param lower_tri If TRUE, only use traits with higher index as exposures for each outcome. This will only estimate the lower triangular part of the matrix representing the graph structure.
#' @param R Nuisance correlation matrix if there is sample overlap (optional).
#' @export
nesmr_complete_mvmr <- function(
    beta_hat, se_beta_hat,
    pval_select = NULL,
    alpha = 5e-8,
    lower_tri = FALSE,
    R = NULL,
    ...
  ) {
  stopifnot(all(dim(beta_hat) == dim(se_beta_hat)))
  d <- ncol(beta_hat)
  stopifnot(d > 1)

  if (is.null(pval_select)) {
    Z_cursed <- beta_hat/se_beta_hat
    pval_cursed <- 2 * pnorm(-abs(Z_cursed))
    pval_select <- pval_cursed
  }

  ivs <- pval_select < alpha
  n_ivs <- colSums(ivs)

  valid_idx <- which(n_ivs > 0)

  MVMR_models <- lapply(seq_len(d), function(i) {
    x_idx <- if (lower_tri) {
      intersect(which(seq_len(d) > i), valid_idx)
    } else {
      intersect(seq_len(d)[-i], valid_idx)
    }

    if (! is.null(R)) {
      R_sub <- R[c(i, x_idx), c(i, x_idx), drop = FALSE]
    } else {
      R_sub <- NULL
    }

    mvmr_minp <- apply(pval_select[,x_idx, drop = FALSE], 1, min)
    mvmr_ix <- which(mvmr_minp < alpha)

    # Estimate G at each step for fair comparison
    tryCatch({
        capture.output({
          mod_res <- esmr(beta_hat_Y = beta_hat[,i],
                      se_Y = se_beta_hat[,i],
                      beta_hat_X = beta_hat[,x_idx],
                      se_X = se_beta_hat[,x_idx],
                      variant_ix = mvmr_ix,
                      R = R_sub,
                      params = list(
                        strict_mode = FALSE
                      ),
                      ...)
          }, file = nullfile())
        rs <- mod_res$remove_suggest
        low_info_flag <- !is.null(rs)
        while(low_info_flag) { # Note: This will be caught if no traits remain
          x_idx <- x_idx[- (rs - 1)]
          # x_idx <- setdiff(x_idx, mod_res$remove_suggest)
          if (length(x_idx) == 0) {
            stop("All traits were suggested for removal due to low information.")
          }
          capture.output({
            mod_res <- esmr(beta_hat_Y = beta_hat[,i],
                          se_Y = se_beta_hat[,i],
                          beta_hat_X = beta_hat[,x_idx],
                          se_X = se_beta_hat[,x_idx],
                          variant_ix = mvmr_ix,
                          R = R_sub,
                          params = list(
                            strict_mode = FALSE
                          ),
                          ...)
          }, file = nullfile())
          rs <- mod_res$remove_suggest
          low_info_flag <- !is.null(rs)
        }
        # Update the indices in the final output to reflect any removed traits
        mod_res$beta$beta_k <- x_idx
        list(beta = mod_res$beta) # Only return the beta component to save memory
      }, error = function(e) {
        warning(e)
        list(beta = data.frame(
          beta_m = rep(0, length(x_idx)),
          beta_s = rep(0, length(x_idx)),
          beta_j = rep(1, length(x_idx)),
          beta_k = seq_along(x_idx)
        ))
      })
  })

  # Combine the results into a single matrix
  # Create matrix from the effects
  mvmr_beta_df <- do.call(
    'rbind.data.frame',
    lapply(seq_along(MVMR_models), function(i) {
      x <- MVMR_models[[i]]

      res <- x$beta[c('beta_m', 'beta_s')]
      beta_to <- as.numeric(x$beta$beta_j)
      beta_from <- as.numeric(x$beta$beta_k)
      res$to <- rep(i, length(beta_to))
      res$from <- beta_from
      res
    })
  )

  mvmr_beta_edgelist <- mvmr_beta_df[
    , c("from", "to", "beta_m", "beta_s")]

  adj_mat_beta <- matrix(0, nrow = d, ncol = d)
  adj_mat_beta[as.matrix(mvmr_beta_edgelist[, 1:2])] <- mvmr_beta_edgelist$beta_m

  mvmr_se <- matrix(0, nrow = d, ncol = d)
  mvmr_se[as.matrix(mvmr_beta_edgelist[, 1:2])] <- mvmr_beta_edgelist$beta_s

  return(list(
    beta_hat = adj_mat_beta,
    se_beta_hat = mvmr_se
  ))
}


#' Initial NESMR graph estimate using n MVMR estimates
#'
#' @export
nesmr_complete <- function(
    beta_hat, se_beta_hat,
    pval_select = NULL,
    alpha = 5e-8,
    R = NULL,
    ...
) {
  stopifnot(all(dim(beta_hat) == dim(se_beta_hat)))
  d <- ncol(beta_hat)
  stopifnot(d > 1)

  if (is.null(pval_select)) {
    Z_cursed <- beta_hat/se_beta_hat
    pval_cursed <- 2 * pnorm(-abs(Z_cursed))
    pval_select <- pval_cursed
  }

  minp <- apply(pval_select, 1, min)
  ix <- which(minp < alpha)

  B_full <- matrix(1, ncol = d, nrow = d) - diag(d)

  nesmr_full <- nesmr(
    beta_hat_X = beta_hat,
    se_X = se_beta_hat,
    variant_ix = ix,
    direct_effect_template = B_full,
    restrict_dag = FALSE,
    R = R,
    ...)

  return(list(
    beta_hat = nesmr_full$beta_mat$beta_hat,
    se_beta_hat = nesmr_full$beta_mat$beta_se
  ))
}
