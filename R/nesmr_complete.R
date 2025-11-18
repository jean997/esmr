#' Initial NESMR graph estimate using n MVMR estimates
#'
#' @export
nesmr_complete_mvmr <- function(
    beta_hat, se_beta_hat,
    pval_select = NULL,
    alpha = 5e-8,
    lower_tri = FALSE,
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

  MVMR_models <- lapply(seq_len(d - 1), function(i) {
    x_idx <- if (lower_tri) {
      which(seq_len(d) > i)
    } else {
      seq_len(d)[-i]
    }

    mvmr_minp <- apply(pval_select[,x_idx, drop = FALSE], 1, min)
    mvmr_ix <- which(mvmr_minp < alpha)

    # Estimate G at each step for fair comparison
    tryCatch({
      esmr(beta_hat_Y = beta_hat[,i],
                    se_Y = se_beta_hat[,i],
                    beta_hat_X = beta_hat[,x_idx],
                    se_X = se_beta_hat[,x_idx],
                    variant_ix = mvmr_ix,
                    ...)
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

      x_idx <- if (lower_tri) {
        which(seq_len(d) > i)
      } else {
        seq_len(d)[-i]
      }

      res <- x$beta[c('beta_m', 'beta_s')]
      beta_to <- as.numeric(x$beta$beta_j)
      beta_from <- as.numeric(x$beta$beta_k)
      res$to <- rep(i, length(beta_to))
      res$from <- c(i, x_idx)[beta_from]
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

  nesmr_full <- esmr(
    beta_hat_X = beta_hat,
    se_X = se_beta_hat,
    variant_ix = ix,
    G = diag(d), # required for network problem
    direct_effect_template = B_full,
    restrict_dag = FALSE,
    ...)

  return(list(
    beta_hat = nesmr_full$beta_mat$beta_hat,
    se_beta_hat = nesmr_full$beta_mat$beta_se
  ))
}
