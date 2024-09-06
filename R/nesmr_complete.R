#' Initial NESMR graph estimate using n MVMR estimates
#'
#' @export
nesmr_complete_mvmr <- function(
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

  MVMR_models <- lapply(seq_len(d), function(i) {
    mvmr_minp <- apply(pval_select[,-i], 1, min)
    mvmr_ix <- which(mvmr_minp < alpha)

    # Estimate G at each step for fair comparison
    esmr(beta_hat_Y = beta_hat[,i],
              se_Y = se_beta_hat[,i],
              beta_hat_X = beta_hat[,-i],
              se_X = se_beta_hat[,-i],
              variant_ix = mvmr_ix,
              beta_joint = TRUE,
              ...)
  })

  # Combine the results into a single matrix
  # Create matrix from the effects
  mvmr_beta_df <- do.call(
    'rbind.data.frame',
    lapply(1:d, function(i) {
      x <- MVMR_models[[i]]

      res <- x$beta[c('beta_m', 'beta_s')]
      res$to <- rep(i, d - 1)
      res$from <- setdiff(1:d, i)
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
    ...)

  nesmr_beta <- t(nesmr_full$f$fgbar) - diag(d)
  nesmr_se <- sqrt(t(nesmr_full$f$fg2bar) - t(nesmr_full$f$fgbar^2))

  return(list(
    beta_hat = nesmr_beta,
    se_beta_hat = nesmr_se
  ))
}
