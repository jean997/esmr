#' Initial NESMR graph estimate using n MVMR estimates
#'
#' @export
nesmr_complete_mvmr <- function(
    beta_hat_X, se_X,
    pval_select = NULL,
    alpha = 5e-8,
    ...
  ) {

  stopifnot(all(dim(beta_hat_X) == dim(se_X)))
  d <- ncol(beta_hat_X)
  stopifnot(d > 1)

  if (is.null(pval_select)) {
    Z_cursed <- with(dat, beta_hat/s_estimate)
    pval_cursed <- 2 * pnorm(-abs(Z_cursed))
    pval_select <- pval_cursed
  }

  MVMR_models <- lapply(seq_len(d), function(i) {
    mvmr_minp <- apply(pval_select[,-i], 1, min)
    mvmr_ix <- which(mvmr_minp < alpha)

    # Estimate G at each step for fair comparison
    with(dat,
         esmr(beta_hat_Y = beta_hat[,i],
              se_Y = s_estimate[,i],
              beta_hat_X = beta_hat[,-i],
              se_X = s_estimate[,-i],
              variant_ix = mvmr_ix,
              G = NULL,
              beta_joint = TRUE,
              ...)
    )
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
    beta = adj_mat_beta,
    se = mvmr_se
  ))
}
