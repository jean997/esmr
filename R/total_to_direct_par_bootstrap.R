#' Parameteric bootstrap for direct effect standard errors and p-values
#'
#' @param beta Estimated total effects vector.
#' @param beta_cov Covariance matrix of estimated total effects.
#' @param beta_ix Indices of total effects in the beta vector.
#' @param d Number of traits.
#' @param bootstrap_samples Number of bootstrap samples to use.
#' @return A list with direct effect standard errors and bootstrap means.
total_to_direct_parameteric_bootstrap <- function(
  beta,
  beta_cov,
  beta_ix,
  d,
  bootstrap_samples = 10000
) {
    # Parameteric bootstrap for direct effect standard errors and p-values
    resample_beta <- MASS::mvrnorm(
      n = bootstrap_samples,
      mu = beta,
      Sigma = beta_cov
    )

    resample_results <- apply(resample_beta, 1, function(x) {
      total_effects <- matrix(0, nrow = d, ncol = d)
      total_effects[beta_ix] <- x
      # Unclear what to do if SR is > 1 ...
      sr <- spectral_radius(total_effects)
      direct_effects <- total_to_direct(total_effects, restrict_dag = FALSE)
      list(
        total_effects = total_effects,
        direct_effects = direct_effects,
        spectral_radius = sr
      )
    })

    direct_effects_mat <- sapply(resample_results, function(x) as.numeric(x$direct_effects))
    direct_effects_se <- matrix(apply(direct_effects_mat, 1, sd), nrow = ncol(G), ncol = ncol(G))
    direct_effects_bootstrap <- matrix(apply(direct_effects_mat, 1, mean), nrow = ncol(G), ncol = ncol(G))
    list(
      direct_effects_se = direct_effects_se,
      direct_effects_bootstrap = direct_effects_bootstrap
    )
}