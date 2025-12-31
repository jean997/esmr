#' Parameteric bootstrap for direct effect standard errors and p-values
#'
#' @param x Either a numeric vector of total effects (beta) or an nesmr object.
#' @param beta_cov Covariance matrix of estimated total effects (required if x is numeric).
#' @param beta_ix Indices of total effects in the beta vector (required if x is numeric).
#' @param d Number of traits (required if x is numeric).
#' @param bootstrap_samples Number of bootstrap samples to use.
#' @return A list with direct effect standard errors and bootstrap means.
#' @export
total_to_direct_parameteric_bootstrap <- function(
  x,
  beta_cov = NULL,
  beta_ix = NULL,
  d = NULL,
  bootstrap_samples = 1000
) {
  UseMethod("total_to_direct_parameteric_bootstrap")
}

#' @rdname total_to_direct_parameteric_bootstrap
#' @export
total_to_direct_parameteric_bootstrap.numeric <- function(
  x,
  beta_cov = NULL,
  beta_ix = NULL,
  d = NULL,
  bootstrap_samples = 10000
) {
  if (is.null(beta_cov) || is.null(beta_ix) || is.null(d)) {
    stop("beta_cov, beta_ix, and d are required when x is numeric")
  }

  beta <- x

  # Parameteric bootstrap for direct effect standard errors and p-values
  resample_beta <- MASS::mvrnorm(
    n = bootstrap_samples,
    mu = beta,
    Sigma = beta_cov
  )

  resample_results <- apply(resample_beta, 1, function(y) {
    total_effects <- matrix(0, nrow = d, ncol = d)
    total_effects[beta_ix] <- y
    # Unclear what to do if SR is > 1 ...
    sr <- spectral_radius(total_effects)
    direct_effects <- total_to_direct(total_effects, restrict_dag = FALSE)
    list(
      total_effects = total_effects,
      direct_effects = direct_effects,
      spectral_radius = sr
    )
  })

  direct_effects_mat <- sapply(resample_results, function(z) as.numeric(z$direct_effects))
  direct_effects_se <- matrix(apply(direct_effects_mat, 1, sd), nrow = d, ncol = d)
  direct_effects_bootstrap <- matrix(apply(direct_effects_mat, 1, mean), nrow = d, ncol = d)
  list(
    direct_effects_se = direct_effects_se,
    direct_effects_bootstrap = direct_effects_bootstrap
  )
}

#' @rdname total_to_direct_parameteric_bootstrap
#' @export
total_to_direct_parameteric_bootstrap.nesmr <- function(
  x,
  beta_cov = NULL,
  beta_ix = NULL,
  d = NULL,
  bootstrap_samples = 1000
) {
  # Extract parameters from nesmr object
  dat <- x

  beta <- dat$beta$beta_m
  beta_cov <- dat$beta$V
  beta_ix <- cbind(dat$beta$beta_k, dat$beta$beta_j)
  d <- ncol(dat$total_effects)

  # Call the numeric method
  total_to_direct_parameteric_bootstrap.numeric(
    x = beta,
    beta_cov = beta_cov,
    beta_ix = beta_ix,
    d = d,
    bootstrap_samples = bootstrap_samples
  )
}