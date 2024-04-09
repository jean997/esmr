#' Ma et al. (2023) Winner's curse adjusted SNP effect size estimate
#'
#' Selected SNP are selected at |beta / se_beta + Z| > lambda; Z ~ N(0, sd = eta)
#'
#' DOI: 10.1214/22-AOS2247
#'
#' @param beta Selected SNP effect size estimates (numeric or matrix)
#' @param se_beta Selected SNP effect size standard error (numeric or matrix)
#' @param eta a prespecified constant that reflects the noise level of the pseudo SNPs
#' @param alpha significance level.
#' @param lambda quantile of the standard normal distribution. If specified, overrides `alpha`.
#'
#' @return list with elements: beta_rb, se_rb if `est_se_beta` = TRUE and only a vector of Rao-Blackwellized estimates if `est_se_beta` = FALSE
#' @export
#'
#' @examples
#' alpha <- 5e-5
#' lambda <- qnorm(1 - alpha / 2)
#' sim_beta1 <- rnorm(1e5, lambda, 1)
#' sim_beta1_rand <- abs(sim_beta1 + rnorm(1e5, 0, eta))
#' res_beta_rb <- t(mapply(
#'  snp_beta_rb, beta = sim_beta1[sim_beta1_rand > lambda], se_beta = 1))
#'
#' hist(as.numeric(res_beta_rb[, 1]))
#' # Add the naive selection
#' hist(sim_beta1[sim_beta1 > lambda], add = TRUE, col = 'red')
#'
snp_beta_rb <- function(
  beta, se_beta,
  eta = 0.5,
  alpha = 5e-8,
  lambda = qnorm(1 - alpha / 2),
  est_se_beta = TRUE
  ) {
    a_term_plus <- - beta / (se_beta * eta) + lambda / eta
    a_term_neg <- - beta / (se_beta * eta) - lambda / eta
    num <- dnorm(a_term_plus) - dnorm(a_term_neg)
    denom <- 1 - pnorm(a_term_plus) + pnorm(a_term_neg)
    ratio_term <- (num / denom)
    beta_rb <- beta - (se_beta / eta) * ratio_term

    if (!est_se_beta) {
      return(beta_rb)
    } else {
      var_term1 <- ratio_term^2
      var_term2 <- (a_term_plus * dnorm(a_term_plus) -
                      a_term_neg * dnorm(a_term_neg)) / denom
      rb_var <- se_beta^2 * (1 + (var_term1 - var_term2) / eta^2)

      return(list(beta_rb = beta_rb, se_rb = sqrt(rb_var)))
    }
}

